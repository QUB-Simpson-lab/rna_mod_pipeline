from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import stats

from rna_mod_pipeline.config import KNOWN_RELATED, canonical_rbp

from .config import REGIONS
from .fdr import (
    add_fdr_columns,
    add_region_fdr_columns,
    bh_adjust_declared,
)
from .influence import leave_one_gene_out_diagnostics
from .opportunities import OpportunityDesign


@dataclass(frozen=True)
class StratifiedComponents:
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    d: np.ndarray
    total: np.ndarray
    mh_numerator_by_stratum: np.ndarray
    mh_denominator_by_stratum: np.ndarray
    score_by_stratum: np.ndarray
    score_variance_by_stratum: np.ndarray


def stratified_components(
    design: OpportunityDesign, exposure: np.ndarray
) -> StratifiedComponents:
    exposure = np.asarray(exposure, dtype=bool)
    if len(exposure) != design.n_observations:
        raise ValueError(
            f"Exposure length {len(exposure)} does not match "
            f"{design.n_observations} opportunities"
        )
    case_exposure = exposure[design.is_case]
    control_exposure = exposure[~design.is_case]
    a = np.bincount(
        design.case_strata[case_exposure], minlength=design.n_strata
    ).astype(float)
    c = np.bincount(
        design.control_strata[control_exposure], minlength=design.n_strata
    ).astype(float)
    b = design.case_counts.astype(float) - a
    d = design.control_counts.astype(float) - c
    total = a + b + c + d
    numerator = np.divide(
        a * d, total, out=np.zeros_like(total), where=total > 0
    )
    denominator = np.divide(
        b * c, total, out=np.zeros_like(total), where=total > 0
    )
    exposed = a + c
    unexposed = b + d
    expected_a = np.divide(
        (a + b) * exposed, total, out=np.zeros_like(total), where=total > 0
    )
    score = a - expected_a
    score_variance = np.divide(
        (a + b) * (c + d) * exposed * unexposed,
        total**2 * (total - 1),
        out=np.zeros_like(total),
        where=total > 1,
    )
    return StratifiedComponents(
        a=a,
        b=b,
        c=c,
        d=d,
        total=total,
        mh_numerator_by_stratum=numerator,
        mh_denominator_by_stratum=denominator,
        score_by_stratum=score,
        score_variance_by_stratum=score_variance,
    )


def _mh_model_log_or_se(components: StratifiedComponents) -> float:
    a, b, c, d, n = (
        components.a,
        components.b,
        components.c,
        components.d,
        components.total,
    )
    ad = a * d
    bc = b * c
    apd = a + d
    adns = float(np.sum(ad / n))
    bcns = float(np.sum(bc / n))
    if adns <= 0 or bcns <= 0:
        return math.nan
    variance = float(np.sum(apd * ad / n**2)) / adns**2
    middle = apd * bc / n**2 + (1 - apd / n) * ad / n
    variance += float(np.sum(middle)) / (adns * bcns)
    variance += float(np.sum((1 - apd / n) * bc / n)) / bcns**2
    return math.sqrt(max(0.0, variance / 2))


def _gene_cluster_log_or_se(
    numerator_by_stratum: np.ndarray,
    denominator_by_stratum: np.ndarray,
    stratum_gene_codes: np.ndarray,
    n_genes: int,
) -> tuple[float, int, np.ndarray, np.ndarray, np.ndarray]:
    numerator = float(numerator_by_stratum.sum())
    denominator = float(denominator_by_stratum.sum())
    gene_numerator = np.bincount(
        stratum_gene_codes,
        weights=numerator_by_stratum,
        minlength=n_genes,
    ).astype(float)
    gene_denominator = np.bincount(
        stratum_gene_codes,
        weights=denominator_by_stratum,
        minlength=n_genes,
    ).astype(float)
    contributing = (gene_numerator > 0) | (gene_denominator > 0)
    count = int(contributing.sum())
    if numerator <= 0 or denominator <= 0 or count <= 1:
        return (
            math.nan,
            count,
            gene_numerator,
            gene_denominator,
            np.full(n_genes, np.nan),
        )
    influence = gene_numerator / numerator - gene_denominator / denominator
    variance = count / (count - 1) * float(
        np.sum(influence[contributing] ** 2)
    )
    return (
        math.sqrt(max(0.0, variance)),
        count,
        gene_numerator,
        gene_denominator,
        influence,
    )


def _subset_design(
    design: OpportunityDesign, exposure: np.ndarray, region: str
) -> tuple[OpportunityDesign, np.ndarray]:
    keep_strata = design.stratum_regions == region
    indexes = np.flatnonzero(keep_strata)
    if not len(indexes):
        raise ValueError(f"No strata exist for region {region}")
    keep_rows = keep_strata[design.stratum_codes]
    old_to_new = np.full(design.n_strata, -1, dtype=np.int32)
    old_to_new[indexes] = np.arange(len(indexes), dtype=np.int32)
    codes = old_to_new[design.stratum_codes[keep_rows]]
    cases = design.is_case[keep_rows]
    old_gene_codes = design.stratum_gene_codes[indexes]
    unique_gene_codes = np.unique(old_gene_codes)
    remap = {
        int(old): new for new, old in enumerate(unique_gene_codes)
    }
    gene_codes = np.asarray(
        [remap[int(code)] for code in old_gene_codes], dtype=np.int32
    )
    subset = OpportunityDesign(
        chroms=design.chroms[keep_rows],
        positions=design.positions[keep_rows],
        strands=design.strands[keep_rows],
        coverages=design.coverages[keep_rows],
        fractions=design.fractions[keep_rows],
        is_case=cases,
        stratum_codes=codes,
        case_strata=codes[cases],
        control_strata=codes[~cases],
        case_counts=design.case_counts[indexes],
        control_counts=design.control_counts[indexes],
        stratum_transcripts=design.stratum_transcripts[indexes],
        stratum_gene_ids=design.stratum_gene_ids[indexes],
        stratum_genes=design.stratum_genes[indexes],
        stratum_regions=design.stratum_regions[indexes],
        stratum_gene_codes=gene_codes,
        genes=design.genes[unique_gene_codes],
        gene_names=design.gene_names[unique_gene_codes],
        case_site_ids=np.array([], dtype=np.int64),
    )
    return subset, np.asarray(exposure, dtype=bool)[keep_rows]


def _effect_direction(odds_ratio: float) -> str:
    if math.isnan(odds_ratio):
        return "not_estimable"
    if odds_ratio > 1:
        return "enriched"
    if odds_ratio < 1:
        return "depleted"
    return "no_effect"


def _crude_fisher(
    design: OpportunityDesign, exposure: np.ndarray
) -> tuple[int, int, float, float]:
    case_overlaps = int(exposure[design.is_case].sum())
    control_overlaps = int(exposure[~design.is_case].sum())
    case_nonoverlaps = design.n_cases - case_overlaps
    control_nonoverlaps = design.n_controls - control_overlaps
    if (
        case_overlaps + control_overlaps == 0
        or case_nonoverlaps + control_nonoverlaps == 0
    ):
        return case_overlaps, control_overlaps, math.nan, math.nan
    odds_ratio, pvalue = stats.fisher_exact(
        [
            [case_overlaps, case_nonoverlaps],
            [control_overlaps, control_nonoverlaps],
        ],
        alternative="two-sided",
    )
    return case_overlaps, control_overlaps, float(odds_ratio), float(pvalue)


def analyse_exposure(
    design: OpportunityDesign,
    exposure: np.ndarray,
    analysis_id: str,
    modification: str,
    database: str,
    rbp: str,
    n_source_intervals: int,
    partial_input: bool,
    region_filter: str | None = None,
) -> tuple[dict, list[dict]]:
    if region_filter is None:
        working_design = design
        working_exposure = np.asarray(exposure, dtype=bool)
    else:
        working_design, working_exposure = _subset_design(
            design, exposure, region_filter
        )

    components = stratified_components(working_design, working_exposure)
    numerator = float(components.mh_numerator_by_stratum.sum())
    denominator = float(components.mh_denominator_by_stratum.sum())
    if numerator > 0 and denominator > 0:
        mh_or = numerator / denominator
        log_or = math.log(mh_or)
        estimable, boundary = True, False
    elif numerator == 0 and denominator > 0:
        mh_or, log_or = 0.0, -math.inf
        estimable, boundary = False, True
    elif numerator > 0 and denominator == 0:
        mh_or, log_or = math.inf, math.inf
        estimable, boundary = False, True
    else:
        mh_or, log_or = math.nan, math.nan
        estimable, boundary = False, False

    score = float(components.score_by_stratum.sum())
    score_variance = float(components.score_variance_by_stratum.sum())
    if score_variance > 0:
        cmh_statistic = score**2 / score_variance
        cmh_pvalue = float(stats.chi2.sf(cmh_statistic, 1))
    else:
        cmh_statistic = cmh_pvalue = math.nan

    model_se = _mh_model_log_or_se(components)
    if estimable and np.isfinite(model_se) and model_se > 0:
        critical = float(stats.norm.ppf(0.975))
        model_lower = math.exp(log_or - critical * model_se)
        model_upper = math.exp(log_or + critical * model_se)
    else:
        model_lower = model_upper = math.nan

    (
        cluster_se,
        n_clusters,
        gene_numerator,
        gene_denominator,
        gene_influence,
    ) = _gene_cluster_log_or_se(
        components.mh_numerator_by_stratum,
        components.mh_denominator_by_stratum,
        working_design.stratum_gene_codes,
        working_design.n_genes,
    )
    cluster_df = max(n_clusters - 1, 0)
    if (
        estimable
        and np.isfinite(cluster_se)
        and cluster_se > 0
        and cluster_df >= 1
    ):
        cluster_statistic = log_or / cluster_se
        cluster_pvalue = float(
            2 * stats.t.sf(abs(cluster_statistic), df=cluster_df)
        )
        critical = float(stats.t.ppf(0.975, df=cluster_df))
        cluster_lower = math.exp(log_or - critical * cluster_se)
        cluster_upper = math.exp(log_or + critical * cluster_se)
    else:
        cluster_statistic = cluster_pvalue = math.nan
        cluster_lower = cluster_upper = math.nan

    case_overlaps, control_overlaps, crude_or, crude_pvalue = _crude_fisher(
        working_design, working_exposure
    )

    informative = components.score_variance_by_stratum > 0
    informative_genes = np.unique(
        working_design.stratum_gene_codes[informative]
    )
    direction = _effect_direction(mh_or)
    (
        maximum_influence,
        stable,
        reversals,
        boundaries,
        influence_rows,
    ) = leave_one_gene_out_diagnostics(
        working_design,
        numerator,
        denominator,
        log_or,
        estimable,
        direction,
        gene_numerator,
        gene_denominator,
        gene_influence,
        analysis_id,
        modification,
        database,
        rbp,
        region_filter,
    )

    biological_modification = "m6a" if modification == "m6a_rep2" else modification
    related = canonical_rbp(rbp) in KNOWN_RELATED.get(
        biological_modification, set()
    )
    row = {
        "analysis_id": analysis_id,
        "modification": modification,
        "database": database,
        "RBP": rbp,
        "region_filter": region_filter or "all",
        "n_source_intervals": int(n_source_intervals),
        "n_case_sites": working_design.n_cases,
        "n_control_sites": working_design.n_controls,
        "n_strata_total": working_design.n_strata,
        "n_transcripts": working_design.n_transcripts,
        "n_genes": working_design.n_genes,
        "n_informative_strata": int(informative.sum()),
        "n_informative_genes": len(informative_genes),
        "n_gene_clusters_for_robust_variance": n_clusters,
        "gene_cluster_df": cluster_df,
        "low_gene_cluster_count": bool(n_clusters < 20),
        "case_overlaps": case_overlaps,
        "case_overlap_rate": case_overlaps / working_design.n_cases,
        "control_overlaps": control_overlaps,
        "control_overlap_rate": control_overlaps / working_design.n_controls,
        "crude_or": crude_or,
        "crude_pvalue": crude_pvalue,
        "mh_numerator": numerator,
        "mh_denominator": denominator,
        "mh_or": mh_or,
        "log2_mh_or": (
            math.log2(mh_or)
            if np.isfinite(mh_or) and mh_or > 0
            else math.nan
        ),
        "cmh_statistic": cmh_statistic,
        "cmh_pvalue": cmh_pvalue,
        "model_log_or_se": model_se,
        "model_ci_lower_or": model_lower,
        "model_ci_upper_or": model_upper,
        "gene_cluster_log_or_se": cluster_se,
        "gene_cluster_test_statistic": cluster_statistic,
        "gene_cluster_pvalue": cluster_pvalue,
        "gene_cluster_ci_lower_or": cluster_lower,
        "gene_cluster_ci_upper_or": cluster_upper,
        "maximum_absolute_gene_influence": maximum_influence,
        "leave_one_gene_out_direction_reversals": reversals,
        "leave_one_gene_out_boundaries": boundaries,
        "direction_stable_after_each_gene_removed": stable,
        "estimable": bool(estimable),
        "boundary_estimate": bool(boundary),
        "direction": direction,
        "is_modification_related": bool(related),
        "partial_input": bool(partial_input),
    }
    return row, influence_rows


def analyse_all_regions(
    design: OpportunityDesign,
    exposure: np.ndarray,
    analysis_id: str,
    modification: str,
    database: str,
    rbp: str,
    n_source_intervals: int,
    partial_input: bool,
) -> tuple[dict, list[dict], list[dict]]:
    row, influence = analyse_exposure(
        design,
        exposure,
        analysis_id,
        modification,
        database,
        rbp,
        n_source_intervals,
        partial_input,
    )
    regional: list[dict] = []
    for region in REGIONS:
        if np.any(design.stratum_regions == region):
            regional.append(
                analyse_exposure(
                    design,
                    exposure,
                    analysis_id,
                    modification,
                    database,
                    rbp,
                    n_source_intervals,
                    partial_input,
                    region_filter=region,
                )[0]
            )
    return row, influence, regional
