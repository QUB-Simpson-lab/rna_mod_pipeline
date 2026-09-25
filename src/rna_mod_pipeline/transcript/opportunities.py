from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .case_assignment import assign_cases, validate_case_table
from .control_stream import (
    seed_case_observations,
    stream_low_fraction_controls,
)
from .transcripts import (
    PrimaryTranscriptIndex,
    TranscriptModel,
    load_complete_basic_transcripts,
)


@dataclass
class OpportunityDesign:
    chroms: np.ndarray
    positions: np.ndarray
    strands: np.ndarray
    coverages: np.ndarray
    fractions: np.ndarray
    is_case: np.ndarray
    stratum_codes: np.ndarray
    case_strata: np.ndarray
    control_strata: np.ndarray
    case_counts: np.ndarray
    control_counts: np.ndarray
    stratum_transcripts: np.ndarray
    stratum_gene_ids: np.ndarray
    stratum_genes: np.ndarray
    stratum_regions: np.ndarray
    stratum_gene_codes: np.ndarray
    genes: np.ndarray
    gene_names: np.ndarray
    case_site_ids: np.ndarray

    @property
    def n_observations(self) -> int:
        return len(self.positions)

    @property
    def n_cases(self) -> int:
        return int(self.is_case.sum())

    @property
    def n_controls(self) -> int:
        return int((~self.is_case).sum())

    @property
    def n_strata(self) -> int:
        return len(self.case_counts)

    @property
    def n_transcripts(self) -> int:
        return len(set(map(str, self.stratum_transcripts)))

    @property
    def n_genes(self) -> int:
        return len(self.genes)


@dataclass(frozen=True)
class OpportunityBuild:
    design: OpportunityDesign
    case_assignments: pd.DataFrame
    strata: pd.DataFrame
    flow: pd.DataFrame
    partial_bed_scan: bool


def _compact_design(
    chroms: list[str],
    positions: list[int],
    strands: list[str],
    coverages: list[int],
    fractions: list[float],
    outcomes: list[bool],
    stratum_codes: list[int],
    stratum_keys: list[tuple[str, str]],
    models: dict[str, TranscriptModel],
    case_site_ids: list[int],
) -> tuple[OpportunityDesign, pd.DataFrame, Counter]:
    chrom_array = np.asarray(chroms, dtype=object)
    position_array = np.asarray(positions, dtype=np.int64)
    strand_array = np.asarray(strands, dtype="U1")
    coverage_array = np.asarray(coverages, dtype=np.int32)
    fraction_array = np.asarray(fractions, dtype=np.float32)
    case_array = np.asarray(outcomes, dtype=bool)
    stratum_array = np.asarray(stratum_codes, dtype=np.int32)

    n_initial = len(stratum_keys)
    case_counts = np.bincount(
        stratum_array[case_array], minlength=n_initial
    ).astype(np.int64)
    control_counts = np.bincount(
        stratum_array[~case_array], minlength=n_initial
    ).astype(np.int64)
    informative = (case_counts > 0) & (control_counts > 0)
    keep = informative[stratum_array]
    remap = np.full(n_initial, -1, dtype=np.int32)
    remap[informative] = np.arange(informative.sum(), dtype=np.int32)
    compact_codes = remap[stratum_array[keep]]
    compact_case = case_array[keep]
    old_codes = np.flatnonzero(informative)
    keys = [stratum_keys[index] for index in old_codes]
    transcripts = np.asarray([key[0] for key in keys], dtype=object)
    regions = np.asarray([key[1] for key in keys], dtype=object)
    gene_ids = np.asarray(
        [models[transcript_id].gene_id for transcript_id in transcripts],
        dtype=object,
    )
    gene_names = np.asarray(
        [models[transcript_id].gene_name for transcript_id in transcripts],
        dtype=object,
    )
    if any(not str(value) for value in gene_ids):
        raise ValueError("A retained transcript has no stable GENCODE gene_id")
    genes = np.asarray(sorted(set(map(str, gene_ids))), dtype=object)
    gene_to_code = {gene: index for index, gene in enumerate(genes)}
    gene_codes = np.asarray(
        [gene_to_code[str(gene)] for gene in gene_ids], dtype=np.int32
    )
    names_by_id: dict[str, str] = {}
    for gene_id, gene_name in zip(gene_ids, gene_names):
        previous = names_by_id.setdefault(str(gene_id), str(gene_name))
        if previous != str(gene_name):
            raise ValueError(f"{gene_id} maps to inconsistent gene names")
    unique_gene_names = np.asarray(
        [names_by_id[str(gene)] for gene in genes], dtype=object
    )

    original_case_rows = np.flatnonzero(case_array)
    retained_case_rows = original_case_rows[keep[original_case_rows]]
    id_by_row = {
        int(row): int(site_id)
        for row, site_id in zip(original_case_rows, case_site_ids)
    }
    retained_case_ids = np.asarray(
        [id_by_row[int(row)] for row in retained_case_rows], dtype=np.int64
    )
    design = OpportunityDesign(
        chroms=chrom_array[keep],
        positions=position_array[keep],
        strands=strand_array[keep],
        coverages=coverage_array[keep],
        fractions=fraction_array[keep],
        is_case=compact_case,
        stratum_codes=compact_codes,
        case_strata=compact_codes[compact_case],
        control_strata=compact_codes[~compact_case],
        case_counts=case_counts[informative],
        control_counts=control_counts[informative],
        stratum_transcripts=transcripts,
        stratum_gene_ids=gene_ids,
        stratum_genes=gene_names,
        stratum_regions=regions,
        stratum_gene_codes=gene_codes,
        genes=genes,
        gene_names=unique_gene_names,
        case_site_ids=retained_case_ids,
    )
    strata = pd.DataFrame(
        {
            "stratum_id": np.arange(design.n_strata, dtype=int),
            "transcript_id": transcripts,
            "gene_id": gene_ids,
            "gene_name": gene_names,
            "region": regions,
            "n_case_sites": design.case_counts,
            "n_control_sites": design.control_counts,
        }
    )
    flow = Counter(
        initial_strata_with_cases=n_initial,
        strata_with_case_and_control=int(informative.sum()),
        strata_without_controls_removed=int(
            ((case_counts > 0) & (control_counts == 0)).sum()
        ),
        case_sites_removed_with_control_free_strata=int(
            case_counts[~informative].sum()
        ),
        control_sites_removed_outside_contributing_strata=int(
            control_counts[~informative].sum()
        ),
        contributing_case_sites=design.n_cases,
        contributing_control_sites=design.n_controls,
        contributing_transcripts=design.n_transcripts,
        contributing_genes=design.n_genes,
    )
    return design, strata, flow


def build_opportunity_design(
    sites_path: Path,
    bedmethyl_path: Path,
    gtf_path: Path,
    modification: str,
    coverage_min: int = 20,
    case_min_fraction: float = 20.0,
    background_max_fraction: float = 20.0,
    chromosomes: Optional[set[str]] = None,
    max_cases: int | None = None,
    max_bed_rows: int | None = None,
    seed: int = 20260728,
) -> OpportunityBuild:
    cases, all_high, case_flow = validate_case_table(
        sites_path,
        modification,
        coverage_min,
        case_min_fraction,
        chromosomes,
    )
    all_expected = {
        (str(row.chrom), int(row.start), str(row.strand)): (
            int(row.Nvalid_cov),
            float(row.fraction_modified),
            str(row.mod_code).lower(),
        )
        for row in cases.itertuples(index=False)
    }
    models = load_complete_basic_transcripts(gtf_path, chromosomes)
    if not models:
        raise RuntimeError("No complete GENCODE transcript models were loaded")
    transcript_index = PrimaryTranscriptIndex(models)
    assigned, assignment_flow = assign_cases(
        cases, transcript_index, models, max_cases, seed
    )
    stratum_keys = sorted(
        {
            (str(row.transcript_id), str(row.region))
            for row in assigned.itertuples(index=False)
        }
    )
    code_by_key = {key: index for index, key in enumerate(stratum_keys)}
    target_strata = set(stratum_keys)

    buffers = seed_case_observations(assigned, code_by_key)
    stream_flow, partial = stream_low_fraction_controls(
        bedmethyl_path,
        buffers,
        all_expected,
        all_high,
        modification,
        coverage_min,
        case_min_fraction,
        background_max_fraction,
        chromosomes,
        max_bed_rows,
        transcript_index,
        target_strata,
        code_by_key,
    )

    design, strata, compact_flow = _compact_design(
        buffers.chroms,
        buffers.positions,
        buffers.strands,
        buffers.coverages,
        buffers.fractions,
        buffers.outcomes,
        buffers.stratum_codes,
        stratum_keys,
        models,
        buffers.case_site_ids,
    )
    if not design.n_cases or not design.n_controls or not design.n_strata:
        raise RuntimeError("No strata contain both cases and controls")
    assignments = pd.DataFrame(buffers.case_records)
    retained_ids = set(map(int, design.case_site_ids))
    assignments["contributes_to_stratified_analysis"] = assignments[
        "case_site_id"
    ].isin(retained_ids)
    flow_rows = []
    for section, values in (
        ("case_loading", case_flow),
        ("primary_assignment", assignment_flow),
        ("bedmethyl_stream", stream_flow),
        ("design_compaction", compact_flow),
    ):
        flow_rows.extend(
            {
                "section": section,
                "metric": metric,
                "value": int(value),
            }
            for metric, value in sorted(values.items())
        )
    return OpportunityBuild(
        design=design,
        case_assignments=assignments,
        strata=strata,
        flow=pd.DataFrame(flow_rows),
        partial_bed_scan=partial,
    )
