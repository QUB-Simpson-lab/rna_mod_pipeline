from __future__ import annotations

import math

import numpy as np
import pandas as pd

from ..crossdb.io import load_database_results
from ..stats import direction, safe_log2
from .profiles import DatasetProfile


DETAIL_COLUMNS = [
    "design",
    "database",
    "RBP",
    "odds_ratio_a",
    "fdr_a",
    "p_value_a",
    "source_file_a",
    "odds_ratio_b",
    "fdr_b",
    "p_value_b",
    "source_file_b",
    "inference_eligible_a",
    "inference_eligible_b",
    "evidence_quality_a",
    "evidence_quality_b",
    "evidence_quality_reasons_a",
    "evidence_quality_reasons_b",
    "gene_cluster_quality_a",
    "gene_cluster_quality_b",
    "leave_one_gene_out_quality_a",
    "leave_one_gene_out_quality_b",
    "single_gene_influence_status_a",
    "single_gene_influence_status_b",
    "membership",
    "present_a",
    "present_b",
    "log2_odds_ratio_a",
    "log2_odds_ratio_b",
    "state_a",
    "state_b",
    "primary_state_a",
    "primary_state_b",
    "significance_call_agreement",
    "primary_significance_call_agreement",
    "both_significant",
    "both_significant_same_direction",
    "both_primary_significant",
    "both_primary_significant_same_direction",
    "eligible_for_primary_effect_comparison",
    "absolute_log2_effect_difference",
]
SUMMARY_COLUMNS = [
    "design",
    "database",
    "n_rbps_a",
    "n_rbps_b",
    "n_shared_rbps",
    "n_finite_shared_effects",
    "n_shared_inference_eligible_effects",
    "log2_effect_pearson",
    "log2_effect_spearman",
    "log2_effect_lins_ccc",
    "median_absolute_log2_effect_difference",
    "finite_effect_direction_agreement_rate",
    "eligible_log2_effect_pearson",
    "eligible_log2_effect_spearman",
    "eligible_log2_effect_lins_ccc",
    "eligible_median_absolute_log2_effect_difference",
    "eligible_effect_direction_agreement_rate",
    "n_significant_a",
    "n_significant_b",
    "n_significant_both",
    "significant_set_jaccard",
    "n_both_significant_same_direction",
    "both_significant_direction_agreement_rate",
    "shared_significance_call_agreement_rate",
    "n_primary_significant_a",
    "n_primary_significant_b",
    "n_primary_significant_both",
    "primary_significant_set_jaccard",
    "n_both_primary_significant_same_direction",
    "both_primary_significant_direction_agreement_rate",
    "shared_primary_significance_call_agreement_rate",
    "comparison_interpretation",
    *[
        f"input_{label}_status_{side}"
        for label in ("modification", "database", "design", "context")
        for side in ("a", "b")
    ],
]


def _finite_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    values = frame[["log2_odds_ratio_a", "log2_odds_ratio_b"]].to_numpy(dtype=float)
    return frame[np.isfinite(values).all(axis=1)].copy()


def _correlation(frame: pd.DataFrame, method: str) -> float:
    if len(frame) < 2:
        return math.nan
    x = frame["log2_odds_ratio_a"].to_numpy(dtype=float)
    y = frame["log2_odds_ratio_b"].to_numpy(dtype=float)
    if method == "spearman":
        x = pd.Series(x).rank(method="average").to_numpy()
        y = pd.Series(y).rank(method="average").to_numpy()
    if np.std(x) == 0 or np.std(y) == 0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def _concordance(frame: pd.DataFrame) -> float:
    if len(frame) < 2:
        return math.nan
    x = frame["log2_odds_ratio_a"].to_numpy(dtype=float)
    y = frame["log2_odds_ratio_b"].to_numpy(dtype=float)
    denominator = np.var(x) + np.var(y) + (np.mean(x) - np.mean(y)) ** 2
    if denominator == 0:
        return math.nan
    covariance = np.mean((x - np.mean(x)) * (y - np.mean(y)))
    return float(2 * covariance / denominator)


def _load(
    profile: DatasetProfile,
    design: str,
    database: str,
    project_root,
) -> pd.DataFrame:
    return load_database_results(
        profile.enrichment[(design, database)],
        database,
        context="all",
        project_root=project_root,
        expected_modification=profile.analysis_id,
        expected_design=design,
    )


def _compare_pair(
    first: pd.DataFrame,
    second: pd.DataFrame,
    design: str,
    database: str,
    fdr_threshold: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    keep = [
        "RBP",
        "odds_ratio",
        "fdr",
        "p_value",
        "source_file",
        "inference_eligible",
        "evidence_quality",
        "evidence_quality_reasons",
        "gene_cluster_quality",
        "leave_one_gene_out_quality",
        "single_gene_influence_status",
    ]
    merged = first[keep].merge(
        second[keep],
        on="RBP",
        how="outer",
        suffixes=("_a", "_b"),
        indicator=True,
        validate="one_to_one",
    )
    merged["membership"] = merged.pop("_merge").map(
        {
            "left_only": "dataset_a_only",
            "right_only": "dataset_b_only",
            "both": "shared",
        }
    )
    merged["present_a"] = merged["membership"].isin(["dataset_a_only", "shared"])
    merged["present_b"] = merged["membership"].isin(["dataset_b_only", "shared"])
    merged["log2_odds_ratio_a"] = merged["odds_ratio_a"].map(safe_log2)
    merged["log2_odds_ratio_b"] = merged["odds_ratio_b"].map(safe_log2)
    merged["state_a"] = [
        direction(float(or_value), float(fdr), fdr_threshold)
        if pd.notna(or_value)
        else "Unavailable"
        for or_value, fdr in zip(merged["odds_ratio_a"], merged["fdr_a"])
    ]
    merged["state_b"] = [
        direction(float(or_value), float(fdr), fdr_threshold)
        if pd.notna(or_value)
        else "Unavailable"
        for or_value, fdr in zip(merged["odds_ratio_b"], merged["fdr_b"])
    ]
    shared = merged["membership"].eq("shared")
    finite = _finite_pairs(merged[shared])
    eligible_finite = finite[
        finite["inference_eligible_a"].astype(bool)
        & finite["inference_eligible_b"].astype(bool)
    ].copy()
    effect_sign_a = np.sign(finite["log2_odds_ratio_a"])
    effect_sign_b = np.sign(finite["log2_odds_ratio_b"])
    finite["effect_direction_agreement"] = effect_sign_a.eq(effect_sign_b)
    eligible_effect_sign_a = np.sign(eligible_finite["log2_odds_ratio_a"])
    eligible_effect_sign_b = np.sign(eligible_finite["log2_odds_ratio_b"])
    eligible_finite["effect_direction_agreement"] = eligible_effect_sign_a.eq(
        eligible_effect_sign_b
    )
    merged["primary_state_a"] = np.where(
        merged["inference_eligible_a"].fillna(False),
        merged["state_a"],
        np.where(
            merged["state_a"].isin(["Enriched", "Depleted"]),
            "Quality excluded",
            merged["state_a"],
        ),
    )
    merged["primary_state_b"] = np.where(
        merged["inference_eligible_b"].fillna(False),
        merged["state_b"],
        np.where(
            merged["state_b"].isin(["Enriched", "Depleted"]),
            "Quality excluded",
            merged["state_b"],
        ),
    )
    merged["significance_call_agreement"] = (
        shared & merged["state_a"].eq(merged["state_b"])
    )
    merged["both_significant"] = (
        shared
        & merged["state_a"].isin(["Enriched", "Depleted"])
        & merged["state_b"].isin(["Enriched", "Depleted"])
    )
    merged["both_significant_same_direction"] = (
        merged["both_significant"] & merged["state_a"].eq(merged["state_b"])
    )
    merged["primary_significance_call_agreement"] = (
        shared & merged["primary_state_a"].eq(merged["primary_state_b"])
    )
    merged["both_primary_significant"] = (
        shared
        & merged["primary_state_a"].isin(["Enriched", "Depleted"])
        & merged["primary_state_b"].isin(["Enriched", "Depleted"])
    )
    merged["both_primary_significant_same_direction"] = (
        merged["both_primary_significant"]
        & merged["primary_state_a"].eq(merged["primary_state_b"])
    )
    finite_a = np.isfinite(merged["log2_odds_ratio_a"])
    finite_b = np.isfinite(merged["log2_odds_ratio_b"])
    merged["eligible_for_primary_effect_comparison"] = (
        shared
        & merged["inference_eligible_a"].fillna(False).astype(bool)
        & merged["inference_eligible_b"].fillna(False).astype(bool)
        & finite_a
        & finite_b
    )
    merged["absolute_log2_effect_difference"] = (
        merged["log2_odds_ratio_b"] - merged["log2_odds_ratio_a"]
    ).abs()

    sig_a = set(merged.loc[merged["state_a"].isin(["Enriched", "Depleted"]), "RBP"])
    sig_b = set(merged.loc[merged["state_b"].isin(["Enriched", "Depleted"]), "RBP"])
    sig_union = sig_a | sig_b
    sig_shared = sig_a & sig_b
    primary_sig_a = set(
        merged.loc[
            merged["primary_state_a"].isin(["Enriched", "Depleted"]), "RBP"
        ]
    )
    primary_sig_b = set(
        merged.loc[
            merged["primary_state_b"].isin(["Enriched", "Depleted"]), "RBP"
        ]
    )
    primary_sig_union = primary_sig_a | primary_sig_b
    primary_sig_shared = primary_sig_a & primary_sig_b
    both_significant = int(merged["both_significant"].sum())
    both_primary_significant = int(merged["both_primary_significant"].sum())
    summary = {
        "design": design,
        "database": database,
        "n_rbps_a": int(merged["present_a"].sum()),
        "n_rbps_b": int(merged["present_b"].sum()),
        "n_shared_rbps": int(shared.sum()),
        "n_finite_shared_effects": len(finite),
        "n_shared_inference_eligible_effects": len(eligible_finite),
        "log2_effect_pearson": _correlation(finite, "pearson"),
        "log2_effect_spearman": _correlation(finite, "spearman"),
        "log2_effect_lins_ccc": _concordance(finite),
        "median_absolute_log2_effect_difference": (
            float(finite["log2_odds_ratio_b"].sub(finite["log2_odds_ratio_a"]).abs().median())
            if len(finite)
            else math.nan
        ),
        "finite_effect_direction_agreement_rate": (
            float(finite["effect_direction_agreement"].mean())
            if len(finite)
            else math.nan
        ),
        "eligible_log2_effect_pearson": _correlation(
            eligible_finite, "pearson"
        ),
        "eligible_log2_effect_spearman": _correlation(
            eligible_finite, "spearman"
        ),
        "eligible_log2_effect_lins_ccc": _concordance(eligible_finite),
        "eligible_median_absolute_log2_effect_difference": (
            float(
                eligible_finite["log2_odds_ratio_b"]
                .sub(eligible_finite["log2_odds_ratio_a"])
                .abs()
                .median()
            )
            if len(eligible_finite)
            else math.nan
        ),
        "eligible_effect_direction_agreement_rate": (
            float(eligible_finite["effect_direction_agreement"].mean())
            if len(eligible_finite)
            else math.nan
        ),
        "n_significant_a": len(sig_a),
        "n_significant_b": len(sig_b),
        "n_significant_both": len(sig_shared),
        "significant_set_jaccard": (
            len(sig_shared) / len(sig_union) if sig_union else math.nan
        ),
        "n_both_significant_same_direction": int(
            merged["both_significant_same_direction"].sum()
        ),
        "both_significant_direction_agreement_rate": (
            int(merged["both_significant_same_direction"].sum()) / both_significant
            if both_significant
            else math.nan
        ),
        "shared_significance_call_agreement_rate": (
            float(merged.loc[shared, "significance_call_agreement"].mean())
            if shared.any()
            else math.nan
        ),
        "n_primary_significant_a": len(primary_sig_a),
        "n_primary_significant_b": len(primary_sig_b),
        "n_primary_significant_both": len(primary_sig_shared),
        "primary_significant_set_jaccard": (
            len(primary_sig_shared) / len(primary_sig_union)
            if primary_sig_union
            else math.nan
        ),
        "n_both_primary_significant_same_direction": int(
            merged["both_primary_significant_same_direction"].sum()
        ),
        "both_primary_significant_direction_agreement_rate": (
            int(merged["both_primary_significant_same_direction"].sum())
            / both_primary_significant
            if both_primary_significant
            else math.nan
        ),
        "shared_primary_significance_call_agreement_rate": (
            float(
                merged.loc[shared, "primary_significance_call_agreement"].mean()
            )
            if shared.any()
            else math.nan
        ),
        "comparison_interpretation": (
            "Dataset/callset robustness using shared RBP resources; not "
            "independent CLIP replication."
        ),
    }
    for label in ("modification", "database", "design", "context"):
        column = f"input_{label}_status"
        summary[f"{column}_a"] = str(first[column].iloc[0])
        summary[f"{column}_b"] = str(second[column].iloc[0])
    merged.insert(0, "database", database)
    merged.insert(0, "design", design)
    return merged, summary


def compare_enrichment(
    first: DatasetProfile,
    second: DatasetProfile,
    project_root,
    fdr_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    details = []
    summaries = []
    common_pairs = sorted(set(first.enrichment) & set(second.enrichment))
    for design, database in common_pairs:
        table_a = _load(first, design, database, project_root)
        table_b = _load(second, design, database, project_root)
        detail, summary = _compare_pair(
            table_a, table_b, design, database, fdr_threshold
        )
        details.append(detail)
        summaries.append(summary)
    return (
        (
            pd.concat(details, ignore_index=True)
            if details
            else pd.DataFrame(columns=DETAIL_COLUMNS)
        ),
        pd.DataFrame(summaries, columns=SUMMARY_COLUMNS),
    )
