from __future__ import annotations

import itertools
import math
from typing import Mapping

import numpy as np
import pandas as pd
from scipy import stats

from ..evidence_quality import QUALITY_COLUMNS, add_evidence_quality
from ..stats import bh_adjust, direction


def _state(
    odds_ratio: float,
    fdr: float,
    threshold: float,
    inference_eligible: bool = True,
    estimate_quality: str = "",
) -> str:
    if not inference_eligible:
        if estimate_quality == "boundary":
            return "boundary_estimate"
        if estimate_quality in {"non_estimable", "finite_estimate_without_fdr"}:
            return "tested_unestimable"
        return "quality_excluded"
    label = direction(odds_ratio, fdr, threshold)
    if label == "Unavailable":
        return "tested_unestimable"
    if label == "Not significant":
        return "tested_not_significant"
    return f"significant_{label.lower()}"


def _classify(
    group: pd.DataFrame,
    requested: tuple[str, ...],
    state_column: str = "state",
) -> str:
    significant = group[group[state_column].str.startswith("significant_")]
    covered = set(group["modification"])
    if significant.empty:
        return "no_significant_result" if covered else "not_covered"
    n = len(significant)
    directions = set(significant["direction"])
    if n == 1:
        return "modification_specific_significant"
    scope = "all_requested" if n == len(requested) else f"shared_{n}"
    return (
        f"{scope}_same_direction"
        if len(directions) == 1
        else f"{scope}_direction_switching"
    )


def _load_long_table(
    tables: Mapping[tuple[str, str], pd.DataFrame],
    modifications: tuple[str, ...],
    databases: tuple[str, ...],
    fdr_threshold: float,
) -> pd.DataFrame:
    rows = []
    for modification in modifications:
        for database in databases:
            table = tables.get((modification, database))
            if table is None:
                continue
            qualified = add_evidence_quality(
                table,
                odds_column="odds_ratio",
                fdr_column="fdr",
            )
            for _, record in qualified.iterrows():
                odds_ratio = float(record["odds_ratio"])
                fdr = float(record["fdr"])
                label = direction(odds_ratio, fdr, fdr_threshold)
                inference_eligible = bool(record["inference_eligible"])
                row = {
                    "database": database,
                    "modification": modification,
                    "RBP": record["RBP"],
                    "covered": True,
                    "odds_ratio": odds_ratio,
                    "log2_odds_ratio": (
                        math.log2(odds_ratio)
                        if math.isfinite(odds_ratio) and odds_ratio > 0
                        else np.nan
                    ),
                    "p_value": record["p_value"],
                    "fdr": fdr,
                    "direction": label,
                    "statistically_significant": label in {"Enriched", "Depleted"},
                    "primary_direction": (
                        label
                        if inference_eligible
                        else "Quality excluded"
                        if label in {"Enriched", "Depleted"}
                        else label
                    ),
                    "statistical_state": _state(odds_ratio, fdr, fdr_threshold),
                    "state": _state(
                        odds_ratio,
                        fdr,
                        fdr_threshold,
                        inference_eligible,
                        str(record["estimate_quality"]),
                    ),
                    "source_file": record["source_file"],
                }
                for column in QUALITY_COLUMNS:
                    row[column] = record[column]
                for column in (
                    "n_gene_clusters_for_robust_variance",
                    "low_gene_cluster_count",
                    "direction_stable_after_each_gene_removed",
                    "leave_one_gene_out_direction_reversals",
                    "leave_one_gene_out_boundaries",
                    "maximum_absolute_gene_influence",
                    "partial_input",
                    "estimable",
                    "boundary_estimate",
                ):
                    if column in record.index:
                        row[column] = record[column]
                rows.append(row)
    long = pd.DataFrame(rows)
    if long.empty:
        raise ValueError("No cross-modification input rows were loaded")
    if long.duplicated(["database", "modification", "RBP"]).any():
        raise ValueError("Duplicate database/modification/RBP result")
    return long


def _complete_matrix(
    long: pd.DataFrame,
    modifications: tuple[str, ...],
    databases: tuple[str, ...],
) -> pd.DataFrame:
    complete_rows = []
    for database in databases:
        names = sorted(
            set(long.loc[long["database"].eq(database), "RBP"].astype(str))
        )
        observed = long[long["database"].eq(database)].set_index(
            ["modification", "RBP"]
        )
        for rbp in names:
            for modification in modifications:
                key = (modification, rbp)
                if key in observed.index:
                    record = observed.loc[key].to_dict()
                    record.update(
                        {
                            "database": database,
                            "modification": modification,
                            "RBP": rbp,
                        }
                    )
                    complete_rows.append(record)
                else:
                    complete_rows.append(
                        {
                            "database": database,
                            "modification": modification,
                            "RBP": rbp,
                            "covered": False,
                            "odds_ratio": np.nan,
                            "log2_odds_ratio": np.nan,
                            "p_value": np.nan,
                            "fdr": np.nan,
                            "direction": "Unavailable",
                            "statistically_significant": False,
                            "primary_direction": "Unavailable",
                            "statistical_state": "not_covered",
                            "state": "not_covered",
                            "source_file": "",
                            "estimate_eligible": False,
                            "inference_eligible": False,
                            "evidence_quality": "not_covered",
                            "evidence_quality_reasons": "not_covered",
                            "inference_exclusion_reasons": "not_covered",
                            "input_completeness": "not_available",
                            "estimate_quality": "not_covered",
                            "gene_cluster_quality": "not_available",
                            "leave_one_gene_out_quality": "not_available",
                            "single_gene_influence_status": "not_available",
                            "quality_diagnostics_available": "none",
                        }
                    )
    return pd.DataFrame(complete_rows)


def _pattern_table(
    matrix: pd.DataFrame,
    modifications: tuple[str, ...],
) -> pd.DataFrame:
    patterns = []
    for (database, rbp), group in matrix.groupby(["database", "RBP"], sort=True):
        significant = group[group["state"].str.startswith("significant_")]
        statistical = group[
            group["statistical_state"].str.startswith("significant_")
        ]
        patterns.append(
            {
                "database": database,
                "RBP": rbp,
                "n_modifications_covered": int(group["covered"].sum()),
                "n_modifications_significant": len(significant),
                "n_modifications_statistically_significant": len(statistical),
                "covered_modifications": ";".join(
                    group.loc[group["covered"], "modification"]
                ),
                "significant_modifications": ";".join(significant["modification"]),
                "significant_directions": ";".join(significant["direction"]),
                "pattern": _classify(group[group["covered"]], modifications),
                "statistical_significant_modifications": ";".join(
                    statistical["modification"]
                ),
                "statistical_significant_directions": ";".join(
                    statistical["direction"]
                ),
                "statistical_pattern": _classify(
                    group[group["covered"]],
                    modifications,
                    state_column="statistical_state",
                ),
            }
        )
    return pd.DataFrame(patterns)


def _correlation_result(
    valid: pd.DataFrame,
    method: str,
) -> tuple[float, float]:
    if len(valid) < 3:
        return np.nan, np.nan
    test = (
        stats.pearsonr(
            valid["log2_odds_ratio_first"],
            valid["log2_odds_ratio_second"],
        )
        if method == "pearson"
        else stats.spearmanr(
            valid["log2_odds_ratio_first"],
            valid["log2_odds_ratio_second"],
        )
    )
    return float(test.statistic), float(test.pvalue)


def _pairwise_tables(
    matrix: pd.DataFrame,
    modifications: tuple[str, ...],
    databases: tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    correlations = []
    overlaps = []
    for database in databases:
        subset = matrix[matrix["database"].eq(database)]
        for first, second in itertools.combinations(modifications, 2):
            left = subset[subset["modification"].eq(first)].set_index("RBP")
            right = subset[subset["modification"].eq(second)].set_index("RBP")
            shared = left[["log2_odds_ratio", "state"]].join(
                right[["log2_odds_ratio", "state"]],
                how="inner",
                lsuffix="_first",
                rsuffix="_second",
            )
            valid = shared.dropna(subset=["log2_odds_ratio_first", "log2_odds_ratio_second"])
            eligible_names = set(
                left.index[left["inference_eligible"].astype(bool)]
            ) & set(right.index[right["inference_eligible"].astype(bool)])
            eligible_valid = valid.loc[valid.index.isin(eligible_names)]
            for method in ("pearson", "spearman"):
                coefficient, p_value = _correlation_result(valid, method)
                eligible_coefficient, eligible_p_value = _correlation_result(
                    eligible_valid, method
                )
                correlations.append(
                    {
                        "database": database,
                        "modification_a": first,
                        "modification_b": second,
                        "method": method,
                        "n_shared_estimable_rbps": len(valid),
                        "coefficient": coefficient,
                        "p_value": p_value,
                        "n_shared_inference_eligible_rbps": len(eligible_valid),
                        "inference_eligible_coefficient": eligible_coefficient,
                        "inference_eligible_p_value": eligible_p_value,
                    }
                )
            sig_a = left.index[left["state"].str.startswith("significant_")]
            sig_b = right.index[right["state"].str.startswith("significant_")]
            common = set(sig_a) & set(sig_b)
            switches = {
                rbp
                for rbp in common
                if left.loc[rbp, "state"] != right.loc[rbp, "state"]
            }
            union = set(sig_a) | set(sig_b)
            raw_sig_a = left.index[
                left["statistical_state"].str.startswith("significant_")
            ]
            raw_sig_b = right.index[
                right["statistical_state"].str.startswith("significant_")
            ]
            raw_common = set(raw_sig_a) & set(raw_sig_b)
            raw_switches = {
                rbp
                for rbp in raw_common
                if left.loc[rbp, "statistical_state"]
                != right.loc[rbp, "statistical_state"]
            }
            raw_union = set(raw_sig_a) | set(raw_sig_b)
            overlaps.append(
                {
                    "database": database,
                    "modification_a": first,
                    "modification_b": second,
                    "n_significant_a": len(sig_a),
                    "n_significant_b": len(sig_b),
                    "n_shared_significant": len(common),
                    "n_direction_switching": len(switches),
                    "jaccard_significant": len(common) / len(union) if union else np.nan,
                    "shared_significant_rbps": ";".join(sorted(common)),
                    "direction_switching_rbps": ";".join(sorted(switches)),
                    "n_statistically_significant_a": len(raw_sig_a),
                    "n_statistically_significant_b": len(raw_sig_b),
                    "n_shared_statistically_significant": len(raw_common),
                    "n_statistical_direction_switching": len(raw_switches),
                    "jaccard_statistically_significant": (
                        len(raw_common) / len(raw_union) if raw_union else np.nan
                    ),
                    "shared_statistically_significant_rbps": ";".join(
                        sorted(raw_common)
                    ),
                    "statistical_direction_switching_rbps": ";".join(
                        sorted(raw_switches)
                    ),
                }
            )
    correlation_table = pd.DataFrame(correlations)
    correlation_table["fdr_within_database"] = correlation_table.groupby(
        "database", sort=False
    )["p_value"].transform(lambda values: bh_adjust(values.to_numpy(dtype=float)))
    correlation_table["fdr_across_run"] = bh_adjust(
        correlation_table["p_value"].to_numpy(dtype=float)
    )
    correlation_table["inference_eligible_fdr_within_database"] = (
        correlation_table.groupby("database", sort=False)[
            "inference_eligible_p_value"
        ].transform(lambda values: bh_adjust(values.to_numpy(dtype=float)))
    )
    correlation_table["inference_eligible_fdr_across_run"] = bh_adjust(
        correlation_table["inference_eligible_p_value"].to_numpy(dtype=float)
    )
    return correlation_table, pd.DataFrame(overlaps)


def compare_modifications(
    tables: Mapping[tuple[str, str], pd.DataFrame],
    *,
    modifications: tuple[str, ...],
    databases: tuple[str, ...],
    fdr_threshold: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if len(modifications) < 2:
        raise ValueError("Cross-modification comparison requires at least two modifications")
    long = _load_long_table(tables, modifications, databases, fdr_threshold)
    matrix = _complete_matrix(long, modifications, databases)
    patterns = _pattern_table(matrix, modifications)
    correlations, overlaps = _pairwise_tables(matrix, modifications, databases)
    return matrix, patterns, correlations, overlaps
