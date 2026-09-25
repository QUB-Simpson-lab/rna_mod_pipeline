from __future__ import annotations

import math
from typing import Mapping

import numpy as np
import pandas as pd
from scipy import stats

from ..crossdb.analysis import recompute_consensus_fields
from ..stats import bh_adjust


FISHER_FIELDS = [
    "status",
    "universe_size",
    "n_targets_in_universe",
    "n_degs_in_universe",
    "a_target_and_deg",
    "b_deg_not_target",
    "c_target_not_deg",
    "d_neither",
    "odds_ratio",
    "p_value",
]
ORTHOGONAL_BASE_COLUMNS = [
    "rbp",
    "dataset_id",
    "cell_line",
    "cell_line_relevance",
    "design",
    "context_region",
    "assay",
    "treated_n",
    "control_n",
    "has_pvalues",
    "n_resolved_degs",
    "n_targets_encori",
    "n_targets_postar3",
    "encori_resource_available",
    "postar3_resource_available",
    "encori_resource_status",
    "postar3_resource_status",
    "encori_availability_basis",
    "postar3_availability_basis",
    "n_targets_union",
    "n_targets_intersection",
    "n_overlap",
    "overlap_genes",
    "overlap_up_genes",
    "overlap_down_genes",
    "small_count_warning",
    "interpretation",
]
REGULATORY_COLUMNS = [
    "knocked_down_rbp",
    "affected_rbp",
    "dataset_id",
    "cell_line",
    "log2fc",
    "padj",
    "deg_direction",
    "relationship_label",
    "consensus_direction",
    "statistical_consensus_direction",
    "n_databases_present",
    "n_databases_significant",
    "n_databases_statistically_significant",
    "n_databases_inference_eligible",
    "triple_significant",
    "triple_statistically_significant",
    "triple_direction_concordant",
    "direction_conflict",
    "statistical_direction_concordant_consensus",
    "statistical_direction_conflict",
    "consensus_evidence_quality",
    "consensus_quality_reasons",
    "arithmetic_mean_significant_or",
    "geometric_mean_significant_or",
    "arithmetic_mean_statistically_significant_or",
    "geometric_mean_statistically_significant_or",
    *[
        f"{database}_{field}"
        for database in ("ornament", "encori", "postar3")
        for field in (
            "present",
            "inference_eligible",
            "evidence_quality",
            "inference_exclusion_reasons",
            "gene_cluster_quality",
            "leave_one_gene_out_quality",
            "single_gene_influence_status",
        )
    ],
]


def _orthogonal_columns() -> list[str]:
    columns = list(ORTHOGONAL_BASE_COLUMNS)
    for prefix in ("primary_modified_gene_universe", "sensitivity_fixed_20000"):
        columns.extend(f"{prefix}_{field}" for field in FISHER_FIELDS)
        columns.append(f"{prefix}_fdr")
    return columns


def _join(values: set[str]) -> str:
    return ";".join(sorted(values))


def _fisher(
    targets: set[str],
    degs: set[str],
    universe: set[str] | None,
    fixed_size: int = 20_000,
) -> dict[str, object]:
    if universe is None:
        active_targets, active_degs, size = set(targets), set(degs), fixed_size
    else:
        active_targets = set(targets) & universe
        active_degs = set(degs) & universe
        size = len(universe)
    if size <= 0 or len(active_targets | active_degs) > size:
        return {"status": "not_testable_invalid_universe", "odds_ratio": np.nan, "p_value": np.nan}
    a = len(active_targets & active_degs)
    b = len(active_degs - active_targets)
    c = len(active_targets - active_degs)
    d = size - len(active_targets | active_degs)
    if not active_targets:
        status = "not_testable_zero_targets"
    elif not active_degs:
        status = "not_testable_zero_degs"
    else:
        status = "tested"
    if status == "tested":
        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    else:
        odds_ratio, p_value = np.nan, np.nan
    return {
        "status": status,
        "universe_size": size,
        "n_targets_in_universe": len(active_targets),
        "n_degs_in_universe": len(active_degs),
        "a_target_and_deg": a,
        "b_deg_not_target": b,
        "c_target_not_deg": c,
        "d_neither": d,
        "odds_ratio": float(odds_ratio),
        "p_value": float(p_value),
    }


def orthogonal_validation(
    deg_data: Mapping[str, Mapping[str, object]],
    targets: Mapping[str, Mapping[str, object]],
    modified_gene_universe: set[str],
    *,
    design: str,
    context: str,
) -> pd.DataFrame:
    rows = []
    for dataset_id, dataset in sorted(deg_data.items()):
        rbp = str(dataset["rbp"])
        if rbp not in targets:
            continue
        genes = dataset["genes"]
        degs = set(genes)
        target = set(targets[rbp]["union"])
        overlap = target & degs
        row = {
            "rbp": rbp,
            "dataset_id": dataset_id,
            "cell_line": dataset["cell_line"],
            "cell_line_relevance": (
                "same-cell-line"
                if dataset["cell_line"] == "MDA-MB-231"
                else "closely-related derivative"
                if dataset["cell_line"] == "MDA-MB-231-LM2"
                else "other context"
            ),
            "design": design,
            "context_region": context,
            "assay": dataset["assay"],
            "treated_n": dataset["treated_n"],
            "control_n": dataset["control_n"],
            "has_pvalues": dataset["has_pvalues"],
            "n_resolved_degs": len(degs),
            "n_targets_encori": len(targets[rbp]["encori"]),
            "n_targets_postar3": len(targets[rbp]["postar3"]),
            "encori_resource_available": bool(
                targets[rbp].get("encori_available", False)
            ),
            "postar3_resource_available": bool(
                targets[rbp].get("postar3_available", False)
            ),
            "encori_resource_status": targets[rbp].get(
                "encori_resource_status", "not_assessed"
            ),
            "postar3_resource_status": targets[rbp].get(
                "postar3_resource_status", "not_assessed"
            ),
            "encori_availability_basis": targets[rbp].get(
                "encori_availability_basis", "not_assessed"
            ),
            "postar3_availability_basis": targets[rbp].get(
                "postar3_availability_basis", "not_assessed"
            ),
            "n_targets_union": len(target),
            "n_targets_intersection": len(targets[rbp]["intersection"]),
            "n_overlap": len(overlap),
            "overlap_genes": _join(overlap),
            "overlap_up_genes": _join(
                {gene for gene in overlap if genes[gene]["direction"] == "up"}
            ),
            "overlap_down_genes": _join(
                {gene for gene in overlap if genes[gene]["direction"] == "down"}
            ),
            "small_count_warning": len(target) < 5 or len(overlap) < 2,
            "interpretation": (
                "association_only; small counts"
                if len(target) < 5 or len(overlap) < 2
                else "association_only; adequate minimum counts"
            ),
        }
        for prefix, result in (
            ("primary_modified_gene_universe", _fisher(target, degs, modified_gene_universe)),
            ("sensitivity_fixed_20000", _fisher(target, degs, None)),
        ):
            row.update({f"{prefix}_{key}": value for key, value in result.items()})
        rows.append(row)
    result = pd.DataFrame(rows)
    if result.empty:
        return pd.DataFrame(columns=_orthogonal_columns())
    for prefix in ("primary_modified_gene_universe", "sensitivity_fixed_20000"):
        result[f"{prefix}_fdr"] = bh_adjust(result[f"{prefix}_p_value"])
    return result.reindex(columns=_orthogonal_columns())


def regulatory_network(
    deg_data: Mapping[str, Mapping[str, object]],
    cross_database: pd.DataFrame,
) -> pd.DataFrame:
    if "RBP" not in cross_database:
        raise ValueError("Cross-database table must contain RBP")
    qualified = recompute_consensus_fields(cross_database)
    annotations = qualified.set_index("RBP", verify_integrity=True).to_dict("index")
    rows = []
    for dataset_id, dataset in sorted(deg_data.items()):
        knocked_down = str(dataset["rbp"])
        for affected, annotation in sorted(annotations.items()):
            change = dataset["genes"].get(affected)
            if change is None:
                continue
            row = {
                "knocked_down_rbp": knocked_down,
                "affected_rbp": affected,
                "dataset_id": dataset_id,
                "cell_line": dataset["cell_line"],
                "log2fc": change["log2fc"],
                "padj": change["padj"],
                "deg_direction": change["direction"],
                "relationship_label": "perturbation-associated; directness not established",
            }
            for column in REGULATORY_COLUMNS[8:]:
                row[column] = annotation.get(column, np.nan)
            rows.append(row)
    return pd.DataFrame(rows, columns=REGULATORY_COLUMNS)


def dataset_summary(
    deg_data: Mapping[str, Mapping[str, object]],
) -> pd.DataFrame:
    rows = []
    for dataset_id, dataset in sorted(deg_data.items()):
        rows.append(
            {
                "dataset_id": dataset_id,
                "rbp": dataset["rbp"],
                "cell_line": dataset["cell_line"],
                "assay": dataset["assay"],
                "perturbation": dataset["perturbation"],
                "treated_n": dataset["treated_n"],
                "control_n": dataset["control_n"],
                "has_pvalues": dataset["has_pvalues"],
                "n_filtered_rows_before_resolution": dataset[
                    "n_filtered_rows_before_resolution"
                ],
                "n_resolved_degs": dataset["n_resolved_degs"],
            }
        )
    return pd.DataFrame(rows)


def target_summary(
    targets: Mapping[str, Mapping[str, object]],
) -> pd.DataFrame:
    rows = []
    for rbp, target in sorted(targets.items()):
        rows.append(
            {
                "rbp": rbp,
                "encori_resource_available": target.get("encori_available", False),
                "postar3_resource_available": target.get("postar3_available", False),
                "encori_resource_status": target.get(
                    "encori_resource_status", "not_assessed"
                ),
                "postar3_resource_status": target.get(
                    "postar3_resource_status", "not_assessed"
                ),
                "encori_availability_basis": target.get(
                    "encori_availability_basis", "not_assessed"
                ),
                "postar3_availability_basis": target.get(
                    "postar3_availability_basis", "not_assessed"
                ),
                "n_targets_encori": len(target["encori"]),
                "n_targets_postar3": len(target["postar3"]),
                "n_targets_union": len(target["union"]),
                "n_targets_intersection": len(target["intersection"]),
                "target_genes_encori": _join(set(target["encori"])),
                "target_genes_postar3": _join(set(target["postar3"])),
                "target_genes_union": _join(set(target["union"])),
                "target_genes_intersection": _join(set(target["intersection"])),
            }
        )
    return pd.DataFrame(rows)


def regulatory_summary(network: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "dataset_id",
        "knocked_down_rbp",
        "cell_line",
        "n_affected_rbps",
        "n_up",
        "n_down",
        "n_triple_direction_concordant",
        "n_sensitivity_only_direction_concordant",
        "n_cross_database_direction_conflicts",
    ]
    if network.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for keys, group in network.groupby(
        ["dataset_id", "knocked_down_rbp", "cell_line"], sort=True
    ):
        rows.append(
            {
                "dataset_id": keys[0],
                "knocked_down_rbp": keys[1],
                "cell_line": keys[2],
                "n_affected_rbps": group["affected_rbp"].nunique(),
                "n_up": int(group["deg_direction"].eq("up").sum()),
                "n_down": int(group["deg_direction"].eq("down").sum()),
                "n_triple_direction_concordant": int(
                    group.get("triple_direction_concordant", False).eq(True).sum()
                    if "triple_direction_concordant" in group
                    else 0
                ),
                "n_sensitivity_only_direction_concordant": int(
                    group.get("consensus_evidence_quality", "").eq(
                        "sensitivity_only_consensus_quality_excluded"
                    ).sum()
                    if "consensus_evidence_quality" in group
                    else 0
                ),
                "n_cross_database_direction_conflicts": int(
                    group.get("direction_conflict", False).eq(True).sum()
                    if "direction_conflict" in group
                    else 0
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)
