from __future__ import annotations

import math

import numpy as np
import pandas as pd

from ..evidence_quality import QUALITY_COLUMNS, add_prefixed_evidence_quality


def _legacy_to_current(table: pd.DataFrame) -> pd.DataFrame:
    required = {
        "RBP_norm",
        "OR_ENCORI",
        "FDR_ENCORI",
        "dir_ENCORI",
        "OR_POSTAR3",
        "FDR_POSTAR3",
        "dir_POSTAR3",
    }
    if not required.issubset(table):
        return table
    result = table.rename(
        columns={
            "RBP_norm": "RBP",
            "OR_ENCORI": "encori_odds_ratio",
            "FDR_ENCORI": "encori_fdr",
            "dir_ENCORI": "encori_direction",
            "OR_POSTAR3": "postar3_odds_ratio",
            "FDR_POSTAR3": "postar3_fdr",
            "dir_POSTAR3": "postar3_direction",
        }
    ).copy()
    result["encori_significant"] = result["encori_fdr"] < 0.05
    result["postar3_significant"] = result["postar3_fdr"] < 0.05
    return result


def _database_quality(table: pd.DataFrame, database: str) -> pd.DataFrame:
    prefix = f"{database}_"
    qualified = add_prefixed_evidence_quality(table, prefix)
    return qualified[[f"{prefix}{column}" for column in QUALITY_COLUMNS]].rename(
        columns={f"{prefix}{column}": column for column in QUALITY_COLUMNS}
    )


def consensus_selection_audit(
    cross_database: pd.DataFrame,
    *,
    top_n: int = 10,
    directions: tuple[str, ...] = ("enriched", "depleted"),
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    table = _legacy_to_current(cross_database).copy()
    required = {
        "RBP",
        "encori_odds_ratio",
        "encori_fdr",
        "encori_direction",
        "postar3_odds_ratio",
        "postar3_fdr",
        "postar3_direction",
    }
    missing = required - set(table)
    if missing:
        raise ValueError(f"Cross-database table is missing: {sorted(missing)}")
    if top_n < 1:
        raise ValueError("top_n must be positive")

    for database in ("encori", "postar3"):
        quality = _database_quality(table, database)
        for column in QUALITY_COLUMNS:
            table[f"{database}_{column}"] = quality[column]

    candidates = []
    for requested_direction in directions:
        label = requested_direction.capitalize()
        subset = table[
            pd.to_numeric(table["encori_fdr"], errors="coerce").lt(fdr_threshold)
            & pd.to_numeric(table["postar3_fdr"], errors="coerce").lt(fdr_threshold)
            & table["encori_direction"].astype(str).str.capitalize().eq(label)
            & table["postar3_direction"].astype(str).str.capitalize().eq(label)
        ].copy()
        if subset.empty:
            continue
        subset["arithmetic_mean_or"] = subset[
            ["encori_odds_ratio", "postar3_odds_ratio"]
        ].astype(float).mean(axis=1)
        positive = (
            subset["encori_odds_ratio"].astype(float).gt(0)
            & subset["postar3_odds_ratio"].astype(float).gt(0)
        )
        subset["geometric_mean_or"] = np.where(
            positive,
            np.exp(
                np.log(
                    subset[["encori_odds_ratio", "postar3_odds_ratio"]].astype(float)
                ).mean(axis=1)
            ),
            math.nan,
        )
        ascending = requested_direction == "depleted"
        subset = subset.sort_values(
            ["arithmetic_mean_or", "RBP"],
            ascending=[ascending, True],
            kind="mergesort",
        )
        subset.insert(0, "direction", requested_direction)
        subset["sensitivity_rank"] = np.arange(1, len(subset) + 1)
        subset["selection_inference_eligible"] = (
            subset["encori_inference_eligible"]
            & subset["postar3_inference_eligible"]
        )
        subset["primary_rank"] = pd.Series(pd.NA, index=subset.index, dtype="Int64")
        eligible_indexes = subset.index[subset["selection_inference_eligible"]]
        subset.loc[eligible_indexes, "primary_rank"] = np.arange(
            1, len(eligible_indexes) + 1
        )
        subset["selected_for_primary_string"] = (
            subset["selection_inference_eligible"]
            & subset["primary_rank"].le(top_n).fillna(False)
        )
        subset["selection_status"] = np.select(
            [
                subset["selected_for_primary_string"],
                ~subset["selection_inference_eligible"],
            ],
            [
                "selected_primary",
                "sensitivity_only_quality_excluded",
            ],
            default="primary_eligible_below_top_n_cutoff",
        )
        subset["selection_evidence_quality"] = np.where(
            subset["selection_inference_eligible"],
            "primary_eligible",
            "sensitivity_only",
        )
        subset["selection_quality_reasons"] = subset.apply(
            lambda row: " | ".join(
                item
                for item in (
                    (
                        f"encori:{row['encori_inference_exclusion_reasons']}"
                        if not bool(row["encori_inference_eligible"])
                        else ""
                    ),
                    (
                        f"postar3:{row['postar3_inference_exclusion_reasons']}"
                        if not bool(row["postar3_inference_eligible"])
                        else ""
                    ),
                )
                if item
            ),
            axis=1,
        )
        subset["selection_rule"] = (
            "ENCORI and POSTAR3 raw FDR below threshold with concordant direction; "
            "primary STRING selection additionally requires both estimates to be "
            "inference eligible and ranks by arithmetic mean OR"
        )
        candidates.append(subset)
    if not candidates:
        return pd.DataFrame()
    return pd.concat(candidates, ignore_index=True)


def select_consensus_rbps(
    cross_database: pd.DataFrame,
    *,
    top_n: int = 10,
    directions: tuple[str, ...] = ("enriched", "depleted"),
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    audit = consensus_selection_audit(
        cross_database,
        top_n=top_n,
        directions=directions,
        fdr_threshold=fdr_threshold,
    )
    if audit.empty:
        return audit
    selected = audit[audit["selected_for_primary_string"]].copy()
    selected.insert(1, "rank", selected["primary_rank"].astype(int))
    return selected.reset_index(drop=True)
