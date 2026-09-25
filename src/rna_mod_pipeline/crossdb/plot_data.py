from __future__ import annotations

import math
from typing import Mapping

import numpy as np
import pandas as pd

from .io import DATABASES


COVERAGE_ORDER = ("both_present", "first_only", "second_only")
INFERENCE_ORDER = (
    "both_significant_same_direction",
    "both_significant_opposite_direction",
    "first_only_significant",
    "second_only_significant",
    "neither_significant",
    "any_quality_excluded",
    "any_boundary_estimate",
    "any_non_estimable",
)


def _as_bool(value: object) -> bool:
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes"}:
            return True
        if lowered in {"false", "0", "no", ""}:
            return False
    return bool(value) if pd.notna(value) else False


def _as_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def database_state(row: pd.Series, database: str, threshold: float) -> str:
    if not _as_bool(row.get(f"{database}_present", False)):
        return "not_covered"

    odds = _as_float(row.get(f"{database}_odds_ratio"))
    fdr = _as_float(row.get(f"{database}_fdr"))
    boundary = _as_bool(row.get(f"{database}_boundary_estimate", False))
    if boundary or odds == 0 or math.isinf(odds):
        return "boundary_estimate"

    explicit_estimable = row.get(f"{database}_estimable")
    if pd.notna(explicit_estimable) and not _as_bool(explicit_estimable):
        return "tested_non_estimable"
    if not math.isfinite(odds) or odds <= 0 or not math.isfinite(fdr):
        return "tested_non_estimable"
    inference_eligible = row.get(f"{database}_inference_eligible")
    if pd.notna(inference_eligible) and not _as_bool(inference_eligible):
        return "quality_excluded"
    if fdr >= threshold or odds == 1:
        return "tested_not_significant"
    return "significant_enriched" if odds > 1 else "significant_depleted"


def add_plot_states(combined: pd.DataFrame, threshold: float) -> pd.DataFrame:
    frame = combined.copy()
    for database in DATABASES:
        if f"{database}_present" not in frame:
            frame[f"{database}_present"] = False
        frame[f"{database}_plot_state"] = frame.apply(
            database_state, axis=1, database=database, threshold=threshold
        )
        odds = pd.to_numeric(
            frame.get(f"{database}_odds_ratio", pd.Series(np.nan, index=frame.index)),
            errors="coerce",
        )
        log_odds = pd.Series(np.nan, index=frame.index, dtype=float)
        finite_positive = np.isfinite(odds) & odds.gt(0)
        log_odds.loc[finite_positive] = np.log2(odds.loc[finite_positive])
        frame[f"{database}_log2_odds_ratio"] = log_odds
    return frame


def pairwise_plot_source(
    combined: pd.DataFrame,
    pairwise: Mapping[str, pd.DataFrame],
    threshold: float,
) -> pd.DataFrame:
    states = add_plot_states(combined, threshold)
    rows: list[dict[str, object]] = []
    for pair in pairwise:
        first, second = pair.split("_vs_", maxsplit=1)
        for row in states.itertuples(index=False):
            values = row._asdict()
            first_present = bool(values[f"{first}_present"])
            second_present = bool(values[f"{second}_present"])
            if first_present and second_present:
                coverage = "both_present"
            elif first_present:
                coverage = "first_only"
            elif second_present:
                coverage = "second_only"
            else:
                continue

            first_state = values[f"{first}_plot_state"]
            second_state = values[f"{second}_plot_state"]
            inference = "not_both_present"
            if coverage == "both_present":
                pair_states = {first_state, second_state}
                if "boundary_estimate" in pair_states:
                    inference = "any_boundary_estimate"
                elif "tested_non_estimable" in pair_states:
                    inference = "any_non_estimable"
                elif "quality_excluded" in pair_states:
                    inference = "any_quality_excluded"
                else:
                    first_sig = str(first_state).startswith("significant_")
                    second_sig = str(second_state).startswith("significant_")
                    if first_sig and second_sig:
                        inference = (
                            "both_significant_same_direction"
                            if first_state == second_state
                            else "both_significant_opposite_direction"
                        )
                    elif first_sig:
                        inference = "first_only_significant"
                    elif second_sig:
                        inference = "second_only_significant"
                    else:
                        inference = "neither_significant"
            rows.append(
                {
                    "pair": pair,
                    "first_database": first,
                    "second_database": second,
                    "RBP": values["RBP"],
                    "coverage_state": coverage,
                    "inference_state": inference,
                    "first_state": first_state,
                    "second_state": second_state,
                    "first_odds_ratio": values.get(f"{first}_odds_ratio", np.nan),
                    "second_odds_ratio": values.get(f"{second}_odds_ratio", np.nan),
                    "first_fdr": values.get(f"{first}_fdr", np.nan),
                    "second_fdr": values.get(f"{second}_fdr", np.nan),
                    "first_log2_odds_ratio": values.get(
                        f"{first}_log2_odds_ratio", np.nan
                    ),
                    "second_log2_odds_ratio": values.get(
                        f"{second}_log2_odds_ratio", np.nan
                    ),
                    "eligible_for_primary_correlation": (
                        coverage == "both_present"
                        and first_state
                        not in {
                            "quality_excluded",
                            "boundary_estimate",
                            "tested_non_estimable",
                        }
                        and second_state
                        not in {
                            "quality_excluded",
                            "boundary_estimate",
                            "tested_non_estimable",
                        }
                    ),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "pair",
            "first_database",
            "second_database",
            "RBP",
            "coverage_state",
            "inference_state",
            "first_state",
            "second_state",
            "first_odds_ratio",
            "second_odds_ratio",
            "first_fdr",
            "second_fdr",
            "first_log2_odds_ratio",
            "second_log2_odds_ratio",
            "eligible_for_primary_correlation",
        ],
    )


def pairwise_count_source(source: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for pair, group in source.groupby("pair", sort=False):
        for category in COVERAGE_ORDER:
            rows.append(
                {
                    "pair": pair,
                    "dimension": "coverage",
                    "category": category,
                    "count": int(group["coverage_state"].eq(category).sum()),
                }
            )
        shared = group[group["coverage_state"].eq("both_present")]
        for category in INFERENCE_ORDER:
            rows.append(
                {
                    "pair": pair,
                    "dimension": "inference_among_shared",
                    "category": category,
                    "count": int(shared["inference_state"].eq(category).sum()),
                }
            )
    return pd.DataFrame(rows, columns=["pair", "dimension", "category", "count"])


def select_pairwise_top(source: pd.DataFrame, top_n: int) -> pd.DataFrame:
    eligible = source[
        source["inference_state"].eq("both_significant_same_direction")
        & np.isfinite(source["first_log2_odds_ratio"])
        & np.isfinite(source["second_log2_odds_ratio"])
    ].copy()
    eligible["mean_log2_odds_ratio"] = eligible[
        ["first_log2_odds_ratio", "second_log2_odds_ratio"]
    ].mean(axis=1)
    eligible["rank_magnitude"] = eligible["mean_log2_odds_ratio"].abs()
    pair_order = {
        pair: index for index, pair in enumerate(source["pair"].drop_duplicates())
    }
    eligible["_pair_order"] = eligible["pair"].map(pair_order)
    eligible = eligible.sort_values(
        ["_pair_order", "rank_magnitude", "RBP"],
        ascending=[True, False, True],
        kind="mergesort",
    )
    selected = eligible.groupby("pair", sort=False).head(top_n).copy()
    selected["rank_within_pair"] = (
        selected.groupby("pair", sort=False).cumcount() + 1
    )
    return selected.drop(columns="_pair_order").reset_index(drop=True)


def triple_plot_source(
    combined: pd.DataFrame,
    threshold: float,
    top_n: int,
) -> pd.DataFrame:
    states = add_plot_states(combined, threshold)
    columns = [
        "RBP",
        "triple_significant",
        "triple_direction_concordant",
        "direction_conflict",
        "consensus_direction",
    ]
    for database in DATABASES:
        columns.extend(
            [
                f"{database}_plot_state",
                f"{database}_odds_ratio",
                f"{database}_fdr",
                f"{database}_log2_odds_ratio",
            ]
        )
    available = [column for column in columns if column in states]
    source = states[states.get("triple_significant", False)].loc[:, available].copy()
    if source.empty:
        source["selected_for_plot"] = pd.Series(dtype=bool)
        source["rank_magnitude"] = pd.Series(dtype=float)
        return source
    log_columns = [f"{database}_log2_odds_ratio" for database in DATABASES]
    source["plot_eligible"] = np.isfinite(source[log_columns]).all(axis=1)
    source["rank_magnitude"] = source[log_columns].abs().mean(axis=1)
    ordered = source.sort_values(
        ["plot_eligible", "rank_magnitude", "RBP"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    selected = set(ordered.loc[ordered["plot_eligible"], "RBP"].head(top_n))
    source["selected_for_plot"] = source["RBP"].isin(selected)
    return source.sort_values(
        ["selected_for_plot", "rank_magnitude", "RBP"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def validation_summary_source(
    combined: pd.DataFrame,
    threshold: float,
    databases: tuple[str, ...],
) -> pd.DataFrame:
    states = add_plot_states(combined, threshold)
    rows: list[dict[str, object]] = []
    maximum = len(databases)
    for count in range(maximum + 1):
        rows.append(
            {
                "dimension": "databases_present",
                "category": str(count),
                "count": int(states["n_databases_present"].eq(count).sum()),
            }
        )
        rows.append(
            {
                "dimension": "databases_significant",
                "category": str(count),
                "count": int(states["n_databases_significant"].eq(count).sum()),
            }
        )
    consensus = {
        "concordant_enriched": (
            states["direction_concordant_consensus"]
            & states["consensus_direction"].eq("Enriched")
        ),
        "concordant_depleted": (
            states["direction_concordant_consensus"]
            & states["consensus_direction"].eq("Depleted")
        ),
        "direction_conflict": states["direction_conflict"],
    }
    for category, mask in consensus.items():
        rows.append(
            {
                "dimension": "multi_database_consensus",
                "category": category,
                "count": int(mask.sum()),
            }
        )
    for state in (
        "significant_enriched",
        "significant_depleted",
        "tested_not_significant",
        "quality_excluded",
        "boundary_estimate",
        "tested_non_estimable",
        "not_covered",
    ):
        count = sum(
            int(states[f"{database}_plot_state"].eq(state).sum())
            for database in databases
        )
        rows.append(
            {
                "dimension": "database_rbp_state",
                "category": state,
                "count": count,
            }
        )
    return pd.DataFrame(rows, columns=["dimension", "category", "count"])


def membership_count_source(
    combined: pd.DataFrame,
    databases: tuple[str, ...],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for dimension, suffix in (
        ("coverage_membership", "present"),
        ("significant_membership", "significant"),
    ):
        memberships = []
        for row in combined.itertuples(index=False):
            values = row._asdict()
            members = [
                database
                for database in databases
                if _as_bool(values.get(f"{database}_{suffix}", False))
            ]
            memberships.append("+".join(members) if members else "none")
        counts = pd.Series(memberships, dtype=str).value_counts()
        for membership, count in counts.items():
            rows.append(
                {
                    "dimension": dimension,
                    "membership": membership,
                    "n_databases": 0 if membership == "none" else membership.count("+") + 1,
                    "count": int(count),
                }
            )
    return pd.DataFrame(
        rows, columns=["dimension", "membership", "n_databases", "count"]
    ).sort_values(
        ["dimension", "n_databases", "membership"],
        ascending=[True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
