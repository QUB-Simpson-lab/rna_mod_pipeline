from __future__ import annotations

import math
from functools import reduce
from typing import Mapping

import numpy as np
import pandas as pd

from ..evidence_quality import add_evidence_quality, add_prefixed_evidence_quality
from ..stats import direction
from .io import DATABASES

def _prefix(frame: pd.DataFrame, database: str, threshold: float) -> pd.DataFrame:
    table = add_evidence_quality(
        frame,
        odds_column="odds_ratio",
        fdr_column="fdr",
    ).drop(columns=["database"]).copy()
    table["present"] = True
    table["statistical_direction"] = [
        direction(or_value, fdr, threshold)
        for or_value, fdr in zip(table["odds_ratio"], table["fdr"])
    ]
    table["statistically_significant"] = table["statistical_direction"].isin(
        ["Enriched", "Depleted"]
    )
    table["significant"] = (
        table["statistically_significant"] & table["inference_eligible"]
    )
    table["direction"] = table["statistical_direction"]
    table["primary_direction"] = table["statistical_direction"]
    table.loc[
        table["statistically_significant"] & ~table["inference_eligible"],
        "primary_direction",
    ] = "Quality excluded"
    return table.rename(
        columns={name: f"{database}_{name}" for name in table if name != "RBP"}
    )


def _mean_effects(row: pd.Series, databases: tuple[str, ...]) -> tuple[float, float, float]:
    values = [
        float(row[f"{database}_odds_ratio"])
        for database in databases
        if bool(row.get(f"{database}_significant", False))
        and pd.notna(row.get(f"{database}_odds_ratio"))
    ]
    if not values:
        return math.nan, math.nan, math.nan
    arithmetic = float(np.mean(values))
    if any(not math.isfinite(value) or value <= 0 for value in values):
        return arithmetic, math.nan, math.nan
    mean_log = float(np.mean(np.log(values)))
    return arithmetic, mean_log, float(math.exp(mean_log))


def _mean_statistical_effects(row: pd.Series, databases: tuple[str, ...]) -> tuple[float, float, float]:
    values = [
        float(row[f"{database}_odds_ratio"])
        for database in databases
        if bool(row.get(f"{database}_statistically_significant", False))
        and pd.notna(row.get(f"{database}_odds_ratio"))
    ]
    if not values:
        return math.nan, math.nan, math.nan
    arithmetic = float(np.mean(values))
    if any(not math.isfinite(value) or value <= 0 for value in values):
        return arithmetic, math.nan, math.nan
    mean_log = float(np.mean(np.log(values)))
    return arithmetic, mean_log, float(math.exp(mean_log))


def _quality_reason(row: pd.Series, databases: tuple[str, ...]) -> str:
    reasons = []
    for database in databases:
        if not bool(row.get(f"{database}_statistically_significant", False)):
            continue
        if bool(row.get(f"{database}_inference_eligible", False)):
            continue
        reason = str(row.get(f"{database}_inference_exclusion_reasons", "")).strip()
        reasons.append(f"{database}:{reason or 'quality_excluded'}")
    return " | ".join(reasons)


def _as_bool(value: object) -> bool:
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"true", "1", "yes"}:
            return True
        if text in {"false", "0", "no", "", "nan"}:
            return False
    return bool(value) if pd.notna(value) else False


def recompute_consensus_fields(
    table: pd.DataFrame,
    *,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Rebuild per-database and aggregate consensus calls from source statistics."""
    if not 0 < fdr_threshold <= 1:
        raise ValueError("fdr_threshold must lie in (0, 1]")
    result = table.copy()
    included = tuple(
        database
        for database in DATABASES
        if {
            f"{database}_odds_ratio",
            f"{database}_fdr",
        }.issubset(result)
    )
    if len(included) < 2:
        raise ValueError(
            "Consensus recomputation requires source odds-ratio and FDR columns "
            "for at least two databases"
        )
    for database in DATABASES:
        prefix = f"{database}_"
        odds_column = f"{prefix}odds_ratio"
        fdr_column = f"{prefix}fdr"
        present_column = f"{prefix}present"
        if odds_column not in result or fdr_column not in result:
            result[present_column] = False
            result[f"{prefix}inference_eligible"] = False
            result[f"{prefix}evidence_quality"] = "not_covered"
            result[f"{prefix}inference_exclusion_reasons"] = "not_covered"
            result[f"{prefix}gene_cluster_quality"] = "not_available"
            result[f"{prefix}leave_one_gene_out_quality"] = "not_available"
            result[f"{prefix}single_gene_influence_status"] = "not_available"
            result[f"{prefix}statistically_significant"] = False
            result[f"{prefix}significant"] = False
            result[f"{prefix}statistical_direction"] = "Unavailable"
            result[f"{prefix}direction"] = "Unavailable"
            result[f"{prefix}primary_direction"] = "Unavailable"
            continue

        result = add_prefixed_evidence_quality(result, prefix)
        if present_column in table:
            result[present_column] = table[present_column].map(_as_bool)
        else:
            source_column = f"{prefix}source_file"
            if source_column in table:
                source_present = table[source_column].notna() & table[
                    source_column
                ].astype(str).str.strip().ne("")
            else:
                source_present = pd.Series(False, index=result.index)
            result[present_column] = (
                source_present
                | pd.to_numeric(result[odds_column], errors="coerce").notna()
                | pd.to_numeric(result[fdr_column], errors="coerce").notna()
            )
        result[f"{prefix}statistical_direction"] = [
            direction(
                float(odds) if pd.notna(odds) else math.nan,
                float(fdr) if pd.notna(fdr) else math.nan,
                fdr_threshold,
            )
            for odds, fdr in zip(result[odds_column], result[fdr_column])
        ]
        result[f"{prefix}statistically_significant"] = result[
            f"{prefix}statistical_direction"
        ].isin(["Enriched", "Depleted"]) & result[present_column]
        result[f"{prefix}inference_eligible"] = (
            result[f"{prefix}inference_eligible"] & result[present_column]
        )
        result[f"{prefix}significant"] = (
            result[f"{prefix}statistically_significant"]
            & result[f"{prefix}inference_eligible"]
            & result[present_column]
        )
        result[f"{prefix}direction"] = result[f"{prefix}statistical_direction"]
        result[f"{prefix}primary_direction"] = result[
            f"{prefix}statistical_direction"
        ]
        result.loc[
            result[f"{prefix}statistically_significant"]
            & ~result[f"{prefix}inference_eligible"],
            f"{prefix}primary_direction",
        ] = "Quality excluded"
        not_present = ~result[present_column]
        result.loc[not_present, f"{prefix}inference_eligible"] = False
        result.loc[not_present, f"{prefix}evidence_quality"] = "not_covered"
        result.loc[
            not_present, f"{prefix}inference_exclusion_reasons"
        ] = "not_covered"
        result.loc[not_present, f"{prefix}gene_cluster_quality"] = (
            "not_available"
        )
        result.loc[
            not_present, f"{prefix}leave_one_gene_out_quality"
        ] = "not_available"
        result.loc[
            not_present, f"{prefix}single_gene_influence_status"
        ] = "not_available"
        for column in (
            "statistical_direction",
            "direction",
            "primary_direction",
        ):
            result.loc[not_present, f"{prefix}{column}"] = "Unavailable"

    result["n_databases_present"] = result[
        [f"{database}_present" for database in DATABASES]
    ].sum(axis=1)
    result["n_databases_significant"] = result[
        [f"{database}_significant" for database in DATABASES]
    ].sum(axis=1)
    result["n_databases_statistically_significant"] = result[
        [f"{database}_statistically_significant" for database in DATABASES]
    ].sum(axis=1)
    result["n_databases_inference_eligible"] = result[
        [f"{database}_inference_eligible" for database in DATABASES]
    ].sum(axis=1)
    result["triple_significant"] = (
        result["n_databases_present"].eq(3)
        & result["n_databases_significant"].eq(3)
    )
    result["triple_statistically_significant"] = (
        result["n_databases_present"].eq(3)
        & result["n_databases_statistically_significant"].eq(3)
    )

    primary_directions = [
        {
            row[f"{database}_direction"]
            for database in DATABASES
            if row[f"{database}_significant"]
        }
        for _, row in result.iterrows()
    ]
    result["direction_concordant_consensus"] = [
        len(values) == 1 and int(count) >= 2
        for values, count in zip(
            primary_directions, result["n_databases_significant"]
        )
    ]
    result["triple_direction_concordant"] = (
        result["triple_significant"] & result["direction_concordant_consensus"]
    )
    result["direction_conflict"] = [
        len(values) > 1 for values in primary_directions
    ]
    result["consensus_direction"] = [
        next(iter(values)) if len(values) == 1 else "Conflicting" if values else "None"
        for values in primary_directions
    ]

    statistical_directions = [
        {
            row[f"{database}_statistical_direction"]
            for database in DATABASES
            if row[f"{database}_statistically_significant"]
            and row[f"{database}_present"]
        }
        for _, row in result.iterrows()
    ]
    result["statistical_direction_concordant_consensus"] = [
        len(values) == 1 and int(count) >= 2
        for values, count in zip(
            statistical_directions,
            result["n_databases_statistically_significant"],
        )
    ]
    result["statistical_direction_conflict"] = [
        len(values) > 1 for values in statistical_directions
    ]
    result["statistical_consensus_direction"] = [
        next(iter(values)) if len(values) == 1 else "Conflicting" if values else "None"
        for values in statistical_directions
    ]
    result["consensus_evidence_quality"] = np.select(
        [
            result["direction_concordant_consensus"],
            result["statistical_direction_concordant_consensus"],
        ],
        [
            "primary_eligible_consensus",
            "sensitivity_only_consensus_quality_excluded",
        ],
        default="no_direction_concordant_consensus",
    )
    result["consensus_quality_reasons"] = result.apply(
        _quality_reason,
        axis=1,
        databases=included,
    )
    effects = result.apply(lambda row: _mean_effects(row, included), axis=1)
    result[
        [
            "arithmetic_mean_significant_or",
            "mean_log_significant_or",
            "geometric_mean_significant_or",
        ]
    ] = pd.DataFrame(effects.tolist(), index=result.index)
    statistical_effects = result.apply(
        lambda row: _mean_statistical_effects(row, included), axis=1
    )
    result[
        [
            "arithmetic_mean_statistically_significant_or",
            "mean_log_statistically_significant_or",
            "geometric_mean_statistically_significant_or",
        ]
    ] = pd.DataFrame(statistical_effects.tolist(), index=result.index)
    return result


def compare_databases(
    tables: Mapping[str, pd.DataFrame],
    *,
    fdr_threshold: float = 0.05,
    allow_missing: bool = False,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    unknown = set(tables) - set(DATABASES)
    if unknown:
        raise ValueError(f"Unknown databases: {sorted(unknown)}")
    missing = [database for database in DATABASES if database not in tables]
    if missing and not allow_missing:
        raise ValueError(
            "Cross-database validation requires oRNAment, ENCORI, and POSTAR3; "
            f"missing: {', '.join(missing)}"
        )
    if len(tables) < 2:
        raise ValueError("At least two database tables are required")

    ordered = tuple(database for database in DATABASES if database in tables)
    prefixed = [_prefix(tables[database], database, fdr_threshold) for database in ordered]
    merged = reduce(lambda left, right: left.merge(right, on="RBP", how="outer"), prefixed)
    for database in DATABASES:
        present_column = f"{database}_present"
        if present_column in merged:
            merged[present_column] = merged[present_column].map(
                lambda value: bool(value) if pd.notna(value) else False
            )
        else:
            merged[present_column] = False
        if f"{database}_significant" not in merged:
            merged[f"{database}_significant"] = False
            merged[f"{database}_statistically_significant"] = False
            merged[f"{database}_inference_eligible"] = False
            merged[f"{database}_direction"] = "Unavailable"
            merged[f"{database}_statistical_direction"] = "Unavailable"
            merged[f"{database}_primary_direction"] = "Unavailable"
        else:
            merged[f"{database}_significant"] = merged[
                f"{database}_significant"
            ].map(lambda value: bool(value) if pd.notna(value) else False)
            merged[f"{database}_statistically_significant"] = merged[
                f"{database}_statistically_significant"
            ].map(lambda value: bool(value) if pd.notna(value) else False)
            merged[f"{database}_inference_eligible"] = merged[
                f"{database}_inference_eligible"
            ].map(lambda value: bool(value) if pd.notna(value) else False)
            merged[f"{database}_direction"] = merged[
                f"{database}_direction"
            ].fillna("Unavailable")
            merged[f"{database}_statistical_direction"] = merged[
                f"{database}_statistical_direction"
            ].fillna("Unavailable")
            merged[f"{database}_primary_direction"] = merged[
                f"{database}_primary_direction"
            ].fillna("Unavailable")

    merged["n_databases_present"] = merged[
        [f"{database}_present" for database in DATABASES]
    ].sum(axis=1)
    merged["n_databases_significant"] = merged[
        [f"{database}_significant" for database in DATABASES]
    ].sum(axis=1)
    merged["n_databases_statistically_significant"] = merged[
        [f"{database}_statistically_significant" for database in DATABASES]
    ].sum(axis=1)
    merged["n_databases_inference_eligible"] = merged[
        [f"{database}_inference_eligible" for database in DATABASES]
    ].sum(axis=1)
    merged["triple_significant"] = (
        merged["n_databases_present"].eq(3)
        & merged["n_databases_significant"].eq(3)
    )
    merged["triple_statistically_significant"] = (
        merged["n_databases_present"].eq(3)
        & merged["n_databases_statistically_significant"].eq(3)
    )

    significant_directions = []
    for _, row in merged.iterrows():
        significant_directions.append(
            {
                row[f"{database}_direction"]
                for database in DATABASES
                if row[f"{database}_significant"]
            }
        )
    merged["direction_concordant_consensus"] = [
        len(values) == 1 and int(n_significant) >= 2
        for values, n_significant in zip(
            significant_directions, merged["n_databases_significant"]
        )
    ]
    merged["triple_direction_concordant"] = (
        merged["triple_significant"] & merged["direction_concordant_consensus"]
    )
    merged["direction_conflict"] = [
        len(values) > 1 for values in significant_directions
    ]
    merged["consensus_direction"] = [
        next(iter(values)) if len(values) == 1 else "Conflicting" if values else "None"
        for values in significant_directions
    ]
    statistical_directions = []
    for _, row in merged.iterrows():
        statistical_directions.append(
            {
                row[f"{database}_statistical_direction"]
                for database in DATABASES
                if row[f"{database}_statistically_significant"]
            }
        )
    merged["statistical_direction_concordant_consensus"] = [
        len(values) == 1 and int(n_significant) >= 2
        for values, n_significant in zip(
            statistical_directions,
            merged["n_databases_statistically_significant"],
        )
    ]
    merged["statistical_direction_conflict"] = [
        len(values) > 1 for values in statistical_directions
    ]
    merged["statistical_consensus_direction"] = [
        next(iter(values)) if len(values) == 1 else "Conflicting" if values else "None"
        for values in statistical_directions
    ]
    merged["consensus_evidence_quality"] = np.select(
        [
            merged["direction_concordant_consensus"],
            merged["statistical_direction_concordant_consensus"],
        ],
        [
            "primary_eligible_consensus",
            "sensitivity_only_consensus_quality_excluded",
        ],
        default="no_direction_concordant_consensus",
    )
    merged["consensus_quality_reasons"] = merged.apply(
        _quality_reason, axis=1, databases=ordered
    )
    effects = merged.apply(lambda row: _mean_effects(row, ordered), axis=1)
    merged[
        ["arithmetic_mean_significant_or", "mean_log_significant_or", "geometric_mean_significant_or"]
    ] = pd.DataFrame(effects.tolist(), index=merged.index)
    statistical_effects = merged.apply(
        lambda row: _mean_statistical_effects(row, ordered), axis=1
    )
    merged[
        [
            "arithmetic_mean_statistically_significant_or",
            "mean_log_statistically_significant_or",
            "geometric_mean_statistically_significant_or",
        ]
    ] = pd.DataFrame(statistical_effects.tolist(), index=merged.index)
    merged = merged.sort_values(
        ["triple_direction_concordant", "n_databases_significant", "RBP"],
        ascending=[False, False, True],
    ).reset_index(drop=True)

    pairwise: dict[str, pd.DataFrame] = {}
    for index, first in enumerate(ordered):
        for second in ordered[index + 1 :]:
            shared = merged[
                merged[f"{first}_present"] & merged[f"{second}_present"]
            ].copy()
            shared["both_significant"] = (
                shared[f"{first}_significant"] & shared[f"{second}_significant"]
            )
            shared["both_statistically_significant"] = (
                shared[f"{first}_statistically_significant"]
                & shared[f"{second}_statistically_significant"]
            )
            shared["same_significant_direction"] = (
                shared["both_significant"]
                & shared[f"{first}_direction"].eq(shared[f"{second}_direction"])
            )
            shared["opposite_significant_direction"] = (
                shared["both_significant"]
                & ~shared[f"{first}_direction"].eq(shared[f"{second}_direction"])
            )
            columns = ["RBP"] + [
                column
                for column in shared
                if column.startswith(f"{first}_") or column.startswith(f"{second}_")
            ] + [
                "both_significant",
                "same_significant_direction",
                "opposite_significant_direction",
                "both_statistically_significant",
            ]
            pairwise[f"{first}_vs_{second}"] = shared[columns]
    return merged, pairwise
