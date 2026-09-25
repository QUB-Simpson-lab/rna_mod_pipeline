from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np
import pandas as pd

from .stats import direction


MIN_GENE_CLUSTERS = 20
EVIDENCE_QUALITY_POLICY_ID = "quality_v1"

QUALITY_COLUMNS = (
    "estimate_eligible",
    "inference_eligible",
    "evidence_quality",
    "evidence_quality_reasons",
    "inference_exclusion_reasons",
    "input_completeness",
    "estimate_quality",
    "gene_cluster_quality",
    "leave_one_gene_out_quality",
    "single_gene_influence_status",
    "quality_diagnostics_available",
)


def _missing(value: object) -> bool:
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _optional_bool(value: object) -> bool | None:
    if _missing(value):
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and value in {0, 1}:
        return bool(value)
    if isinstance(value, (float, np.floating)) and value in {0.0, 1.0}:
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n", ""}:
        return False
    return None


def _optional_float(value: object) -> float | None:
    if _missing(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _value(record: Mapping[str, object], name: str) -> object:
    return record.get(name, np.nan)


def classify_evidence_quality(
    record: Mapping[str, object],
    *,
    odds_column: str = "odds_ratio",
    fdr_column: str = "fdr",
    minimum_gene_clusters: int = MIN_GENE_CLUSTERS,
) -> dict[str, object]:
    if minimum_gene_clusters < 2:
        raise ValueError("minimum_gene_clusters must be at least 2")

    odds = _optional_float(_value(record, odds_column))
    fdr = _optional_float(_value(record, fdr_column))
    partial = _optional_bool(_value(record, "partial_input"))
    explicit_estimable = _optional_bool(_value(record, "estimable"))
    explicit_boundary = _optional_bool(_value(record, "boundary_estimate"))
    estimate_status = str(_value(record, "estimate_status")).strip().lower()

    boundary = bool(explicit_boundary) or (
        odds is not None and (odds == 0 or math.isinf(odds))
    ) or estimate_status in {"boundary", "boundary_estimate"}
    finite_positive = odds is not None and math.isfinite(odds) and odds > 0
    estimable = explicit_estimable is not False and finite_positive and not boundary

    core_reasons: list[str] = []
    if partial is True:
        core_reasons.append("partial_input")
    if boundary:
        core_reasons.append("boundary_estimate")
    elif explicit_estimable is False or not estimable:
        core_reasons.append("non_estimable_estimate")
    if odds is None or not math.isfinite(odds):
        if "boundary_estimate" not in core_reasons:
            core_reasons.append("non_finite_odds_ratio")
    elif odds <= 0 and "boundary_estimate" not in core_reasons:
        core_reasons.append("non_positive_odds_ratio")
    if fdr is None or not math.isfinite(fdr):
        core_reasons.append("fdr_unavailable")
    elif not 0 <= fdr <= 1:
        core_reasons.append("invalid_fdr")

    estimate_eligible = not core_reasons

    low_clusters = _optional_bool(_value(record, "low_gene_cluster_count"))
    n_clusters = _optional_float(
        _value(record, "n_gene_clusters_for_robust_variance")
    )
    if n_clusters is not None and not math.isfinite(n_clusters):
        n_clusters = None
    cluster_available = low_clusters is not None or n_clusters is not None
    cluster_failure = bool(low_clusters) or (
        n_clusters is not None
        and math.isfinite(n_clusters)
        and n_clusters < minimum_gene_clusters
    )
    if not cluster_available:
        cluster_quality = "not_available"
    elif cluster_failure:
        cluster_quality = "low_gene_cluster_count"
    else:
        cluster_quality = "adequate"

    stable = _optional_bool(
        _value(record, "direction_stable_after_each_gene_removed")
    )
    reversals = _optional_float(
        _value(record, "leave_one_gene_out_direction_reversals")
    )
    loo_boundaries = _optional_float(
        _value(record, "leave_one_gene_out_boundaries")
    )
    if reversals is not None and not math.isfinite(reversals):
        reversals = None
    if loo_boundaries is not None and not math.isfinite(loo_boundaries):
        loo_boundaries = None
    loo_available = any(value is not None for value in (stable, reversals, loo_boundaries))
    loo_failures: list[str] = []
    if stable is False:
        loo_failures.append("leave_one_gene_out_direction_unstable")
    if reversals is not None and reversals > 0:
        loo_failures.append("leave_one_gene_out_direction_reversal")
    if loo_boundaries is not None and loo_boundaries > 0:
        loo_failures.append("leave_one_gene_out_boundary")
    if not loo_available:
        loo_quality = "not_available"
    elif loo_failures:
        loo_quality = "unstable"
    elif stable is True and (reversals in {None, 0.0}) and (
        loo_boundaries in {None, 0.0}
    ):
        loo_quality = "stable"
    else:
        loo_quality = "partially_available_no_instability_observed"

    maximum_influence = _optional_float(
        _value(record, "maximum_absolute_gene_influence")
    )
    influence_invalid = bool(
        maximum_influence is not None
        and math.isfinite(maximum_influence)
        and maximum_influence < 0
    )
    influence_available = bool(
        maximum_influence is not None
        and math.isfinite(maximum_influence)
        and maximum_influence >= 0
    )
    if influence_invalid:
        influence_status = "invalid_negative_value"
    elif not influence_available:
        influence_status = "not_available"
    elif loo_failures:
        influence_status = "available_with_leave_one_gene_out_instability"
    elif loo_quality == "stable":
        influence_status = "available_direction_stable"
    else:
        influence_status = "available_without_complete_leave_one_out_status"

    robustness_reasons: list[str] = []
    if cluster_failure:
        robustness_reasons.append("low_gene_cluster_count")
    robustness_reasons.extend(loo_failures)

    inference_eligible = estimate_eligible and not robustness_reasons
    diagnostics = []
    if cluster_available:
        diagnostics.append("gene_clusters")
    if loo_available:
        diagnostics.append("leave_one_gene_out")
    if influence_available:
        diagnostics.append("single_gene_influence")

    missing_diagnostics = []
    if not cluster_available:
        missing_diagnostics.append("gene_cluster_diagnostics_unavailable")
    if not loo_available:
        missing_diagnostics.append("leave_one_gene_out_diagnostics_unavailable")
    if not influence_available:
        missing_diagnostics.append("single_gene_influence_unavailable")
    if influence_invalid:
        missing_diagnostics.append("invalid_single_gene_influence")

    if not estimate_eligible:
        evidence_quality = "not_inference_eligible"
    elif robustness_reasons:
        evidence_quality = "sensitivity_only_robustness_warning"
    elif cluster_available and loo_quality == "stable" and influence_available:
        evidence_quality = "primary_eligible_robust"
    elif diagnostics:
        evidence_quality = "primary_eligible_partial_diagnostics"
    else:
        evidence_quality = "primary_eligible_core_only"

    if partial is True:
        input_completeness = "partial"
    elif partial is False:
        input_completeness = "complete"
    else:
        input_completeness = "not_reported_assumed_complete"

    if boundary:
        estimate_quality = "boundary"
    elif not estimable:
        estimate_quality = "non_estimable"
    elif fdr is None or not math.isfinite(fdr):
        estimate_quality = "finite_estimate_without_fdr"
    elif not 0 <= fdr <= 1:
        estimate_quality = "finite_estimate_with_invalid_fdr"
    else:
        estimate_quality = "finite_estimable"

    exclusion_reasons = [*core_reasons, *robustness_reasons]
    quality_reasons = [*exclusion_reasons, *missing_diagnostics]
    return {
        "estimate_eligible": bool(estimate_eligible),
        "inference_eligible": bool(inference_eligible),
        "evidence_quality": evidence_quality,
        "evidence_quality_reasons": ";".join(dict.fromkeys(quality_reasons)),
        "inference_exclusion_reasons": ";".join(
            dict.fromkeys(exclusion_reasons)
        ),
        "input_completeness": input_completeness,
        "estimate_quality": estimate_quality,
        "gene_cluster_quality": cluster_quality,
        "leave_one_gene_out_quality": loo_quality,
        "single_gene_influence_status": influence_status,
        "quality_diagnostics_available": (
            ";".join(diagnostics) if diagnostics else "none"
        ),
    }


def add_evidence_quality(
    frame: pd.DataFrame,
    *,
    odds_column: str | None = None,
    fdr_column: str | None = None,
    minimum_gene_clusters: int = MIN_GENE_CLUSTERS,
) -> pd.DataFrame:
    result = frame.copy()
    if odds_column is None:
        odds_column = "mh_or" if "mh_or" in result else "odds_ratio"
    if fdr_column is None:
        fdr_column = next(
            (
                name
                for name in (
                    "fdr_within_database_region",
                    "fdr_within_database",
                    "fdr_within_analysis",
                    "fdr",
                )
                if name in result
            ),
            "fdr",
        )
    if odds_column not in result:
        raise ValueError(f"Evidence-quality input lacks {odds_column}")
    if fdr_column not in result:
        raise ValueError(f"Evidence-quality input lacks {fdr_column}")
    if result.empty:
        for column in QUALITY_COLUMNS:
            result[column] = pd.Series(
                dtype=bool if column in {"estimate_eligible", "inference_eligible"} else object
            )
        return result
    records = [
        classify_evidence_quality(
            record,
            odds_column=odds_column,
            fdr_column=fdr_column,
            minimum_gene_clusters=minimum_gene_clusters,
        )
        for record in result.to_dict("records")
    ]
    quality = pd.DataFrame(records, index=result.index)
    for column in QUALITY_COLUMNS:
        result[column] = quality[column]
    return result


def add_prefixed_evidence_quality(
    frame: pd.DataFrame,
    prefix: str,
    *,
    odds_suffix: str = "odds_ratio",
    fdr_suffix: str = "fdr",
    minimum_gene_clusters: int = MIN_GENE_CLUSTERS,
) -> pd.DataFrame:
    result = frame.copy()
    odds_column = f"{prefix}{odds_suffix}"
    fdr_column = f"{prefix}{fdr_suffix}"
    if odds_column not in result or fdr_column not in result:
        raise ValueError(
            f"Prefixed evidence-quality input lacks {odds_column} or {fdr_column}"
        )
    source = pd.DataFrame(
        {
            "odds_ratio": pd.to_numeric(result[odds_column], errors="coerce"),
            "fdr": pd.to_numeric(result[fdr_column], errors="coerce"),
        },
        index=result.index,
    )
    for column in (
        "partial_input",
        "estimable",
        "boundary_estimate",
        "n_gene_clusters_for_robust_variance",
        "low_gene_cluster_count",
        "direction_stable_after_each_gene_removed",
        "leave_one_gene_out_direction_reversals",
        "leave_one_gene_out_boundaries",
        "maximum_absolute_gene_influence",
        "estimate_status",
    ):
        candidate = f"{prefix}{column}"
        if candidate in result:
            source[column] = result[candidate]
    quality = add_evidence_quality(
        source,
        odds_column="odds_ratio",
        fdr_column="fdr",
        minimum_gene_clusters=minimum_gene_clusters,
    )
    for column in QUALITY_COLUMNS:
        result[f"{prefix}{column}"] = quality[column]
    return result


def add_significance_calls(
    frame: pd.DataFrame,
    *,
    odds_column: str | None = None,
    fdr_column: str | None = None,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Add parallel raw-statistical and quality-qualified significance calls."""
    if not 0 < fdr_threshold <= 1:
        raise ValueError("fdr_threshold must lie in (0, 1]")
    result = add_evidence_quality(
        frame,
        odds_column=odds_column,
        fdr_column=fdr_column,
    )
    active_odds = odds_column or ("mh_or" if "mh_or" in result else "odds_ratio")
    active_fdr = fdr_column or next(
        (
            name
            for name in (
                "fdr_within_database_region",
                "fdr_within_database",
                "fdr_within_analysis",
                "fdr",
            )
            if name in result
        ),
        "fdr",
    )
    directions = []
    for odds, fdr in zip(result[active_odds], result[active_fdr]):
        odds_value = _optional_float(odds)
        fdr_value = _optional_float(fdr)
        directions.append(
            direction(
                odds_value if odds_value is not None else math.nan,
                fdr_value if fdr_value is not None else math.nan,
                fdr_threshold,
            )
        )
    result["statistical_direction"] = directions
    result["statistically_significant"] = result[
        "statistical_direction"
    ].isin(["Enriched", "Depleted"])
    result["significant"] = (
        result["statistically_significant"] & result["inference_eligible"]
    )
    result["primary_direction"] = result["statistical_direction"]
    result.loc[
        result["statistically_significant"] & ~result["inference_eligible"],
        "primary_direction",
    ] = "Quality excluded"
    return result


def evidence_quality_manifest(
    frame: pd.DataFrame,
    *,
    fdr_threshold: float = 0.05,
) -> dict[str, object]:
    """Return compact, JSON-safe policy and row-count provenance."""
    if not 0 < fdr_threshold <= 1:
        raise ValueError("fdr_threshold must lie in (0, 1]")
    inference = (
        frame["inference_eligible"].map(_optional_bool).fillna(False).astype(bool)
        if "inference_eligible" in frame
        else pd.Series(False, index=frame.index, dtype=bool)
    )
    raw = (
        frame["statistically_significant"]
        .map(_optional_bool)
        .fillna(False)
        .astype(bool)
        if "statistically_significant" in frame
        else pd.Series(False, index=frame.index, dtype=bool)
    )
    primary = (
        frame["significant"].map(_optional_bool).fillna(False).astype(bool)
        if "significant" in frame
        else raw & inference
    )

    def counts(column: str) -> dict[str, int]:
        if column not in frame:
            return {"not_reported": int(len(frame))}
        values = frame[column].fillna("missing").astype(str)
        return {str(key): int(value) for key, value in values.value_counts().items()}

    return {
        "policy_id": EVIDENCE_QUALITY_POLICY_ID,
        "minimum_gene_clusters": MIN_GENE_CLUSTERS,
        "requires_complete_input_when_reported": True,
        "requires_finite_positive_nonboundary_odds_ratio": True,
        "requires_valid_fdr": True,
        "requires_no_observed_leave_one_gene_out_instability": True,
        "missing_diagnostics_policy": (
            "primary_eligible_core_only; missing completeness is labelled "
            "not_reported_assumed_complete"
        ),
        "single_gene_influence_policy": (
            "reported without an arbitrary magnitude cutoff; observed leave-one-gene-out "
            "reversal or boundary determines exclusion"
        ),
        "fdr_threshold": fdr_threshold,
        "n_rows_total": int(len(frame)),
        "n_inference_eligible": int(inference.sum()),
        "n_raw_statistically_significant": int(raw.sum()),
        "n_primary_significant": int(primary.sum()),
        "n_raw_significant_quality_excluded": int((raw & ~inference).sum()),
        "input_completeness_counts": counts("input_completeness"),
        "diagnostics_availability_counts": counts(
            "quality_diagnostics_available"
        ),
    }
