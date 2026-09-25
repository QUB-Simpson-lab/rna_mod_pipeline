from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class FisherResult:
    odds_ratio: float
    p_value: float


def bh_adjust(values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjustment while preserving missing values."""
    pvalues = np.asarray(values, dtype=float)
    adjusted = np.full(len(pvalues), np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(pvalues))
    if len(finite) == 0:
        return adjusted

    observed = pvalues[finite]
    if np.any((observed < 0) | (observed > 1)):
        raise ValueError("p-values must lie between 0 and 1")

    order = np.argsort(observed, kind="mergesort")
    ranked = observed[order]
    scaled = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    monotone = np.minimum.accumulate(scaled[::-1])[::-1]
    restored = np.empty(len(ranked), dtype=float)
    restored[order] = np.clip(monotone, 0, 1)
    adjusted[finite] = restored
    return adjusted


def fisher_test(
    case_bound: int,
    case_total: int,
    control_bound: int,
    control_total: int,
    alternative: str = "two-sided",
) -> FisherResult:
    counts = (case_bound, case_total, control_bound, control_total)
    if any(int(value) != value or value < 0 for value in counts):
        raise ValueError("Fisher counts must be non-negative integers")
    if case_bound > case_total or control_bound > control_total:
        raise ValueError("Bound counts cannot exceed their totals")
    if case_total == 0 or control_total == 0:
        return FisherResult(math.nan, math.nan)

    table = [
        [case_bound, case_total - case_bound],
        [control_bound, control_total - control_bound],
    ]
    odds_ratio, p_value = stats.fisher_exact(table, alternative=alternative)
    return FisherResult(float(odds_ratio), float(p_value))


def direction(odds_ratio: float, fdr: float, threshold: float = 0.05) -> str:
    if math.isnan(odds_ratio) or not math.isfinite(fdr) or odds_ratio < 0:
        return "Unavailable"
    if fdr >= threshold or odds_ratio == 1:
        return "Not significant"
    return "Enriched" if odds_ratio > 1 else "Depleted"


def safe_log2(value: float) -> float:
    if not math.isfinite(value) or value <= 0:
        return math.nan
    return math.log2(value)
