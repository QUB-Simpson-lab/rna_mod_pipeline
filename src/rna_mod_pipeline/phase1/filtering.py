from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import FILTERED_COLUMNS
from ..io import bedmethyl_chunks, write_filtered_sites


COVERAGE_BINS = np.linspace(0, 500, 51)
COVERAGE_PLOT_MAX = 200
FRACTION_BINS = np.linspace(0, 100, 51)
COVERAGE_THRESHOLDS = (10, 20, 30, 50)
FRACTION_THRESHOLDS = (0.0, 10.0, 20.0, 50.0, 90.0)


@dataclass(frozen=True)
class FilterSummary:
    total_sites: int
    retained_sites: int
    coverage_min: int
    fraction_min: float
    chromosome_counts: dict[str, int]
    coverage_histogram: np.ndarray
    coverage_at_or_above_plot_max: int
    coverage_above_histogram_max: int
    fraction_histogram: np.ndarray
    coverage_threshold_counts: np.ndarray
    fraction_threshold_counts: np.ndarray
    threshold_counts: np.ndarray


def filter_bedmethyl(
    input_path: str | Path,
    output_path: str | Path,
    coverage_min: int = 20,
    fraction_min: float = 20.0,
    chunk_size: int = 2_000_000,
    expected_mod_codes: frozenset[str] | None = None,
) -> tuple[pd.DataFrame, FilterSummary]:
    """Retain sites meeting the validated coverage and modification thresholds."""
    if coverage_min < 0:
        raise ValueError("coverage_min must be non-negative")
    if not 0 <= fraction_min <= 100:
        raise ValueError("fraction_min must lie between 0 and 100")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")

    retained: list[pd.DataFrame] = []
    total_sites = 0
    chromosome_counts: dict[str, int] = {}
    coverage_histogram = np.zeros(len(COVERAGE_BINS) - 1, dtype=np.int64)
    coverage_at_or_above_plot_max = 0
    coverage_above_histogram_max = 0
    fraction_histogram = np.zeros(len(FRACTION_BINS) - 1, dtype=np.int64)
    coverage_threshold_counts = np.zeros(
        len(COVERAGE_THRESHOLDS), dtype=np.int64
    )
    fraction_threshold_counts = np.zeros(
        len(FRACTION_THRESHOLDS), dtype=np.int64
    )
    threshold_counts = np.zeros(
        (len(COVERAGE_THRESHOLDS), len(FRACTION_THRESHOLDS)),
        dtype=np.int64,
    )
    for chunk in bedmethyl_chunks(input_path, chunk_size):
        if expected_mod_codes is not None:
            observed_codes = {
                str(value).strip().lower()
                for value in chunk["mod_code"].dropna().unique()
            }
            expected_codes = {code.lower() for code in expected_mod_codes}
            unexpected = observed_codes.difference(expected_codes)
            if unexpected:
                raise ValueError(
                    "bedMethyl modification code does not match the selected "
                    f"analysis: observed {sorted(unexpected)}, expected one of "
                    f"{sorted(expected_codes)}"
                )
        total_sites += len(chunk)
        for chromosome, count in chunk["chrom"].astype(str).value_counts().items():
            chromosome_counts[chromosome] = (
                chromosome_counts.get(chromosome, 0) + int(count)
            )
        coverage = chunk["Nvalid_cov"].to_numpy()
        fraction = chunk["fraction_modified"].to_numpy()
        coverage_histogram += np.histogram(coverage, bins=COVERAGE_BINS)[0]
        coverage_at_or_above_plot_max += int(
            np.sum(coverage >= COVERAGE_PLOT_MAX)
        )
        coverage_above_histogram_max += int(
            np.sum(coverage > COVERAGE_BINS[-1])
        )
        fraction_histogram += np.histogram(
            fraction[fraction > 0],
            bins=FRACTION_BINS,
        )[0]
        fraction_masks = []
        for column, threshold_fraction in enumerate(FRACTION_THRESHOLDS):
            fraction_mask = (
                fraction > 0
                if threshold_fraction == 0
                else fraction >= threshold_fraction
            )
            fraction_masks.append(fraction_mask)
            fraction_threshold_counts[column] += int(np.sum(fraction_mask))
        for row, threshold_coverage in enumerate(COVERAGE_THRESHOLDS):
            coverage_mask = coverage >= threshold_coverage
            coverage_threshold_counts[row] += int(np.sum(coverage_mask))
            for column, fraction_mask in enumerate(fraction_masks):
                threshold_counts[row, column] += int(
                    np.sum(coverage_mask & fraction_mask)
                )
        mask = (
            chunk["Nvalid_cov"].ge(coverage_min)
            & chunk["fraction_modified"].ge(fraction_min)
        )
        if mask.any():
            retained.append(chunk.loc[mask, FILTERED_COLUMNS])

    if retained:
        filtered = pd.concat(retained, ignore_index=True)
    else:
        filtered = pd.DataFrame(columns=FILTERED_COLUMNS)

    for column in ("start", "end", "Nvalid_cov", "Nmod", "Ncanonical"):
        if column in filtered:
            filtered[column] = pd.to_numeric(filtered[column], downcast="integer")

    write_filtered_sites(filtered, output_path)
    summary = FilterSummary(
        total_sites=total_sites,
        retained_sites=len(filtered),
        coverage_min=coverage_min,
        fraction_min=fraction_min,
        chromosome_counts=chromosome_counts,
        coverage_histogram=coverage_histogram,
        coverage_at_or_above_plot_max=coverage_at_or_above_plot_max,
        coverage_above_histogram_max=coverage_above_histogram_max,
        fraction_histogram=fraction_histogram,
        coverage_threshold_counts=coverage_threshold_counts,
        fraction_threshold_counts=fraction_threshold_counts,
        threshold_counts=threshold_counts,
    )
    return filtered, summary


def threshold_grid(
    input_path: str | Path,
    coverage_thresholds: tuple[int, ...] = COVERAGE_THRESHOLDS,
    fraction_thresholds: tuple[float, ...] = FRACTION_THRESHOLDS,
    chunk_size: int = 2_000_000,
) -> tuple[np.ndarray, int]:
    """Count sites passing each coverage-by-fraction threshold pair."""
    counts = np.zeros(
        (len(coverage_thresholds), len(fraction_thresholds)),
        dtype=np.int64,
    )
    total = 0
    for chunk in bedmethyl_chunks(input_path, chunk_size):
        total += len(chunk)
        coverage = chunk["Nvalid_cov"].to_numpy()
        fraction = chunk["fraction_modified"].to_numpy()
        for row, threshold_coverage in enumerate(coverage_thresholds):
            coverage_mask = coverage >= threshold_coverage
            for column, threshold_fraction in enumerate(fraction_thresholds):
                fraction_mask = (
                    fraction > 0
                    if threshold_fraction == 0
                    else fraction >= threshold_fraction
                )
                counts[row, column] += int(np.sum(coverage_mask & fraction_mask))
    return counts, total
