from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class ChromIntervals:
    starts: np.ndarray
    ends: np.ndarray
    maximum_end: np.ndarray


IntervalIndex = Mapping[str, ChromIntervals]


def build_interval_index(
    intervals: Mapping[str, Iterable[tuple[int, int]]],
    merge: bool = False,
) -> dict[str, ChromIntervals]:
    result: dict[str, ChromIntervals] = {}
    for chrom, values in intervals.items():
        ordered = sorted((int(start), int(end)) for start, end in values if end > start)
        if merge:
            compact: list[list[int]] = []
            for start, end in ordered:
                if not compact or start > compact[-1][1]:
                    compact.append([start, end])
                else:
                    compact[-1][1] = max(compact[-1][1], end)
            ordered = [(start, end) for start, end in compact]
        if not ordered:
            continue
        starts = np.asarray([item[0] for item in ordered], dtype=np.int64)
        ends = np.asarray([item[1] for item in ordered], dtype=np.int64)
        result[chrom] = ChromIntervals(starts, ends, np.maximum.accumulate(ends))
    return result


def group_positions(chromosomes: Sequence[object]) -> dict[str, np.ndarray]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, chrom in enumerate(chromosomes):
        grouped[str(chrom)].append(index)
    return {
        chrom: np.asarray(indexes, dtype=np.int64)
        for chrom, indexes in grouped.items()
    }


def overlap_mask(
    chromosomes: Sequence[object],
    positions: Sequence[int],
    intervals: IntervalIndex,
    window: int,
) -> np.ndarray:
    if window < 0:
        raise ValueError("window must be non-negative")
    positions_array = np.asarray(positions, dtype=np.int64)
    mask = np.zeros(len(positions_array), dtype=bool)

    for chrom, indexes in group_positions(chromosomes).items():
        indexed = intervals.get(chrom)
        if indexed is None or len(indexed.starts) == 0:
            continue
        site_positions = positions_array[indexes]
        right = np.searchsorted(
            indexed.starts,
            site_positions + window + 1,
            side="left",
        )
        candidates = right > 0
        if not np.any(candidates):
            continue
        furthest_end = indexed.maximum_end[right[candidates] - 1]
        mask[indexes[candidates]] = furthest_end > site_positions[candidates] - window
    return mask
