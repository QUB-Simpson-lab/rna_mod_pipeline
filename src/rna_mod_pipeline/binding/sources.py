from __future__ import annotations

import re
from array import array
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

from ..intervals import ChromIntervals, IntervalIndex
from ..io import open_text, require_directory, require_file
from .catalog import CatalogEntry


@dataclass(frozen=True)
class BindingTrack:
    rbp_name: str
    canonical_rbp: str
    intervals: IntervalIndex
    source_names: tuple[str, ...]
    metadata: dict[str, object] = field(default_factory=dict)

    @property
    def interval_count(self) -> int:
        return sum(len(chromosome.starts) for chromosome in self.intervals.values())


def _bed_path(directory: Path, source_name: str) -> Path:
    compressed = directory / f"{source_name}.bed.gz"
    plain = directory / f"{source_name}.bed"
    if compressed.is_file():
        return compressed
    if plain.is_file():
        return plain
    raise FileNotFoundError(
        f"Binding track not found for {source_name}: expected {compressed} or {plain}"
    )


def _load_bed_compact(
    paths: list[Path],
    minimum_score: int,
    merge: bool,
    allow_zero_length: bool,
) -> dict[str, ChromIntervals]:
    buffers: dict[str, tuple[array, array]] = defaultdict(
        lambda: (array("i"), array("i"))
    )
    maximum_coordinate = np.iinfo(np.int32).max
    for path in paths:
        with open_text(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                if minimum_score > 1 and len(fields) > 4:
                    try:
                        if int(float(fields[4])) < minimum_score:
                            continue
                    except ValueError:
                        continue
                start, end = int(fields[1]), int(fields[2])
                invalid_length = (
                    end < start if allow_zero_length else end <= start
                )
                if start < 0 or invalid_length or end > maximum_coordinate:
                    continue
                starts, ends = buffers[fields[0]]
                starts.append(start)
                ends.append(end)

    result: dict[str, ChromIntervals] = {}
    for chrom, (start_buffer, end_buffer) in buffers.items():
        starts = np.frombuffer(start_buffer, dtype=np.int32)
        ends = np.frombuffer(end_buffer, dtype=np.int32)
        if len(starts) > 1 and np.any(starts[1:] < starts[:-1]):
            order = np.argsort(starts, kind="stable")
            starts = starts[order]
            ends = ends[order]
        if merge and len(starts):
            prefix = np.maximum.accumulate(ends)
            group_starts = np.r_[
                0,
                np.flatnonzero(starts[1:] > prefix[:-1]) + 1,
            ]
            starts = starts[group_starts]
            ends = np.maximum.reduceat(ends, group_starts)
        result[chrom] = ChromIntervals(
            starts=starts,
            ends=ends,
            maximum_end=np.maximum.accumulate(ends),
        )
    return result


def load_bed_directory(
    directory: str | Path,
    entries: Iterable[CatalogEntry],
    minimum_score: int = 1,
    allow_empty_resources: bool = False,
) -> list[BindingTrack]:
    return list(
        iter_bed_directory(
            directory,
            entries,
            minimum_score=minimum_score,
            allow_empty_resources=allow_empty_resources,
        )
    )


def iter_bed_directory(
    directory: str | Path,
    entries: Iterable[CatalogEntry],
    minimum_score: int = 1,
    allow_empty_resources: bool = False,
) -> Iterable[BindingTrack]:
    root = require_directory(directory, "binding-track directory")
    for entry in entries:
        paths = [_bed_path(root, source) for source in entry.source_names]
        intervals = _load_bed_compact(
            paths,
            minimum_score=minimum_score,
            merge=len(paths) > 1,
            allow_zero_length=entry.database == "ornament",
        )
        interval_count = sum(
            len(chromosome.starts) for chromosome in intervals.values()
        )
        if interval_count == 0 and not allow_empty_resources:
            raise RuntimeError(
                f"{entry.database}/{entry.rbp_name} has zero usable intervals"
            )
        track = BindingTrack(
            rbp_name=entry.rbp_name,
            canonical_rbp=entry.canonical_rbp,
            intervals=intervals,
            source_names=entry.source_names,
            metadata={
                "resource_status": (
                    "loaded" if interval_count else "resource_empty"
                )
            },
        )
        yield track
        del track, intervals


def load_postar3(
    path: str | Path,
    entries: Iterable[CatalogEntry],
    cell_types: set[str] | None = None,
    methods: set[str] | None = None,
    allow_empty_resources: bool = False,
) -> list[BindingTrack]:
    """Load the validated POSTAR3 panel in one streaming pass."""
    entries = list(entries)
    source_to_entry: dict[str, CatalogEntry] = {}
    for entry in entries:
        for source_name in entry.source_names:
            existing = source_to_entry.get(source_name)
            if existing is not None and existing.rbp_name != entry.rbp_name:
                raise ValueError(f"POSTAR3 source alias is ambiguous: {source_name}")
            source_to_entry[source_name] = entry

    intervals: dict[str, dict[str, tuple[array, array]]] = defaultdict(
        lambda: defaultdict(lambda: (array("q"), array("q")))
    )
    metadata: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {
            "methods": set(),
            "cell_types": set(),
            "accessions": set(),
        }
    )
    source_path = require_file(path, "POSTAR3 human.txt")
    with open_text(source_path) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            entry = source_to_entry.get(fields[5])
            if entry is None:
                continue
            cell_type = fields[7]
            method_field = fields[6]
            clip_method = method_field.split(",", 1)[0]
            if cell_types and cell_type not in cell_types:
                continue
            if methods and clip_method not in methods:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid POSTAR3 coordinates at row {line_number}"
                ) from exc
            if end <= start:
                continue
            starts, ends = intervals[entry.rbp_name][fields[0]]
            starts.append(start)
            ends.append(end)
            item = metadata[entry.rbp_name]
            item["methods"].add(method_field)
            item["cell_types"].add(cell_type)
            item["accessions"].update(
                value.strip()
                for value in re.split(r"[;,]", fields[8])
                if value.strip()
            )

    tracks: list[BindingTrack] = []
    for entry in entries:
        index: dict[str, ChromIntervals] = {}
        raw = intervals.pop(entry.rbp_name, {})
        for chrom, (start_buffer, end_buffer) in raw.items():
            starts = np.frombuffer(start_buffer, dtype=np.int64)
            ends = np.frombuffer(end_buffer, dtype=np.int64)
            order = np.argsort(starts, kind="stable")
            starts = starts[order]
            ends = ends[order]
            index[chrom] = ChromIntervals(
                starts=starts,
                ends=ends,
                maximum_end=np.maximum.accumulate(ends),
            )
        interval_count = sum(
            len(chromosome.starts) for chromosome in index.values()
        )
        if interval_count == 0 and not allow_empty_resources:
            raise RuntimeError(
                f"postar3/{entry.rbp_name} has zero usable intervals"
            )
        item = metadata.get(
            entry.rbp_name,
            {"methods": set(), "cell_types": set(), "accessions": set()},
        )
        tracks.append(
            BindingTrack(
                rbp_name=entry.rbp_name,
                canonical_rbp=entry.canonical_rbp,
                intervals=index,
                source_names=entry.source_names,
                metadata={
                    "resource_status": (
                        "loaded" if interval_count else "resource_empty"
                    ),
                    "n_experiments": len(item["accessions"]),
                    "n_cell_types": len(item["cell_types"]),
                    "methods": "; ".join(sorted(item["methods"])),
                },
            )
        )
    return tracks


def binding_source_files(
    database: str,
    source: str | Path,
    entries: Iterable[CatalogEntry],
) -> list[Path]:
    """Resolve the files represented by one binding-source invocation."""
    if database == "postar3":
        return [require_file(source, "POSTAR3 human.txt")]
    if database not in {"ornament", "encori"}:
        raise ValueError(f"Unsupported binding database: {database}")
    root = require_directory(source, "binding-track directory")
    paths = {
        _bed_path(root, source_name)
        for entry in entries
        for source_name in entry.source_names
    }
    return sorted(paths)


def load_binding_source(
    database: str,
    source: str | Path,
    entries: Iterable[CatalogEntry],
    minimum_clip_experiments: int = 1,
    cell_types: set[str] | None = None,
    methods: set[str] | None = None,
    allow_empty_resources: bool = False,
) -> list[BindingTrack]:
    if database == "ornament":
        return load_bed_directory(
            source,
            entries,
            minimum_score=1,
            allow_empty_resources=allow_empty_resources,
        )
    if database == "encori":
        filtered = minimum_clip_experiments > 1
        return load_bed_directory(
            source,
            entries,
            minimum_score=minimum_clip_experiments,
            allow_empty_resources=allow_empty_resources or filtered,
        )
    if database == "postar3":
        filtered = bool(cell_types or methods)
        return load_postar3(
            source,
            entries,
            cell_types=cell_types,
            methods=methods,
            allow_empty_resources=allow_empty_resources or filtered,
        )
    raise ValueError(f"Unsupported binding database: {database}")


def iter_binding_source(
    database: str,
    source: str | Path,
    entries: Iterable[CatalogEntry],
    minimum_clip_experiments: int = 1,
    cell_types: set[str] | None = None,
    methods: set[str] | None = None,
    allow_empty_resources: bool = False,
) -> Iterable[BindingTrack]:
    """Yield BED-backed tracks one at a time to bound memory use."""
    entries = list(entries)
    if database in {"ornament", "encori"}:
        minimum_score = 1 if database == "ornament" else minimum_clip_experiments
        filtered = database == "encori" and minimum_clip_experiments > 1
        yield from iter_bed_directory(
            source,
            entries,
            minimum_score=minimum_score,
            allow_empty_resources=allow_empty_resources or filtered,
        )
        return
    if database == "postar3":
        filtered = bool(cell_types or methods)
        yield from load_postar3(
            source,
            entries,
            cell_types=cell_types,
            methods=methods,
            allow_empty_resources=allow_empty_resources or filtered,
        )
        return
    raise ValueError(f"Unsupported binding database: {database}")
