from __future__ import annotations

import gzip
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from rna_mod_pipeline.binding.catalog import load_catalog
from rna_mod_pipeline.config import RBP_ALIASES
from rna_mod_pipeline.provenance import portable_path

from .config import DATABASE_LABELS
from .opportunities import OpportunityDesign


POSTAR_EXCLUDED = {"RBP_OCCUPANCY"}


def canonical_rbp(name: object) -> str:
    cleaned = str(name).strip().upper()
    aliases = {str(key).upper(): str(value).upper() for key, value in RBP_ALIASES.items()}
    return aliases.get(cleaned, cleaned)


def load_rbp_catalog(
    path: Path,
    databases: Sequence[str],
    requested_rbps: set[str] | None = None,
    maximum_per_database: int | None = None,
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    requested = (
        {canonical_rbp(rbp) for rbp in requested_rbps}
        if requested_rbps
        else None
    )
    lists: dict[str, list[str]] = {}
    catalog_rows = []
    for database in databases:
        entries = load_catalog(path, database, validate_count=False)
        names = {canonical_rbp(entry.canonical_rbp) for entry in entries}
        catalog_rows.extend(
            {
                "catalog_version": entry.version,
                "database_source": entry.database,
                "RBP_source": entry.canonical_rbp,
                "database_key": entry.database,
                "RBP": canonical_rbp(entry.canonical_rbp),
                "database": DATABASE_LABELS[entry.database],
                "rbp_name": entry.rbp_name,
                "canonical_rbp": entry.canonical_rbp,
                "source_names": ";".join(entry.source_names),
                "duplicate_catalog_entry": False,
            }
            for entry in entries
        )
        if requested is not None:
            names &= requested
        ordered = sorted(names)
        if maximum_per_database is not None:
            ordered = ordered[:maximum_per_database]
        if not ordered:
            raise RuntimeError(
                f"No RBPs remain for {DATABASE_LABELS[database]}"
            )
        lists[database] = ordered
    return lists, pd.DataFrame(catalog_rows)


def catalog_alias_lookup(
    catalog: pd.DataFrame, database_key: str
) -> dict[str, str]:
    aliases: dict[str, str] = {}
    subset = catalog[catalog["database_key"].eq(database_key)]
    for row in subset.itertuples(index=False):
        canonical = canonical_rbp(row.RBP)
        values = [getattr(row, "RBP_source", ""), getattr(row, "rbp_name", "")]
        values.extend(
            re.split(r"[,;|]", str(getattr(row, "source_names", "")))
        )
        for value in values:
            cleaned = str(value).strip().upper()
            if not cleaned or cleaned == "NAN":
                continue
            previous = aliases.setdefault(cleaned, canonical)
            if previous != canonical:
                raise ValueError(
                    f"Catalogue source name {cleaned} maps to both "
                    f"{previous} and {canonical}"
                )
    return aliases


def iter_bed_lines(path: Path) -> Iterable[list[str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                yield line.rstrip("\n").split("\t")


def resource_file_map(
    directory: Path, catalog_aliases: dict[str, str] | None = None
) -> dict[str, list[Path]]:
    mapping: dict[str, list[Path]] = defaultdict(list)
    for path in sorted(directory.iterdir()):
        if not path.is_file() or path.name.endswith(".tbi"):
            continue
        if not (path.name.endswith(".bed") or path.name.endswith(".bed.gz")):
            continue
        suffix = ".bed.gz" if path.name.endswith(".bed.gz") else ".bed"
        source_name = path.name[: -len(suffix)]
        canonical = (
            catalog_aliases.get(source_name.strip().upper())
            if catalog_aliases
            else None
        ) or canonical_rbp(source_name)
        mapping[canonical].append(path)
    return dict(mapping)


class IntervalCollection:
    def __init__(self, intervals: dict[object, list[tuple[int, int]]]):
        self.data: dict[object, tuple[np.ndarray, np.ndarray]] = {}
        self.n_invalid_intervals = 0
        for key, values in intervals.items():
            valid = []
            for start, end in values:
                if start < 0 or end <= start:
                    self.n_invalid_intervals += 1
                else:
                    valid.append((int(start), int(end)))
            if not valid:
                continue
            valid.sort()
            starts = np.asarray([start for start, _ in valid], dtype=np.int64)
            ends = np.asarray([end for _, end in valid], dtype=np.int64)
            self.data[key] = (starts, np.maximum.accumulate(ends))
        self.n_intervals = sum(len(starts) for starts, _ in self.data.values())

    def overlap_mask(
        self,
        chroms: np.ndarray,
        positions: np.ndarray,
        strands: np.ndarray,
        window: int,
        strand_mode: str,
    ) -> np.ndarray:
        mask = np.zeros(len(positions), dtype=bool)
        groups: dict[object, list[int]] = defaultdict(list)
        for index, (chrom, strand) in enumerate(zip(chroms, strands)):
            key = (
                (str(chrom), str(strand))
                if strand_mode == "same"
                else str(chrom)
            )
            groups[key].append(index)
        for key, indexes in groups.items():
            indexed = self.data.get(key)
            if indexed is None:
                continue
            starts, prefix_max_end = indexed
            rows = np.asarray(indexes, dtype=np.int64)
            positions_here = positions[rows]
            right = np.searchsorted(
                starts, positions_here + window + 1, side="left"
            )
            possible = right > 0
            if not np.any(possible):
                continue
            candidate = right[possible] - 1
            mask[rows[possible]] = (
                prefix_max_end[candidate] > positions_here[possible] - window
            )
        return mask


@dataclass(frozen=True)
class LoadedExposure:
    exposure: np.ndarray
    status: str
    n_source_intervals: int
    resource_row: dict


def load_interval_exposure(
    database_key: str,
    rbp: str,
    mapping: dict[str, list[Path]],
    design: OpportunityDesign,
    window: int,
    strand_mode: str,
    minimum_encori_support: int,
    project_root: Path,
) -> LoadedExposure:
    if database_key not in {"ornament", "encori"}:
        raise ValueError(database_key)
    if database_key == "ornament" and strand_mode == "same":
        raise ValueError("oRNAment resources do not contain strand information")
    paths = mapping.get(rbp, [])
    zero = np.zeros(design.n_observations, dtype=bool)
    if not paths:
        return LoadedExposure(
            zero,
            "missing",
            0,
            {
                "database": DATABASE_LABELS[database_key],
                "RBP": rbp,
                "status": "missing",
                "n_source_files": 0,
                "n_source_intervals": 0,
                "source_files": "",
                "malformed_rows": 0,
                "invalid_rows": 0,
                "filtered_rows": 0,
                "source_identity_present": False,
                "note": "No canonical or alias BED resource was found.",
            },
        )

    intervals: dict[object, list[tuple[int, int]]] = defaultdict(list)
    malformed = invalid = filtered = 0
    for path in paths:
        for fields in iter_bed_lines(path):
            try:
                if database_key == "ornament":
                    if len(fields) < 3:
                        raise ValueError
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    key: object = chrom
                else:
                    if len(fields) < 6:
                        raise ValueError
                    support = int(float(fields[4]))
                    if support < minimum_encori_support:
                        filtered += 1
                        continue
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                    if strand_mode == "same" and fields[5] not in {"+", "-"}:
                        invalid += 1
                        continue
                    key = (
                        (chrom, fields[5])
                        if strand_mode == "same"
                        else chrom
                    )
                if start < 0 or end <= start:
                    invalid += 1
                    continue
            except (ValueError, IndexError):
                malformed += 1
                continue
            intervals[key].append((start, end))

    track = IntervalCollection(intervals)
    invalid += track.n_invalid_intervals
    status = "loaded" if track.n_intervals else "empty_after_filter"
    exposure = (
        track.overlap_mask(
            design.chroms,
            design.positions,
            design.strands,
            window,
            strand_mode,
        )
        if track.n_intervals
        else zero
    )
    return LoadedExposure(
        exposure,
        status,
        track.n_intervals,
        {
            "database": DATABASE_LABELS[database_key],
            "RBP": rbp,
            "status": status,
            "n_source_files": len(paths),
            "n_source_intervals": track.n_intervals,
            "source_files": ",".join(
                portable_path(path, project_root) for path in paths
            ),
            "malformed_rows": malformed,
            "invalid_rows": invalid,
            "filtered_rows": filtered,
            "source_identity_present": True,
            "note": "",
        },
    )


class OpportunityPointIndex:
    def __init__(self, design: OpportunityDesign, strand_mode: str):
        groups: dict[object, list[tuple[int, int]]] = defaultdict(list)
        for row, (chrom, position, strand) in enumerate(
            zip(design.chroms, design.positions, design.strands)
        ):
            key = (
                (str(chrom), str(strand))
                if strand_mode == "same"
                else str(chrom)
            )
            groups[key].append((int(position), row))
        self.groups: dict[object, tuple[np.ndarray, np.ndarray]] = {}
        for key, values in groups.items():
            values.sort()
            self.groups[key] = (
                np.asarray([position for position, _ in values], dtype=np.int64),
                np.asarray([row for _, row in values], dtype=np.int64),
            )

    def mark_interval(
        self,
        mask: np.ndarray,
        chrom: str,
        strand: str,
        start: int,
        end: int,
        window: int,
        strand_mode: str,
    ) -> None:
        key: object = (chrom, strand) if strand_mode == "same" else chrom
        indexed = self.groups.get(key)
        if indexed is None:
            return
        positions, rows = indexed
        left = int(np.searchsorted(positions, start - window, side="left"))
        right = int(np.searchsorted(positions, end + window, side="left"))
        if right > left:
            mask[rows[left:right]] = True


@dataclass
class Postar3Matrix:
    path: Path
    shape: tuple[int, int]
    interval_counts: Counter
    seen_counts: Counter
    filtered_counts: Counter
    invalid_counts: Counter
    metadata: dict[str, dict[str, set[str]]]
    physical_rows: int
    partial_scan: bool

    def open(self, mode: str = "r") -> np.memmap:
        return np.memmap(self.path, dtype=np.uint8, mode=mode, shape=self.shape)


def stream_postar3(
    path: Path,
    requested_rbps: Sequence[str],
    design: OpportunityDesign,
    window: int,
    strand_mode: str,
    cell_types: set[str] | None,
    methods: set[str] | None,
    max_rows: int | None,
    temporary_directory: Path,
    catalog_aliases: dict[str, str] | None = None,
) -> Postar3Matrix:
    requested = list(requested_rbps)
    rbp_to_row = {rbp: index for index, rbp in enumerate(requested)}
    required_bytes = len(requested) * design.n_observations
    if required_bytes + 512 * 1024**2 > shutil.disk_usage(temporary_directory).free:
        raise RuntimeError(
            "Insufficient temporary space for the POSTAR3 opportunity matrix"
        )
    matrix_path = temporary_directory / "postar3_opportunity_masks.uint8"
    masks = np.memmap(
        matrix_path,
        dtype=np.uint8,
        mode="w+",
        shape=(len(requested), design.n_observations),
    )
    masks[:] = 0
    points = OpportunityPointIndex(design, strand_mode)
    intervals: Counter = Counter()
    seen: Counter = Counter()
    filtered: Counter = Counter()
    invalid: Counter = Counter()
    metadata = {
        rbp: {"cell_types": set(), "methods": set(), "accessions": set()}
        for rbp in requested
    }
    partial = False
    physical_rows = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if max_rows is not None and line_number > max_rows:
                partial = True
                break
            physical_rows = line_number
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                invalid["__MALFORMED__"] += 1
                continue
            source_rbp = fields[5].strip().upper()
            rbp = (
                catalog_aliases.get(source_rbp)
                if catalog_aliases
                else None
            ) or canonical_rbp(fields[5])
            if rbp in POSTAR_EXCLUDED or rbp not in rbp_to_row:
                continue
            seen[rbp] += 1
            method_field = fields[6]
            method = method_field.split(",")[0]
            cell_type = fields[7]
            if cell_types and cell_type not in cell_types:
                filtered[rbp] += 1
                continue
            if methods and method not in methods:
                filtered[rbp] += 1
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                invalid[rbp] += 1
                continue
            if start < 0 or end <= start or (
                strand_mode == "same" and fields[4] not in {"+", "-"}
            ):
                invalid[rbp] += 1
                continue
            points.mark_interval(
                masks[rbp_to_row[rbp]],
                fields[0],
                fields[4],
                start,
                end,
                window,
                strand_mode,
            )
            intervals[rbp] += 1
            metadata[rbp]["cell_types"].add(cell_type)
            metadata[rbp]["methods"].add(method_field)
            metadata[rbp]["accessions"].update(
                value.strip()
                for value in re.split(r"[;,]", fields[8])
                if value.strip()
            )
    masks.flush()
    del masks
    return Postar3Matrix(
        path=matrix_path,
        shape=(len(requested), design.n_observations),
        interval_counts=intervals,
        seen_counts=seen,
        filtered_counts=filtered,
        invalid_counts=invalid,
        metadata=metadata,
        physical_rows=physical_rows,
        partial_scan=partial,
    )
