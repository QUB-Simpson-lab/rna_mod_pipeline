from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

from .binding.catalog import EXPECTED_CATALOG_COUNTS, load_catalog
from .binding.sources import binding_source_files
from .config import canonical_rbp
from .gui.resources import (
    RESOURCE_KINDS,
    RESOURCE_LABELS,
    ResourceProfile,
    ResourceProfileStore,
    fasta_index_checks,
)
from .io import open_text


CORE_RESOURCES = ("reference_fasta", "annotation_gtf", "rbp_catalog")
WORKFLOW_SCRIPTS = (
    "phase1_filter.py",
    "phase1_drach.py",
    "phase1_metagene.py",
    "loose_overlap.py",
    "transcript_region_overlap.py",
    "plot_enrichment.py",
    "cross_database.py",
    "integrate_expression.py",
    "knockrbp_validation.py",
    "string_analysis.py",
    "compare_modifications.py",
    "compare_datasets.py",
)


@dataclass(frozen=True)
class DoctorCheck:
    status: str
    check: str
    details: str

    @property
    def failed(self) -> bool:
        return self.status == "FAIL"


def _check(label: str, operation) -> DoctorCheck:
    try:
        details = operation()
        return DoctorCheck("PASS", label, str(details or "OK"))
    except Exception as exc:
        return DoctorCheck("FAIL", label, f"{type(exc).__name__}: {exc}")


def _resource_status(
    profile: ResourceProfile,
    key: str,
    required: set[str],
) -> DoctorCheck:
    value = profile.path(key)
    label = RESOURCE_LABELS[key]
    if value is None:
        if key in required:
            return DoctorCheck("FAIL", label, "not configured")
        return DoctorCheck("SKIP", label, "not configured; not required for this check")
    kind = RESOURCE_KINDS[key]
    exists = value.is_dir() if kind == "directory" else value.is_file()
    if not exists:
        return DoctorCheck("FAIL", label, f"configured path does not exist: {value}")
    if kind == "file" and value.stat().st_size == 0:
        return DoctorCheck("FAIL", label, f"configured file is empty: {value}")
    return DoctorCheck("PASS", label, str(value))


def _catalog_check(profile: ResourceProfile) -> str:
    catalog = profile.require("rbp_catalog")
    counts = {
        database: len(load_catalog(catalog, database))
        for database in EXPECTED_CATALOG_COUNTS
    }
    return ", ".join(f"{database}={count}" for database, count in counts.items())


def _binding_check(profile: ResourceProfile, database: str) -> str:
    catalog = profile.require("rbp_catalog")
    entries = load_catalog(catalog, database)
    key = {
        "ornament": "ornament_dir",
        "encori": "encori_dir",
        "postar3": "postar3",
    }[database]
    source = profile.require(key)
    files = binding_source_files(database, source, entries)
    empty = [path for path in files if path.stat().st_size == 0]
    if empty:
        raise ValueError(f"{len(empty)} resolved binding files are empty")
    return f"{len(files)} source file(s) resolve for {len(entries)} RBPs"


def _knockrbp_check(profile: ResourceProfile) -> str:
    directory = profile.require("knockrbp_dir")
    metadata_path = profile.require("knockrbp_metadata")
    metadata = pd.read_csv(metadata_path, sep="\t")
    required = {"dataset_id", "rbp", "cell_line", "has_pvalues"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"metadata is missing columns: {sorted(missing)}")
    if metadata["dataset_id"].duplicated().any():
        raise ValueError("metadata contains duplicate dataset_id values")
    files = tuple(sorted(directory.glob("*_degs.json")))
    if not files:
        raise ValueError("no *_degs.json files were found")
    indexed = metadata.set_index("dataset_id")
    filename_pattern = re.compile(r"(.+)_DataSet_(\d+)_degs\.json$")
    observed: set[str] = set()
    for path in files:
        match = filename_pattern.fullmatch(path.name)
        if not match:
            raise ValueError(f"unexpected DEG filename: {path.name}")
        dataset_id = f"DataSet_{match.group(2)}"
        if dataset_id not in indexed.index:
            raise ValueError(f"{dataset_id} has no metadata row")
        file_rbp = canonical_rbp(match.group(1)).upper()
        metadata_rbp = canonical_rbp(indexed.loc[dataset_id, "rbp"]).upper()
        if file_rbp != metadata_rbp:
            raise ValueError(f"{dataset_id} RBP differs between metadata and filename")
        observed.add(dataset_id)
    if set(indexed.index.astype(str)) != observed:
        raise ValueError("metadata dataset IDs and DEG filenames do not match")
    return f"{len(metadata)} metadata rows; {len(files)} DEG files"


def _postar3_schema_check(profile: ResourceProfile) -> str:
    path = profile.require("postar3")
    inspected = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            inspected += 1
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                if inspected == 1:
                    continue
                raise ValueError(f"non-numeric coordinates at row {line_number}")
            if start < 0 or end <= start:
                raise ValueError(f"invalid interval at row {line_number}")
            if not fields[0] or not fields[5]:
                raise ValueError(f"missing chromosome or RBP at row {line_number}")
            return f"usable 10-column-or-wider record found at row {line_number}"

    raise ValueError("no usable 10-column POSTAR3 binding record was found")


def _code_check(code_root: Path) -> str:
    missing = [name for name in WORKFLOW_SCRIPTS if not (code_root / "scripts" / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing workflow scripts: {missing}")
    return f"{len(WORKFLOW_SCRIPTS)} workflow entry points present"


def run_resource_doctor(
    project_root: str | Path,
    *,
    profile_path: str | Path | None = None,
    required_resources: Iterable[str] = CORE_RESOURCES,
) -> list[DoctorCheck]:
    """Validate configured handoff resources without reading legacy results."""
    root = Path(project_root).expanduser().resolve()
    code_root = Path(__file__).resolve().parents[2]
    store = ResourceProfileStore(
        root,
        path=profile_path,
        code_root=code_root,
    )
    try:
        profile = store.load().resolved(root)
    except Exception as exc:
        return [
            DoctorCheck(
                "FAIL",
                "Resource profile",
                f"{type(exc).__name__}: {exc}",
            )
        ]

    required = set(required_resources)
    unknown = required.difference(RESOURCE_KINDS)
    if unknown:
        raise KeyError(f"Unknown required resources: {sorted(unknown)}")
    profile_source = (
        str(store.path)
        if store.path.is_file()
        else "automatic discovery (no saved profile yet)"
    )
    checks = [
        DoctorCheck("PASS", "Resource profile", profile_source),
        _check("Workflow code", lambda: _code_check(code_root)),
    ]
    if root.is_dir() and os.access(root, os.W_OK):
        checks.append(DoctorCheck("PASS", "Project workspace", f"writable: {root}"))
    else:
        checks.append(DoctorCheck("FAIL", "Project workspace", f"not writable: {root}"))
    checks.extend(
        _resource_status(profile, key, required)
        for key in RESOURCE_KINDS
    )
    if profile.reference_fasta and profile.reference_fasta.is_file():
        index_errors, index_warnings = fasta_index_checks(profile.reference_fasta)
        checks.extend(
            DoctorCheck("FAIL", "Reference FASTA index", message)
            for message in index_errors
        )
        checks.extend(
            DoctorCheck("WARN", "Reference FASTA index", message)
            for message in index_warnings
        )
        if not index_errors and not index_warnings:
            checks.append(
                DoctorCheck(
                    "PASS",
                    "Reference FASTA index",
                    str(Path(f"{profile.reference_fasta}.fai")),
                )
            )
    if profile.rbp_catalog and profile.rbp_catalog.is_file():
        checks.append(_check("RBP catalogue panels", lambda: _catalog_check(profile)))
    for database, key in (
        ("ornament", "ornament_dir"),
        ("encori", "encori_dir"),
        ("postar3", "postar3"),
    ):
        if profile.path(key) is not None and profile.rbp_catalog is not None:
            checks.append(
                _check(
                    f"{database} catalogue/source mapping",
                    lambda database=database: _binding_check(profile, database),
                )
            )
    if profile.postar3 is not None and profile.postar3.is_file():
        checks.append(_check("POSTAR3 table schema", lambda: _postar3_schema_check(profile)))
    if profile.knockrbp_dir is not None and profile.knockrbp_metadata is not None:
        checks.append(_check("KnockRBP mapping", lambda: _knockrbp_check(profile)))
    if not profile.genome_build:
        checks.append(
            DoctorCheck(
                "WARN",
                "Genome-build label",
                "not recorded in the resource profile",
            )
        )
    if not profile.annotation_release:
        checks.append(
            DoctorCheck(
                "WARN",
                "Annotation-release label",
                "not recorded in the resource profile",
            )
        )
    return checks
