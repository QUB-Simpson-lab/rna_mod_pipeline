from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from ..config import canonical_rbp
from ..io import require_file

DATABASES = ("encori", "postar3")


def _valid_genes(values: pd.Series) -> pd.Series:
    genes = values.astype("string").str.strip()
    return genes.where(
        genes.notna() & ~genes.str.lower().isin({"", "nan", "none", "na", "<na>"})
    )


def _resource_availability(
    path: str | Path | None,
    database: str,
) -> dict[str, tuple[bool, str, str]]:
    if path is None:
        return {}
    frame = pd.read_csv(require_file(path, f"{database} resource table"), sep="\t")
    if "database" in frame:
        frame = frame[
            frame["database"].astype(str).str.lower() == database
        ].copy()
    if "RBP" not in frame:
        raise ValueError(f"{path} must contain RBP")
    frame["canonical_rbp"] = frame["RBP"].map(
        lambda value: canonical_rbp(value).upper()
    )
    if frame["canonical_rbp"].duplicated().any():
        raise ValueError(f"{path} contains duplicate canonical RBP rows")
    count_column = next(
        (
            column
            for column in (
                "n_source_intervals",
                "total_binding_sites",
            )
            if column in frame
        ),
        None,
    )
    result = {}
    for row in frame.itertuples(index=False):
        record = row._asdict()
        name = str(record["canonical_rbp"])
        count = (
            int(pd.to_numeric(record[count_column], errors="coerce"))
            if count_column and pd.notna(record[count_column])
            else 0
        )
        raw_status = str(
            record.get("resource_status", record.get("status", "not_reported"))
        )
        status = raw_status.strip().lower()
        usable = count > 0 and status not in {
            "resource_empty",
            "empty",
            "empty_after_filter",
            "missing",
        }
        basis = (
            f"{Path(path).name}:{count_column or 'no_interval_count'}"
        )
        result[name] = usable, raw_status, basis
    return result


def load_loose_targets(
    annotated_paths: dict[str, str | Path],
    rbps: Iterable[str],
    *,
    context: str = "all",
    chunk_size: int = 25_000,
    resource_paths: dict[str, str | Path] | None = None,
) -> tuple[dict[str, dict[str, object]], set[str]]:
    requested = sorted({canonical_rbp(rbp).upper() for rbp in rbps})
    targets = {
        rbp: {
            "encori": set(),
            "postar3": set(),
            "union": set(),
            "intersection": set(),
            "encori_available": False,
            "postar3_available": False,
            "encori_resource_status": "not_assessed",
            "postar3_resource_status": "not_assessed",
            "encori_availability_basis": "not_assessed",
            "postar3_availability_basis": "not_assessed",
        }
        for rbp in requested
    }
    universe: set[str] = set()
    resource_maps = {
        database: _resource_availability(
            (resource_paths or {}).get(database), database
        )
        for database in DATABASES
    }
    for database, raw_path in annotated_paths.items():
        if database not in DATABASES:
            raise ValueError(f"Unknown target database: {database}")
        path = require_file(raw_path, f"{database} annotated sites")
        header = pd.read_csv(path, sep="\t", nrows=0)
        prefix = f"{database}_"
        source_columns: dict[str, list[str]] = {}
        for column in header.columns:
            if column.startswith(prefix):
                name = canonical_rbp(column[len(prefix) :]).upper()
                source_columns.setdefault(name, []).append(column)
        use = ["gene_name"] + (["region"] if "region" in header else [])
        use += [
            column
            for rbp in requested
            for column in source_columns.get(rbp, [])
        ]
        use = list(dict.fromkeys(use))
        for rbp in requested:
            columns_present = bool(source_columns.get(rbp))
            resource = resource_maps[database].get(rbp)
            if resource is None:
                targets[rbp][f"{database}_available"] = columns_present
                targets[rbp][f"{database}_resource_status"] = (
                    "annotation_column_present"
                    if columns_present
                    else "annotation_column_absent"
                )
                targets[rbp][f"{database}_availability_basis"] = (
                    "annotation_only_unverified"
                )
            else:
                usable, status, basis = resource
                targets[rbp][f"{database}_available"] = (
                    columns_present and usable
                )
                targets[rbp][f"{database}_resource_status"] = status
                targets[rbp][f"{database}_availability_basis"] = basis
        for chunk in pd.read_csv(path, sep="\t", usecols=use, chunksize=chunk_size):
            if context != "all":
                if "region" not in chunk:
                    raise ValueError(f"{path} cannot be filtered by transcript region")
                chunk = chunk[chunk["region"].astype(str).str.upper() == context.upper()]
            genes = _valid_genes(chunk["gene_name"])
            universe.update(genes.dropna().astype(str))
            for rbp in requested:
                columns = source_columns.get(rbp, [])
                if not columns or not targets[rbp][f"{database}_available"]:
                    continue
                bound = pd.concat(
                    [
                        pd.to_numeric(chunk[column], errors="coerce").fillna(0).eq(1)
                        for column in columns
                    ],
                    axis=1,
                ).any(axis=1)
                targets[rbp][database].update(genes[bound].dropna().astype(str))
    for rbp in requested:
        first, second = targets[rbp]["encori"], targets[rbp]["postar3"]
        targets[rbp]["union"] = first | second
        targets[rbp]["intersection"] = first & second
    return targets, universe


def _load_assignments(run_dir: Path, context: str) -> pd.DataFrame:
    path = require_file(run_dir / "case_site_assignments.tsv.gz", "case assignments")
    frame = pd.read_csv(path, sep="\t")
    required = {
        "case_site_id",
        "gene_name",
        "region",
        "contributes_to_stratified_analysis",
    }
    if not required.issubset(frame):
        raise ValueError(f"Case assignments missing: {sorted(required - set(frame))}")
    contributes = frame["contributes_to_stratified_analysis"].astype(str).str.lower().isin(
        {"true", "1", "yes"}
    )
    frame = frame[contributes].copy()
    frame["case_site_id"] = pd.to_numeric(frame["case_site_id"], errors="raise").astype(int)
    if frame["case_site_id"].duplicated().any():
        raise ValueError(f"Duplicate contributing case_site_id values in {path}")
    if context != "all":
        frame = frame[frame["region"].astype(str).str.upper() == context.upper()].copy()
    frame["gene_name"] = _valid_genes(frame["gene_name"])
    return frame.dropna(subset=["gene_name"])


def _matrix_masks(
    run_dir: Path,
    database: str,
    requested: set[str],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    rows_path = require_file(
        run_dir / f"case_binding_matrix_{database}_rows.tsv", f"{database} row table"
    )
    matrix_path = require_file(
        run_dir / f"case_binding_matrix_{database}.npz", f"{database} binding matrix"
    )
    rows = pd.read_csv(rows_path, sep="\t")
    with np.load(matrix_path, allow_pickle=False) as archive:
        packed = archive["packed_matrix"]
        rbps = archive["rbps"].astype(str)
        site_ids = archive["case_site_ids"].astype(np.int64)
        n_sites = int(archive["n_case_sites"][0])
        bitorder = str(archive["bitorder"][0])
    if bitorder != "little" or packed.shape != (len(rbps), math.ceil(n_sites / 8)):
        raise ValueError(f"Invalid packed matrix: {matrix_path}")
    if list(rows.columns) != ["matrix_row", "RBP"]:
        raise ValueError(f"Invalid matrix row table: {rows_path}")
    row_numbers = pd.to_numeric(rows["matrix_row"], errors="raise").astype(int)
    if not np.array_equal(row_numbers.to_numpy(), np.arange(len(rows))):
        raise ValueError(f"Non-contiguous matrix rows: {rows_path}")
    canonical = [canonical_rbp(value).upper() for value in rbps]
    row_rbps = [canonical_rbp(value).upper() for value in rows["RBP"]]
    if row_rbps != canonical:
        raise ValueError(f"Matrix row table disagrees with archive: {matrix_path}")
    if len(canonical) != len(set(canonical)):
        raise ValueError(f"Duplicate canonical RBP in {matrix_path}")
    if len(site_ids) != n_sites or len(set(site_ids.tolist())) != n_sites:
        raise ValueError(f"Invalid or duplicate case-site IDs: {matrix_path}")
    masks = {}
    for index, rbp in enumerate(canonical):
        if rbp in requested:
            masks[rbp] = np.unpackbits(
                packed[index], bitorder="little", count=n_sites
            ).astype(bool)
    return site_ids, masks


def load_transcript_targets(
    run_dir: str | Path,
    rbps: Iterable[str],
    *,
    context: str = "all",
) -> tuple[dict[str, dict[str, object]], set[str]]:
    run = Path(run_dir).resolve()
    requested = {canonical_rbp(rbp).upper() for rbp in rbps}
    all_assignments = _load_assignments(run, "all")
    selected_assignments = _load_assignments(run, context)
    selected_ids = set(selected_assignments["case_site_id"].astype(int))
    targets = {
        rbp: {
            "encori": set(),
            "postar3": set(),
            "union": set(),
            "intersection": set(),
            "encori_available": False,
            "postar3_available": False,
            "encori_resource_status": "not_assessed",
            "postar3_resource_status": "not_assessed",
            "encori_availability_basis": "not_assessed",
            "postar3_availability_basis": "not_assessed",
        }
        for rbp in sorted(requested)
    }
    report_path = run / "resource_loading_report.tsv"
    resource_maps = {
        database: _resource_availability(
            report_path if report_path.is_file() else None,
            database,
        )
        for database in DATABASES
    }
    for database in DATABASES:
        site_ids, masks = _matrix_masks(run, database, requested)
        assignment_ids = set(all_assignments["case_site_id"].astype(int))
        if assignment_ids != set(site_ids.tolist()):
            raise ValueError(
                f"{database} matrix case IDs do not exactly match assignments"
            )
        aligned = all_assignments.set_index("case_site_id").reindex(site_ids)
        if aligned["gene_name"].isna().any():
            raise ValueError(f"{database} matrix and case assignments are not aligned")
        context_mask = np.asarray([int(site_id) in selected_ids for site_id in site_ids])
        for rbp, mask in masks.items():
            resource = resource_maps[database].get(rbp)
            if resource is None:
                targets[rbp][f"{database}_available"] = True
                targets[rbp][f"{database}_resource_status"] = (
                    "matrix_row_present"
                )
                targets[rbp][f"{database}_availability_basis"] = (
                    "matrix_only_unverified"
                )
            else:
                usable, status, basis = resource
                targets[rbp][f"{database}_available"] = usable
                targets[rbp][f"{database}_resource_status"] = status
                targets[rbp][f"{database}_availability_basis"] = basis
            if targets[rbp][f"{database}_available"]:
                targets[rbp][database] = set(
                    aligned.loc[mask & context_mask, "gene_name"].astype(str)
                )
    for rbp in targets:
        first, second = targets[rbp]["encori"], targets[rbp]["postar3"]
        targets[rbp]["union"] = first | second
        targets[rbp]["intersection"] = first & second
    return targets, set(selected_assignments["gene_name"].astype(str))
