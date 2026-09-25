from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
from typing import TYPE_CHECKING, Sequence

import numpy as np
import pandas as pd

from .association import analyse_all_regions
from .config import DATABASE_LABELS
from .outputs import save_case_binding_matrix
from .resources import (
    catalog_alias_lookup,
    load_interval_exposure,
    resource_file_map,
    stream_postar3,
)

if TYPE_CHECKING:
    from .opportunities import OpportunityDesign
    from .workflow import RunOptions


@dataclass(frozen=True)
class DatabaseRunResult:
    result_rows: list[dict]
    influence_rows: list[dict]
    region_rows: list[dict]
    resource_rows: list[dict]
    resource_inputs: list[Path]
    partial_postar: bool


def _database_filename(database_key: str) -> str:
    return {
        "ornament": "rbp_enrichment_results.tsv",
        "encori": "encori_enrichment_results.tsv",
        "postar3": "postar3_enrichment_results.tsv",
    }[database_key]


def _process_interval_database(
    database_key: str,
    rbps: Sequence[str],
    directory: Path,
    design: OpportunityDesign,
    options: RunOptions,
    output_dir: Path,
    partial_run: bool,
    catalog: pd.DataFrame,
) -> tuple[list[dict], list[dict], list[dict], list[dict], list[Path]]:
    mapping = resource_file_map(
        directory, catalog_alias_lookup(catalog, database_key)
    )
    results: list[dict] = []
    influence: list[dict] = []
    regional: list[dict] = []
    resources: list[dict] = []
    input_files: list[Path] = []
    matrix_rbps: list[str] = []
    packed_rows: list[np.ndarray] = []
    database = DATABASE_LABELS[database_key]
    for rbp in rbps:
        loaded = load_interval_exposure(
            database_key,
            rbp,
            mapping,
            design,
            options.window,
            options.strand_mode,
            options.minimum_encori_support,
            options.project_root,
        )
        input_files.extend(mapping.get(rbp, []))
        row, gene_rows, region_rows = analyse_all_regions(
            design,
            loaded.exposure,
            options.analysis_id,
            options.modification,
            database,
            rbp,
            loaded.n_source_intervals,
            partial_run
            or loaded.status == "missing"
            or (
                loaded.status == "empty_after_filter"
                and not (
                    database_key == "encori"
                    and options.minimum_encori_support != 1
                )
            ),
        )
        results.append(row)
        influence.extend(gene_rows)
        regional.extend(region_rows)
        resources.append(loaded.resource_row)
        matrix_rbps.append(rbp)
        packed_rows.append(
            np.packbits(loaded.exposure[design.is_case], bitorder="little")
        )
    save_case_binding_matrix(
        output_dir,
        database_key,
        matrix_rbps,
        packed_rows,
        design.case_site_ids,
        design.n_cases,
    )
    return results, influence, regional, resources, input_files


def _process_postar3(
    rbps: Sequence[str],
    design: OpportunityDesign,
    options: RunOptions,
    output_dir: Path,
    temporary_directory: Path,
    partial_run: bool,
    catalog: pd.DataFrame,
) -> tuple[list[dict], list[dict], list[dict], list[dict], bool]:
    assert options.postar3 is not None
    streamed = stream_postar3(
        options.postar3,
        rbps,
        design,
        options.window,
        options.strand_mode,
        set(options.postar_cell_types) if options.postar_cell_types else None,
        set(options.postar_methods) if options.postar_methods else None,
        options.max_postar_rows,
        temporary_directory,
        catalog_alias_lookup(catalog, "postar3"),
    )
    results: list[dict] = []
    influence: list[dict] = []
    regional: list[dict] = []
    resources: list[dict] = []
    packed_rows: list[np.ndarray] = []
    masks = streamed.open()
    try:
        for index, rbp in enumerate(rbps):
            count = int(streamed.interval_counts[rbp])
            if count:
                status = "partial_loaded" if streamed.partial_scan else "loaded"
            elif streamed.partial_scan:
                status = "partial_missing"
            elif streamed.seen_counts[rbp]:
                status = "empty_after_filter"
            else:
                status = "missing"
            exposure = np.asarray(masks[index], dtype=bool)
            row, gene_rows, region_rows = analyse_all_regions(
                design,
                exposure,
                options.analysis_id,
                options.modification,
                "POSTAR3",
                rbp,
                count,
                partial_run
                or status in {"missing", "partial_missing", "partial_loaded"}
                or (
                    status == "empty_after_filter"
                    and not (
                        options.postar_cell_types or options.postar_methods
                    )
                ),
            )
            results.append(row)
            influence.extend(gene_rows)
            regional.extend(region_rows)
            packed_rows.append(
                np.packbits(exposure[design.is_case], bitorder="little")
            )
            meta = streamed.metadata[rbp]
            resources.append(
                {
                    "database": "POSTAR3",
                    "RBP": rbp,
                    "status": status,
                    "n_source_files": 1,
                    "n_source_intervals": count,
                    "source_files": options.postar3.relative_to(
                        options.project_root
                    ).as_posix()
                    if options.postar3.is_relative_to(options.project_root)
                    else options.postar3.as_posix(),
                    "malformed_rows": 0,
                    "invalid_rows": int(streamed.invalid_counts[rbp]),
                    "filtered_rows": int(streamed.filtered_counts[rbp]),
                    "source_identity_present": bool(streamed.seen_counts[rbp]),
                    "note": (
                        f"{len(meta['cell_types'])} cell types; "
                        f"{len(meta['methods'])} method labels; "
                        f"{len(meta['accessions'])} accessions."
                    ),
                }
            )
    finally:
        del masks
        streamed.path.unlink(missing_ok=True)
    save_case_binding_matrix(
        output_dir,
        "postar3",
        rbps,
        packed_rows,
        design.case_site_ids,
        design.n_cases,
    )
    resources.append(
        {
            "database": "POSTAR3",
            "RBP": "__STREAM_SUMMARY__",
            "status": "partial" if streamed.partial_scan else "complete",
            "n_source_files": 1,
            "n_source_intervals": int(sum(streamed.interval_counts.values())),
            "source_files": options.postar3.relative_to(
                options.project_root
            ).as_posix()
            if options.postar3.is_relative_to(options.project_root)
            else options.postar3.as_posix(),
            "malformed_rows": int(
                streamed.invalid_counts["__MALFORMED__"]
            ),
            "invalid_rows": int(sum(streamed.invalid_counts.values())),
            "filtered_rows": int(sum(streamed.filtered_counts.values())),
            "source_identity_present": True,
            "note": f"{streamed.physical_rows} physical rows scanned.",
        }
    )
    return results, influence, regional, resources, streamed.partial_scan


def process_databases(
    lists: dict[str, list[str]],
    catalog: pd.DataFrame,
    design: OpportunityDesign,
    options: RunOptions,
    output: Path,
    partial_run: bool,
) -> DatabaseRunResult:
    result_rows: list[dict] = []
    influence_rows: list[dict] = []
    region_rows: list[dict] = []
    resource_rows: list[dict] = []
    resource_inputs: list[Path] = []
    partial_postar = False
    for database_key in options.databases:
        if database_key in {"ornament", "encori"}:
            directory = (
                options.ornament_dir
                if database_key == "ornament"
                else options.encori_dir
            )
            assert directory is not None
            rows, influences, regions, resources, source_files = (
                _process_interval_database(
                    database_key,
                    lists[database_key],
                    directory,
                    design,
                    options,
                    output,
                    partial_run,
                    catalog,
                )
            )
            resource_inputs.extend(source_files)
        else:
            rows, influences, regions, resources, partial_postar = (
                _process_postar3(
                    lists[database_key],
                    design,
                    options,
                    output,
                    output,
                    partial_run,
                    catalog,
                )
            )
        result_rows.extend(rows)
        influence_rows.extend(influences)
        region_rows.extend(regions)
        resource_rows.extend(resources)
    return DatabaseRunResult(
        result_rows,
        influence_rows,
        region_rows,
        resource_rows,
        resource_inputs,
        partial_postar,
    )
