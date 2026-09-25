from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from rna_mod_pipeline.provenance import (
    portable_path,
    sha256_file,
    validate_file_record,
)

from .association import bh_adjust_declared
from .config import DATABASE_KEYS, REGIONS


MAIN_COLUMNS = {
    "analysis_id",
    "modification",
    "database",
    "RBP",
    "region_filter",
    "n_case_sites",
    "n_control_sites",
    "n_strata_total",
    "case_overlaps",
    "control_overlaps",
    "mh_numerator",
    "mh_denominator",
    "mh_or",
    "gene_cluster_pvalue",
    "cmh_pvalue",
    "estimable",
    "boundary_estimate",
    "partial_input",
}


def _truthy(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )


def _check_fdr(
    frame: pd.DataFrame,
    group_columns: list[str],
    p_column: str,
    fdr_column: str,
) -> list[str]:
    errors = []
    for key, group in frame.groupby(group_columns, sort=False, dropna=False):
        expected = bh_adjust_declared(
            pd.to_numeric(group[p_column], errors="coerce")
        )
        observed = pd.to_numeric(group[fdr_column], errors="coerce").to_numpy()
        if not np.allclose(expected, observed, rtol=1e-10, atol=1e-12, equal_nan=True):
            errors.append(f"{fdr_column} is inconsistent for group {key}")
    return errors


def _validate_results(frame: pd.DataFrame, regional: bool) -> list[str]:
    errors = []
    required = MAIN_COLUMNS | (
        {"fdr_within_database_region", "cmh_fdr_within_database_region"}
        if regional
        else {"fdr_within_database", "cmh_fdr_within_database"}
    )
    missing = required - set(frame.columns)
    if missing:
        return [f"result table is missing columns: {sorted(missing)}"]
    duplicate_keys = ["database", "RBP", "region_filter"] if regional else [
        "database",
        "RBP",
    ]
    if frame.duplicated(duplicate_keys).any():
        errors.append(f"duplicate result keys: {duplicate_keys}")
    if regional:
        invalid_region = ~frame["region_filter"].isin(REGIONS)
    else:
        invalid_region = ~frame["region_filter"].eq("all")
    if invalid_region.any():
        errors.append("unexpected region_filter values")

    numeric = [
        "n_case_sites",
        "n_control_sites",
        "n_strata_total",
        "case_overlaps",
        "control_overlaps",
        "mh_numerator",
        "mh_denominator",
        "mh_or",
    ]
    converted = {column: pd.to_numeric(frame[column], errors="coerce") for column in numeric}
    if (
        (converted["n_case_sites"] <= 0).any()
        or (converted["n_control_sites"] <= 0).any()
        or (converted["n_strata_total"] <= 0).any()
    ):
        errors.append("non-positive analysis counts")
    if (
        (converted["case_overlaps"] < 0).any()
        or (converted["control_overlaps"] < 0).any()
        or (converted["case_overlaps"] > converted["n_case_sites"]).any()
        or (converted["control_overlaps"] > converted["n_control_sites"]).any()
    ):
        errors.append("overlap counts are outside their totals")

    estimable = _truthy(frame["estimable"])
    boundary = _truthy(frame["boundary_estimate"])
    if estimable.isna().any() or boundary.isna().any():
        errors.append("invalid boolean flags")
    if (estimable.fillna(False) & boundary.fillna(False)).any():
        errors.append("an estimate cannot be both estimable and boundary")
    expected_estimable = (
        converted["mh_numerator"].gt(0) & converted["mh_denominator"].gt(0)
    )
    expected_boundary = (
        converted["mh_numerator"].gt(0)
        ^ converted["mh_denominator"].gt(0)
    )
    if not expected_estimable.equals(estimable.astype(bool)):
        errors.append("estimable flag is inconsistent with MH components")
    if not expected_boundary.equals(boundary.astype(bool)):
        errors.append("boundary flag is inconsistent with MH components")
    finite = expected_estimable
    expected_or = (
        converted["mh_numerator"][finite] / converted["mh_denominator"][finite]
    )
    if not np.allclose(
        expected_or,
        converted["mh_or"][finite],
        rtol=1e-10,
        atol=1e-12,
    ):
        errors.append("MH odds ratios do not match their components")

    if regional:
        errors.extend(
            _check_fdr(
                frame,
                ["database", "region_filter"],
                "gene_cluster_pvalue",
                "fdr_within_database_region",
            )
        )
        errors.extend(
            _check_fdr(
                frame,
                ["database", "region_filter"],
                "cmh_pvalue",
                "cmh_fdr_within_database_region",
            )
        )
    else:
        errors.extend(
            _check_fdr(
                frame,
                ["database"],
                "gene_cluster_pvalue",
                "fdr_within_database",
            )
        )
        errors.extend(
            _check_fdr(
                frame,
                ["database"],
                "cmh_pvalue",
                "cmh_fdr_within_database",
            )
        )
    return errors


def _validate_matrix(
    run_dir: Path,
    database: str,
    results: pd.DataFrame,
    n_cases: int,
    case_site_ids: np.ndarray,
) -> list[str]:
    errors = []
    key = DATABASE_KEYS.get(database)
    if key is None:
        return [f"unknown database label {database}"]
    matrix_path = run_dir / f"case_binding_matrix_{key}.npz"
    rows_path = run_dir / f"case_binding_matrix_{key}_rows.tsv"
    if not matrix_path.is_file() or not rows_path.is_file():
        return [f"missing case-binding matrix for {database}"]
    rows = pd.read_csv(rows_path, sep="\t")
    subset = results[results["database"].eq(database)].sort_values(
        "RBP", kind="mergesort"
    )
    try:
        with np.load(matrix_path) as archive:
            packed = archive["packed_matrix"]
            rbps = archive["rbps"].astype(str)
            stored_ids = archive["case_site_ids"]
            stored_n_cases = int(archive["n_case_sites"][0])
            bitorder = str(archive["bitorder"][0])
    except Exception as exc:
        return [f"cannot read {matrix_path.name}: {exc}"]
    if bitorder != "little":
        errors.append(f"{database} matrix uses unexpected bit order")
    if stored_n_cases != n_cases or not np.array_equal(stored_ids, case_site_ids):
        errors.append(f"{database} matrix case-site identity mismatch")
    if packed.shape != (len(rbps), math.ceil(n_cases / 8)):
        errors.append(f"{database} matrix has wrong shape")
        return errors
    if list(rbps) != rows["RBP"].astype(str).tolist():
        errors.append(f"{database} matrix row labels disagree")
    if set(rbps) != set(subset["RBP"].astype(str)):
        errors.append(f"{database} matrix RBP universe disagrees with results")
    unpacked = np.unpackbits(packed, axis=1, bitorder="little")[:, :n_cases]
    overlap_by_rbp = dict(zip(rbps, unpacked.sum(axis=1)))
    observed = subset.set_index("RBP")["case_overlaps"]
    for rbp, count in observed.items():
        if int(count) != int(overlap_by_rbp[str(rbp)]):
            errors.append(f"{database}/{rbp} matrix overlap count mismatch")
            break
    return errors


def _validate_manifest_outputs(
    manifest: dict,
    run_dir: Path,
    project_root: Path,
    intended_destination: Path | None,
) -> list[str]:
    errors = []
    destination = (
        portable_path(intended_destination, project_root)
        if intended_destination is not None
        else portable_path(run_dir, project_root)
    )
    prefix = "" if destination == "." else destination.rstrip("/") + "/"
    for record in manifest.get("outputs", []):
        recorded = str(record.get("path", ""))
        if prefix and not recorded.startswith(prefix):
            errors.append(f"manifest output is outside run directory: {recorded}")
            continue
        relative = recorded[len(prefix) :] if prefix else recorded
        candidate = run_dir / relative
        if not candidate.is_file():
            errors.append(f"manifest output is missing: {relative}")
            continue
        if candidate.stat().st_size != int(record.get("bytes", -1)):
            errors.append(f"manifest output size mismatch: {relative}")
        elif sha256_file(candidate) != record.get("sha256"):
            errors.append(f"manifest output checksum mismatch: {relative}")
    return errors


def validate_run(
    run_dir: Path,
    project_root: Path,
    check_input_hashes: bool = True,
    intended_destination: Path | None = None,
) -> list[str]:
    run_dir = run_dir.expanduser().resolve()
    project_root = project_root.expanduser().resolve()
    required = [
        "README.md",
        "run_manifest.json",
        "stratified_enrichment_results.tsv",
        "region_specific_results.tsv",
        "resource_loading_report.tsv",
        "rbp_list_provenance.tsv",
        "opportunity_design.npz",
        "case_site_assignments.tsv.gz",
        "transcript_region_strata.tsv",
        "opportunity_flow.tsv",
        "coverage_balance.tsv",
    ]
    missing = [name for name in required if not (run_dir / name).is_file()]
    if missing:
        return [f"missing required files: {missing}"]
    errors: list[str] = []
    try:
        results = pd.read_csv(
            run_dir / "stratified_enrichment_results.tsv", sep="\t"
        )
        regional = pd.read_csv(run_dir / "region_specific_results.tsv", sep="\t")
        strata = pd.read_csv(run_dir / "transcript_region_strata.tsv", sep="\t")
        resources = pd.read_csv(
            run_dir / "resource_loading_report.tsv", sep="\t"
        )
    except Exception as exc:
        return [f"cannot read output tables: {exc}"]
    errors.extend(_validate_results(results, regional=False))
    errors.extend(_validate_results(regional, regional=True))
    if set(regional["database"]) != set(results["database"]):
        errors.append("regional and overall database sets differ")
    expected_regions = set(strata["region"].astype(str))
    for database, rbps in results.groupby("database"):
        present = regional[regional["database"].eq(database)]
        expected = {
            (str(rbp), region)
            for rbp in rbps["RBP"]
            for region in expected_regions
        }
        observed = set(zip(present["RBP"].astype(str), present["region_filter"]))
        if observed != expected:
            errors.append(f"{database} regional RBP/region grid is incomplete")
    if not {"database", "RBP", "status", "n_source_intervals"} <= set(
        resources.columns
    ):
        errors.append("resource report schema is incomplete")

    try:
        with np.load(run_dir / "opportunity_design.npz") as archive:
            positions = archive["positions"]
            is_case = archive["is_case"].astype(bool)
            stratum_codes = archive["stratum_codes"]
            case_site_ids = archive["case_site_ids"]
    except Exception as exc:
        return errors + [f"cannot read opportunity design: {exc}"]
    if not (
        len(positions) == len(is_case) == len(stratum_codes)
        and len(case_site_ids) == int(is_case.sum())
    ):
        errors.append("opportunity-design array lengths disagree")
    if len(stratum_codes) and (
        stratum_codes.min() < 0 or stratum_codes.max() >= len(strata)
    ):
        errors.append("opportunity-design stratum code is out of range")
    if len(results):
        expected_cases = int(is_case.sum())
        if not results["n_case_sites"].eq(expected_cases).all():
            errors.append("overall case totals disagree with opportunity design")
        for database in sorted(set(results["database"])):
            errors.extend(
                _validate_matrix(
                    run_dir,
                    database,
                    results,
                    expected_cases,
                    case_site_ids,
                )
            )

    try:
        manifest = json.loads((run_dir / "run_manifest.json").read_text())
    except Exception as exc:
        return errors + [f"cannot read run manifest: {exc}"]
    if manifest.get("workflow") != "transcript-region":
        errors.append("manifest workflow is not transcript-region")
    if manifest.get("parameters", {}).get("design") != "transcript-region":
        errors.append("manifest design does not identify transcript-region")
    if check_input_hashes:
        for record in manifest.get("inputs", []):
            errors.extend(validate_file_record(record, project_root))
    errors.extend(
        _validate_manifest_outputs(
            manifest, run_dir, project_root, intended_destination
        )
    )

    for manifest_name in ("plot_manifest.tsv",):
        candidate = run_dir / manifest_name
        if candidate.is_file():
            plot_manifest = pd.read_csv(candidate, sep="\t")
            if "path" not in plot_manifest:
                errors.append(f"{manifest_name} has no path column")
            else:
                for path in plot_manifest["path"].dropna():
                    if not (run_dir / str(path)).is_file():
                        errors.append(f"missing plotted file: {path}")
    return errors
