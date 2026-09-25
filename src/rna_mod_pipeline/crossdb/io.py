from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import RBP_ALIASES, canonical_rbp
from ..evidence_quality import add_evidence_quality
from ..io import require_file
from ..provenance import portable_path
from ..schemas import validate_result_identity

DATABASES = ("ornament", "encori", "postar3")


def _canonical_case_insensitive(value: object) -> str:
    cleaned = str(value).strip()
    direct = canonical_rbp(cleaned)
    if direct != cleaned:
        return direct
    aliases = {key.casefold(): target for key, target in RBP_ALIASES.items()}
    return aliases.get(cleaned.casefold(), cleaned.upper())


def load_database_results(
    path: str | Path,
    database: str,
    *,
    context: str = "all",
    project_root: str | Path | None = None,
    expected_modification: str | None = None,
    expected_design: str | None = None,
) -> pd.DataFrame:
    database = database.lower()
    if database not in DATABASES:
        raise ValueError(f"Unknown database: {database}")
    frame = pd.read_csv(require_file(path, f"{database} enrichment table"), sep="\t")
    if "database" in frame.columns:
        frame = frame[
            frame["database"].astype(str).str.lower() == database
        ].copy()
        if frame.empty:
            raise ValueError(f"No {database} rows were found in {path}")
    identity = validate_result_identity(
        frame,
        expected_modification=expected_modification,
        expected_database=database,
        expected_design=expected_design,
        context=context,
        table_name=f"{database} enrichment table",
    )
    if "region_filter" in frame.columns:
        frame = frame[frame["region_filter"].astype(str).str.lower() == context.lower()]
    if frame.empty:
        raise ValueError(
            f"No {database} rows were found for context {context!r} in {path}"
        )

    or_column = "mh_or" if "mh_or" in frame.columns else "odds_ratio"
    fdr_candidates = (
        "fdr_within_database_region",
        "fdr_within_database",
        "fdr_within_analysis",
        "fdr",
    )
    p_candidates = ("gene_cluster_pvalue", "p_value", "pvalue", "cmh_pvalue")
    fdr_column = next((name for name in fdr_candidates if name in frame), None)
    p_column = next((name for name in p_candidates if name in frame), None)
    required = {"RBP", or_column}
    if fdr_column is None:
        required.add("fdr")
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")

    frame = frame.copy()
    frame["RBP"] = frame["RBP"].map(_canonical_case_insensitive)
    if frame["RBP"].duplicated().any():
        names = sorted(frame.loc[frame["RBP"].duplicated(False), "RBP"].unique())
        raise ValueError(
            f"{database} contains duplicate canonical RBP names: {names[:10]}"
        )

    standard = pd.DataFrame(
        {
            "RBP": frame["RBP"],
            "odds_ratio": pd.to_numeric(frame[or_column], errors="coerce"),
            "fdr": pd.to_numeric(frame[fdr_column], errors="coerce"),
            "p_value": (
                pd.to_numeric(frame[p_column], errors="coerce")
                if p_column
                else np.nan
            ),
        }
    )
    retained = [
        name
        for name in (
            "total_binding_sites",
            "n_source_intervals",
            "mod_overlaps",
            "case_overlaps",
            "n_case_sites",
            "bg_overlaps",
            "control_overlaps",
            "n_control_sites",
            "n_source_tracks",
            "n_experiments",
            "n_cell_types",
            "n_informative_strata",
            "n_informative_genes",
            "n_gene_clusters_for_robust_variance",
            "low_gene_cluster_count",
            "direction_stable_after_each_gene_removed",
            "leave_one_gene_out_direction_reversals",
            "leave_one_gene_out_boundaries",
            "maximum_absolute_gene_influence",
            "estimable",
            "boundary_estimate",
            "partial_input",
        )
        if name in frame.columns
    ]
    for name in retained:
        standard[name] = frame[name].to_numpy()
    standard = add_evidence_quality(
        standard,
        odds_column="odds_ratio",
        fdr_column="fdr",
    )
    standard["database"] = database
    standard["source_file"] = (
        portable_path(path, project_root)
        if project_root is not None
        else Path(path).as_posix()
    )
    standard["source_or_column"] = or_column
    standard["source_fdr_column"] = fdr_column
    standard["source_p_column"] = p_column or "not_available"
    for label, value in identity.items():
        standard[f"input_{label}_status"] = value
    return standard.reset_index(drop=True)


def default_input_paths(
    project_root: str | Path,
    analysis: str,
    design: str,
    context: str = "all",
) -> dict[str, Path]:
    root = Path(project_root).resolve()
    if design == "loose":
        base = root / "refactored_outputs" / analysis / "loose"
        return {
            "ornament": base / "ornament" / "rbp_enrichment_results.tsv",
            "encori": base / "encori" / "encori_enrichment_results.tsv",
            "postar3": base / "postar3" / "postar3_enrichment_results.tsv",
        }
    if design in {"transcript", "transcript-region"}:
        run = root / "refactored_outputs" / analysis / "transcript_region"
        table = run / (
            "stratified_enrichment_results.tsv"
            if context.lower() == "all"
            else "region_specific_results.tsv"
        )
        return {database: table for database in DATABASES}
    raise ValueError("design must be 'loose' or 'transcript-region'")
