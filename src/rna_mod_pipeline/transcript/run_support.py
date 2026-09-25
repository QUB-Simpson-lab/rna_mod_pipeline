from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from rna_mod_pipeline.evidence_quality import add_evidence_quality

from .config import DATABASE_LABELS, EXPECTED_DATABASE_COUNTS
from .fdr import add_fdr_columns, add_region_fdr_columns
from .outputs import write_run_manifest, write_run_readme
from .processing import _database_filename
from .resources import load_rbp_catalog

if TYPE_CHECKING:
    from .opportunities import OpportunityBuild, OpportunityDesign
    from .workflow import RunOptions


def prepare_catalog(
    options: RunOptions, output: Path
) -> tuple[dict[str, list[str]], pd.DataFrame, bool]:
    lists, catalog = load_rbp_catalog(
        options.rbp_catalog,
        options.databases,
        set(options.rbps) if options.rbps else None,
        options.maximum_rbps_per_database,
    )
    selected_pairs = {
        (database, rbp)
        for database, rbps in lists.items()
        for rbp in rbps
    }
    catalog["selected_for_run"] = [
        (database, rbp) in selected_pairs
        for database, rbp in zip(catalog["database_key"], catalog["RBP"])
    ]
    catalog.to_csv(output / "rbp_list_provenance.tsv", sep="\t", index=False)
    complete = (
        options.rbps is None and options.maximum_rbps_per_database is None
    )
    if complete:
        wrong = {
            database: len(lists[database])
            for database in options.databases
            if len(lists[database]) != EXPECTED_DATABASE_COUNTS[database]
        }
        if wrong:
            raise ValueError(
                "Canonical RBP catalogue count mismatch for a complete "
                f"run: {wrong}"
            )
    return lists, catalog, complete


def validate_resource_table(
    resources: pd.DataFrame,
    complete_rbp_universe: bool,
    partial_postar: bool,
    resource_filtered: bool,
) -> None:
    if not complete_rbp_universe:
        return
    checked = resources[~resources["RBP"].eq("__STREAM_SUMMARY__")].copy()
    if partial_postar:
        checked = checked[~checked["database"].eq("POSTAR3")]
    unavailable = checked[
        ~checked["source_identity_present"].astype(bool)
    ]
    if not unavailable.empty:
        summary = unavailable.groupby("database").size().to_dict()
        raise RuntimeError(
            f"A complete run has unavailable RBP resources: {summary}"
        )
    if not resource_filtered:
        unusable = checked[~checked["status"].eq("loaded")]
        if not unusable.empty:
            summary = unusable.groupby("database").size().to_dict()
            raise RuntimeError(
                "Default resource loading produced empty or invalid "
                f"RBP tracks: {summary}"
            )


def write_analysis_tables(
    output: Path,
    options: RunOptions,
    result_rows: list[dict],
    region_rows: list[dict],
    influence_rows: list[dict],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    results = add_evidence_quality(
        add_fdr_columns(pd.DataFrame(result_rows)),
        odds_column="mh_or",
        fdr_column="fdr_within_database",
    ).sort_values(
        ["database", "RBP"], kind="mergesort"
    )
    regional = add_evidence_quality(
        add_region_fdr_columns(pd.DataFrame(region_rows)),
        odds_column="mh_or",
        fdr_column="fdr_within_database_region",
    ).sort_values(
        ["database", "region_filter", "RBP"], kind="mergesort"
    )
    results.to_csv(
        output / "stratified_enrichment_results.tsv", sep="\t", index=False
    )
    regional.to_csv(
        output / "region_specific_results.tsv", sep="\t", index=False
    )
    pd.DataFrame(influence_rows).to_csv(
        output / "top_gene_influence.tsv", sep="\t", index=False
    )
    for database_key in options.databases:
        database = DATABASE_LABELS[database_key]
        folder = output / database_key
        folder.mkdir()
        results[results["database"].eq(database)].to_csv(
            folder / _database_filename(database_key), sep="\t", index=False
        )
    return results, regional


def generate_plots(
    output: Path,
    options: RunOptions,
    results: pd.DataFrame,
    regional: pd.DataFrame,
) -> None:
    from rna_mod_pipeline.plotting.enrichment import generate_all_plots
    from rna_mod_pipeline.plotting.regional import generate_regional_plots

    generate_all_plots(
        results,
        output,
        options.analysis_id,
        options.fdr_threshold,
        options.top_n,
        options.dpi,
    )
    generate_regional_plots(
        regional,
        output / "regional_top_plots",
        options.analysis_id,
        options.fdr_threshold,
        min(options.top_n, 30),
        options.dpi,
    )


def write_metadata(
    output: Path,
    options: RunOptions,
    build: OpportunityBuild,
    design: OpportunityDesign,
    resource_inputs: list[Path],
    partial_postar: bool,
    actual_partial: bool,
    database_subset: bool,
    custom_parameters: bool,
    resource_filtered: bool,
) -> None:
    write_run_readme(
        output / "README.md",
        options.analysis_id,
        options.modification,
        [DATABASE_LABELS[key] for key in options.databases],
        design.n_cases,
        design.n_controls,
        actual_partial,
        database_subset,
        custom_parameters,
        resource_filtered,
    )
    parameters = {
        "analysis_id": options.analysis_id,
        "modification": options.modification,
        "databases": list(options.databases),
        "coverage_min": options.coverage_min,
        "case_min_fraction": options.case_min_fraction,
        "background_max_fraction": options.background_max_fraction,
        "window": options.window,
        "strand_mode": options.strand_mode,
        "minimum_encori_support": options.minimum_encori_support,
        "fdr_threshold": options.fdr_threshold,
        "top_n": options.top_n,
        "seed": options.seed,
        "chromosomes": sorted(options.chromosomes or ()),
        "requested_rbps": sorted(options.rbps or ()),
        "maximum_rbps_per_database": options.maximum_rbps_per_database,
        "max_cases": options.max_cases,
        "max_bed_rows": options.max_bed_rows,
        "max_postar_rows": options.max_postar_rows,
        "postar_cell_types": sorted(options.postar_cell_types or ()),
        "postar_methods": sorted(options.postar_methods or ()),
        "partial_input": actual_partial,
        "database_subset": database_subset,
        "custom_parameters": custom_parameters,
        "resource_filtered": resource_filtered,
        "partial_bedmethyl_scan": build.partial_bed_scan,
        "partial_postar3_scan": partial_postar,
        "design": "transcript-region",
        "all_rbps_overview_policy": (
            "finite positive non-boundary estimates on a descending log2 "
            "Mantel-Haenszel odds-ratio scale; OR=0 boundaries represented "
            "by individual negative-infinity bars; other boundary and "
            "non-estimable tests counted on the figure and retained in the "
            "source TSV"
        ),
        "ornament_coordinate_policy": (
            "historical transcript interpretation; end<=start rows excluded"
            if "ornament" in options.databases
            else "not_applicable"
        ),
    }
    inputs = [
        options.sites,
        options.bedmethyl,
        options.gtf,
        options.rbp_catalog,
        *resource_inputs,
    ]
    if options.postar3 is not None and "postar3" in options.databases:
        inputs.append(options.postar3)
    write_run_manifest(
        output / "run_manifest.json",
        output,
        options.output_dir,
        options.project_root,
        inputs,
        parameters,
    )
