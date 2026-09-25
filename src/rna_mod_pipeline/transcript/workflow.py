from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from rna_mod_pipeline.cli import (
    StagedOutputDirectory,
    project_protected_trees,
    validate_output_location,
)

from .config import (
    DATABASE_LABELS,
)
from .opportunities import build_opportunity_design
from .outputs import (
    save_opportunity_design,
)
from .processing import (
    _database_filename,
    _process_interval_database,
    _process_postar3,
    process_databases,
)
from .run_support import (
    generate_plots,
    prepare_catalog,
    validate_resource_table,
    write_analysis_tables,
    write_metadata,
)


@dataclass(frozen=True)
class RunOptions:
    project_root: Path
    analysis_id: str
    modification: str
    sites: Path
    bedmethyl: Path
    gtf: Path
    rbp_catalog: Path
    output_dir: Path
    databases: tuple[str, ...] = ("ornament", "encori", "postar3")
    ornament_dir: Path | None = None
    encori_dir: Path | None = None
    postar3: Path | None = None
    coverage_min: int = 20
    case_min_fraction: float = 20.0
    background_max_fraction: float = 20.0
    window: int = 10
    strand_mode: str = "ignore"
    minimum_encori_support: int = 1
    fdr_threshold: float = 0.05
    top_n: int = 30
    dpi: int = 180
    seed: int = 20260728
    rbps: frozenset[str] | None = None
    chromosomes: frozenset[str] | None = None
    postar_cell_types: frozenset[str] | None = None
    postar_methods: frozenset[str] | None = None
    maximum_rbps_per_database: int | None = None
    max_cases: int | None = None
    max_bed_rows: int | None = None
    max_postar_rows: int | None = None
    make_plots: bool = True
    overwrite: bool = False


def _validate_options(options: RunOptions) -> None:
    allowed = set(DATABASE_LABELS)
    unknown = set(options.databases) - allowed
    if unknown or not options.databases:
        raise ValueError(f"Invalid database selection: {sorted(unknown)}")
    if len(set(options.databases)) != len(options.databases):
        raise ValueError("Each database may be selected only once")
    if options.strand_mode not in {"ignore", "same"}:
        raise ValueError("strand_mode must be ignore or same")
    if "ornament" in options.databases and options.strand_mode == "same":
        raise ValueError("oRNAment cannot be used with same-strand overlap")
    if options.coverage_min < 1 or options.window < 0:
        raise ValueError("coverage_min must be positive and window non-negative")
    if options.minimum_encori_support < 1:
        raise ValueError("minimum_encori_support must be positive")
    if (
        "encori" not in options.databases
        and options.minimum_encori_support != 1
    ):
        raise ValueError(
            "minimum_encori_support applies only when ENCORI is selected"
        )
    if "postar3" not in options.databases and (
        options.postar_cell_types
        or options.postar_methods
        or options.max_postar_rows
    ):
        raise ValueError(
            "POSTAR3 filters and max_postar_rows require POSTAR3 selection"
        )
    if options.top_n < 1 or options.dpi < 1:
        raise ValueError("top_n and dpi must be positive")
    for name, value in (
        ("maximum_rbps_per_database", options.maximum_rbps_per_database),
        ("max_cases", options.max_cases),
        ("max_bed_rows", options.max_bed_rows),
        ("max_postar_rows", options.max_postar_rows),
    ):
        if value is not None and value < 1:
            raise ValueError(f"{name} must be positive")
    if not 0 <= options.background_max_fraction <= 100:
        raise ValueError("background_max_fraction must be between 0 and 100")
    if not 0 <= options.case_min_fraction <= 100:
        raise ValueError("case_min_fraction must be between 0 and 100")
    if options.background_max_fraction > options.case_min_fraction:
        raise ValueError(
            "background_max_fraction cannot exceed case_min_fraction"
        )
    if not 0 < options.fdr_threshold <= 1:
        raise ValueError("fdr_threshold must be in (0, 1]")
    required_files = {
        "sites": options.sites,
        "bedMethyl": options.bedmethyl,
        "GTF": options.gtf,
        "RBP catalogue": options.rbp_catalog,
    }
    if "postar3" in options.databases:
        required_files["POSTAR3"] = options.postar3
    for label, path in required_files.items():
        if path is None or not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
    required_directories = {}
    if "ornament" in options.databases:
        required_directories["oRNAment"] = options.ornament_dir
    if "encori" in options.databases:
        required_directories["ENCORI"] = options.encori_dir
    for label, path in required_directories.items():
        if path is None or not path.is_dir():
            raise FileNotFoundError(f"{label} directory not found: {path}")
    protected_files = [
        path.expanduser().resolve()
        for path in required_files.values()
        if path is not None
    ]
    protected_directories = [
        path.expanduser().resolve()
        for path in required_directories.values()
        if path is not None
    ]
    output = validate_output_location(
        options.output_dir,
        [*protected_files, *protected_directories],
        protected_roots=(options.project_root,),
        protected_trees=project_protected_trees(options.project_root),
    )
    if any(
        output == path or output.is_relative_to(path)
        for path in protected_directories
    ):
        raise ValueError(
            "Output directory must be separate from binding-resource directories"
        )


def _has_truncated_input(options: RunOptions) -> bool:
    return bool(
        options.rbps
        or options.maximum_rbps_per_database
        or options.max_cases
        or options.max_bed_rows
        or options.max_postar_rows
        or options.chromosomes
    )


def _uses_custom_parameters(options: RunOptions) -> bool:
    return bool(
        options.coverage_min != 20
        or not math.isclose(options.case_min_fraction, 20.0)
        or not math.isclose(options.background_max_fraction, 20.0)
        or options.window != 10
        or options.strand_mode != "ignore"
        or (
            "encori" in options.databases
            and options.minimum_encori_support != 1
        )
        or (
            "postar3" in options.databases
            and (options.postar_cell_types or options.postar_methods)
        )
    )


def _uses_resource_filters(options: RunOptions) -> bool:
    return bool(
        (
            "encori" in options.databases
            and options.minimum_encori_support != 1
        )
        or (
            "postar3" in options.databases
            and (options.postar_cell_types or options.postar_methods)
        )
    )


def run_analysis(options: RunOptions) -> Path:
    _validate_options(options)
    truncated_input = _has_truncated_input(options)
    custom_parameters = _uses_custom_parameters(options)
    database_subset = set(options.databases) != set(DATABASE_LABELS)
    publisher = StagedOutputDirectory(
        options.output_dir,
        options.overwrite,
        protected_roots=(options.project_root,),
        protected_trees=project_protected_trees(options.project_root),
        expected_workflow="transcript-region",
    )
    with publisher as output:
        build = build_opportunity_design(
            options.sites,
            options.bedmethyl,
            options.gtf,
            options.modification,
            options.coverage_min,
            options.case_min_fraction,
            options.background_max_fraction,
            set(options.chromosomes) if options.chromosomes else None,
            options.max_cases,
            options.max_bed_rows,
            options.seed,
        )
        design = build.design
        save_opportunity_design(build, output)
        lists, catalog, complete_rbp_universe = prepare_catalog(options, output)
        processed = process_databases(
            lists,
            catalog,
            design,
            options,
            output,
            truncated_input or build.partial_bed_scan,
        )
        resources = pd.DataFrame(processed.resource_rows)
        resources.to_csv(
            output / "resource_loading_report.tsv", sep="\t", index=False
        )
        resource_filtered = _uses_resource_filters(options)
        validate_resource_table(
            resources,
            complete_rbp_universe,
            processed.partial_postar,
            resource_filtered,
        )
        results, regional = write_analysis_tables(
            output,
            options,
            processed.result_rows,
            processed.region_rows,
            processed.influence_rows,
        )
        if options.make_plots:
            generate_plots(output, options, results, regional)

        actual_partial = (
            truncated_input
            or build.partial_bed_scan
            or processed.partial_postar
        )
        write_metadata(
            output,
            options,
            build,
            design,
            processed.resource_inputs,
            processed.partial_postar,
            actual_partial,
            database_subset,
            custom_parameters,
            resource_filtered,
        )

        from .validation import validate_run

        errors = validate_run(
            output,
            options.project_root,
            check_input_hashes=False,
            intended_destination=options.output_dir,
        )
        if errors:
            raise RuntimeError(
                "Staged transcript-region output failed validation:\n- "
                + "\n- ".join(errors)
            )
        publisher.publish()
    return options.output_dir
