#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.cli import (
    StagedOutputFiles,
    csv_values,
    percentage,
    positive_int,
    project_protected_trees,
    resolve_path,
)
from rna_mod_pipeline.config import MODIFICATIONS, modification, project_path
from rna_mod_pipeline.binding.catalog import load_catalog
from rna_mod_pipeline.binding.loose import run_loose_overlap
from rna_mod_pipeline.binding.manifest import loose_manifest_parameters
from rna_mod_pipeline.binding.plots import plot_loose_results
from rna_mod_pipeline.binding.paths import (
    default_phase1_sites,
    loose_defaults,
    loose_output_paths,
    validate_loose_output_ownership,
)
from rna_mod_pipeline.binding.sources import (
    binding_source_files,
    iter_binding_source,
)
from rna_mod_pipeline.io import write_tsv
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_files,
    write_manifest,
)


DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CATALOG = Path(__file__).resolve().parents[1] / "config" / "rbp_catalog.tsv"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Run one loose, all-context RBP overlap analysis. Cases and controls "
            "are not matched by transcript, region, motif, or coverage bin."
        )
    )
    result.add_argument(
        "--modification",
        required=True,
        choices=MODIFICATIONS,
        help="Run exactly one modification or m6A replicate.",
    )
    result.add_argument(
        "--database",
        required=True,
        choices=["ornament", "encori", "postar3"],
        help="Run exactly one binding database.",
    )
    result.add_argument("--project-root", default=DEFAULT_PROJECT_ROOT)
    result.add_argument("--sites", help="Phase 1 metagene site table.")
    result.add_argument("--bedmethyl", help="Raw bedMethyl used for controls.")
    result.add_argument(
        "--binding-source",
        help="oRNAment/ENCORI directory or POSTAR3 human.txt.",
    )
    result.add_argument("--catalog", default=DEFAULT_CATALOG)
    result.add_argument("--output-dir", help="Database-specific result directory.")
    result.add_argument("--plot-dir", help="Database-specific plot directory.")
    result.add_argument("--window", type=int, default=10)
    result.add_argument("--background-min-coverage", type=positive_int, default=20)
    result.add_argument("--background-max-fraction", type=percentage, default=20.0)
    result.add_argument("--case-min-coverage", type=positive_int, default=20)
    result.add_argument("--case-min-fraction", type=percentage, default=20.0)
    result.add_argument("--fdr-threshold", type=float, default=0.05)
    result.add_argument("--minimum-clip-experiments", type=positive_int, default=1)
    result.add_argument(
        "--cell-types",
        type=csv_values,
        help="POSTAR3 exact cell-type names, comma-separated.",
    )
    result.add_argument(
        "--methods",
        type=csv_values,
        help="POSTAR3 CLIP method names, comma-separated.",
    )
    result.add_argument(
        "--site-annotation",
        choices=["all", "significant", "none"],
        default="all",
        help=(
            "all: binary column per tested RBP plus significant list; "
            "significant: list only; none: no per-site binding annotation."
        ),
    )
    result.add_argument("--top-n", type=positive_int, default=30)
    result.add_argument(
        "--allow-site-subset",
        action="store_true",
        help=(
            "Allow an intentional subset of otherwise qualifying calls, such "
            "as a DRACH-only site table. All supplied sites must still be "
            "valid calls from the selected raw bedMethyl."
        ),
    )
    result.add_argument(
        "--allow-empty-resources",
        action="store_true",
        help=(
            "Retain declared resources with no usable intervals as "
            "non-estimable rows. Intended for exploratory subsets."
        ),
    )
    result.add_argument("--skip-plots", action="store_true")
    result.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace only this run's known output files.",
    )
    return result


def main() -> None:
    args = parser().parse_args()
    if args.window < 0:
        raise ValueError("--window must be non-negative")
    if not 0 < args.fdr_threshold <= 1:
        raise ValueError("--fdr-threshold must lie in (0, 1]")
    if args.database != "postar3" and (args.cell_types or args.methods):
        raise ValueError("--cell-types and --methods apply only to POSTAR3")
    if args.database != "encori" and args.minimum_clip_experiments != 1:
        raise ValueError(
            "--minimum-clip-experiments applies only to ENCORI"
        )

    config = modification(args.modification)
    root = Path(args.project_root).expanduser().resolve()
    default_source, default_output, default_plots = loose_defaults(
        root, args.modification, args.database
    )
    sites_path = (
        resolve_path(root, args.sites)
        if args.sites
        else default_phase1_sites(root, args.modification, config.short_name)
    )
    bedmethyl_path = (
        resolve_path(root, args.bedmethyl)
        if args.bedmethyl
        else project_path(root, config.bedmethyl)
    )
    source = (
        resolve_path(root, args.binding_source)
        if args.binding_source
        else default_source
    )
    catalog_path = resolve_path(root, args.catalog)
    output_dir = (
        resolve_path(root, args.output_dir)
        if args.output_dir
        else default_output
    )
    plot_dir = (
        resolve_path(root, args.plot_dir)
        if args.plot_dir
        else default_plots
    )
    enrichment_path, annotation_path, manifest_path = loose_output_paths(
        output_dir, config.short_name, args.database
    )
    validate_loose_output_ownership(
        output_dir,
        args.modification,
        args.database,
    )
    plot_paths = [
        plot_dir / f"{args.database}_top_enriched.png",
        plot_dir / f"{args.database}_top_depleted.png",
        plot_dir / f"{args.database}_all_rbps_overview.png",
        plot_dir / f"{args.database}_machinery_heatmap.png",
    ]
    entries = load_catalog(catalog_path, args.database)
    source_files = binding_source_files(args.database, source, entries)
    inputs = [sites_path, bedmethyl_path, catalog_path, *source_files]
    output_paths = [
        enrichment_path,
        annotation_path,
        *([] if args.skip_plots else plot_paths),
    ]
    destinations = [*output_paths, manifest_path]
    if plot_dir in {enrichment_path, annotation_path, manifest_path} or len(
        set(destinations)
    ) != len(destinations):
        raise ValueError("Enrichment, annotation, manifest, and plot paths must differ")
    publisher = StagedOutputFiles(
        destinations,
        overwrite=args.overwrite,
        expected_workflow="loose_overlap",
        ownership_manifest=manifest_path,
        project_root=root,
        expected_parameters={
            "modification": args.modification,
            "database": args.database,
        },
        inputs=inputs,
        protected_roots=[root],
        protected_trees=project_protected_trees(root),
    )
    with publisher as staged:
        tracks = iter_binding_source(
            args.database,
            source,
            entries,
            minimum_clip_experiments=args.minimum_clip_experiments,
            cell_types=set(args.cell_types) if args.cell_types else None,
            methods=set(args.methods) if args.methods else None,
            allow_empty_resources=args.allow_empty_resources,
        )
        result = run_loose_overlap(
            sites_path,
            bedmethyl_path,
            tracks,
            modification=args.modification,
            database=args.database,
            window=args.window,
            background_min_coverage=args.background_min_coverage,
            background_max_fraction=args.background_max_fraction,
            fdr_threshold=args.fdr_threshold,
            site_annotation=args.site_annotation,
            case_min_coverage=args.case_min_coverage,
            case_min_fraction=args.case_min_fraction,
            allow_site_subset=args.allow_site_subset,
        )
        write_tsv(result.enrichment, staged.path(enrichment_path))
        write_tsv(result.annotated_sites, staged.path(annotation_path))
        if not args.skip_plots:
            plot_loose_results(
                result.enrichment,
                config.display,
                args.database,
                staged.directory(plot_dir),
                args.top_n,
                args.fdr_threshold,
            )
        staged_outputs = [staged.path(path) for path in output_paths]
        manifest = make_manifest(
            "loose_overlap",
            root,
            inputs,
            loose_manifest_parameters(
                result,
                modification=args.modification,
                database=args.database,
                catalog_version=entries[0].version,
                options=vars(args),
            ),
            staged_outputs,
        )
        remap_manifest_output_files(
            manifest, staged_outputs, output_paths, root
        )
        write_manifest(manifest, staged.path(manifest_path))
        publisher.publish()
    significant = int(result.enrichment["significant"].sum())
    raw_significant = int(result.enrichment["statistically_significant"].sum())
    print(
        f"Tested {len(result.enrichment)} RBPs; {significant} primary significant "
        f"({raw_significant} raw statistical calls). "
        f"Wrote {enrichment_path} and {annotation_path}"
    )


if __name__ == "__main__":
    main()
