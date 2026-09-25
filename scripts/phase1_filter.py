#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.cli import (
    StagedOutputFiles,
    percentage,
    positive_int,
    project_protected_trees,
    resolve_path,
)
from rna_mod_pipeline.config import MODIFICATIONS, modification, project_path
from rna_mod_pipeline.phase1.filtering import (
    COVERAGE_BINS,
    COVERAGE_PLOT_MAX,
    filter_bedmethyl,
)
from rna_mod_pipeline.phase1.plots import plot_filter_qc
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_files,
    write_manifest,
)


DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parents[2]


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Filter one bedMethyl dataset by coverage and modification percentage."
    )
    result.add_argument(
        "--modification",
        required=True,
        choices=MODIFICATIONS,
        help="Run exactly one modification or m6A replicate.",
    )
    result.add_argument("--project-root", default=DEFAULT_PROJECT_ROOT)
    result.add_argument("--input", help="Raw bedMethyl; defaults from the modification.")
    result.add_argument("--output", help="Filtered TSV; defaults from the modification.")
    result.add_argument("--plot-dir", help="QC plot directory.")
    result.add_argument("--manifest", help="Run-manifest JSON path.")
    result.add_argument("--min-coverage", type=positive_int, default=20)
    result.add_argument("--min-fraction", type=percentage, default=20.0)
    result.add_argument("--chunk-size", type=positive_int, default=2_000_000)
    result.add_argument("--skip-plots", action="store_true")
    result.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace the selected output file; unrelated files are never removed.",
    )
    return result


def main() -> None:
    args = parser().parse_args()
    config = modification(args.modification)
    root = Path(args.project_root).expanduser().resolve()
    phase_dir = root / "refactored_outputs" / args.modification / "phase1"
    input_path = (
        resolve_path(root, args.input)
        if args.input
        else project_path(root, config.bedmethyl)
    )
    output_path = (
        resolve_path(root, args.output)
        if args.output
        else phase_dir / f"filtered_{config.short_name}.tsv"
    )
    plot_dir = (
        resolve_path(root, args.plot_dir)
        if args.plot_dir
        else phase_dir / "plots"
    )
    manifest_path = (
        resolve_path(root, args.manifest)
        if args.manifest
        else phase_dir / "phase1_filter_manifest.json"
    )
    plot_paths = [
        plot_dir / "01_coverage_distribution.png",
        plot_dir / "02_fraction_modified_distribution.png",
        plot_dir / "03_chromosome_distribution.png",
        plot_dir / "04_threshold_sensitivity.png",
        plot_dir / "04a_joint_threshold_heatmap.png",
        plot_dir / "05_filtered_chromosome_distribution.png",
        plot_dir / "06_filtered_fraction_distribution.png",
        plot_dir / "07_coverage_vs_fraction.png",
    ]
    destinations = [
        output_path,
        manifest_path,
        *([] if args.skip_plots else plot_paths),
    ]
    if plot_dir in {output_path, manifest_path} or len(set(destinations)) != len(
        destinations
    ):
        raise ValueError("Output, manifest, and plot paths must be different")
    publisher = StagedOutputFiles(
        destinations,
        overwrite=args.overwrite,
        expected_workflow="phase1_filter",
        ownership_manifest=manifest_path,
        project_root=root,
        expected_parameters={"modification": args.modification},
        inputs=[input_path],
        protected_roots=[root],
        protected_trees=project_protected_trees(root),
    )
    with publisher as staged:
        staged_output = staged.path(output_path)
        filtered, summary = filter_bedmethyl(
            input_path,
            staged_output,
            coverage_min=args.min_coverage,
            fraction_min=args.min_fraction,
            chunk_size=args.chunk_size,
            expected_mod_codes=config.expected_mod_codes,
        )
        if not args.skip_plots:
            plot_filter_qc(
                summary,
                filtered,
                staged.directory(plot_dir),
                config.display,
            )
        staged_outputs = [staged_output]
        final_outputs = [output_path]
        for final_path in plot_paths if not args.skip_plots else []:
            staged_outputs.append(staged.path(final_path))
            final_outputs.append(final_path)
        manifest = make_manifest(
            "phase1_filter",
            root,
            [input_path],
            {
                "modification": args.modification,
                "minimum_coverage": args.min_coverage,
                "minimum_fraction_percent": args.min_fraction,
                "chunk_size": args.chunk_size,
                "expected_mod_codes": sorted(config.expected_mod_codes),
                "total_sites": summary.total_sites,
                "retained_sites": summary.retained_sites,
                "coverage_plot_max": COVERAGE_PLOT_MAX,
                "coverage_at_or_above_plot_max": (
                    summary.coverage_at_or_above_plot_max
                ),
                "coverage_histogram_max": float(COVERAGE_BINS[-1]),
                "coverage_above_histogram_max": (
                    summary.coverage_above_histogram_max
                ),
            },
            staged_outputs,
        )
        remap_manifest_output_files(
            manifest, staged_outputs, final_outputs, root
        )
        write_manifest(manifest, staged.path(manifest_path))
        publisher.publish()
    print(
        f"Retained {summary.retained_sites:,}/{summary.total_sites:,} sites; "
        f"wrote {output_path}"
    )


if __name__ == "__main__":
    main()
