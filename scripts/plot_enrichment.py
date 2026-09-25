#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.cli import (
    StagedOutputDirectory,
    project_protected_trees,
    resolve_path,
)
from rna_mod_pipeline.config import MODIFICATIONS
from rna_mod_pipeline.plotting import generate_all_plots, generate_regional_plots
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_paths,
    write_manifest,
)
from rna_mod_pipeline.schemas import validate_result_identity


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description="Plot refactored RBP enrichment results without confidence intervals"
    )
    command.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[2]),
    )
    command.add_argument("--modification", required=True, choices=tuple(MODIFICATIONS))
    command.add_argument(
        "--design",
        default="transcript-region",
        choices=("transcript-region",),
    )
    command.add_argument("--context", default="all", choices=("all",))
    command.add_argument("--results", required=True)
    command.add_argument("--regional-results")
    command.add_argument("--output-dir", required=True)
    command.add_argument("--analysis-id")
    command.add_argument("--fdr", type=float, default=0.05)
    command.add_argument("--top-n", type=int, default=30)
    command.add_argument("--dpi", type=int, default=180)
    command.add_argument("--overwrite", action="store_true")
    return command


def main() -> None:
    args = parser().parse_args()
    if not 0 < args.fdr <= 1:
        raise ValueError("--fdr must lie in (0, 1]")
    if args.top_n <= 0:
        raise ValueError("--top-n must be positive")
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive")
    root = Path(args.project_root).expanduser().resolve()
    results_path = resolve_path(root, args.results)
    regional_path = (
        resolve_path(root, args.regional_results)
        if args.regional_results
        else None
    )
    results = pd.read_csv(results_path, sep="\t")
    identity_status = {
        "results": validate_result_identity(
            results,
            expected_modification=args.modification,
            expected_design=args.design,
            context=args.context,
            table_name="primary plot input",
        )
    }
    embedded_id = None
    if "analysis_id" in results and results["analysis_id"].nunique() == 1:
        embedded_id = str(results["analysis_id"].iloc[0])
    expected_id = embedded_id or args.modification
    if args.analysis_id and args.analysis_id.strip().lower() != expected_id.lower():
        raise ValueError(
            f"--analysis-id {args.analysis_id!r} does not match input identity "
            f"{expected_id!r}"
        )
    analysis_id = expected_id
    inputs = [results_path, *([regional_path] if regional_path is not None else [])]
    target = resolve_path(root, args.output_dir)
    publisher = StagedOutputDirectory(
        target,
        args.overwrite,
        "plot_enrichment",
        inputs=inputs,
        protected_roots=(root,),
        protected_trees=project_protected_trees(root),
    )
    with publisher as output:
        generate_all_plots(
            results, output, analysis_id, args.fdr, args.top_n, args.dpi
        )
        if regional_path is not None:
            regional = pd.read_csv(regional_path, sep="\t")
            if "region_filter" not in regional:
                raise ValueError(
                    "regional plot input is missing column: 'region_filter'"
                )
            regions = {
                str(value).strip().upper()
                for value in regional["region_filter"].dropna()
            }
            allowed_regions = {"5UTR", "CDS", "3UTR"}
            if not regions or not regions.issubset(allowed_regions):
                raise ValueError(
                    "regional plot input has unexpected region_filter values: "
                    f"{sorted(regions)}"
                )
            identity_status["regional_results"] = validate_result_identity(
                regional,
                expected_modification=args.modification,
                expected_design=args.design,
                context=sorted(regions)[0],
                table_name="regional plot input",
            )
            generate_regional_plots(
                regional,
                output / "regional_top_plots",
                analysis_id,
                args.fdr,
                args.top_n,
                args.dpi,
            )
        output_files = [
            path
            for path in output.rglob("*")
            if path.is_file() and path.name != "run_manifest.json"
        ]
        manifest = make_manifest(
            "plot_enrichment",
            root,
            inputs,
            {
                "analysis_id": analysis_id,
                "modification": args.modification,
                "design": args.design,
                "context": args.context,
                "input_identity_status": identity_status,
                "fdr_threshold": args.fdr,
                "top_n": args.top_n,
                "dpi": args.dpi,
                "regional_results_included": regional_path is not None,
                "confidence_intervals_drawn": False,
                "all_rbps_overview_policy": (
                    "finite positive non-boundary estimates on a descending "
                    "log2 Mantel-Haenszel odds-ratio scale; OR=0 boundaries "
                    "represented by individual negative-infinity bars; "
                    "other boundary and non-estimable tests counted on the "
                    "figure and retained in the source TSV"
                ),
            },
            output_files,
        )
        remap_manifest_output_paths(
            manifest, output_files, output, target, root
        )
        write_manifest(manifest, output / "run_manifest.json")
        publisher.publish()
    print(f"Plots written to {target}")


if __name__ == "__main__":
    main()
