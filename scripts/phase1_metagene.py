#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.config import MODIFICATIONS, modification
from rna_mod_pipeline.cli import (
    StagedOutputFiles,
    project_protected_trees,
    resolve_path,
)
from rna_mod_pipeline.phase1.metagene import annotate_metagene
from rna_mod_pipeline.phase1.plots import plot_metagene_summary
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_files,
    write_manifest,
)


DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parents[2]


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description="Map one modification dataset to complete GENCODE Basic transcripts."
    )
    result.add_argument(
        "--modification",
        required=True,
        choices=MODIFICATIONS,
    )
    result.add_argument("--project-root", default=DEFAULT_PROJECT_ROOT)
    result.add_argument(
        "--input",
        help="Site TSV. m6A defaults to the DRACH-annotated all-context table.",
    )
    result.add_argument("--gtf", help="GENCODE GTF; default data/gencode.v44.annotation.gtf.")
    result.add_argument("--output", help="Metagene-annotated TSV.")
    result.add_argument("--plot-dir", help="QC plot directory.")
    result.add_argument("--manifest", help="Run-manifest JSON path.")
    result.add_argument("--skip-plots", action="store_true")
    result.add_argument("--overwrite", action="store_true")
    return result


def main() -> None:
    args = parser().parse_args()
    config = modification(args.modification)
    root = Path(args.project_root).expanduser().resolve()
    phase_dir = root / "refactored_outputs" / args.modification / "phase1"
    input_suffix = "_drach" if config.drach_sites else ""
    input_path = (
        resolve_path(root, args.input)
        if args.input
        else phase_dir / f"filtered_{config.short_name}{input_suffix}.tsv"
    )
    gtf_path = (
        resolve_path(root, args.gtf)
        if args.gtf
        else root / "data" / "gencode.v44.annotation.gtf"
    )
    output_path = (
        resolve_path(root, args.output)
        if args.output
        else phase_dir / f"filtered_{config.short_name}_metagene.tsv"
    )
    plot_dir = (
        resolve_path(root, args.plot_dir)
        if args.plot_dir
        else phase_dir / "plots"
    )
    manifest_path = (
        resolve_path(root, args.manifest)
        if args.manifest
        else phase_dir / "phase1_metagene_manifest.json"
    )
    inputs = [input_path, gtf_path]
    plot_paths = [
        plot_dir / "12_metagene_profile.png",
        plot_dir / "14_region_distribution.png",
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
        expected_workflow="phase1_metagene",
        ownership_manifest=manifest_path,
        project_root=root,
        expected_parameters={"modification": args.modification},
        inputs=inputs,
        protected_roots=[root],
        protected_trees=project_protected_trees(root),
    )
    with publisher as staged:
        staged_output = staged.path(output_path)
        annotated = annotate_metagene(
            input_path,
            gtf_path,
            staged_output,
            expected_modification=args.modification,
        )
        if not args.skip_plots:
            plot_metagene_summary(
                annotated, staged.directory(plot_dir), config.display
            )
        mapped = int(annotated["metagene_pos"].notna().sum())
        staged_outputs = [staged_output]
        final_outputs = [output_path]
        for final_path in plot_paths if not args.skip_plots else []:
            staged_outputs.append(staged.path(final_path))
            final_outputs.append(final_path)
        manifest = make_manifest(
            "phase1_metagene",
            root,
            inputs,
            {
                "modification": args.modification,
                "input_modification_code_validation": (
                    "verified; m6A replicate identity is not encoded by mod_code"
                ),
                "transcript_policy": (
                    "complete GENCODE Basic protein-coding transcripts"
                ),
                "transcript_selection": (
                    "longest complete transcript, then transcript_id"
                ),
                "sites": len(annotated),
                "mapped_sites": mapped,
            },
            staged_outputs,
        )
        remap_manifest_output_files(
            manifest, staged_outputs, final_outputs, root
        )
        write_manifest(manifest, staged.path(manifest_path))
        publisher.publish()
    print(f"Mapped {mapped:,}/{len(annotated):,} sites; wrote {output_path}")


if __name__ == "__main__":
    main()
