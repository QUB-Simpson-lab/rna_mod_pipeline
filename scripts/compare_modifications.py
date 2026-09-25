#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.cli import (
    StagedOutputDirectory,
    project_protected_trees,
    resolve_path,
)
from rna_mod_pipeline.crossdb.io import DATABASES, default_input_paths, load_database_results
import pandas as pd

from rna_mod_pipeline.config import modification as modification_config
from rna_mod_pipeline.evidence_quality import evidence_quality_manifest
from rna_mod_pipeline.crossmod import (
    compare_modifications,
    create_crossmod_plots,
    phase1_summary,
    plot_metagene_overlay,
    plot_site_summary,
)
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_paths,
    write_manifest,
)
from rna_mod_pipeline.schemas import validate_modification_codes


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Compare RBP effects across modifications")
    command.add_argument(
        "--modifications",
        nargs="+",
        choices=("m6a", "m5c", "pseu", "m6a_rep2"),
        default=["m6a", "m5c", "pseu"],
    )
    command.add_argument("--databases", nargs="+", choices=DATABASES, default=list(DATABASES))
    command.add_argument(
        "--design",
        choices=("loose", "transcript-region"),
        default="loose",
    )
    command.add_argument("--context", default="all", choices=("all", "5UTR", "CDS", "3UTR"))
    command.add_argument(
        "--input",
        action="append",
        nargs=3,
        metavar=("MODIFICATION", "DATABASE", "PATH"),
        help="Override one default input; may be repeated",
    )
    command.add_argument(
        "--metagene",
        action="append",
        nargs=2,
        metavar=("MODIFICATION", "PATH"),
        help="Override a Phase-1 metagene table; may be repeated",
    )
    command.add_argument("--skip-phase1-context", action="store_true")
    command.add_argument("--allow-missing", action="store_true")
    command.add_argument("--fdr", type=float, default=0.05)
    command.add_argument(
        "--project-root", default=str(Path(__file__).resolve().parents[2])
    )
    command.add_argument("--output-dir", required=True)
    command.add_argument("--overwrite", action="store_true")
    return command


def main() -> None:
    args = parser().parse_args()
    root = Path(args.project_root).expanduser().resolve()
    if not 0 < args.fdr <= 1:
        raise ValueError("--fdr must lie in (0, 1]")
    if len(args.modifications) != len(set(args.modifications)):
        raise ValueError("--modifications contains duplicates")
    if len(args.databases) != len(set(args.databases)):
        raise ValueError("--databases contains duplicates")
    design = args.design
    input_rows = args.input or []
    input_keys = [(modification, database) for modification, database, _ in input_rows]
    if len(input_keys) != len(set(input_keys)):
        raise ValueError("--input contains duplicate modification/database overrides")
    invalid_inputs = [
        key
        for key in input_keys
        if key[0] not in args.modifications or key[1] not in args.databases
    ]
    if invalid_inputs:
        raise ValueError(f"--input override is outside the selected scope: {invalid_inputs}")
    overrides = {
        (modification, database): resolve_path(root, path)
        for modification, database, path in input_rows
    }
    tables = {}
    paths = []
    missing = []
    for modification in args.modifications:
        defaults = default_input_paths(
            args.project_root, modification, design, context=args.context
        )
        for database in args.databases:
            path = overrides.get((modification, database), defaults[database])
            if not path.is_file():
                missing.append((modification, database, str(path)))
                continue
            tables[(modification, database)] = load_database_results(
                path,
                database,
                context=args.context,
                project_root=args.project_root,
                expected_modification=modification,
                expected_design=design,
            )
            paths.append(path)
    if missing and not args.allow_missing:
        details = "\n".join(f"{mod}/{db}: {path}" for mod, db, path in missing)
        raise FileNotFoundError(f"Cross-modification inputs are missing:\n{details}")
    matrix, patterns, correlations, overlaps = compare_modifications(
        tables,
        modifications=tuple(args.modifications),
        databases=tuple(args.databases),
        fdr_threshold=args.fdr,
    )
    quality_summary_input = matrix.copy()
    quality_summary_input["significant"] = quality_summary_input[
        "state"
    ].str.startswith("significant_")
    enrichment_identity = {
        f"{modification}/{database}": {
            column.removeprefix("input_").removesuffix("_status"): str(
                table[column].iloc[0]
            )
            for column in table
            if column.startswith("input_") and column.endswith("_status")
        }
        for (modification, database), table in tables.items()
    }
    metagene_inputs = {}
    missing_metagene = []
    metagene_identity = {}
    if not args.skip_phase1_context:
        metagene_rows = args.metagene or []
        metagene_keys = [modification for modification, _ in metagene_rows]
        if len(metagene_keys) != len(set(metagene_keys)):
            raise ValueError("--metagene contains duplicate modification overrides")
        invalid_metagene = [
            key for key in metagene_keys if key not in args.modifications
        ]
        if invalid_metagene:
            raise ValueError(
                f"--metagene override is outside the selected scope: {invalid_metagene}"
            )
        overrides_metagene = {
            modification: resolve_path(root, path)
            for modification, path in metagene_rows
        }
        for modification in args.modifications:
            short_name = modification_config(modification).short_name
            default = (
                root
                / "refactored_outputs"
                / modification
                / "phase1"
                / f"filtered_{short_name}_metagene.tsv"
            )
            path = overrides_metagene.get(modification, default)
            if path.is_file():
                table = pd.read_csv(path, sep="\t")
                if not {"region", "metagene_pos"}.issubset(table):
                    raise ValueError(
                        f"Metagene table lacks region/metagene_pos columns: {path}"
                    )
                metagene_identity[modification] = validate_modification_codes(
                    table,
                    modification,
                    f"{modification} metagene table",
                )
                metagene_inputs[modification] = table
                paths.append(path)
            else:
                missing_metagene.append((modification, str(path)))
        if missing_metagene and not args.allow_missing:
            details = "\n".join(
                f"{modification}: {path}"
                for modification, path in missing_metagene
            )
            raise FileNotFoundError(
                f"Phase-1 metagene context inputs are missing:\n{details}"
            )
    target = resolve_path(root, args.output_dir)
    publisher = StagedOutputDirectory(
        target,
        args.overwrite,
        "cross_modification_comparison",
        inputs=paths,
        protected_roots=(root,),
        protected_trees=project_protected_trees(root),
    )
    with publisher as output:
        output_files = []
        for name, table in (
            ("rbp_modification_matrix.tsv", matrix),
            ("rbp_crossmod_patterns.tsv", patterns),
            ("pairwise_effect_correlations.tsv", correlations),
            ("pairwise_significant_overlaps.tsv", overlaps),
        ):
            path = output / name
            table.to_csv(path, sep="\t", index=False)
            output_files.append(path)
        output_files.extend(
            create_crossmod_plots(matrix, patterns, correlations, output)
        )
        if metagene_inputs:
            summary = phase1_summary(metagene_inputs)
            summary_path = output / "phase1_site_summary.tsv"
            overlay_path = output / "phase1_metagene_overlay.png"
            site_path = output / "phase1_site_count_summary.png"
            summary.to_csv(summary_path, sep="\t", index=False)
            plot_metagene_overlay(metagene_inputs, overlay_path)
            plot_site_summary(summary, site_path)
            output_files.extend([summary_path, overlay_path, site_path])
        manifest = make_manifest(
            "cross_modification_comparison",
            args.project_root,
            paths,
            {
                "modifications": args.modifications,
                "databases": args.databases,
                "design": design,
                "context": args.context,
                "fdr_threshold": args.fdr,
                "missing_inputs": missing,
                "input_identity_status": enrichment_identity,
                "absence_definition": "not_covered means no row/resource, not a null result",
                "phase1_context_included": bool(metagene_inputs),
                "phase1_context_included_modifications": [
                    modification
                    for modification in args.modifications
                    if modification in metagene_inputs
                ],
                "phase1_context_missing_modifications": [
                    modification for modification, _ in missing_metagene
                ],
                "phase1_context_missing_inputs": [
                    {"modification": modification, "path": path}
                    for modification, path in missing_metagene
                ],
                "phase1_context_identity_status": metagene_identity,
                "omitted_legacy_plot": (
                    "Manual functional-category pie charts are intentionally omitted "
                    "because that annotation layer is outside this refactor."
                ),
                "plot_inference_scope": {
                    "crossmod_effect_correlations.png": (
                        "quality-qualified finite effects"
                    ),
                    "crossmod_effect_correlations_raw_sensitivity.png": (
                        "all finite effects; sensitivity only"
                    ),
                    "crossmod_effect_heatmap.png": (
                        "quality-qualified finite effects; Q marks excluded cells"
                    ),
                    "crossmod_effect_heatmap_raw_sensitivity.png": (
                        "all finite effects; Q marks quality-excluded cells"
                    ),
                },
                "evidence_quality": evidence_quality_manifest(
                    quality_summary_input,
                    fdr_threshold=args.fdr,
                ),
            },
            output_files,
        )
        remap_manifest_output_paths(
            manifest, output_files, output, target, args.project_root
        )
        write_manifest(manifest, output / "run_manifest.json")
        publisher.publish()
    print(f"Wrote cross-modification comparison to {target}")


if __name__ == "__main__":
    main()
