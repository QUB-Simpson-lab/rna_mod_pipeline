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
from rna_mod_pipeline.crossdb import compare_databases, load_database_results
from rna_mod_pipeline.crossdb.io import DATABASES, default_input_paths
from rna_mod_pipeline.crossdb.plots import create_crossdb_plots
from rna_mod_pipeline.evidence_quality import evidence_quality_manifest
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_paths,
    write_manifest,
)


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Compare RBP enrichment databases")
    command.add_argument("--project-root", default=str(Path(__file__).resolve().parents[2]))
    command.add_argument("--modification", required=True, choices=("m6a", "m5c", "pseu", "m6a_rep2"))
    command.add_argument(
        "--design",
        choices=("loose", "transcript-region"),
        default="loose",
    )
    command.add_argument("--context", choices=("all", "5UTR", "CDS", "3UTR"), default="all")
    command.add_argument("--ornament")
    command.add_argument("--encori")
    command.add_argument("--postar3")
    command.add_argument("--databases", nargs="+", choices=DATABASES, default=list(DATABASES))
    command.add_argument("--fdr", type=float, default=0.05)
    command.add_argument("--allow-missing", action="store_true")
    command.add_argument("--skip-plots", action="store_true")
    command.add_argument("--plot-top-n", type=int, default=20)
    command.add_argument("--plot-dpi", type=int, default=180)
    command.add_argument("--output-dir", required=True)
    command.add_argument("--overwrite", action="store_true")
    return command


def main() -> None:
    args = parser().parse_args()
    root = Path(args.project_root).expanduser().resolve()
    if not 0 < args.fdr <= 1:
        raise ValueError("--fdr must lie in (0, 1]")
    if args.plot_top_n < 1:
        raise ValueError("--plot-top-n must be at least 1")
    if args.plot_dpi < 72:
        raise ValueError("--plot-dpi must be at least 72")
    if len(args.databases) != len(set(args.databases)):
        raise ValueError("--databases contains duplicates")
    design = args.design
    defaults = default_input_paths(
        args.project_root, args.modification, design, context=args.context
    )
    explicit = {database: getattr(args, database) for database in DATABASES}
    ignored_overrides = [
        database
        for database, value in explicit.items()
        if value and database not in args.databases
    ]
    if ignored_overrides:
        raise ValueError(
            "Explicit input paths were supplied for databases not selected by "
            f"--databases: {', '.join(ignored_overrides)}"
        )
    tables = {}
    inputs = {}
    for database in args.databases:
        path = resolve_path(root, explicit[database]) if explicit[database] else defaults[database]
        if not path.is_file():
            if args.allow_missing:
                continue
            raise FileNotFoundError(f"{database} enrichment table not found: {path}")
        tables[database] = load_database_results(
            path,
            database,
            context=args.context,
            project_root=args.project_root,
            expected_modification=args.modification,
            expected_design=design,
        )
        inputs[database] = str(path)

    combined, pairwise = compare_databases(
        tables, fdr_threshold=args.fdr, allow_missing=args.allow_missing
    )
    quality_by_database = {}
    for database in tables:
        quality_frame = pd.DataFrame(
            {
                suffix: combined[f"{database}_{suffix}"]
                for suffix in (
                    "inference_eligible",
                    "statistically_significant",
                    "significant",
                    "input_completeness",
                    "quality_diagnostics_available",
                )
                if f"{database}_{suffix}" in combined
            }
        )
        quality_by_database[database] = evidence_quality_manifest(
            quality_frame,
            fdr_threshold=args.fdr,
        )
    for name, value in (
        ("modification", args.modification),
        ("design", design),
        ("context_region", args.context),
    ):
        combined.insert(0, name, value)
        for table in pairwise.values():
            table.insert(0, name, value)
    target = resolve_path(root, args.output_dir)
    publisher = StagedOutputDirectory(
        target,
        args.overwrite,
        "cross_database",
        inputs=inputs.values(),
        protected_roots=(root,),
        protected_trees=project_protected_trees(root),
    )
    with publisher as output:
        output_files = []
        for name, table in (
            ("cross_database_results.tsv", combined),
            ("triple_significant.tsv", combined[combined["triple_significant"]]),
            (
                "triple_statistically_significant.tsv",
                combined[combined["triple_statistically_significant"]],
            ),
            (
                "triple_direction_concordant.tsv",
                combined[combined["triple_direction_concordant"]],
            ),
            (
                "sensitivity_only_direction_concordant.tsv",
                combined[
                    combined["consensus_evidence_quality"].eq(
                        "sensitivity_only_consensus_quality_excluded"
                    )
                ],
            ),
            ("direction_conflicts.tsv", combined[combined["direction_conflict"]]),
        ):
            destination = output / name
            table.to_csv(destination, sep="\t", index=False)
            output_files.append(destination)
        for name, table in pairwise.items():
            destination = output / f"pairwise_{name}.tsv"
            table.to_csv(destination, sep="\t", index=False)
            output_files.append(destination)
        if not args.skip_plots:
            output_files.extend(
                create_crossdb_plots(
                    combined,
                    pairwise,
                    output,
                    modification=args.modification,
                    design=design,
                    context=args.context,
                    fdr_threshold=args.fdr,
                    top_n=args.plot_top_n,
                    dpi=args.plot_dpi,
                )
            )
        manifest = make_manifest(
            "cross_database",
            args.project_root,
            inputs.values(),
            {
                "modification": args.modification,
                "design": design,
                "context": args.context,
                "fdr_threshold": args.fdr,
                "allow_missing": args.allow_missing,
                "plots_generated": not args.skip_plots,
                "plot_top_n": args.plot_top_n,
                "plot_dpi": args.plot_dpi,
                "plot_confidence_intervals_drawn": False,
                "input_identity_status": {
                    database: {
                        column.removeprefix("input_").removesuffix("_status"): str(
                            table[column].iloc[0]
                        )
                        for column in table
                        if column.startswith("input_") and column.endswith("_status")
                    }
                    for database, table in tables.items()
                },
                "definitions": {
                    "triple_significant": (
                        "FDR below threshold and inference eligible in all three databases"
                    ),
                    "triple_statistically_significant": (
                        "Raw FDR below threshold in all three databases, retained "
                        "for sensitivity auditing regardless of evidence quality"
                    ),
                    "triple_direction_concordant": (
                        "quality-qualified triple significant with the same OR direction"
                    ),
                    "direction_conflict": "at least two significant databases with opposite directions",
                },
                "plot_correlation_scope": (
                    "primary Spearman rho uses shared finite inference-eligible "
                    "estimates; raw all-finite rho is labelled sensitivity"
                ),
                "evidence_quality_by_database": quality_by_database,
                "consensus_counts": {
                    "n_triple_primary_significant": int(
                        combined["triple_significant"].sum()
                    ),
                    "n_triple_raw_statistically_significant": int(
                        combined["triple_statistically_significant"].sum()
                    ),
                    "n_primary_direction_concordant": int(
                        combined["direction_concordant_consensus"].sum()
                    ),
                    "n_sensitivity_only_direction_concordant": int(
                        combined["consensus_evidence_quality"].eq(
                            "sensitivity_only_consensus_quality_excluded"
                        ).sum()
                    ),
                },
            },
            output_files,
        )
        remap_manifest_output_paths(
            manifest, output_files, output, target, args.project_root
        )
        write_manifest(manifest, output / "run_manifest.json")
        publisher.publish()
    print(f"Wrote {len(combined)} canonical RBPs to {target}")


if __name__ == "__main__":
    main()
