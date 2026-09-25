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
from rna_mod_pipeline.knockrbp import (
    load_and_resolve_degs,
    load_loose_targets,
    load_transcript_targets,
    dataset_summary,
    orthogonal_validation,
    regulatory_network,
    regulatory_summary,
    target_summary,
)
from rna_mod_pipeline.knockrbp.plots import plot_orthogonal, plot_regulatory_network
from rna_mod_pipeline.knockrbp.provenance import (
    annotated_identity,
    cross_database_schema_status,
    enrichment_identity,
    knockrbp_manifest_parameters,
    transcript_run_identity,
    validate_design_inputs,
)
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_paths,
    write_manifest,
)
from rna_mod_pipeline.schemas import validate_result_identity


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Validate RBP targets with KnockRBP")
    command.add_argument(
        "--modification", required=True, choices=("m6a", "m5c", "pseu", "m6a_rep2")
    )
    command.add_argument("--project-root", default=str(Path(__file__).resolve().parents[2]))
    command.add_argument(
        "--design",
        required=True,
        choices=("loose", "transcript-region"),
    )
    command.add_argument("--context", default="all", choices=("all", "5UTR", "CDS", "3UTR"))
    command.add_argument("--knockrbp-dir", required=True)
    command.add_argument("--metadata", required=True)
    command.add_argument("--cell-lines", nargs="+", default=["MDA-MB-231", "MDA-MB-231-LM2"])
    command.add_argument("--dataset-ids", nargs="*")
    command.add_argument("--encori-annotated")
    command.add_argument("--postar3-annotated")
    command.add_argument("--encori-enrichment")
    command.add_argument("--postar3-enrichment")
    command.add_argument("--transcript-run-dir")
    command.add_argument("--cross-database")
    command.add_argument("--log2fc-cutoff", type=float, default=0.5)
    command.add_argument("--padj-cutoff", type=float, default=0.05)
    command.add_argument("--regulatory-plot-top-n", type=int, default=40)
    command.add_argument("--output-dir", required=True)
    command.add_argument("--overwrite", action="store_true")
    return command


def main() -> None:
    args = parser().parse_args()
    root = Path(args.project_root).expanduser().resolve()
    if args.log2fc_cutoff < 0:
        raise ValueError("--log2fc-cutoff must be non-negative")
    if not 0 < args.padj_cutoff <= 1:
        raise ValueError("--padj-cutoff must lie in (0, 1]")
    if args.regulatory_plot_top_n < 1:
        raise ValueError("--regulatory-plot-top-n must be at least 1")
    if args.dataset_ids and len(args.dataset_ids) != len(set(args.dataset_ids)):
        raise ValueError("--dataset-ids contains duplicates")
    validate_design_inputs(
        args.design,
        encori_annotated=args.encori_annotated,
        postar3_annotated=args.postar3_annotated,
        encori_enrichment=args.encori_enrichment,
        postar3_enrichment=args.postar3_enrichment,
        transcript_run_dir=args.transcript_run_dir,
    )
    metadata_path = resolve_path(root, args.metadata)
    knockrbp_dir = resolve_path(root, args.knockrbp_dir)
    design = args.design
    degs, audit = load_and_resolve_degs(
        knockrbp_dir,
        metadata_path,
        allowed_cell_lines=args.cell_lines,
        dataset_ids=args.dataset_ids,
        log2fc_cutoff=args.log2fc_cutoff,
        padj_cutoff=args.padj_cutoff,
    )
    if not degs:
        raise ValueError("No eligible KnockRBP datasets were loaded")
    input_paths = [metadata_path]
    input_paths.extend(entry["source_file"] for entry in degs.values())
    requested_rbps = sorted({str(entry["rbp"]) for entry in degs.values()})
    identity_status = {}
    if design == "loose":
        if not args.encori_annotated or not args.postar3_annotated:
            raise ValueError("Loose design requires both annotated site tables")
        encori_path = resolve_path(root, args.encori_annotated)
        postar3_path = resolve_path(root, args.postar3_annotated)
        default_resource_paths = {
            "encori": encori_path.parent / "encori_enrichment_results.tsv",
            "postar3": postar3_path.parent / "postar3_enrichment_results.tsv",
        }
        resource_paths = {
            "encori": (
                resolve_path(root, args.encori_enrichment)
                if args.encori_enrichment
                else default_resource_paths["encori"]
            ),
            "postar3": (
                resolve_path(root, args.postar3_enrichment)
                if args.postar3_enrichment
                else default_resource_paths["postar3"]
            ),
        }
        for database, explicit_value in (
            ("encori", args.encori_enrichment),
            ("postar3", args.postar3_enrichment),
        ):
            if explicit_value and not resource_paths[database].is_file():
                raise FileNotFoundError(
                    f"{database} enrichment table not found: "
                    f"{resource_paths[database]}"
                )
        resource_paths = {
            database: path
            for database, path in resource_paths.items()
            if path.is_file()
        }
        identity_status = {
            "encori": annotated_identity(
                encori_path, "encori", args.modification
            ),
            "postar3": annotated_identity(
                postar3_path, "postar3", args.modification
            ),
        }
        identity_status.update(
            {
                f"{database}_resource": enrichment_identity(
                    path, database, args.modification
                )
                for database, path in resource_paths.items()
            }
        )
        targets, universe = load_loose_targets(
            {
                "encori": encori_path,
                "postar3": postar3_path,
            },
            requested_rbps,
            context=args.context,
            resource_paths=resource_paths,
        )
        input_paths.extend([encori_path, postar3_path, *resource_paths.values()])
    else:
        if not args.transcript_run_dir:
            raise ValueError("Transcript design requires --transcript-run-dir")
        run = resolve_path(root, args.transcript_run_dir)
        identity_status = {
            "transcript_run": transcript_run_identity(
                run, args.modification
            )
        }
        targets, universe = load_transcript_targets(
            run, requested_rbps, context=args.context
        )
        input_paths.extend(
            [
                run / "case_site_assignments.tsv.gz",
                run / "case_binding_matrix_encori.npz",
                run / "case_binding_matrix_encori_rows.tsv",
                run / "case_binding_matrix_postar3.npz",
                run / "case_binding_matrix_postar3_rows.tsv",
                *(
                    [run / "resource_loading_report.tsv"]
                    if (run / "resource_loading_report.tsv").is_file()
                    else []
                ),
                *(
                    [run / "run_manifest.json"]
                    if (run / "run_manifest.json").is_file()
                    else []
                ),
            ]
        )
    orthogonal = orthogonal_validation(
        degs, targets, universe, design=design, context=args.context
    )
    if not orthogonal.empty:
        orthogonal.insert(0, "modification", args.modification)
    datasets = dataset_summary(degs)
    targets_table = target_summary(targets)
    for table in (datasets, targets_table):
        table.insert(0, "context_region", args.context)
        table.insert(0, "design", design)
        table.insert(0, "modification", args.modification)
    network = pd.DataFrame()
    cross_database_quality_schema = "not_supplied"
    if args.cross_database:
        cross_path = resolve_path(root, args.cross_database)
        cross = pd.read_csv(cross_path, sep="\t")
        identity_status["cross_database"] = validate_result_identity(
            cross,
            expected_modification=args.modification,
            expected_design=design,
            context=args.context,
            table_name="cross-database table",
        )
        input_paths.append(cross_path)
        cross_database_quality_schema = cross_database_schema_status(cross)
        network = regulatory_network(degs, cross)
        if not network.empty:
            network.insert(0, "context_region", args.context)
            network.insert(0, "design", design)
            network.insert(0, "modification", args.modification)
    target = resolve_path(root, args.output_dir)
    publisher = StagedOutputDirectory(
        target,
        args.overwrite,
        "knockrbp_validation",
        inputs=input_paths,
        protected_roots=(root,),
        protected_trees=project_protected_trees(root),
    )
    with publisher as output:
        output_files = [
            output / "orthogonal_validation_results.tsv",
            output / "duplicate_deg_resolution_audit.tsv",
            output / "orthogonal_validation.png",
            output / "eligible_knockrbp_datasets.tsv",
            output / "binding_target_summary.tsv",
        ]
        orthogonal.to_csv(output_files[0], sep="\t", index=False)
        audit.to_csv(output_files[1], sep="\t", index=False)
        plot_orthogonal(orthogonal, output_files[2])
        datasets.to_csv(output_files[3], sep="\t", index=False)
        targets_table.to_csv(output_files[4], sep="\t", index=False)

        if args.cross_database:
            network_path = output / "rbp_regulatory_network.tsv"
            summary_path = output / "rbp_regulatory_summary.tsv"
            plot_source_path = output / "rbp_regulatory_network_plot_source.tsv"
            plot_path = output / "rbp_regulatory_network.png"
            output_files.extend(
                [network_path, summary_path, plot_source_path, plot_path]
            )
            network.to_csv(network_path, sep="\t", index=False)
            summary = regulatory_summary(network)
            summary.insert(0, "context_region", args.context)
            summary.insert(0, "design", design)
            summary.insert(0, "modification", args.modification)
            summary.to_csv(summary_path, sep="\t", index=False)
            plot_source = plot_regulatory_network(
                network,
                plot_path,
                top_n=args.regulatory_plot_top_n,
            )
            plot_source.to_csv(plot_source_path, sep="\t", index=False)

        manifest = make_manifest(
            "knockrbp_validation",
            args.project_root,
            input_paths,
            knockrbp_manifest_parameters(
                options=vars(args),
                n_datasets=len(degs),
                network=network,
                schema_status=cross_database_quality_schema,
                identity_status=identity_status,
            ),
            output_files,
        )
        remap_manifest_output_paths(
            manifest, output_files, output, target, args.project_root
        )
        write_manifest(manifest, output / "run_manifest.json")
        publisher.publish()
    print(f"Wrote KnockRBP results to {target}")


if __name__ == "__main__":
    main()
