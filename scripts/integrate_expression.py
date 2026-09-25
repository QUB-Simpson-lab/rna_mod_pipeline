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
from rna_mod_pipeline.expression import (
    annotate_rbp_expression,
    load_depmap,
    load_nanopore,
    modification_expression_association,
)
from rna_mod_pipeline.expression.plots import plot_gene_association, plot_rbp_expression
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_paths,
    write_manifest,
)
from rna_mod_pipeline.schemas import (
    validate_modification_codes,
    validate_result_identity,
)


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Integrate RBP and transcript expression")
    command.add_argument(
        "--modification", required=True, choices=("m6a", "m5c", "pseu", "m6a_rep2")
    )
    command.add_argument("--project-root", default=str(Path(__file__).resolve().parents[2]))
    command.add_argument("--source", required=True, choices=("nanopore", "depmap"))
    command.add_argument(
        "--design",
        choices=("loose", "transcript-region"),
        default="loose",
    )
    command.add_argument("--context", choices=("all", "5UTR", "CDS", "3UTR"), default="all")
    command.add_argument("--database", choices=("ornament", "encori", "postar3"))
    command.add_argument("--expression", required=True)
    command.add_argument("--enrichment", required=True)
    command.add_argument("--metagene", required=True)
    command.add_argument("--output-dir", required=True)
    command.add_argument("--model-id", default="ACH-000768")
    command.add_argument("--tpm-threshold", type=float, default=1.0)
    command.add_argument("--nanopore-value-column", default="cpm_0h")
    command.add_argument("--nanopore-count-column", default="count_0h")
    command.add_argument("--gene-lengths")
    command.add_argument("--overwrite", action="store_true")
    return command


def main() -> None:
    args = parser().parse_args()
    root = Path(args.project_root).expanduser().resolve()
    if args.tpm_threshold < 0:
        raise ValueError("--tpm-threshold must be non-negative")
    expression_path = resolve_path(root, args.expression)
    enrichment_path = resolve_path(root, args.enrichment)
    metagene_path = resolve_path(root, args.metagene)
    gene_lengths = resolve_path(root, args.gene_lengths) if args.gene_lengths else None
    enrichment_labels = pd.read_csv(enrichment_path, sep="\t")
    database_column = next(
        (
            column
            for column in ("database", "analysis_database")
            if column in enrichment_labels
        ),
        None,
    )
    if args.database and database_column:
        enrichment_labels = enrichment_labels[
            enrichment_labels[database_column].astype(str).str.lower()
            == args.database
        ]
    identity = validate_result_identity(
        enrichment_labels,
        expected_modification=args.modification,
        expected_database=args.database,
        expected_design=args.design,
        context=args.context,
        table_name="RBP enrichment table",
    )
    metagene_header = pd.read_csv(metagene_path, sep="\t", nrows=0)
    metagene_codes = (
        pd.read_csv(metagene_path, sep="\t", usecols=["mod_code"])
        if "mod_code" in metagene_header
        else metagene_header
    )
    metagene_identity = validate_modification_codes(
        metagene_codes,
        args.modification,
        "metagene table",
    )
    if args.source == "nanopore":
        profile = load_nanopore(
            expression_path,
            value_column=args.nanopore_value_column,
            count_column=args.nanopore_count_column,
        )
    else:
        profile = load_depmap(
            expression_path, model_id=args.model_id, threshold=args.tpm_threshold
        )
    rbps = annotate_rbp_expression(
        enrichment_path,
        profile,
        database=args.database,
        context=args.context,
    )
    genes, correlations = modification_expression_association(
        metagene_path,
        profile,
        gene_lengths_path=gene_lengths,
        context=args.context,
    )
    for table in (rbps, genes, correlations):
        table["modification"] = args.modification
        table["design"] = args.design
        table["context_region"] = args.context
        table["database_scope"] = args.database or "all_rows_in_input"
    inputs = [expression_path, enrichment_path, metagene_path]
    if gene_lengths:
        inputs.append(gene_lengths)
    target = resolve_path(root, args.output_dir)
    publisher = StagedOutputDirectory(
        target,
        args.overwrite,
        "expression_integration",
        inputs=inputs,
        protected_roots=(root,),
        protected_trees=project_protected_trees(root),
    )
    with publisher as output:
        output_files = [
            output / "rbp_expression.tsv",
            output / "gene_modification_expression.tsv",
            output / "expression_associations.tsv",
            output / "rbp_expression.png",
            output / "expression_modification_count.png",
        ]
        rbps.to_csv(output_files[0], sep="\t", index=False)
        genes.to_csv(output_files[1], sep="\t", index=False)
        correlations.to_csv(output_files[2], sep="\t", index=False)
        plot_rbp_expression(rbps, output_files[3])
        plot_gene_association(genes, output_files[4])
        manifest = make_manifest(
            "expression_integration",
            args.project_root,
            inputs,
            {
                "modification": args.modification,
                "design": args.design,
                "context": args.context,
                "database": args.database,
                "source": profile.source,
                "cell_line": profile.cell_line,
                "threshold": profile.threshold,
                "caveat": profile.absence_interpretation,
                "input_identity_status": {
                    "enrichment": identity,
                    "metagene": metagene_identity,
                },
            },
            output_files,
        )
        remap_manifest_output_paths(
            manifest, output_files, output, target, args.project_root
        )
        write_manifest(manifest, output / "run_manifest.json")
        publisher.publish()
    print(f"Wrote expression integration to {target}")


if __name__ == "__main__":
    main()
