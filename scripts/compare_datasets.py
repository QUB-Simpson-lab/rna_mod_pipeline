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
from rna_mod_pipeline.dataset_compare import run_dataset_comparison
from rna_mod_pipeline.evidence_quality import EVIDENCE_QUALITY_POLICY_ID
from rna_mod_pipeline.dataset_compare.profiles import (
    enrichment_pair_inventory,
    load_profile,
    validate_compatibility,
)
from rna_mod_pipeline.provenance import (
    make_manifest,
    remap_manifest_output_paths,
    write_manifest,
)


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description=(
            "Compare two independently processed same-modification datasets as "
            "a robustness analysis; do not pool raw biological replicates."
        )
    )
    command.add_argument(
        "--project-root", default=str(Path(__file__).resolve().parents[2])
    )
    command.add_argument(
        "--left-profile",
        "--profile-a",
        dest="left_profile",
        required=True,
        help="DatasetProfile JSON for the left/dataset-A workspace.",
    )
    command.add_argument(
        "--right-profile",
        "--profile-b",
        dest="right_profile",
        required=True,
        help="DatasetProfile JSON for the right/dataset-B workspace.",
    )
    command.add_argument("--output-dir", required=True)
    command.add_argument("--fdr", type=float, default=0.05)
    command.add_argument("--dpi", type=int, default=180)
    command.add_argument("--skip-plots", action="store_true")
    command.add_argument("--overwrite", action="store_true")
    return command


def _pair_records(pairs):
    return [
        {"design": design, "database": database}
        for design, database in pairs
    ]


def main() -> None:
    args = parser().parse_args()
    if not 0 < args.fdr <= 1:
        raise ValueError("--fdr must lie in (0, 1]")
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive")
    root = Path(args.project_root).expanduser().resolve()
    first = load_profile(resolve_path(root, args.left_profile))
    second = load_profile(resolve_path(root, args.right_profile))
    compatibility = validate_compatibility(first, second).iloc[0]
    pair_inventory = enrichment_pair_inventory(first, second)
    inputs = list(dict.fromkeys([*first.input_paths, *second.input_paths]))
    target = resolve_path(root, args.output_dir)
    publisher = StagedOutputDirectory(
        target,
        args.overwrite,
        "dataset_robustness_comparison",
        inputs=inputs,
        protected_roots=(root,),
        protected_trees=project_protected_trees(root),
    )
    with publisher as output:
        output_files = run_dataset_comparison(
            first,
            second,
            output,
            project_root=root,
            fdr_threshold=args.fdr,
            make_plots=not args.skip_plots,
            dpi=args.dpi,
        )
        manifest = make_manifest(
            "dataset_robustness_comparison",
            root,
            inputs,
            {
                "comparison_type": "same_modification_dataset_robustness",
                "dataset_a": first.dataset_id,
                "dataset_b": second.dataset_id,
                "analysis_id_a": first.analysis_id,
                "analysis_id_b": second.analysis_id,
                "modification": first.modification,
                "genome_build": compatibility["genome_build"],
                "annotation_release": compatibility["annotation_release"],
                "fasta_identity_status": compatibility["fasta_identity_status"],
                "gtf_identity_status": compatibility["gtf_identity_status"],
                "enrichment_pairs_compared": _pair_records(
                    pair_inventory["common"]
                ),
                "enrichment_pairs_a_only": _pair_records(
                    pair_inventory["dataset_a_only"]
                ),
                "enrichment_pairs_b_only": _pair_records(
                    pair_inventory["dataset_b_only"]
                ),
                "fdr_threshold": args.fdr,
                "plots_created": not args.skip_plots,
                "independent_clip_replication": False,
                "evidence_quality_policy_id": EVIDENCE_QUALITY_POLICY_ID,
                "primary_effect_comparison_scope": (
                    "shared finite effects that are inference eligible in both datasets"
                ),
                "raw_sensitivity_scope": "all shared finite effects",
                "quality_count_source": "enrichment_concordance_summary.tsv",
                "interpretation": (
                    "Robustness to dataset/callset variation; no raw replicate "
                    "pooling and no independent CLIP replication."
                ),
            },
            output_files,
        )
        remap_manifest_output_paths(manifest, output_files, output, target, root)
        write_manifest(manifest, output / "run_manifest.json")
        publisher.publish()
    print(
        f"Wrote dataset robustness comparison to {target}. "
        "This is not independent CLIP replication."
    )


if __name__ == "__main__":
    main()
