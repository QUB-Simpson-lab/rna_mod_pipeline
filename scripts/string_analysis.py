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
    validate_output_location,
)
from rna_mod_pipeline.evidence_quality import EVIDENCE_QUALITY_POLICY_ID
from rna_mod_pipeline.provenance import make_manifest, portable_path, write_manifest
from rna_mod_pipeline.schemas import validate_result_identity
from rna_mod_pipeline.stringdb import (
    CachedStringClient,
    consensus_selection_audit,
    run_string_groups,
)
from rna_mod_pipeline.stringdb.plots import plot_enrichment, plot_network


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description="Run STRING on consensus RBP groups")
    command.add_argument(
        "--modification", required=True, choices=("m6a", "m5c", "pseu", "m6a_rep2")
    )
    command.add_argument(
        "--design",
        required=True,
        choices=("loose", "transcript-region"),
    )
    command.add_argument("--context", default="all", choices=("all", "5UTR", "CDS", "3UTR"))
    command.add_argument("--cross-database", required=True)
    command.add_argument("--directions", nargs="+", choices=("enriched", "depleted"), default=["enriched", "depleted"])
    command.add_argument("--top-n", type=int, default=10)
    command.add_argument("--fdr", type=float, default=0.05)
    command.add_argument("--cache-dir", required=True)
    command.add_argument("--offline", action="store_true")
    command.add_argument("--retries", type=int, default=3)
    command.add_argument("--api-delay", type=float, default=1.0)
    command.add_argument("--project-root", default=str(Path(__file__).resolve().parents[2]))
    command.add_argument("--output-dir", required=True)
    command.add_argument("--overwrite", action="store_true")
    return command


def _set_metadata(
    table: pd.DataFrame,
    values: tuple[tuple[str, str], ...],
) -> pd.DataFrame:
    result = table.copy()
    for name, value in values:
        if name in result:
            observed = set(result[name].dropna().astype(str))
            if observed and observed != {value}:
                raise ValueError(
                    f"Cross-database {name} metadata {sorted(observed)} "
                    f"does not match requested value {value!r}"
                )
        result[name] = value
    return result


def main() -> None:
    args = parser().parse_args()
    root = Path(args.project_root).expanduser().resolve()
    if args.top_n <= 0:
        raise ValueError("--top-n must be positive")
    if not 0 < args.fdr <= 1:
        raise ValueError("--fdr must lie in (0, 1]")
    if args.retries <= 0 or args.api_delay < 0:
        raise ValueError("--retries must be positive and --api-delay non-negative")
    if len(args.directions) != len(set(args.directions)):
        raise ValueError("--directions contains duplicates")
    design = args.design
    cross_path = resolve_path(root, args.cross_database)
    cross = pd.read_csv(cross_path, sep="\t")
    identity = validate_result_identity(
        cross,
        expected_modification=args.modification,
        expected_design=design,
        context=args.context,
        table_name="cross-database table",
    )
    selection_audit = consensus_selection_audit(
        cross,
        top_n=args.top_n,
        directions=tuple(args.directions),
        fdr_threshold=args.fdr,
    )
    selection = (
        selection_audit[
            selection_audit["selected_for_primary_string"]
        ].copy()
        if not selection_audit.empty
        else selection_audit.copy()
    )
    if not selection.empty:
        selection.insert(1, "rank", selection["primary_rank"].astype(int))
        selection = selection.reset_index(drop=True)
    metadata = (
        ("modification", args.modification),
        ("design", design),
        ("context_region", args.context),
    )
    selection = _set_metadata(selection, metadata)
    selection_audit = _set_metadata(selection_audit, metadata)
    cache_dir = resolve_path(root, args.cache_dir)
    if args.offline:
        if not cache_dir.is_dir():
            raise FileNotFoundError(
                f"Offline STRING cache directory not found: {cache_dir}"
            )
    else:
        validate_output_location(
            cache_dir,
            (cross_path,),
            protected_roots=(root,),
            protected_trees=project_protected_trees(root),
        )
    target = resolve_path(root, args.output_dir)
    if (
        target == cache_dir
        or cache_dir in target.parents
        or target in cache_dir.parents
    ):
        raise ValueError(
            "STRING cache and output directories must be separate, "
            "non-nested paths"
        )
    publisher = StagedOutputDirectory(
        target,
        args.overwrite,
        "string_analysis",
        inputs=(cross_path, cache_dir),
        protected_roots=(root,),
        protected_trees=project_protected_trees(root),
    )
    with publisher as output:
        selection_path = output / "rbp_selection.tsv"
        selection_audit_path = output / "rbp_selection_candidates.tsv"
        selection.to_csv(selection_path, sep="\t", index=False)
        selection_audit.to_csv(selection_audit_path, sep="\t", index=False)
        client = CachedStringClient(
            cache_dir,
            offline=args.offline,
            retries=args.retries,
            delay_seconds=args.api_delay,
        )
        output_files, summary = run_string_groups(
            selection,
            client,
            output,
            requested_directions=tuple(args.directions),
            fdr_threshold=args.fdr,
        )
        summary = _set_metadata(summary, metadata)
        summary_path = output / "group_summary.tsv"
        summary.to_csv(summary_path, sep="\t", index=False)
        output_files.extend(
            [selection_path, selection_audit_path, summary_path]
        )
        for direction in args.directions:
            resolved_path = output / f"string_ids_{direction}.tsv"
            resolved = pd.read_csv(resolved_path, sep="\t")
            resolved_nodes = list(
                resolved.loc[
                    resolved["resolved"].astype(str).str.lower().isin(
                        {"true", "1", "yes"}
                    ),
                    "preferred_name",
                ].astype(str)
            )
            for confidence in ("medium", "high"):
                path = output / f"interactions_{direction}_{confidence}.tsv"
                plot_path = output / f"network_{direction}_{confidence}.png"
                if path.is_file():
                    plot_network(
                        pd.read_csv(path, sep="\t"),
                        plot_path,
                        f"STRING {direction} network ({confidence} confidence)",
                        resolved_nodes,
                    )
                    output_files.append(plot_path)
            path = output / f"functional_enrichment_{direction}.tsv"
            plot_path = output / f"functional_enrichment_{direction}.png"
            if path.is_file():
                plot_enrichment(
                    pd.read_csv(path, sep="\t"),
                    plot_path,
                    f"STRING functional enrichment: {direction} RBPs",
                )
                output_files.append(plot_path)
        requests = pd.DataFrame(
            client.request_log,
            columns=[
                "endpoint",
                "response_format",
                "cache_key",
                "response_source",
                "response_bytes",
                "response_sha256",
                "cache_file",
            ],
        )
        request_path = output / "string_api_requests.tsv"
        requests.to_csv(request_path, sep="\t", index=False)
        output_files.append(request_path)
        cache_inputs = [
            cache_dir / name for name in requests["cache_file"] if str(name)
        ]
        manifest = make_manifest(
            "string_analysis",
            root,
            [cross_path, *cache_inputs],
            {
                "modification": args.modification,
                "design": design,
                "context": args.context,
                "directions": args.directions,
                "top_n": args.top_n,
                "selection_rule": (
                    "ENCORI and POSTAR3 raw FDR below threshold with concordant "
                    "direction; primary selection additionally requires both "
                    "estimates to be inference eligible and ranks by arithmetic "
                    "mean OR. All raw candidates remain in rbp_selection_candidates.tsv."
                ),
                "evidence_quality_policy_id": EVIDENCE_QUALITY_POLICY_ID,
                "n_raw_consensus_candidates": int(len(selection_audit)),
                "n_primary_selected": int(len(selection)),
                "n_quality_excluded_candidates": int(
                    selection_audit.get(
                        "selection_status", pd.Series(dtype=str)
                    )
                    .eq("sensitivity_only_quality_excluded")
                    .sum()
                ),
                "confidence_thresholds": {"medium": 400, "high": 700},
                "species": 9606,
                "offline": args.offline,
                "input_identity_status": identity,
                "string_release": "not pinned by the STRING API",
                "interpretation_caveat": (
                    "STRING connectivity does not demonstrate interaction in MDA-MB-231."
                ),
            },
            output_files,
        )
        for record, path in zip(manifest["outputs"], output_files):
            final_path = target / path.relative_to(output)
            record["path"] = portable_path(final_path, root)
        write_manifest(manifest, output / "run_manifest.json")
        publisher.publish()
    print(f"Wrote STRING analysis to {target}")


if __name__ == "__main__":
    main()
