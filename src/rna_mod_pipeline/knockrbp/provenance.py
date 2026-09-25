from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from ..evidence_quality import EVIDENCE_QUALITY_POLICY_ID
from ..schemas import validate_result_identity


def validate_design_inputs(
    design: str,
    *,
    encori_annotated: str | None,
    postar3_annotated: str | None,
    encori_enrichment: str | None,
    postar3_enrichment: str | None,
    transcript_run_dir: str | None,
) -> None:
    loose_inputs = {
        "--encori-annotated": encori_annotated,
        "--postar3-annotated": postar3_annotated,
        "--encori-enrichment": encori_enrichment,
        "--postar3-enrichment": postar3_enrichment,
    }
    if design == "loose" and transcript_run_dir:
        raise ValueError(
            "--transcript-run-dir applies only to --design transcript-region"
        )
    if design == "transcript-region":
        supplied = [option for option, value in loose_inputs.items() if value]
        if supplied:
            raise ValueError(
                "Loose-design inputs cannot be used with "
                f"--design transcript-region: {', '.join(supplied)}"
            )


def annotated_identity(
    path: Path,
    database: str,
    modification: str,
) -> dict[str, str]:
    header = pd.read_csv(path, sep="\t", nrows=0)
    labels = [
        column
        for column in (
            "modification",
            "analysis_modification",
            "database",
            "analysis_database",
            "design",
            "analysis_design",
        )
        if column in header
    ]
    frame = pd.read_csv(path, sep="\t", usecols=labels) if labels else header
    return validate_result_identity(
        frame,
        expected_modification=modification,
        expected_database=database,
        expected_design="loose",
        context="all",
        table_name=f"{database} annotated site table",
    )


def enrichment_identity(
    path: Path,
    database: str,
    modification: str,
) -> dict[str, str]:
    frame = pd.read_csv(path, sep="\t")
    return validate_result_identity(
        frame,
        expected_modification=modification,
        expected_database=database,
        expected_design="loose",
        context="all",
        table_name=f"{database} enrichment resource table",
    )


def transcript_run_identity(
    run: Path,
    modification: str,
) -> dict[str, str]:
    manifest_path = run / "run_manifest.json"
    if not manifest_path.is_file():
        return {
            "modification": "unverified_no_run_manifest",
            "database": "not_applicable",
            "design": "unverified_no_run_manifest",
            "context": "internal_region_filter",
        }
    manifest = json.loads(manifest_path.read_text())
    parameters = manifest.get("parameters", {})
    analysis_id = manifest.get("analysis_id") or parameters.get("analysis_id")
    biological_modification = (
        manifest.get("modification") or parameters.get("modification")
    )
    observed_design = (
        manifest.get("workflow")
        or parameters.get("design")
        or "transcript-region"
    )
    record = {
        "analysis_design": observed_design,
        "modification": biological_modification,
    }
    if analysis_id is not None:
        record["analysis_id"] = analysis_id
    frame = pd.DataFrame([record])
    status = validate_result_identity(
        frame,
        expected_modification=modification,
        expected_design="transcript-region",
        context="all",
        table_name="transcript-region run manifest",
    )
    status["context"] = "internal_region_filter"
    return status


def cross_database_schema_status(table: pd.DataFrame) -> str:
    quality_columns = {
        f"{database}_{column}"
        for database in ("ornament", "encori", "postar3")
        for column in (
            "inference_eligible",
            "evidence_quality",
            "inference_exclusion_reasons",
        )
    }
    if quality_columns.issubset(table.columns):
        return "current_quality_columns_present_and_recomputed"
    return "legacy_or_partial_quality_schema_recomputed"


def knockrbp_manifest_parameters(
    *,
    options: Mapping[str, object],
    n_datasets: int,
    network: pd.DataFrame,
    schema_status: str,
    identity_status: Mapping[str, object],
) -> dict[str, object]:
    """Build the filtering, consensus, and interpretation metadata."""
    consensus = network.get(
        "consensus_evidence_quality", pd.Series(dtype=str)
    )
    return {
        "modification": options["modification"],
        "design": options["design"],
        "context": options["context"],
        "cell_lines": options["cell_lines"],
        "dataset_ids": options["dataset_ids"],
        "log2fc_filter": f"absolute log2FC > {options['log2fc_cutoff']}",
        "padj_filter_when_available": (
            f"adjusted p-value < {options['padj_cutoff']}"
        ),
        "regulatory_network_plot_top_n": options["regulatory_plot_top_n"],
        "regulatory_network_plot_selection": (
            "affected RBPs ranked by maximum absolute retained DEG log2FC; "
            "the complete regulatory network remains in its TSV"
        ),
        "primary_background": (
            "genes with contributing modified sites in the selected context"
        ),
        "sensitivity_background": 20_000,
        "n_datasets": n_datasets,
        "n_regulatory_edges": len(network),
        "n_regulatory_edges_primary_consensus": int(
            consensus.eq("primary_eligible_consensus").sum()
        ),
        "n_regulatory_edges_sensitivity_only_consensus": int(
            consensus.eq("sensitivity_only_consensus_quality_excluded").sum()
        ),
        "evidence_quality_policy_id": EVIDENCE_QUALITY_POLICY_ID,
        "cross_database_schema_status": schema_status,
        "cross_database_consensus_recomputed": bool(options["cross_database"]),
        "cross_database_recompute_fdr_threshold": 0.05,
        "input_identity_status": dict(identity_status),
        "causal_caveat": "Knockdown-associated changes may be indirect.",
    }
