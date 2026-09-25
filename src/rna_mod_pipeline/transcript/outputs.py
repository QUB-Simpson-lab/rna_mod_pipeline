from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from rna_mod_pipeline.provenance import file_record, portable_path, software_record

from .config import REGIONS
from .opportunities import OpportunityBuild, OpportunityDesign


def save_opportunity_design(build: OpportunityBuild, output_dir: Path) -> None:
    design = build.design
    chromosomes = np.asarray(sorted(set(map(str, design.chroms))), dtype="U")
    chromosome_to_code = {
        chrom: index for index, chrom in enumerate(chromosomes)
    }
    chromosome_codes = np.asarray(
        [chromosome_to_code[str(value)] for value in design.chroms],
        dtype=np.int16,
    )
    strand_codes = np.where(design.strands == "+", 1, -1).astype(np.int8)
    np.savez_compressed(
        output_dir / "opportunity_design.npz",
        positions=design.positions,
        chrom_codes=chromosome_codes,
        chrom_vocabulary=chromosomes,
        strand_codes=strand_codes,
        coverages=design.coverages,
        fractions=design.fractions,
        is_case=design.is_case,
        stratum_codes=design.stratum_codes,
        case_site_ids=design.case_site_ids,
    )
    build.case_assignments.to_csv(
        output_dir / "case_site_assignments.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    build.strata.to_csv(
        output_dir / "transcript_region_strata.tsv", sep="\t", index=False
    )
    build.flow.to_csv(output_dir / "opportunity_flow.tsv", sep="\t", index=False)
    _coverage_balance(design).to_csv(
        output_dir / "coverage_balance.tsv", sep="\t", index=False
    )


def _coverage_balance(design: OpportunityDesign) -> pd.DataFrame:
    observation_regions = design.stratum_regions[design.stratum_codes]
    rows = []
    for region in ("all", *REGIONS):
        subset = (
            np.ones(design.n_observations, dtype=bool)
            if region == "all"
            else observation_regions == region
        )
        case = design.coverages[subset & design.is_case].astype(float)
        control = design.coverages[subset & ~design.is_case].astype(float)
        if not len(case) or not len(control):
            rows.append(
                {
                    "region": region,
                    "n_case_sites": len(case),
                    "n_control_sites": len(control),
                    "case_coverage_mean": math.nan,
                    "control_coverage_mean": math.nan,
                    "case_coverage_median": math.nan,
                    "control_coverage_median": math.nan,
                    "case_coverage_q25": math.nan,
                    "case_coverage_q75": math.nan,
                    "control_coverage_q25": math.nan,
                    "control_coverage_q75": math.nan,
                    "log1p_coverage_standardised_mean_difference": math.nan,
                }
            )
            continue
        case_log, control_log = np.log1p(case), np.log1p(control)
        denominator = len(case_log) + len(control_log) - 2
        pooled = (
            (
                (len(case_log) - 1) * np.var(case_log, ddof=1)
                + (len(control_log) - 1) * np.var(control_log, ddof=1)
            )
            / denominator
            if len(case_log) > 1 and len(control_log) > 1
            else math.nan
        )
        smd = (
            (float(case_log.mean()) - float(control_log.mean()))
            / math.sqrt(pooled)
            if np.isfinite(pooled) and pooled > 0
            else math.nan
        )
        rows.append(
            {
                "region": region,
                "n_case_sites": len(case),
                "n_control_sites": len(control),
                "case_coverage_mean": float(case.mean()),
                "control_coverage_mean": float(control.mean()),
                "case_coverage_median": float(np.median(case)),
                "control_coverage_median": float(np.median(control)),
                "case_coverage_q25": float(np.quantile(case, 0.25)),
                "case_coverage_q75": float(np.quantile(case, 0.75)),
                "control_coverage_q25": float(np.quantile(control, 0.25)),
                "control_coverage_q75": float(np.quantile(control, 0.75)),
                "log1p_coverage_standardised_mean_difference": smd,
            }
        )
    return pd.DataFrame(rows)


def save_case_binding_matrix(
    output_dir: Path,
    database_key: str,
    rbps: Sequence[str],
    packed_rows: Sequence[np.ndarray],
    case_site_ids: np.ndarray,
    n_case_sites: int,
) -> None:
    packed = (
        np.vstack(packed_rows).astype(np.uint8, copy=False)
        if packed_rows
        else np.zeros((0, math.ceil(n_case_sites / 8)), dtype=np.uint8)
    )
    np.savez_compressed(
        output_dir / f"case_binding_matrix_{database_key}.npz",
        packed_matrix=packed,
        rbps=np.asarray(rbps, dtype="U"),
        case_site_ids=np.asarray(case_site_ids, dtype=np.int64),
        n_case_sites=np.asarray([n_case_sites], dtype=np.int64),
        bitorder=np.asarray(["little"], dtype="U"),
    )
    pd.DataFrame(
        {"matrix_row": np.arange(len(rbps), dtype=int), "RBP": rbps}
    ).to_csv(
        output_dir / f"case_binding_matrix_{database_key}_rows.tsv",
        sep="\t",
        index=False,
    )


def write_run_readme(
    path: Path,
    analysis_id: str,
    modification: str,
    databases: Sequence[str],
    n_cases: int,
    n_controls: int,
    partial: bool,
    database_subset: bool = False,
    custom_parameters: bool = False,
    resource_filtered: bool = False,
) -> None:
    path.write_text(
        "\n".join(
            [
                f"# Transcript-region RBP analysis: {analysis_id}",
                "",
                "Cases and low-modification controls are compared within the "
                "same primary transcript and 5′UTR/CDS/3′UTR region.",
                "",
                f"- Biological modification: `{modification}`",
                f"- Databases: {', '.join(databases)}",
                f"- Contributing cases: {n_cases:,}",
                f"- Contributing controls: {n_controls:,}",
                f"- Partial/test run: `{str(partial).lower()}`",
                f"- Independently selected database subset: "
                f"`{str(database_subset).lower()}`",
                f"- Non-default analysis parameters: "
                f"`{str(custom_parameters).lower()}`",
                f"- RBP source records filtered: "
                f"`{str(resource_filtered).lower()}`",
                "",
                "`stratified_enrichment_results.tsv` is the primary table. "
                "Its Mantel–Haenszel odds ratio is accompanied by a "
                "gene-cluster-robust p-value and database-wise BH FDR.",
                "",
                "The crude Fisher result is descriptive. Regional results "
                "are in `region_specific_results.tsv`. Plots do not display "
                "confidence intervals.",
            ]
        )
        + "\n"
    )


def write_run_manifest(
    path: Path,
    stage_dir: Path,
    destination: Path,
    project_root: Path,
    input_files: Iterable[Path],
    parameters: dict,
) -> None:
    inputs = [file_record(item, project_root) for item in dict.fromkeys(input_files)]
    outputs = []
    destination_relative = portable_path(destination, project_root)
    for output in sorted(stage_dir.rglob("*")):
        if not output.is_file() or output == path:
            continue
        relative = output.relative_to(stage_dir)
        record = file_record(output, stage_dir)
        record["path"] = (
            f"{destination_relative}/{relative.as_posix()}"
            if destination_relative != "."
            else relative.as_posix()
        )
        outputs.append(record)
    manifest = {
        "schema_version": 1,
        "workflow": "transcript-region",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": ".",
        "parameters": parameters,
        "inputs": inputs,
        "outputs": outputs,
        "software": software_record(project_root),
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
