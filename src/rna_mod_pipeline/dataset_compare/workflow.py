from __future__ import annotations

from pathlib import Path

import pandas as pd

from .enrichment import compare_enrichment
from .profiles import (
    DatasetProfile,
    enrichment_pair_inventory,
    validate_compatibility,
)
from .sites import compare_site_tables


INTERPRETATION = """Same-modification dataset robustness comparison

This workflow compares independently processed biological dataset/callset
workspaces. It does not pool raw replicates and does not estimate a pooled
biological effect.

Any enrichment comparison reuses the RBP-binding resources declared by the
two analyses. Agreement therefore measures robustness to the modification
dataset/callset and is not independent replication of CLIP experiments.
"""


def _summary(
    compatibility: pd.DataFrame,
    site_summary: pd.DataFrame,
    enrichment_summary: pd.DataFrame,
) -> pd.DataFrame:
    rows = [
        {
            "section": "scope",
            "group": "all",
            "metric": "comparison_type",
            "value": "same_modification_dataset_robustness",
        },
        {
            "section": "scope",
            "group": "all",
            "metric": "independent_clip_replication",
            "value": "false",
        },
    ]
    identity = compatibility.iloc[0]
    for key in (
        "dataset_a",
        "dataset_b",
        "modification",
        "genome_build",
        "annotation_release",
        "fasta_identity_status",
        "gtf_identity_status",
        "n_enrichment_pairs_a",
        "n_enrichment_pairs_b",
        "n_enrichment_pairs_compared",
        "enrichment_pairs_a",
        "enrichment_pairs_b",
        "enrichment_pairs_compared",
        "enrichment_pairs_a_only",
        "enrichment_pairs_b_only",
    ):
        rows.append(
            {
                "section": "identity",
                "group": "all",
                "metric": key,
                "value": identity[key],
            }
        )
    for record in site_summary.to_dict("records"):
        for key, value in record.items():
            if key != "stage":
                rows.append(
                    {
                        "section": "sites",
                        "group": record["stage"],
                        "metric": key,
                        "value": value,
                    }
                )
    for record in enrichment_summary.to_dict("records"):
        group = f"{record['design']}/{record['database']}"
        for key, value in record.items():
            if key not in {"design", "database", "comparison_interpretation"}:
                rows.append(
                    {
                        "section": "enrichment",
                        "group": group,
                        "metric": key,
                        "value": value,
                    }
                )
    return pd.DataFrame(rows)


def _pair_inventory(
    first: DatasetProfile,
    second: DatasetProfile,
) -> pd.DataFrame:
    inventory = enrichment_pair_inventory(first, second)
    all_pairs = sorted(set(inventory["dataset_a"]) | set(inventory["dataset_b"]))
    rows = []
    for design, database in all_pairs:
        in_a = (design, database) in first.enrichment
        in_b = (design, database) in second.enrichment
        rows.append(
            {
                "design": design,
                "database": database,
                "declared_by_dataset_a": in_a,
                "declared_by_dataset_b": in_b,
                "comparison_status": (
                    "compared"
                    if in_a and in_b
                    else "dataset_a_only"
                    if in_a
                    else "dataset_b_only"
                ),
                "absence_interpretation": (
                    "undeclared/not compared; not a null enrichment result"
                    if not (in_a and in_b)
                    else "both inputs declared"
                ),
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "design",
            "database",
            "declared_by_dataset_a",
            "declared_by_dataset_b",
            "comparison_status",
            "absence_interpretation",
        ],
    )


def run_dataset_comparison(
    first: DatasetProfile,
    second: DatasetProfile,
    output_dir: str | Path,
    *,
    project_root: str | Path,
    fdr_threshold: float = 0.05,
    make_plots: bool = True,
    dpi: int = 180,
) -> list[Path]:
    if not 0 < fdr_threshold <= 1:
        raise ValueError("fdr_threshold must lie in (0, 1]")
    if dpi <= 0:
        raise ValueError("dpi must be positive")
    output = Path(output_dir)
    compatibility = validate_compatibility(first, second)
    membership, site_summary = compare_site_tables(first, second)
    enrichment, enrichment_summary = compare_enrichment(
        first, second, project_root, fdr_threshold
    )

    outputs = []
    tables = (
        ("dataset_compatibility.tsv", compatibility),
        ("enrichment_pair_inventory.tsv", _pair_inventory(first, second)),
        ("site_membership.tsv.gz", membership),
        ("site_overlap_summary.tsv", site_summary),
        ("enrichment_concordance.tsv", enrichment),
        ("enrichment_concordance_summary.tsv", enrichment_summary),
        (
            "comparison_summary.tsv",
            _summary(compatibility, site_summary, enrichment_summary),
        ),
    )
    for name, table in tables:
        path = output / name
        table.to_csv(path, sep="\t", index=False)
        outputs.append(path)
    interpretation = output / "INTERPRETATION.txt"
    interpretation.write_text(INTERPRETATION)
    outputs.append(interpretation)
    if make_plots:
        from .plots import create_plots

        manifest = create_plots(
            membership,
            site_summary,
            enrichment,
            output,
            first.dataset_id,
            second.dataset_id,
            dpi,
        )
        outputs.append(output / "plot_manifest.tsv")
        outputs.extend(
            output / path for path in manifest["path"].astype(str)
        )
    return outputs
