from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .enrichment import COLOURS, PREFIXES, _empty_plot
from .schema import prepare_results


REGION_PREFIX = {"5UTR": "5utr", "CDS": "cds", "3UTR": "3utr"}


def _regional_bar(
    data: pd.DataFrame,
    path: Path,
    title: str,
    direction: str,
    top_n: int,
    dpi: int,
) -> tuple[int, set[str]]:
    if direction == "enriched":
        selected = data[
            data["plot_significant"] & data["mh_or"].gt(1)
        ].nlargest(top_n, "log2_mh_or")
    else:
        selected = data[
            data["plot_significant"] & data["mh_or"].lt(1)
        ].nsmallest(top_n, "log2_mh_or")
    selected = selected.sort_values("log2_mh_or", kind="mergesort")
    if selected.empty:
        _empty_plot(path, title, "No FDR-significant estimable RBPs", dpi)
        return 0, set()
    fig, axis = plt.subplots(figsize=(9, max(4.5, len(selected) * 0.3)))
    axis.barh(
        np.arange(len(selected)),
        selected["log2_mh_or"],
        color=COLOURS[direction],
    )
    axis.set_yticks(np.arange(len(selected)), selected["RBP"], fontsize=8)
    axis.axvline(0, color="#333333", linewidth=0.8)
    axis.set_xlabel("log₂ Mantel–Haenszel odds ratio")
    axis.set_title(title, fontweight="bold")
    axis.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(selected), set(selected["RBP"])


def generate_regional_plots(
    regional_results: pd.DataFrame,
    output_dir: Path,
    analysis_id: str,
    fdr_threshold: float = 0.05,
    top_n: int = 30,
    dpi: int = 180,
) -> pd.DataFrame:
    output_dir.mkdir(parents=True, exist_ok=True)
    data = prepare_results(
        regional_results, "fdr_within_database_region", fdr_threshold
    )
    data["selected_for_regional_top_plot"] = False
    data["regional_rank"] = np.nan
    manifest_rows = []
    for (database, region), indexes in data.groupby(
        ["database", "region_filter"], sort=True
    ).groups.items():
        group = data.loc[list(indexes)]
        prefix = PREFIXES.get(database, database.lower())
        region_prefix = REGION_PREFIX.get(str(region), str(region).lower())
        for direction in ("enriched", "depleted"):
            path = output_dir / (
                f"{prefix}_{region_prefix}_top_{direction}.png"
            )
            count, rbps = _regional_bar(
                group,
                path,
                f"{analysis_id}: {database} {region} {direction} RBPs",
                direction,
                top_n,
                dpi,
            )
            selected = (
                data.index.isin(indexes)
                & data["RBP"].isin(rbps)
                & data["plot_direction"].eq(direction)
            )
            ranked = data.loc[selected].sort_values(
                "log2_mh_or",
                ascending=direction == "depleted",
                kind="mergesort",
            )
            data.loc[ranked.index, "selected_for_regional_top_plot"] = True
            data.loc[ranked.index, "regional_rank"] = np.arange(
                1, len(ranked) + 1
            )
            manifest_rows.append(
                {
                    "database": database,
                    "region_filter": region,
                    "direction": direction,
                    "path": path.name,
                    "n_selected": count,
                }
            )
    data.to_csv(
        output_dir / "regional_top_plot_source_data.tsv", sep="\t", index=False
    )
    manifest = pd.DataFrame(manifest_rows)
    manifest["analysis_id"] = analysis_id
    manifest["fdr_threshold"] = fdr_threshold
    manifest["top_n"] = top_n
    manifest["confidence_intervals_drawn"] = False
    manifest.to_csv(
        output_dir / "regional_plot_manifest.tsv", sep="\t", index=False
    )
    (output_dir / "regional_plot_run.json").write_text(
        json.dumps(
            {
                "analysis_id": analysis_id,
                "fdr_threshold": fdr_threshold,
                "top_n": top_n,
                "dpi": dpi,
                "confidence_intervals_drawn": False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    (output_dir / "README.md").write_text(
        "# Regional top-RBP plots\n\n"
        "Each panel ranks FDR-significant transcript-region-stratified "
        "Mantel–Haenszel effects. Confidence intervals are not drawn.\n"
    )
    return manifest
