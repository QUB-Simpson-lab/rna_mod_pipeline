from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCOPE_NOTE = (
    "Dataset/callset robustness comparison; shared RBP resources are not "
    "independent CLIP replication."
)


def _save(figure: plt.Figure, path: Path, dpi: int) -> None:
    figure.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_site_overlap(summary: pd.DataFrame, path: Path, dpi: int) -> None:
    labels = summary["stage"].tolist()
    x = np.arange(len(labels))
    a_only = summary["n_sites_a_only"].to_numpy()
    shared = summary["n_shared_sites"].to_numpy()
    b_only = summary["n_sites_b_only"].to_numpy()
    figure, axis = plt.subplots(figsize=(8, 5))
    axis.bar(x, a_only, label="Dataset A only", color="#457B9D")
    axis.bar(x, shared, bottom=a_only, label="Shared", color="#65A765")
    axis.bar(x, b_only, bottom=a_only + shared, label="Dataset B only", color="#E76F51")
    axis.set_xticks(x, labels)
    axis.set_ylabel("Genomic sites in union")
    axis.set_title("Exact site-set overlap")
    axis.legend(frameon=False)
    axis.text(0.01, -0.16, SCOPE_NOTE, transform=axis.transAxes, fontsize=8)
    figure.tight_layout()
    _save(figure, path, dpi)


def _scatter_or_hexbin(
    axis: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    xlabel: str,
    ylabel: str,
) -> None:
    if len(x) > 2_000:
        axis.hexbin(x, y, gridsize=45, mincnt=1, cmap="viridis")
    else:
        axis.scatter(x, y, s=12, alpha=0.5, color="#457B9D", edgecolors="none")
    if len(x):
        low = min(float(np.min(x)), float(np.min(y)))
        high = max(float(np.max(x)), float(np.max(y)))
        axis.plot([low, high], [low, high], "--", color="#555555", linewidth=0.8)
    axis.set_xlabel(xlabel)
    axis.set_ylabel(ylabel)
    axis.grid(color="#E6E6E6", linewidth=0.5)


def plot_shared_site_concordance(
    membership: pd.DataFrame,
    path: Path,
    first_label: str,
    second_label: str,
    dpi: int,
) -> None:
    shared = membership[
        membership["stage"].eq("filtered") & membership["membership"].eq("shared")
    ]
    figure, axes = plt.subplots(1, 2, figsize=(11, 5))
    _scatter_or_hexbin(
        axes[0],
        shared["fraction_modified_a"].to_numpy(dtype=float),
        shared["fraction_modified_b"].to_numpy(dtype=float),
        f"{first_label} modification fraction",
        f"{second_label} modification fraction",
    )
    _scatter_or_hexbin(
        axes[1],
        np.log1p(shared["Nvalid_cov_a"].to_numpy(dtype=float)),
        np.log1p(shared["Nvalid_cov_b"].to_numpy(dtype=float)),
        f"{first_label} log(coverage + 1)",
        f"{second_label} log(coverage + 1)",
    )
    figure.suptitle(f"Shared filtered-site concordance (n={len(shared):,})")
    figure.text(0.01, 0.01, SCOPE_NOTE, fontsize=8)
    figure.tight_layout(rect=(0, 0.04, 1, 0.95))
    _save(figure, path, dpi)


def plot_enrichment_concordance(
    details: pd.DataFrame,
    path: Path,
    first_label: str,
    second_label: str,
    dpi: int,
) -> None:
    pairs = list(details[["design", "database"]].drop_duplicates().itertuples(index=False))
    columns = min(3, max(1, len(pairs)))
    rows = math.ceil(len(pairs) / columns)
    figure, axes = plt.subplots(
        rows, columns, figsize=(5 * columns, 4.2 * rows), squeeze=False
    )
    for axis, pair in zip(axes.flat, pairs):
        subset = details[
            details["design"].eq(pair.design)
            & details["database"].eq(pair.database)
        ].dropna(subset=["log2_odds_ratio_a", "log2_odds_ratio_b"])
        x = subset["log2_odds_ratio_a"].to_numpy(dtype=float)
        y = subset["log2_odds_ratio_b"].to_numpy(dtype=float)
        primary = subset[
            subset["eligible_for_primary_effect_comparison"]
        ]
        axis.scatter(
            x,
            y,
            s=12,
            alpha=0.35,
            color="#B8B8B8",
            edgecolors="none",
            label="Raw all-finite sensitivity",
        )
        axis.scatter(
            primary["log2_odds_ratio_a"],
            primary["log2_odds_ratio_b"],
            s=14,
            alpha=0.65,
            color="#457B9D",
            edgecolors="none",
            label="Primary quality-qualified",
        )
        if len(subset):
            low = min(float(np.min(x)), float(np.min(y)))
            high = max(float(np.max(x)), float(np.max(y)))
            axis.plot([low, high], [low, high], "--", color="#555555", linewidth=0.8)
        primary_rho = (
            primary["log2_odds_ratio_a"].corr(
                primary["log2_odds_ratio_b"], method="spearman"
            )
            if len(primary) >= 3
            else np.nan
        )
        raw_rho = (
            subset["log2_odds_ratio_a"].corr(
                subset["log2_odds_ratio_b"], method="spearman"
            )
            if len(subset) >= 3
            else np.nan
        )
        axis.text(
            0.03,
            0.97,
            f"Primary eligible n={len(primary):,}; ρ={primary_rho:.3f}\n"
            f"Raw sensitivity n={len(subset):,}; ρ={raw_rho:.3f}",
            transform=axis.transAxes,
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "#DDDDDD", "alpha": 0.9},
        )
        axis.set_xlabel(f"{first_label} log₂ OR")
        axis.set_ylabel(f"{second_label} log₂ OR")
        axis.grid(color="#E6E6E6", linewidth=0.5)
        axis.set_title(f"{pair.design} | {pair.database}")
        axis.legend(frameon=False, fontsize=7, loc="lower right")
    for axis in axes.flat[len(pairs) :]:
        axis.axis("off")
    figure.suptitle(
        "RBP enrichment-effect concordance "
        "(primary quality-qualified; raw all-finite sensitivity in grey)"
    )
    figure.text(0.01, 0.01, SCOPE_NOTE, fontsize=8)
    figure.tight_layout(rect=(0, 0.035, 1, 0.96))
    _save(figure, path, dpi)


def create_plots(
    membership: pd.DataFrame,
    site_summary: pd.DataFrame,
    enrichment: pd.DataFrame,
    output_dir: Path,
    first_label: str,
    second_label: str,
    dpi: int,
) -> pd.DataFrame:
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    site_overlap = plot_dir / "site_set_overlap.png"
    plot_site_overlap(site_summary, site_overlap, dpi)
    rows.append(
        {
            "plot_type": "site_set_overlap",
            "path": site_overlap.relative_to(output_dir).as_posix(),
            "source_table": "site_overlap_summary.tsv",
            "n_records": len(site_summary),
        }
    )
    shared = plot_dir / "shared_site_concordance.png"
    plot_shared_site_concordance(
        membership, shared, first_label, second_label, dpi
    )
    rows.append(
        {
            "plot_type": "shared_site_concordance",
            "path": shared.relative_to(output_dir).as_posix(),
            "source_table": "site_membership.tsv.gz",
            "n_records": int(
                (
                    membership["stage"].eq("filtered")
                    & membership["membership"].eq("shared")
                ).sum()
            ),
        }
    )
    if not enrichment.empty:
        effect = plot_dir / "enrichment_effect_concordance.png"
        plot_enrichment_concordance(
            enrichment, effect, first_label, second_label, dpi
        )
        rows.append(
            {
                "plot_type": "enrichment_effect_concordance",
                "path": effect.relative_to(output_dir).as_posix(),
                "source_table": "enrichment_concordance.tsv",
                "n_records": len(enrichment),
                "primary_correlation_scope": (
                    "both effects finite and inference eligible"
                ),
                "raw_sensitivity_scope": "all shared finite effects",
            }
        )
    manifest = pd.DataFrame(rows)
    manifest["comparison_scope"] = "same_modification_dataset_robustness"
    manifest["independent_clip_replication"] = False
    manifest.to_csv(output_dir / "plot_manifest.tsv", sep="\t", index=False)
    return manifest
