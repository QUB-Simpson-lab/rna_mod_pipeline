from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _pattern_plot(patterns: pd.DataFrame, destination: Path) -> None:
    counts = patterns.groupby(["database", "pattern"]).size().unstack(fill_value=0)
    ax = counts.plot.bar(figsize=(11, 6))
    ax.set_ylabel("Number of RBPs")
    ax.set_title("Cross-modification RBP patterns")
    ax.legend(title="Pattern", fontsize=7)
    ax.figure.tight_layout()
    ax.figure.savefig(destination, dpi=180)
    plt.close(ax.figure)


def _correlation_plot(
    correlations: pd.DataFrame,
    destination: Path,
    *,
    coefficient_column: str,
    title: str,
) -> None:
    frame = correlations[correlations["method"].eq("spearman")].copy()
    frame["pair"] = frame["modification_a"] + " vs " + frame["modification_b"]
    pivot = frame.pivot(
        index="database", columns="pair", values=coefficient_column
    )
    fig, ax = plt.subplots(figsize=(8, max(3, 0.8 * len(pivot))))
    colours = plt.get_cmap("RdBu_r").copy()
    colours.set_bad("#D9D9D9")
    values = np.ma.masked_invalid(pivot.to_numpy(dtype=float))
    image = ax.imshow(values, cmap=colours, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="Spearman ρ")
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def _effect_heatmap(
    matrix: pd.DataFrame,
    destination: Path,
    *,
    primary: bool,
) -> set[str]:
    source = matrix.copy()
    source["_finite_effect"] = np.isfinite(
        pd.to_numeric(source["log2_odds_ratio"], errors="coerce")
    )
    source["_eligible_effect"] = (
        source["_finite_effect"]
        & source["inference_eligible"].fillna(False).astype(bool)
    )
    present = source[
        source["_eligible_effect"] if primary else source["_finite_effect"]
    ].copy()
    if present.empty:
        fig, ax = plt.subplots(figsize=(9, 4))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No quality-qualified finite effects" if primary else "No finite effects",
            ha="center",
            va="center",
        )
        ax.set_title(
            "Primary cross-modification RBP effect-size overview"
            if primary
            else "Raw all-finite effect-size sensitivity overview"
        )
        fig.tight_layout()
        fig.savefig(destination, dpi=180)
        plt.close(fig)
        return set()
    source["column"] = source["database"] + " | " + source["modification"]
    present["column"] = present["database"] + " | " + present["modification"]
    ranking = present.pivot(
        index="RBP", columns="column", values="log2_odds_ratio"
    )
    rank_table = pd.DataFrame(
        {
            "variance": ranking.var(axis=1, skipna=True).fillna(-np.inf),
            "mean_absolute_effect": ranking.abs().mean(axis=1, skipna=True),
        }
    ).sort_values(
        ["variance", "mean_absolute_effect"],
        ascending=[False, False],
        kind="mergesort",
    )
    selected_rbps = list(rank_table.head(60).index)
    display = source[source["RBP"].isin(selected_rbps)].copy()
    if primary:
        display.loc[~display["_eligible_effect"], "log2_odds_ratio"] = np.nan
    pivot = display.pivot(
        index="RBP", columns="column", values="log2_odds_ratio"
    ).reindex(index=selected_rbps)
    fig, ax = plt.subplots(figsize=(10, max(8, len(pivot) * 0.18)))
    finite_values = np.abs(pivot.to_numpy(dtype=float))
    finite_values = finite_values[np.isfinite(finite_values)]
    limit = np.nanpercentile(finite_values, 95) if len(finite_values) else 0.5
    limit = max(float(limit), 0.5)
    colours = plt.get_cmap("RdBu_r").copy()
    colours.set_bad("#D9D9D9")
    image = ax.imshow(
        np.ma.masked_invalid(pivot.to_numpy(dtype=float)),
        cmap=colours,
        vmin=-limit,
        vmax=limit,
        aspect="auto",
    )
    state = source.copy()
    state_pivot = state.pivot(index="RBP", columns="column", values="state").reindex(
        index=pivot.index, columns=pivot.columns
    )
    for row_index in range(len(pivot.index)):
        for column_index in range(len(pivot.columns)):
            state_value = state_pivot.iat[row_index, column_index]
            quality_excluded = state_value == "quality_excluded"
            if primary and quality_excluded:
                label = "Q"
            elif not primary and quality_excluded and pd.notna(
                pivot.iat[row_index, column_index]
            ):
                label = "Q"
            elif pd.isna(pivot.iat[row_index, column_index]):
                label = (
                    "×"
                    if state_value == "not_covered"
                    else "B"
                    if state_value == "boundary_estimate"
                    else "NA"
                )
            else:
                continue
            ax.text(
                column_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=5,
                color="#444444",
            )
    ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)), pivot.index, fontsize=6)
    ax.set_title(
        "Primary cross-modification RBP effects (Q = quality excluded)"
        if primary
        else "Raw all-finite effect sensitivity (Q = quality excluded)"
    )
    fig.colorbar(image, ax=ax, label="log₂ odds ratio")
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)
    return set(selected_rbps)


def create_crossmod_plots(
    matrix: pd.DataFrame,
    patterns: pd.DataFrame,
    correlations: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths = [
        output / "crossmod_pattern_counts.png",
        output / "crossmod_effect_correlations.png",
        output / "crossmod_effect_correlations_raw_sensitivity.png",
        output / "crossmod_effect_heatmap.png",
        output / "crossmod_effect_heatmap_raw_sensitivity.png",
        output / "crossmod_plot_correlation_source.tsv",
        output / "crossmod_plot_effect_source.tsv",
    ]
    _pattern_plot(patterns, paths[0])
    _correlation_plot(
        correlations,
        paths[1],
        coefficient_column="inference_eligible_coefficient",
        title="Primary Spearman correlation of quality-qualified log₂ odds ratios",
    )
    _correlation_plot(
        correlations,
        paths[2],
        coefficient_column="coefficient",
        title="Raw all-finite Spearman correlation (sensitivity)",
    )
    primary_rbps = _effect_heatmap(matrix, paths[3], primary=True)
    raw_rbps = _effect_heatmap(matrix, paths[4], primary=False)
    correlation_source = correlations.copy()
    correlation_source["primary_plot_metric"] = (
        "inference_eligible_coefficient"
    )
    correlation_source["raw_sensitivity_metric"] = "coefficient"
    correlation_source.to_csv(paths[5], sep="\t", index=False)
    effect_source = matrix.copy()
    effect_source["shown_in_primary_effect_heatmap"] = effect_source["RBP"].isin(
        primary_rbps
    )
    effect_source["shown_in_raw_sensitivity_heatmap"] = effect_source["RBP"].isin(
        raw_rbps
    )
    effect_source["primary_plot_scope"] = (
        "finite_effect_and_inference_eligible"
    )
    effect_source["raw_sensitivity_scope"] = "all_finite_effects"
    effect_source.to_csv(paths[6], sep="\t", index=False)
    return [path for path in paths if path.is_file()]


def phase1_summary(metagene_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for modification, table in metagene_tables.items():
        regions = table.get("region", pd.Series(dtype=str)).astype(str)
        mapped = regions.isin(["5UTR", "CDS", "3UTR"])
        row = {
            "modification": modification,
            "n_filtered_sites": len(table),
            "n_transcript_region_mapped": int(mapped.sum()),
        }
        for region in ("5UTR", "CDS", "3UTR"):
            count = int(regions.eq(region).sum())
            row[f"n_{region}"] = count
            row[f"pct_mapped_{region}"] = (
                100 * count / int(mapped.sum()) if mapped.any() else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def plot_metagene_overlay(
    metagene_tables: dict[str, pd.DataFrame], path: str | Path
) -> None:
    destination = Path(path)
    fig, ax = plt.subplots(figsize=(9, 5))
    for modification, table in metagene_tables.items():
        values = pd.to_numeric(table.get("metagene_pos"), errors="coerce").dropna()
        if len(values):
            density, edges = np.histogram(values, bins=90, range=(0, 3), density=True)
            centres = (edges[:-1] + edges[1:]) / 2
            ax.plot(centres, density, label=f"{modification} (n={len(values):,})")
    ax.axvline(1, color="grey", linestyle="--", linewidth=0.8)
    ax.axvline(2, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xlim(0, 3)
    ax.set_xticks([0.5, 1.5, 2.5], ["5′UTR", "CDS", "3′UTR"])
    ax.set_ylabel("Site density")
    ax.set_title("Modification-site metagene profiles")
    ax.legend()
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def plot_site_summary(summary: pd.DataFrame, path: str | Path) -> None:
    destination = Path(path)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].bar(summary["modification"], summary["n_filtered_sites"], color="#457B9D")
    axes[0].set_ylabel("Filtered sites")
    axes[0].set_title("Modification-site counts")
    bottom = np.zeros(len(summary))
    for region, colour in (("5UTR", "#A8DADC"), ("CDS", "#457B9D"), ("3UTR", "#E76F51")):
        values = summary[f"pct_mapped_{region}"].fillna(0).to_numpy()
        axes[1].bar(summary["modification"], values, bottom=bottom, label=region, color=colour)
        bottom += values
    axes[1].set_ylabel("% of transcript-region-mapped sites")
    axes[1].set_title("Transcript-region distribution")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)
