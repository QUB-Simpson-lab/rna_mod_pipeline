from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _dataset_labels(table: pd.DataFrame) -> pd.Series:
    return (
        table["rbp"].astype(str)
        + " | "
        + table["dataset_id"].astype(str)
        + " | "
        + table["cell_line"].astype(str)
    )


def _regulatory_matrix(table: pd.DataFrame) -> pd.DataFrame:
    source = table.copy()
    source["perturbation_dataset"] = (
        source["knocked_down_rbp"].astype(str)
        + " | "
        + source["dataset_id"].astype(str)
        + " | "
        + source["cell_line"].astype(str)
    )
    return source.pivot(
        index="affected_rbp",
        columns="perturbation_dataset",
        values="log2fc",
    )


def regulatory_plot_source(
    table: pd.DataFrame,
    top_n: int = 40,
) -> pd.DataFrame:
    """Select a readable, deterministic subset while retaining full edge rows."""
    if top_n < 1:
        raise ValueError("top_n must be at least 1")
    if table.empty:
        result = table.copy()
        result["regulatory_plot_rank"] = pd.Series(dtype="Int64")
        result["maximum_absolute_log2fc"] = pd.Series(dtype=float)
        return result
    source = table.copy()
    source["log2fc"] = pd.to_numeric(source["log2fc"], errors="coerce")
    source = source.dropna(subset=["affected_rbp", "log2fc"])
    maxima = (
        source.assign(_absolute_log2fc=source["log2fc"].abs())
        .groupby("affected_rbp", as_index=False)["_absolute_log2fc"]
        .max()
        .rename(columns={"_absolute_log2fc": "maximum_absolute_log2fc"})
        .sort_values(
            ["maximum_absolute_log2fc", "affected_rbp"],
            ascending=[False, True],
            kind="mergesort",
        )
        .head(top_n)
        .reset_index(drop=True)
    )
    maxima["regulatory_plot_rank"] = maxima.index + 1
    selected = source.merge(maxima, on="affected_rbp", how="inner")
    return selected.sort_values(
        ["regulatory_plot_rank", "dataset_id"], kind="mergesort"
    ).reset_index(drop=True)


def plot_orthogonal(table: pd.DataFrame, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, max(3, 0.45 * max(len(table), 1))))
    if table.empty:
        ax.text(0.5, 0.5, "No testable KnockRBP datasets", ha="center", va="center")
        ax.axis("off")
    else:
        values = -np.log10(
            pd.to_numeric(
                table["primary_modified_gene_universe_fdr"], errors="coerce"
            ).clip(lower=np.finfo(float).tiny)
        )
        colours = np.where(table["small_count_warning"], "#B8B8B8", "#2A9D8F")
        ax.barh(_dataset_labels(table), values, color=colours)
        ax.axvline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8)
        ax.set_xlabel("−log₁₀ primary-background FDR")
        ax.set_title("KnockRBP target–DEG overlap")
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def plot_regulatory_network(
    table: pd.DataFrame,
    path: str | Path,
    *,
    top_n: int = 40,
) -> pd.DataFrame:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = regulatory_plot_source(table, top_n=top_n)
    n_rbps = source["affected_rbp"].nunique() if not source.empty else 0
    fig, ax = plt.subplots(figsize=(12, max(6, 0.24 * n_rbps)))
    if source.empty:
        ax.text(0.5, 0.5, "No RBP–RBP perturbation relationships", ha="center", va="center")
        ax.axis("off")
    else:
        pivot = _regulatory_matrix(source)
        order = (
            source[["affected_rbp", "regulatory_plot_rank"]]
            .drop_duplicates()
            .sort_values("regulatory_plot_rank")["affected_rbp"]
        )
        pivot = pivot.reindex(order)
        values = pivot.to_numpy(dtype=float)
        colour_limit = float(np.nanmax(np.abs(values)))
        if not np.isfinite(colour_limit) or colour_limit == 0:
            colour_limit = 1.0
        colour_map = plt.get_cmap("RdBu_r").copy()
        colour_map.set_bad("#E6E6E6")
        image = ax.imshow(
            np.ma.masked_invalid(values),
            cmap=colour_map,
            vmin=-colour_limit,
            vmax=colour_limit,
            aspect="auto",
        )
        ax.set_xticks(
            range(len(pivot.columns)), pivot.columns, rotation=45, ha="right"
        )
        ax.set_yticks(range(len(pivot.index)), pivot.index)
        ax.tick_params(axis="y", labelsize=8)
        ax.set_title(
            f"Top {n_rbps} affected RBPs by maximum |log₂FC|\n"
            "Filtered KnockRBP DEGs only; grey = no retained DEG"
        )
        fig.colorbar(image, ax=ax, label="log₂ fold change")
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)
    return source
