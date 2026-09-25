from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_rbp_expression(table: pd.DataFrame, path: str | Path) -> None:
    or_column = "mh_or" if "mh_or" in table else "odds_ratio"
    odds = pd.to_numeric(table[or_column], errors="coerce")
    frame = table[odds.gt(0) & np.isfinite(odds)].copy()
    frame["log2_or"] = np.log2(pd.to_numeric(frame[or_column], errors="coerce"))
    colours = np.where(
        frame["expression_status"].isin(["detected", "expressed"]), "#2A9D8F", "#B8B8B8"
    )
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(frame["log2_or"], frame["expression_value"], c=colours, alpha=0.75)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("log₂ odds ratio")
    ax.set_ylabel("Expression value")
    ax.set_title("RBP enrichment and MDA-MB-231 expression")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def plot_gene_association(table: pd.DataFrame, path: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(
        np.log2(table["expression_value"].astype(float) + 1),
        table["modification_site_count"],
        s=10,
        alpha=0.3,
    )
    ax.set_xlabel("log₂(expression + 1)")
    ax.set_ylabel("Modification-site count per gene")
    ax.set_title("Expression versus modification-site count")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)
