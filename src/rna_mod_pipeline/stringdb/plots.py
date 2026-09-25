from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


def plot_network(
    interactions: pd.DataFrame,
    path: str | Path,
    title: str,
    nodes: list[str] | tuple[str, ...] = (),
) -> None:
    graph = nx.Graph()
    graph.add_nodes_from(str(node) for node in nodes if str(node))
    for _, edge in interactions.iterrows():
        graph.add_edge(
            str(edge["protein_A"]),
            str(edge["protein_B"]),
            weight=float(edge.get("combined_score", 0)),
        )
    fig, ax = plt.subplots(figsize=(8, 7))
    if not graph.nodes:
        ax.text(0.5, 0.5, "No STRING edges at this threshold", ha="center", va="center")
        ax.axis("off")
    else:
        positions = nx.spring_layout(graph, seed=42, weight="weight")
        isolates = set(nx.isolates(graph))
        nx.draw_networkx(
            graph,
            positions,
            ax=ax,
            node_size=700,
            node_color=[
                "#B8B8B8" if node in isolates else "#6BAED6"
                for node in graph.nodes
            ],
            font_size=8,
            width=[max(graph[u][v]["weight"] * 2, 0.5) for u, v in graph.edges],
        )
        if isolates:
            ax.text(
                0.01,
                0.01,
                f"Grey nodes: {len(isolates)} resolved isolate(s)",
                transform=ax.transAxes,
                fontsize=8,
            )
        ax.axis("off")
    ax.set_title(title)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def plot_enrichment(table: pd.DataFrame, path: str | Path, title: str) -> None:
    frame = table.sort_values("fdr").head(20).copy()
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * max(len(frame), 1))))
    if frame.empty:
        ax.text(0.5, 0.5, "No STRING terms passed the FDR threshold", ha="center", va="center")
        ax.axis("off")
    else:
        values = -np.log10(pd.to_numeric(frame["fdr"], errors="coerce").clip(1e-300))
        ax.barh(frame["description"].astype(str), values, color="#457B9D")
        ax.invert_yaxis()
        ax.set_xlabel("−log₁₀ STRING FDR")
    ax.set_title(title)
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(destination, dpi=180, bbox_inches="tight")
    plt.close(fig)
