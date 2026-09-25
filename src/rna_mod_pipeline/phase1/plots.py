from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .filtering import (
    COVERAGE_BINS,
    COVERAGE_PLOT_MAX,
    COVERAGE_THRESHOLDS,
    FRACTION_BINS,
    FRACTION_THRESHOLDS,
    FilterSummary,
)


def _save(figure: plt.Figure, directory: Path, filename: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(directory / filename, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _chromosome_sort_key(chromosome: str) -> tuple[int, int, str]:
    label = chromosome.removeprefix("chr")
    if label.isdigit():
        return (0, int(label), "")
    main = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if label in main:
        return (0, main[label], "")
    return (1, 0, label)


def _plot_chromosomes(
    counts: dict[str, int],
    directory: Path,
    filename: str,
    title: str,
    ylabel: str,
    colour: str,
) -> None:
    ordered = sorted(counts.items(), key=lambda item: _chromosome_sort_key(item[0]))
    width = max(12, min(24, 0.35 * max(len(ordered), 1)))
    figure, axis = plt.subplots(figsize=(width, 5))
    if ordered:
        labels, values = zip(*ordered)
        positions = np.arange(len(labels))
        axis.bar(positions, values, color=colour, edgecolor="white", linewidth=0.3)
        axis.set_xticks(positions, labels, rotation=45, ha="right", fontsize=8)
        axis.set_ylabel(ylabel)
    else:
        axis.axis("off")
        axis.text(0.5, 0.5, "No sites available", ha="center", va="center")
    axis.set_title(title)
    _save(figure, directory, filename)


def _threshold_labels(
    thresholds: tuple[int | float, ...],
    prefix: str,
) -> list[str]:
    return [
        f"{prefix}>0" if threshold == 0 else f"{prefix}≥{threshold:g}"
        for threshold in thresholds
    ]


def _threshold_panel(
    axis: plt.Axes,
    labels: list[str],
    counts: np.ndarray,
    total: int,
    colour: str,
    title: str,
) -> None:
    positions = np.arange(len(labels))
    axis.barh(positions, counts, color=colour)
    axis.set_yticks(positions, labels)
    limit = max(float(np.max(counts)) if len(counts) else 0.0, 1.0)
    axis.set_xlim(0, limit * 1.28)
    for position, count in enumerate(counts):
        percentage = 100 * int(count) / total if total else 0.0
        axis.text(
            int(count) + limit * 0.015,
            position,
            f"{int(count):,} ({percentage:.1f}%)",
            va="center",
            fontsize=8,
        )
    axis.set_xlabel("Sites passing")
    axis.set_title(title)


def _plot_threshold_sensitivity(
    summary: FilterSummary,
    directory: Path,
    display: str,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(14, 5))
    _threshold_panel(
        axes[0],
        _threshold_labels(COVERAGE_THRESHOLDS, "coverage "),
        summary.coverage_threshold_counts,
        summary.total_sites,
        "#4682B4",
        "Sites passing coverage thresholds",
    )
    _threshold_panel(
        axes[1],
        _threshold_labels(FRACTION_THRESHOLDS, "fraction "),
        summary.fraction_threshold_counts,
        summary.total_sites,
        "#E67E22",
        "Sites passing modification-fraction thresholds",
    )
    figure.suptitle(f"{display} threshold sensitivity")
    _save(figure, directory, "04_threshold_sensitivity.png")


def _plot_filtered_qc(
    filtered: pd.DataFrame,
    summary: FilterSummary,
    directory: Path,
    display: str,
) -> None:
    chromosome_counts = (
        filtered["chrom"].astype(str).value_counts().astype(int).to_dict()
        if not filtered.empty
        else {}
    )
    threshold_text = (
        f"coverage ≥{summary.coverage_min}, "
        f"fraction ≥{summary.fraction_min:g}%"
    )
    _plot_chromosomes(
        chromosome_counts,
        directory,
        "05_filtered_chromosome_distribution.png",
        f"High-confidence {display} sites per chromosome ({threshold_text})",
        f"High-confidence {display} sites",
        "#228B22",
    )

    figure, axis = plt.subplots(figsize=(10, 5))
    if filtered.empty:
        axis.axis("off")
        axis.text(0.5, 0.5, "No sites passed filtering", ha="center", va="center")
    else:
        axis.hist(
            pd.to_numeric(filtered["fraction_modified"], errors="coerce").dropna(),
            bins=50,
            color="#228B22",
            edgecolor="white",
            linewidth=0.3,
        )
        axis.set_xlabel("Modified reads (%)")
        axis.set_ylabel("High-confidence sites")
    axis.set_title(
        f"{display} modification-fraction distribution after filtering "
        f"(n={len(filtered):,})"
    )
    _save(figure, directory, "06_filtered_fraction_distribution.png")

    figure, axis = plt.subplots(figsize=(8, 6))
    if filtered.empty:
        axis.axis("off")
        axis.text(0.5, 0.5, "No sites passed filtering", ha="center", va="center")
    else:
        plotted = (
            filtered
            if len(filtered) <= 50_000
            else filtered.sample(50_000, random_state=42)
        )
        axis.scatter(
            plotted["Nvalid_cov"],
            plotted["fraction_modified"],
            alpha=0.3,
            s=5,
            color="#228B22",
        )
        axis.set_xlabel("Valid read coverage")
        axis.set_ylabel("Modified reads (%)")
    axis.set_title(f"{display} coverage versus modification fraction")
    _save(figure, directory, "07_coverage_vs_fraction.png")


def plot_filter_qc(
    summary: FilterSummary,
    filtered: pd.DataFrame,
    plot_dir: str | Path,
    display: str,
) -> None:
    directory = Path(plot_dir)

    figure, axis = plt.subplots(figsize=(10, 5))
    centers = (COVERAGE_BINS[:-1] + COVERAGE_BINS[1:]) / 2
    axis.bar(
        centers,
        summary.coverage_histogram,
        width=np.diff(COVERAGE_BINS),
        color="#4682B4",
    )
    axis.axvline(summary.coverage_min, color="#C0392B", linestyle="--")
    axis.set_xlim(0, COVERAGE_PLOT_MAX)
    axis.text(
        0.98,
        0.96,
        (
            f"Coverage ≥{COVERAGE_PLOT_MAX}: "
            f"{summary.coverage_at_or_above_plot_max:,}\n"
            f"Coverage >{COVERAGE_BINS[-1]:g} "
            f"(outside histogram bins): "
            f"{summary.coverage_above_histogram_max:,}"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "#CCCCCC", "alpha": 0.9},
    )
    axis.set_xlabel("Valid read coverage")
    axis.set_ylabel("Sites")
    axis.set_title(f"{display} coverage distribution")
    _save(figure, directory, "01_coverage_distribution.png")

    figure, axis = plt.subplots(figsize=(10, 5))
    centers = (FRACTION_BINS[:-1] + FRACTION_BINS[1:]) / 2
    axis.bar(
        centers,
        summary.fraction_histogram,
        width=np.diff(FRACTION_BINS),
        color="#E67E22",
    )
    axis.axvline(summary.fraction_min, color="#C0392B", linestyle="--")
    axis.set_xlabel("Modified reads (%)")
    axis.set_ylabel("Sites with modification signal")
    axis.set_title(f"{display} modification-fraction distribution")
    _save(figure, directory, "02_fraction_modified_distribution.png")

    _plot_chromosomes(
        summary.chromosome_counts,
        directory,
        "03_chromosome_distribution.png",
        f"{display} assessed sites per chromosome before filtering",
        "Assessed sites",
        "#4682B4",
    )
    _plot_threshold_sensitivity(summary, directory, display)

    grid = summary.threshold_counts
    figure, axis = plt.subplots(figsize=(9, 5))
    image = axis.imshow(grid, cmap="YlOrRd", aspect="auto")
    for row in range(grid.shape[0]):
        for column in range(grid.shape[1]):
            percentage = (
                100 * int(grid[row, column]) / summary.total_sites
                if summary.total_sites
                else 0.0
            )
            axis.text(
                column,
                row,
                f"{grid[row, column]:,}\n({percentage:.1f}%)",
                ha="center",
                va="center",
                fontsize=8,
            )
    axis.set_xticks(
        range(len(FRACTION_THRESHOLDS)),
        [
            ">0%" if threshold == 0 else f"≥{threshold:g}%"
            for threshold in FRACTION_THRESHOLDS
        ],
    )
    axis.set_yticks(
        range(len(COVERAGE_THRESHOLDS)),
        [f"≥{threshold:g}" for threshold in COVERAGE_THRESHOLDS],
    )
    coverage_index = (
        COVERAGE_THRESHOLDS.index(summary.coverage_min)
        if summary.coverage_min in COVERAGE_THRESHOLDS
        else None
    )
    fraction_index = (
        FRACTION_THRESHOLDS.index(summary.fraction_min)
        if summary.fraction_min in FRACTION_THRESHOLDS
        else None
    )
    if coverage_index is not None and fraction_index is not None:
        axis.add_patch(
            plt.Rectangle(
                (fraction_index - 0.5, coverage_index - 0.5),
                1,
                1,
                fill=False,
                edgecolor="#1F4E79",
                linewidth=2.5,
            )
        )
        axis.text(
            0.99,
            -0.16,
            "Outlined cell = active filter",
            transform=axis.transAxes,
            ha="right",
            fontsize=8,
            color="#1F4E79",
        )
    axis.set_xlabel("Modification fraction")
    axis.set_ylabel("Coverage")
    axis.set_title(f"{display} combined threshold counts")
    figure.colorbar(image, ax=axis, label="Sites")
    _save(figure, directory, "04a_joint_threshold_heatmap.png")
    _plot_filtered_qc(filtered, summary, directory, display)


def plot_drach_summary(
    annotated: pd.DataFrame,
    plot_dir: str | Path,
) -> None:
    directory = Path(plot_dir)
    assessed = annotated["is_DRACH"].notna()
    drach = int(annotated.loc[assessed, "is_DRACH"].sum())
    non_drach = int(assessed.sum()) - drach
    unavailable = int((~assessed).sum())
    values = [drach, non_drach]
    labels = ["DRACH", "Non-DRACH"]
    colors = ["#2ECC71", "#E74C3C"]
    if unavailable:
        values.append(unavailable)
        labels.append("Context unavailable")
        colors.append("#95A5A6")
    figure, axis = plt.subplots(figsize=(7, 7))
    axis.pie(
        values,
        labels=labels,
        autopct="%1.1f%%",
        colors=colors,
        startangle=90,
    )
    axis.set_title("DRACH annotation of filtered m⁶A sites")
    _save(figure, directory, "08_drach_pie.png")

    motifs = (
        annotated.loc[annotated["motif_5mer"].ne("N/A"), "motif_5mer"]
        .value_counts()
        .head(20)
    )
    figure, axis = plt.subplots(figsize=(12, 6))
    colors = [
        "#2ECC71"
        if bool(
            annotated.loc[
                annotated["motif_5mer"].eq(motif), "is_DRACH"
            ].iloc[0]
        )
        else "#E74C3C"
        for motif in motifs.index
    ]
    axis.bar(np.arange(len(motifs)), motifs.values, color=colors)
    axis.set_xticks(np.arange(len(motifs)), motifs.index, rotation=45, ha="right")
    axis.set_ylabel("Sites")
    axis.set_title("Most frequent m⁶A-centred five-mers")
    _save(figure, directory, "09_top_motifs.png")


def plot_metagene_summary(
    annotated: pd.DataFrame,
    plot_dir: str | Path,
    display: str,
) -> None:
    directory = Path(plot_dir)
    mapped = annotated.loc[annotated["metagene_pos"].notna()]
    figure, axis = plt.subplots(figsize=(12, 5))
    axis.hist(mapped["metagene_pos"], bins=150, density=True, color="#4682B4")
    axis.axvline(1, color="black", linestyle="--", linewidth=0.8)
    axis.axvline(2, color="#C0392B", linestyle="--", linewidth=0.8)
    axis.set_xlim(0, 3)
    axis.set_xticks([0.5, 1.5, 2.5], ["5′UTR", "CDS", "3′UTR"])
    axis.set_ylabel("Density")
    axis.set_title(f"{display} metagene profile")
    _save(figure, directory, "12_metagene_profile.png")

    order = ["5UTR", "CDS", "3UTR", "intergenic/intronic"]
    counts = annotated["region"].value_counts().reindex(order, fill_value=0)
    figure, axis = plt.subplots(figsize=(8, 5))
    axis.bar(
        np.arange(len(order)),
        counts.values,
        color=["#3498DB", "#2ECC71", "#E74C3C", "#95A5A6"],
    )
    axis.set_xticks(np.arange(len(order)), ["5′UTR", "CDS", "3′UTR", "Unmapped"])
    axis.set_ylabel("Sites")
    axis.set_title(f"{display} sites by transcript region")
    _save(figure, directory, "14_region_distribution.png")
