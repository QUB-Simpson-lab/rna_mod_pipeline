from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .plot_data import (
    COVERAGE_ORDER,
    INFERENCE_ORDER,
    membership_count_source,
    pairwise_count_source,
    pairwise_plot_source,
    select_pairwise_top,
    triple_plot_source,
    validation_summary_source,
)

DISPLAY = {"ornament": "oRNAment", "encori": "ENCORI", "postar3": "POSTAR3"}
STATE_COLOURS = {
    "both_significant_same_direction": "#277A4B",
    "both_significant_opposite_direction": "#E07A28",
    "first_only_significant": "#C44E52",
    "second_only_significant": "#4C72B0",
    "neither_significant": "#B7BBC2",
    "any_quality_excluded": "#D6A84B",
    "any_boundary_estimate": "#7A3E9D",
    "any_non_estimable": "#555555",
}

def _label(value: str) -> str:
    return value.replace("_", " ").capitalize()


def _placeholder(axis: plt.Axes, message: str) -> None:
    axis.axis("off")
    axis.text(0.5, 0.5, message, ha="center", va="center")


def _pairwise_scatter(source: pd.DataFrame, path: Path, title: str, dpi: int) -> int:
    pairs = list(source["pair"].drop_duplicates())
    fig, axes = plt.subplots(
        1, max(1, len(pairs)), figsize=(6 * max(1, len(pairs)), 5.4), squeeze=False
    )
    plotted_total = 0
    if not pairs:
        _placeholder(axes[0, 0], "No database pairs were available")
    for axis, pair in zip(axes[0], pairs):
        group = source[
            source["pair"].eq(pair)
            & source["coverage_state"].eq("both_present")
            & np.isfinite(source["first_log2_odds_ratio"])
            & np.isfinite(source["second_log2_odds_ratio"])
        ]
        first = source.loc[source["pair"].eq(pair), "first_database"].iloc[0]
        second = source.loc[source["pair"].eq(pair), "second_database"].iloc[0]
        if group.empty:
            _placeholder(axis, "No shared finite odds-ratio estimates")
            axis.set_title(f"{DISPLAY[first]} vs {DISPLAY[second]}")
            continue
        plotted_total += len(group)
        for state in INFERENCE_ORDER:
            selected = group[group["inference_state"].eq(state)]
            if selected.empty:
                continue
            axis.scatter(
                selected["first_log2_odds_ratio"],
                selected["second_log2_odds_ratio"],
                s=28,
                alpha=0.8,
                color=STATE_COLOURS[state],
                label=_label(state),
            )
        primary = group[group["eligible_for_primary_correlation"]]
        primary_rho = np.nan
        if len(primary) >= 3:
            primary_rho = primary["first_log2_odds_ratio"].corr(
                primary["second_log2_odds_ratio"], method="spearman"
            )
        raw_rho = np.nan
        if len(group) >= 3:
            raw_rho = group["first_log2_odds_ratio"].corr(
                group["second_log2_odds_ratio"], method="spearman"
            )
        axis.text(
            0.03,
            0.97,
            f"Primary eligible n={len(primary)}; ρ={primary_rho:.3f}\n"
            f"Raw all-finite sensitivity n={len(group)}; ρ={raw_rho:.3f}",
            transform=axis.transAxes,
            va="top",
            fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "#DDDDDD", "alpha": 0.9},
        )
        axis.axhline(0, color="#555555", linestyle="--", linewidth=0.7)
        axis.axvline(0, color="#555555", linestyle="--", linewidth=0.7)
        axis.set_xlabel(f"{DISPLAY[first]} log₂ odds ratio")
        axis.set_ylabel(f"{DISPLAY[second]} log₂ odds ratio")
        axis.set_title(f"{DISPLAY[first]} vs {DISPLAY[second]}")
        axis.grid(color="#ECECEC", linewidth=0.6)
        axis.legend(frameon=False, fontsize=6, loc="lower right")
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return plotted_total


def _stacked(
    axis: plt.Axes,
    table: pd.DataFrame,
    pairs: list[str],
    order: tuple[str, ...],
    colours: dict[str, str],
    title: str,
) -> None:
    if not pairs:
        _placeholder(axis, "No database pairs were available")
        return
    bottom = np.zeros(len(pairs), dtype=float)
    for category in order:
        values = [
            int(
                table.loc[
                    table["pair"].eq(pair) & table["category"].eq(category), "count"
                ].sum()
            )
            for pair in pairs
        ]
        axis.bar(
            range(len(pairs)),
            values,
            bottom=bottom,
            color=colours[category],
            label=_label(category),
        )
        bottom += values
    axis.set_xticks(
        range(len(pairs)),
        [
            " vs ".join(DISPLAY.get(item, item) for item in pair.split("_vs_"))
            for pair in pairs
        ],
        rotation=20,
        ha="right",
    )
    axis.set_ylabel("Number of RBPs")
    axis.set_title(title)
    axis.legend(frameon=False, fontsize=7)
    axis.grid(axis="y", color="#ECECEC", linewidth=0.6)


def _pairwise_agreement(
    counts: pd.DataFrame, path: Path, title: str, dpi: int
) -> int:
    pairs = list(counts["pair"].drop_duplicates())
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    coverage = counts[counts["dimension"].eq("coverage")]
    inference = counts[counts["dimension"].eq("inference_among_shared")]
    _stacked(
        axes[0],
        coverage,
        pairs,
        COVERAGE_ORDER,
        {
            "both_present": "#4C72B0",
            "first_only": "#DD8452",
            "second_only": "#55A868",
        },
        "Database coverage",
    )
    _stacked(
        axes[1],
        inference,
        pairs,
        INFERENCE_ORDER,
        STATE_COLOURS,
        "Inference among RBPs covered by both databases",
    )
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(counts)


def _pairwise_top(source: pd.DataFrame, path: Path, title: str, dpi: int) -> int:
    pairs = list(source["pair"].drop_duplicates())
    fig, axes = plt.subplots(
        1, max(1, len(pairs)), figsize=(6 * max(1, len(pairs)), 7), squeeze=False
    )
    if not pairs:
        _placeholder(axes[0, 0], "No concordant significant pairwise RBPs")
    for axis, pair in zip(axes[0], pairs):
        group = source[source["pair"].eq(pair)].sort_values(
            ["mean_log2_odds_ratio", "RBP"], kind="mergesort"
        )
        if group.empty:
            _placeholder(axis, "No concordant significant finite estimates")
            axis.set_title(pair.replace("_vs_", " vs "))
            continue
        y = np.arange(len(group))
        first = group["first_database"].iloc[0]
        second = group["second_database"].iloc[0]
        axis.hlines(
            y,
            group["first_log2_odds_ratio"],
            group["second_log2_odds_ratio"],
            color="#B7BBC2",
            linewidth=1,
        )
        axis.scatter(
            group["first_log2_odds_ratio"],
            y,
            color="#C44E52",
            s=30,
            label=DISPLAY[first],
        )
        axis.scatter(
            group["second_log2_odds_ratio"],
            y,
            color="#4C72B0",
            s=30,
            label=DISPLAY[second],
        )
        axis.axvline(0, color="#555555", linewidth=0.7)
        axis.set_yticks(y, group["RBP"], fontsize=7)
        axis.set_xlabel("log₂ odds ratio")
        axis.set_title(f"{DISPLAY[first]} vs {DISPLAY[second]}")
        axis.grid(axis="x", color="#ECECEC", linewidth=0.6)
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(source)


def _triple(source: pd.DataFrame, path: Path, title: str, dpi: int) -> int:
    selected = source[source.get("selected_for_plot", False)].copy()
    fig, axis = plt.subplots(figsize=(9, max(4.5, 0.34 * max(1, len(selected)))))
    if selected.empty:
        _placeholder(axis, "No triple-significant RBPs with three finite estimates")
        plotted = 0
    else:
        selected = selected.sort_values(
            ["rank_magnitude", "RBP"], ascending=[True, True], kind="mergesort"
        )
        y = np.arange(len(selected))
        colours = {"ornament": "#C44E52", "encori": "#4C72B0", "postar3": "#55A868"}
        for database in ("ornament", "encori", "postar3"):
            axis.scatter(
                selected[f"{database}_log2_odds_ratio"],
                y,
                s=32,
                color=colours[database],
                label=DISPLAY[database],
            )
        labels = [
            f"{row.RBP}{' †' if row.direction_conflict else ''}"
            for row in selected.itertuples()
        ]
        axis.set_yticks(y, labels, fontsize=8)
        axis.axvline(0, color="#555555", linewidth=0.7)
        axis.set_xlabel("log₂ odds ratio")
        axis.set_ylabel("† direction conflict")
        axis.grid(axis="x", color="#ECECEC", linewidth=0.6)
        axis.legend(frameon=False)
        plotted = len(selected)
    axis.set_title(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return plotted


def _validation_summary(
    source: pd.DataFrame, path: Path, title: str, dpi: int
) -> int:
    dimensions = [
        ("databases_present", "Number of databases covering each RBP"),
        ("databases_significant", "Number of significant databases"),
        ("multi_database_consensus", "Consensus among ≥2 significant databases"),
        ("database_rbp_state", "All database × RBP result states"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for axis, (dimension, subtitle) in zip(axes.flat, dimensions):
        group = source[source["dimension"].eq(dimension)]
        if group.empty:
            _placeholder(axis, "No results")
            continue
        axis.bar(
            np.arange(len(group)),
            group["count"],
            color=plt.get_cmap("tab10").colors[: len(group)],
        )
        axis.set_xticks(
            np.arange(len(group)),
            [_label(value) for value in group["category"]],
            rotation=25,
            ha="right",
            fontsize=8,
        )
        axis.set_ylabel("Number of RBPs")
        axis.set_title(subtitle)
        axis.grid(axis="y", color="#ECECEC", linewidth=0.6)
        for index, count in enumerate(group["count"]):
            axis.text(index, count, str(int(count)), ha="center", va="bottom", fontsize=8)
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(source)


def _membership_counts(
    source: pd.DataFrame, path: Path, title: str, dpi: int
) -> int:
    dimensions = (
        ("coverage_membership", "Exact database coverage membership"),
        ("significant_membership", "Exact significant-result membership"),
    )
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for axis, (dimension, subtitle) in zip(axes, dimensions):
        group = source[source["dimension"].eq(dimension)].copy()
        if group.empty:
            _placeholder(axis, "No membership results")
            continue
        group["display"] = group["membership"].map(
            lambda value: " + ".join(
                DISPLAY.get(item, item) for item in value.split("+")
            )
            if value != "none"
            else "None"
        )
        axis.barh(
            np.arange(len(group)),
            group["count"],
            color="#4C72B0" if dimension == "coverage_membership" else "#277A4B",
        )
        axis.set_yticks(np.arange(len(group)), group["display"], fontsize=8)
        axis.invert_yaxis()
        axis.set_xlabel("Number of RBPs")
        axis.set_title(subtitle)
        axis.grid(axis="x", color="#ECECEC", linewidth=0.6)
        for index, count in enumerate(group["count"]):
            axis.text(count, index, f" {int(count)}", va="center", fontsize=8)
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(source)


def create_crossdb_plots(
    combined: pd.DataFrame,
    pairwise: dict[str, pd.DataFrame],
    output_dir: str | Path,
    *,
    modification: str,
    design: str,
    context: str,
    fdr_threshold: float,
    top_n: int = 20,
    dpi: int = 180,
) -> list[Path]:
    output = Path(output_dir)
    plot_dir = output / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    pair_source = pairwise_plot_source(combined, pairwise, fdr_threshold)
    count_source = pairwise_count_source(pair_source)
    top_source = select_pairwise_top(pair_source, top_n)
    triple_source = triple_plot_source(combined, fdr_threshold, top_n)
    databases = tuple(
        database
        for database in ("ornament", "encori", "postar3")
        if f"{database}_odds_ratio" in combined
    )
    summary_source = validation_summary_source(
        combined, fdr_threshold, databases
    )
    membership_source = membership_count_source(combined, databases)

    sources = {
        "pairwise_points": pair_source,
        "pairwise_counts": count_source,
        "pairwise_shared_top": top_source,
        "triple_validation": triple_source,
        "validation_summary": summary_source,
        "database_membership_counts": membership_source,
    }
    source_paths: dict[str, Path] = {}
    for name, table in sources.items():
        path = plot_dir / f"{name}_source.tsv"
        table.to_csv(path, sep="\t", index=False)
        source_paths[name] = path

    title_stub = f"{modification} | {design} | {context}"
    plot_specs = [
        (
            "pairwise_scatter",
            _pairwise_scatter,
            pair_source,
            source_paths["pairwise_points"],
            "Pairwise effect-size concordance "
            "(primary ρ quality-qualified; raw sensitivity shown)",
        ),
        (
            "pairwise_agreement",
            _pairwise_agreement,
            count_source,
            source_paths["pairwise_counts"],
            "Pairwise coverage and agreement",
        ),
        (
            "pairwise_shared_top",
            _pairwise_top,
            top_source,
            source_paths["pairwise_shared_top"],
            f"Top {top_n} concordant significant RBPs per pair",
        ),
        (
            "triple_validation",
            _triple,
            triple_source,
            source_paths["triple_validation"],
            f"Top {top_n} triple-significant RBPs",
        ),
        (
            "validation_summary",
            _validation_summary,
            summary_source,
            source_paths["validation_summary"],
            "Cross-database validation summary",
        ),
        (
            "database_membership_counts",
            _membership_counts,
            membership_source,
            source_paths["database_membership_counts"],
            "Exact database memberships (area-free Venn alternative)",
        ),
    ]
    manifest_rows = []
    plot_paths = []
    render_units = {
        "pairwise_scatter": "RBP-pair point",
        "pairwise_agreement": "summary-category bar segment",
        "pairwise_shared_top": "selected RBP-pair row",
        "triple_validation": "selected RBP",
        "validation_summary": "summary-category bar",
        "database_membership_counts": "membership-category bar",
    }
    nonexclusive_counts = {
        "pairwise_agreement",
        "validation_summary",
        "database_membership_counts",
    }
    for name, function, source, source_path, subtitle in plot_specs:
        path = plot_dir / f"{name}.png"
        n_rendered = function(source, path, f"{title_stub}: {subtitle}", dpi)
        plot_paths.append(path)
        manifest_rows.append(
            {
                "plot_type": name,
                "path": path.relative_to(output).as_posix(),
                "source_data": source_path.relative_to(output).as_posix(),
                "n_source_rows": len(source),
                "n_unique_rbps": (
                    int(source["RBP"].nunique())
                    if "RBP" in source
                    else int(combined["RBP"].nunique())
                ),
                "n_rendered_records": n_rendered,
                "rendered_record_unit": render_units[name],
                "cross_panel_count_values_nonexclusive": (
                    name in nonexclusive_counts
                ),
                "fdr_threshold": fdr_threshold,
                "top_n": top_n,
                "confidence_intervals_drawn": False,
                "primary_correlation_scope": (
                    "both database estimates finite and inference eligible"
                    if name == "pairwise_scatter"
                    else "not_applicable"
                ),
                "raw_sensitivity_scope": (
                    "all shared finite effects, displayed separately in the annotation"
                    if name == "pairwise_scatter"
                    else "not_applicable"
                ),
            }
        )
    manifest_path = plot_dir / "crossdb_plot_manifest.tsv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)
    return [*source_paths.values(), *plot_paths, manifest_path]
