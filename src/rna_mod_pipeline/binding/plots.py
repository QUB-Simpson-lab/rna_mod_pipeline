from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from ..evidence_quality import add_significance_calls


DATABASE_LABELS = {
    "ornament": "oRNAment",
    "encori": "ENCORI",
    "postar3": "POSTAR3",
}


def _save(figure: plt.Figure, directory: Path, filename: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(directory / filename, dpi=180, bbox_inches="tight")
    plt.close(figure)


def _placeholder(directory: Path, filename: str, title: str, message: str) -> None:
    figure, axis = plt.subplots(figsize=(9, 4))
    axis.axis("off")
    axis.set_title(title)
    axis.text(0.5, 0.5, message, ha="center", va="center")
    _save(figure, directory, filename)


def _display_effects(
    odds_ratios: pd.Series,
) -> tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    values = pd.to_numeric(odds_ratios, errors="coerce").to_numpy(dtype=float)
    finite_positive = np.isfinite(values) & (values > 0)
    log_values = np.full(len(values), np.nan)
    log_values[finite_positive] = np.log2(values[finite_positive])
    observed = np.abs(log_values[np.isfinite(log_values)])
    cap = max(1.0, float(np.ceil(observed.max()))) if len(observed) else 1.0
    boundary = (values == 0) | np.isposinf(values)
    log_values[values == 0] = -cap
    log_values[np.isposinf(values)] = cap
    non_estimable = np.isnan(values) | np.isneginf(values)
    return log_values, cap, boundary, non_estimable


def _top_effects(
    results: pd.DataFrame,
    display: str,
    database: str,
    directory: Path,
    enriched: bool,
    top_n: int = 30,
) -> None:
    label = "enriched" if enriched else "depleted"
    filename = f"{database}_top_{label}.png"
    selected = _top_effect_rows(results, enriched=enriched, top_n=top_n)
    if selected.empty:
        _placeholder(
            directory,
            filename,
            f"Top {label} {DATABASE_LABELS[database]} RBPs at {display} sites",
            "No RBPs pass the quality-qualified FDR and direction criteria.",
        )
        return
    selected = selected.sort_values("odds_ratio", ascending=enriched)
    values, cap, boundary, _ = _display_effects(selected["odds_ratio"])
    figure, axis = plt.subplots(figsize=(10, max(5, 0.32 * len(selected))))
    colors = [
        "#E74C3C" if known else ("#2ECC71" if enriched else "#3498DB")
        for known in selected["is_known_related"]
    ]
    axis.barh(np.arange(len(selected)), values, color=colors)
    if boundary.any():
        marker = ">" if enriched else "<"
        axis.scatter(
            values[boundary],
            np.flatnonzero(boundary),
            marker=marker,
            color="black",
            s=25,
            zorder=3,
            label="Boundary estimate",
        )
    axis.set_yticks(np.arange(len(selected)))
    axis.set_yticklabels(selected["RBP"], fontsize=8)
    axis.axvline(0, color="black", linewidth=0.8)
    axis.set_xlabel("log₂(odds ratio)")
    axis.set_title(
        f"Top {label} {DATABASE_LABELS[database]} RBPs at {display} sites"
    )
    handles = [
        Patch(facecolor="#E74C3C", label="Modification-related"),
        Patch(
            facecolor="#2ECC71" if enriched else "#3498DB",
            label=label.title(),
        ),
    ]
    if boundary.any():
        handles.append(
            plt.Line2D(
                [0],
                [0],
                marker=">" if enriched else "<",
                color="black",
                linestyle="none",
                label=f"Boundary shown at ±{cap:g}",
            )
        )
    axis.legend(handles=handles, fontsize=8)
    _save(figure, directory, filename)


def _top_effect_rows(
    results: pd.DataFrame,
    *,
    enriched: bool,
    top_n: int,
) -> pd.DataFrame:
    selected = results.loc[
        results["primary_direction"].eq(
            "Enriched" if enriched else "Depleted"
        )
    ].copy()
    return (
        selected.nlargest(top_n, "odds_ratio")
        if enriched
        else selected.nsmallest(top_n, "odds_ratio")
    )


def _loose_overview_data(
    results: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, int, pd.DataFrame]:
    selected = results.copy()
    selected["_odds_ratio"] = pd.to_numeric(
        selected["odds_ratio"], errors="coerce"
    )
    values = selected["_odds_ratio"].to_numpy(dtype=float)
    case_overlaps = pd.to_numeric(
        selected["mod_overlaps"], errors="coerce"
    )
    control_overlaps = pd.to_numeric(
        selected["bg_overlaps"], errors="coerce"
    )
    zero_zero = case_overlaps.eq(0) & control_overlaps.eq(0)
    shown = np.isfinite(values) & (values > 0) & ~zero_zero
    zero_boundary = (values == 0) & ~zero_zero.to_numpy()
    positive_infinity = np.isposinf(values) & ~zero_zero.to_numpy()
    non_estimable = ~(shown | zero_boundary | positive_infinity)
    ordered = (
        selected.loc[shown]
        .assign(_log2_or=lambda frame: np.log2(frame["_odds_ratio"]))
        .sort_values(["_log2_or", "RBP"], ascending=[False, True])
        .reset_index(drop=True)
    )
    boundary_rows = selected.loc[zero_boundary].copy()
    boundary_rows["_fdr"] = pd.to_numeric(
        boundary_rows["fdr"], errors="coerce"
    )
    boundary_rows = boundary_rows.sort_values(
        ["_fdr", "RBP"],
        ascending=[True, True],
        na_position="last",
        kind="mergesort",
    ).reset_index(drop=True)
    return (
        ordered,
        boundary_rows,
        int(positive_infinity.sum()),
        selected.loc[non_estimable].copy(),
    )


def plot_all_rbps(
    results: pd.DataFrame,
    display: str,
    database: str,
    directory: Path,
    fdr_threshold: float = 0.05,
) -> None:
    selected, zero_boundary, n_positive_infinity, non_estimable = (
        _loose_overview_data(results)
    )
    colors = []
    for row in selected.itertuples(index=False):
        if row.statistically_significant and not row.inference_eligible:
            colors.append("#D6A84B")
        elif row.is_known_related:
            colors.append("#E74C3C")
        elif row.primary_direction == "Enriched":
            colors.append("#2ECC71")
        elif row.primary_direction == "Depleted":
            colors.append("#3498DB")
        else:
            colors.append("#BDC3C7")
    figure, axis = plt.subplots(figsize=(16, 5.5))
    if selected.empty and zero_boundary.empty:
        axis.axis("off")
        axis.text(
            0.5,
            0.5,
            "No finite or OR=0 boundary estimates to display",
            ha="center",
            va="center",
        )
    else:
        axis.bar(
            np.arange(len(selected)),
            selected["_log2_or"],
            color=colors,
            width=1,
            edgecolor="none",
            zorder=2,
        )
        for index, row in selected.iterrows():
            if bool(row["is_known_related"]):
                offset = 4 if row["_log2_or"] >= 0 else -4
                axis.annotate(
                    str(row["RBP"]),
                    (index, row["_log2_or"]),
                    xytext=(1.5, offset),
                    textcoords="offset points",
                    rotation=90,
                    ha="center",
                    va="bottom" if row["_log2_or"] >= 0 else "top",
                    fontsize=5.2,
                    clip_on=False,
                )
        maximum = (
            float(selected["_log2_or"].abs().max())
            if not selected.empty else 1.0
        )
        limit = max(2.0, float(np.ceil(maximum * 4) / 4 + 0.5))
        zero_count = len(zero_boundary)
        if zero_count:
            boundary_x = np.arange(
                len(selected), len(selected) + zero_count
            )
            axis.bar(
                boundary_x,
                np.full(zero_count, -limit),
                width=1.0,
                color="#6C5CE7",
                edgecolor="white",
                linewidth=0.45,
                zorder=2,
            )
            label_size = 3.6 if len(results) > 250 else 4.2
            for x_value, row in zip(
                boundary_x, zero_boundary.itertuples(index=False)
            ):
                axis.text(
                    x_value,
                    -0.56 * limit,
                    str(row.RBP),
                    color="white",
                    fontsize=label_size,
                    rotation=90,
                    ha="center",
                    va="center",
                    clip_on=True,
                )
            axis.set_xlim(
                -0.5, len(selected) + zero_count - 0.5
            )
        elif len(selected):
            axis.set_xlim(-0.5, len(selected) - 0.5)
        axis.axhline(0, color="black", linewidth=0.8, zorder=3)
        axis.set_ylim(-limit, limit)
        axis.set_xlabel("Finite RBP estimates ranked by log₂(odds ratio)")
        axis.set_ylabel("log₂(odds ratio)")
        axis.grid(axis="y", color="#E6E6E6", linewidth=0.7)

    zero_count = len(zero_boundary)
    zero_significant = int(
        pd.to_numeric(zero_boundary.get("fdr"), errors="coerce")
        .lt(fdr_threshold)
        .sum()
    ) if zero_count else 0
    note = (
        f"Finite log₂(OR) bars shown: {len(selected)} of {len(results)} RBPs."
    )
    if zero_count:
        note += (
            f" {zero_count} observed OR=0 result"
            f"{'s' if zero_count != 1 else ''} "
            f"({zero_significant} raw FDR<{fdr_threshold:g}; sensitivity-only) "
            f"{'are' if zero_count != 1 else 'is'} shown as individual "
            "−∞ boundary bars."
        )
    if n_positive_infinity:
        note += (
            f" {n_positive_infinity} OR=∞ boundary result(s) are omitted."
        )
    if not non_estimable.empty:
        names = ", ".join(non_estimable["RBP"].astype(str))
        reasons = (
            set(non_estimable["estimate_reason"].dropna().astype(str))
            if "estimate_reason" in non_estimable
            else set()
        )
        zero_overlaps = (
            pd.to_numeric(non_estimable["mod_overlaps"], errors="coerce").eq(0)
            & pd.to_numeric(non_estimable["bg_overlaps"], errors="coerce").eq(0)
        ).all()
        if reasons == {"resource_empty"}:
            reason_text = "empty binding resource"
        elif zero_overlaps:
            reason_text = "zero overlaps in both groups"
        elif reasons == {"all_case_and_control_sites_bound"}:
            reason_text = "all case and comparison sites bound; no unbound sites"
        else:
            reason_text = "see estimate_reason in the enrichment TSV"
        note += (
            f" {len(non_estimable)} non-estimable RBP"
            f"{'s' if len(non_estimable) != 1 else ''} "
            f"({names}; {reason_text}) "
            f"{'are' if len(non_estimable) != 1 else 'is'} omitted."
        )
    figure.text(0.06, 0.012, note, fontsize=7.5, color="#555555", va="bottom")
    axis.set_title(
        f"{DATABASE_LABELS[database]} all-RBP overview at {display} sites "
        "— loose overlap"
    )
    handles = [
        Patch(facecolor="#E74C3C", label="Modification-related"),
        Patch(facecolor="#2ECC71", label="Primary significant enriched"),
        Patch(facecolor="#3498DB", label="Primary significant depleted"),
        Patch(
            facecolor="#D6A84B",
            label="Raw FDR call excluded from primary inference",
        ),
        Patch(facecolor="#BDC3C7", label="Not significant"),
    ]
    if zero_count:
        handles.append(
            Patch(
                facecolor="#6C5CE7",
                edgecolor="white",
                label=(
                    "Observed OR = 0 (−∞ boundary; sensitivity-only)"
                ),
            )
        )
    axis.legend(handles=handles, fontsize=8)
    directory.mkdir(parents=True, exist_ok=True)
    figure.tight_layout(rect=(0, 0.12, 1, 1))
    figure.savefig(
        directory / f"{database}_all_rbps_overview.png",
        dpi=180,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_machinery(
    results: pd.DataFrame,
    display: str,
    database: str,
    directory: Path,
) -> None:
    machinery = results.loc[results["is_known_related"]].copy()
    filename = f"{database}_machinery_heatmap.png"
    if machinery.empty:
        _placeholder(
            directory,
            filename,
            f"{display}-related RBPs in {DATABASE_LABELS[database]}",
            "No modification-related RBPs are represented in this panel.",
        )
        return
    machinery = machinery.sort_values("odds_ratio")
    values, cap, boundary, non_estimable = _display_effects(
        machinery["odds_ratio"]
    )
    labels = _machinery_labels(machinery)
    matrix = np.ma.masked_invalid(values[:, None])
    colormap = plt.get_cmap("RdBu_r").copy()
    colormap.set_bad("#BDBDBD")
    figure, axis = plt.subplots(figsize=(6, max(4, 0.35 * len(machinery))))
    image = axis.imshow(
        matrix,
        aspect="auto",
        cmap=colormap,
        vmin=-cap,
        vmax=cap,
    )
    axis.set_yticks(np.arange(len(machinery)))
    axis.set_yticklabels(labels, fontsize=8)
    axis.set_xticks([0], ["log₂(OR)"])
    for row, value in enumerate(values):
        if non_estimable[row]:
            label = "NA"
        elif boundary[row]:
            label = "−∞" if value < 0 else "+∞"
        else:
            label = f"{value:.2f}"
        axis.text(0, row, label, ha="center", va="center", fontsize=7)
    axis.set_title(
        f"{display}-related RBPs in {DATABASE_LABELS[database]} "
        "(* quality-qualified FDR-significant)"
    )
    figure.colorbar(image, ax=axis, label="log₂(odds ratio)", shrink=0.8)
    _save(figure, directory, filename)


def _machinery_labels(machinery: pd.DataFrame) -> list[str]:
    return [
        f"{rbp}{'*' if significant else ''}"
        for rbp, significant in zip(machinery["RBP"], machinery["significant"])
    ]


def plot_loose_results(
    results: pd.DataFrame,
    display: str,
    database: str,
    plot_dir: str | Path,
    top_n: int = 30,
    fdr_threshold: float = 0.05,
) -> None:
    directory = Path(plot_dir)
    qualified = add_significance_calls(
        results,
        odds_column="odds_ratio",
        fdr_column="fdr",
        fdr_threshold=fdr_threshold,
    )
    if "direction" not in qualified:
        qualified["direction"] = qualified["statistical_direction"]
    _top_effects(qualified, display, database, directory, True, top_n)
    _top_effects(qualified, display, database, directory, False, top_n)
    plot_all_rbps(
        qualified,
        display,
        database,
        directory,
        fdr_threshold=fdr_threshold,
    )
    plot_machinery(qualified, display, database, directory)
