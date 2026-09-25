from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from .schema import prepare_results


COLOURS = {
    "enriched": "#D94F45",
    "depleted": "#3478B7",
    "not_significant": "#A9ADB3",
}
PREFIXES = {"oRNAment": "ornament", "ENCORI": "encori", "POSTAR3": "postar3"}


def _empty_plot(path: Path, title: str, message: str, dpi: int) -> None:
    fig, axis = plt.subplots(figsize=(8, 4))
    axis.axis("off")
    axis.set_title(title, fontweight="bold")
    axis.text(0.5, 0.5, message, ha="center", va="center")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _overview_mask(data: pd.DataFrame) -> pd.Series:
    return (
        data["estimable"]
        & ~data["boundary_estimate"]
        & data["mh_or"].gt(0)
        & np.isfinite(data["mh_or"])
        & np.isfinite(data["log2_mh_or"])
    )


def _overview_exclusion_reason(data: pd.DataFrame) -> pd.Series:
    shown = _overview_mask(data)
    return pd.Series(
        np.select(
            [
                shown,
                data["boundary_estimate"] & data["mh_or"].eq(0),
                data["boundary_estimate"] & np.isposinf(data["mh_or"]),
                data["boundary_estimate"],
            ],
            [
                "shown",
                "boundary_zero",
                "boundary_infinite",
                "boundary_other",
            ],
            default="not_estimable",
        ),
        index=data.index,
        dtype="object",
    )


def _overview_log_limit(data: pd.DataFrame, shown: pd.Series) -> float:
    effects = data.loc[shown, "log2_mh_or"].abs()
    if effects.empty:
        return 2.0
    return max(2.0, float(np.ceil(float(effects.max()) * 4) / 4 + 0.5))


def _overview_note(
    data: pd.DataFrame,
    shown: pd.Series,
) -> tuple[str, int, int]:
    boundary = data["boundary_estimate"]
    not_estimable = ~shown & ~boundary
    n_boundary = int(boundary.sum())
    n_missing = int(not_estimable.sum())
    zero = boundary & data["mh_or"].eq(0)
    infinite = boundary & np.isposinf(data["mh_or"])
    other_boundary = boundary & ~zero & ~infinite
    note = (
        f"Finite log₂(OR) bars shown: {int(shown.sum())} of "
        f"{len(data)} RBPs."
    )
    if zero.any():
        count = int(zero.sum())
        note += (
            f" {count} OR=0 boundary estimate"
            f"{'s are' if count != 1 else ' is'} shown as individual "
            "−∞ bars (descriptive; no finite robust FDR estimate)."
        )
    if infinite.any():
        count = int(infinite.sum())
        note += (
            f" {count} OR=∞ boundary estimate"
            f"{'s are' if count != 1 else ' is'} omitted."
        )
    if other_boundary.any():
        note += f" {int(other_boundary.sum())} other boundary test(s) are omitted."
    if n_missing:
        names = ", ".join(data.loc[not_estimable, "RBP"].astype(str))
        note += (
            f" {n_missing} non-estimable RBP"
            f"{'s' if n_missing != 1 else ''} ({names}) "
            f"{'are' if n_missing != 1 else 'is'} omitted."
        )
    if n_boundary or n_missing:
        note += " Exact rows and statuses remain in the source TSV."
    if n_boundary or n_missing:
        related = data[(~shown) & data["is_modification_related"]]
        if not related.empty:
            reasons = _overview_exclusion_reason(data)
            labels = [
                f"{data.loc[index, 'RBP']} "
                f"({reasons.loc[index].replace('_', ' ')})"
                for index in related.index
            ]
            note += (
                "\nPredefined modification-related RBP not shown among the "
                "finite estimates: "
                f"{', '.join(labels)}."
            )
    return note, n_boundary, n_missing


def _all_rbps(
    data: pd.DataFrame, path: Path, title: str, dpi: int
) -> tuple[int, int, int, float, int]:
    shown = _overview_mask(data)
    ordered = data.loc[shown].sort_values(
        ["log2_mh_or", "RBP"],
        ascending=[False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    limit = _overview_log_limit(data, shown)
    note, n_boundary, n_missing = _overview_note(data, shown)
    zero_boundary = data["boundary_estimate"] & data["mh_or"].eq(0)
    boundary_rows = data.loc[zero_boundary].sort_values(
        "RBP", kind="mergesort"
    ).reset_index(drop=True)
    zero_count = len(boundary_rows)
    width = max(11, min(22, (len(ordered) + zero_count) / 14))
    fig, axis = plt.subplots(figsize=(width, 6))
    if ordered.empty and zero_count == 0:
        axis.axis("off")
        axis.text(
            0.5,
            0.55,
            "No finite or OR=0 boundary estimates to display",
            ha="center",
            va="center",
        )
    else:
        axis.bar(
            np.arange(len(ordered)),
            ordered["log2_mh_or"],
            color=ordered["plot_direction"].map(COLOURS),
            edgecolor="none",
            width=1.0,
        )
        if zero_count:
            boundary_x = np.arange(
                len(ordered), len(ordered) + zero_count
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
            label_size = 3.6 if len(data) > 250 else 4.2
            for x_value, row in zip(
                boundary_x, boundary_rows.itertuples(index=False)
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
            axis.set_xlim(-0.5, len(ordered) + zero_count - 0.5)
        elif len(ordered):
            axis.set_xlim(-0.5, len(ordered) - 0.5)
        axis.axhline(0, color="#333333", linewidth=0.8)
        axis.set_ylim(-limit, limit)
    axis.set_ylabel("log₂(Mantel–Haenszel common odds ratio)")
    axis.set_xlabel(
        "Finite RBP estimates ranked by log₂(Mantel–Haenszel odds ratio)"
    )
    axis.set_title(
        f"{title}\n{len(ordered)} finite estimates shown "
        f"of {len(data)} tested",
        fontweight="bold",
    )
    axis.set_xticks([])
    axis.grid(axis="y", color="#E6E6E6", linewidth=0.7)

    candidates = ordered[ordered["plot_significant"]].copy()
    candidates["magnitude"] = candidates["log2_mh_or"].abs()
    for index in candidates.nlargest(12, "magnitude").index:
        y_value = float(ordered.loc[index, "log2_mh_or"])
        offset = 4 if y_value >= 0 else -4
        axis.annotate(
            ordered.loc[index, "RBP"],
            (index, y_value),
            xytext=(1.5, offset),
            textcoords="offset points",
            ha="center",
            va="bottom" if y_value >= 0 else "top",
            rotation=90,
            fontsize=5.2,
            clip_on=False,
        )
    legend = [
        Patch(facecolor=COLOURS[key], edgecolor="none", label=label)
        for key, label in (
            ("enriched", "FDR-significant enriched"),
            ("depleted", "FDR-significant depleted"),
            ("not_significant", "Not significant"),
        )
    ]
    if zero_count:
        legend.append(
            Patch(
                facecolor="#6C5CE7",
                edgecolor="white",
                label="OR = 0 (individual log₂(OR) = −∞ bars)",
            )
        )
    axis.legend(handles=legend, frameon=False, fontsize=8, ncol=2)
    fig.text(0.06, 0.015, note, fontsize=7.4, color="#555555", va="bottom")
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(ordered), n_boundary, n_missing, limit, 0


def _top_direction(
    data: pd.DataFrame,
    direction: str,
    path: Path,
    title: str,
    top_n: int,
    dpi: int,
) -> tuple[int, set[str]]:
    if direction == "enriched":
        eligible = data[
            data["plot_significant"] & data["mh_or"].gt(1)
        ].nlargest(top_n, "log2_mh_or")
    else:
        eligible = data[
            data["plot_significant"] & data["mh_or"].lt(1)
        ].nsmallest(top_n, "log2_mh_or")
    eligible = eligible.sort_values("log2_mh_or", kind="mergesort")
    if eligible.empty:
        _empty_plot(path, title, "No FDR-significant estimable RBPs", dpi)
        return 0, set()
    height = max(4.5, len(eligible) * 0.3)
    fig, axis = plt.subplots(figsize=(9, height))
    axis.barh(
        np.arange(len(eligible)),
        eligible["log2_mh_or"],
        color=COLOURS[direction],
        alpha=0.9,
    )
    axis.set_yticks(np.arange(len(eligible)), eligible["RBP"], fontsize=8)
    axis.axvline(0, color="#333333", linewidth=0.8)
    axis.set_xlabel("log₂ Mantel–Haenszel odds ratio")
    axis.set_title(title, fontweight="bold")
    axis.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(eligible), set(eligible["RBP"])


def _machinery(
    data: pd.DataFrame, path: Path, title: str, dpi: int
) -> tuple[int, set[str]]:
    selected = data[data["is_modification_related"]].sort_values(
        "log2_mh_or", kind="mergesort"
    )
    if selected.empty:
        _empty_plot(path, title, "No predefined modification-related RBPs", dpi)
        return 0, set()
    values = selected["log2_mh_or"].to_numpy(dtype=float)
    finite = np.isfinite(values)
    limit = max(1.0, float(np.nanmax(np.abs(values[finite]))) if finite.any() else 1.0)
    matrix = np.where(finite, values, np.nan)[:, None]
    fig, axis = plt.subplots(figsize=(5.5, max(4.5, len(selected) * 0.3)))
    image = axis.imshow(
        matrix,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
    )
    axis.set_xticks([0], ["Effect"])
    labels = [
        f"{row.RBP}{'*' if row.plot_significant else ''}"
        for row in selected.itertuples()
    ]
    axis.set_yticks(np.arange(len(selected)), labels, fontsize=8)
    axis.set_title(title, fontweight="bold")
    colourbar = fig.colorbar(image, ax=axis, fraction=0.08)
    colourbar.set_label("log₂ MH odds ratio")
    for row_index, is_finite in enumerate(finite):
        if not is_finite:
            axis.text(0, row_index, "NA", ha="center", va="center", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return len(selected), set(selected["RBP"])


def generate_all_plots(
    results: pd.DataFrame,
    output_dir: Path,
    analysis_id: str,
    fdr_threshold: float = 0.05,
    top_n: int = 30,
    dpi: int = 180,
) -> pd.DataFrame:
    data = prepare_results(results, "fdr_within_database", fdr_threshold)
    manifest_rows = []
    for database, group in data.groupby("database", sort=True):
        prefix = PREFIXES.get(database, database.lower())
        plot_dir = output_dir / f"plots_{prefix}"
        plot_dir.mkdir(parents=True, exist_ok=True)
        source = group.copy()
        source["shown_in_all_rbps_overview"] = _overview_mask(source)
        source["all_rbps_overview_exclusion_reason"] = (
            _overview_exclusion_reason(source)
        )
        overview_path = plot_dir / f"{prefix}_all_rbps_overview.png"
        (
            count,
            n_boundary,
            n_not_estimable,
            overview_cap,
            n_overview_capped,
        ) = _all_rbps(
            source,
            overview_path,
            f"{analysis_id}: {database} transcript-region associations",
            dpi,
        )
        overview_counts = {
            "n_rbps_tested": len(source),
            "n_boundary_excluded_from_overview": (
                n_boundary
                - int(
                    (
                        source["boundary_estimate"]
                        & source["mh_or"].eq(0)
                    ).sum()
                )
            ),
            "n_not_estimable_excluded_from_overview": n_not_estimable,
            "all_rbps_overview_effect_scale": (
                "log2_mantel_haenszel_odds_ratio"
            ),
            "all_rbps_overview_display_limit": overview_cap,
            "n_zero_boundary_estimates_displayed_individually": int(
                (
                    source["boundary_estimate"]
                    & source["mh_or"].eq(0)
                ).sum()
            ),
            "n_positive_infinity_estimates_excluded": int(
                (
                    source["boundary_estimate"]
                    & np.isposinf(source["mh_or"])
                ).sum()
            ),
            "n_finite_overview_estimates_capped": n_overview_capped,
        }
        manifest_rows.append(
            {
                "database": database,
                "plot_type": "all_rbps_overview",
                "path": overview_path.relative_to(output_dir).as_posix(),
                "n_selected": count,
                **overview_counts,
            }
        )
        for direction in ("enriched", "depleted"):
            path = plot_dir / f"{prefix}_top_{direction}.png"
            count, rbps = _top_direction(
                source,
                direction,
                path,
                f"{analysis_id}: top {database} {direction} RBPs",
                top_n,
                dpi,
            )
            source[f"selected_top_{direction}"] = source["RBP"].isin(rbps)
            manifest_rows.append(
                {
                    "database": database,
                    "plot_type": f"top_{direction}",
                    "path": path.relative_to(output_dir).as_posix(),
                    "n_selected": count,
                    **overview_counts,
                }
            )
        machinery_path = plot_dir / f"{prefix}_machinery_heatmap.png"
        count, rbps = _machinery(
            source,
            machinery_path,
            f"{analysis_id}: modification-related RBP panel ({database})",
            dpi,
        )
        source["selected_machinery"] = source["RBP"].isin(rbps)
        manifest_rows.append(
            {
                "database": database,
                "plot_type": "machinery_heatmap",
                "path": machinery_path.relative_to(output_dir).as_posix(),
                "n_selected": count,
                **overview_counts,
            }
        )
        source.to_csv(
            plot_dir / "stratified_plot_source_data.tsv", sep="\t", index=False
        )
    manifest = pd.DataFrame(manifest_rows)
    manifest["analysis_id"] = analysis_id
    manifest["fdr_threshold"] = fdr_threshold
    manifest["top_n"] = top_n
    manifest["confidence_intervals_drawn"] = False
    manifest.to_csv(output_dir / "plot_manifest.tsv", sep="\t", index=False)
    return manifest
