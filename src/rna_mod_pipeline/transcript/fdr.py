from __future__ import annotations

import numpy as np
import pandas as pd


def bh_adjust_declared(
    values: pd.Series | np.ndarray | list[float],
) -> np.ndarray:
    """BH correction where non-estimable declared tests count as p=1."""
    pvalues = np.asarray(values, dtype=float)
    missing = ~np.isfinite(pvalues)
    observed = np.where(missing, 1.0, pvalues)
    if np.any((observed < 0) | (observed > 1)):
        raise ValueError("p-values must lie between zero and one")
    order = np.argsort(observed, kind="mergesort")
    ranked = observed[order]
    scaled = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    monotone = np.minimum.accumulate(scaled[::-1])[::-1]
    restored = np.empty(len(ranked), dtype=float)
    restored[order] = np.clip(monotone, 0, 1)
    restored[missing] = np.nan
    return restored


def add_fdr_columns(results: pd.DataFrame) -> pd.DataFrame:
    frame = results.copy()
    frame["fdr_within_analysis"] = bh_adjust_declared(
        frame["gene_cluster_pvalue"]
    )
    frame["cmh_fdr_within_analysis"] = bh_adjust_declared(frame["cmh_pvalue"])
    frame["fdr_within_database"] = np.nan
    frame["cmh_fdr_within_database"] = np.nan
    for indexes in frame.groupby("database", sort=False).groups.values():
        rows = list(indexes)
        frame.loc[rows, "fdr_within_database"] = bh_adjust_declared(
            frame.loc[rows, "gene_cluster_pvalue"]
        )
        frame.loc[rows, "cmh_fdr_within_database"] = bh_adjust_declared(
            frame.loc[rows, "cmh_pvalue"]
        )
    return frame


def add_region_fdr_columns(results: pd.DataFrame) -> pd.DataFrame:
    frame = results.copy()
    frame["fdr_within_database_region"] = np.nan
    frame["cmh_fdr_within_database_region"] = np.nan
    groups = frame.groupby(
        ["database", "region_filter"], sort=False
    ).groups.values()
    for indexes in groups:
        rows = list(indexes)
        frame.loc[rows, "fdr_within_database_region"] = bh_adjust_declared(
            frame.loc[rows, "gene_cluster_pvalue"]
        )
        frame.loc[rows, "cmh_fdr_within_database_region"] = (
            bh_adjust_declared(frame.loc[rows, "cmh_pvalue"])
        )
    return frame
