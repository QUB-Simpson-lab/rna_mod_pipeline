from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..schemas import require_columns, validate_modification_codes, validate_unique
from .profiles import DatasetProfile


SITE_KEYS = ["chrom", "start", "end", "strand"]
SITE_REQUIRED = {*SITE_KEYS, "mod_code", "Nvalid_cov", "fraction_modified"}


def _finite_pairs(first: pd.Series, second: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    x = pd.to_numeric(first, errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(second, errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    return x[keep], y[keep]


def _correlation(first: pd.Series, second: pd.Series, method: str) -> float:
    x, y = _finite_pairs(first, second)
    if len(x) < 2:
        return math.nan
    if method == "spearman":
        x = pd.Series(x).rank(method="average").to_numpy()
        y = pd.Series(y).rank(method="average").to_numpy()
    if np.std(x) == 0 or np.std(y) == 0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def _concordance(first: pd.Series, second: pd.Series) -> float:
    x, y = _finite_pairs(first, second)
    if len(x) < 2:
        return math.nan
    denominator = np.var(x) + np.var(y) + (np.mean(x) - np.mean(y)) ** 2
    if denominator == 0:
        return math.nan
    covariance = np.mean((x - np.mean(x)) * (y - np.mean(y)))
    return float(2 * covariance / denominator)


def _load_sites(path: Path, profile: DatasetProfile, stage: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = set(SITE_REQUIRED)
    if stage == "metagene":
        required.update({"region", "metagene_pos"})
    require_columns(table, required, f"{profile.dataset_id} {stage} sites")
    validate_unique(table, SITE_KEYS, f"{profile.dataset_id} {stage} sites")
    validate_modification_codes(
        table,
        profile.modification,
        f"{profile.dataset_id} {stage} sites",
    )
    for column in ("start", "end", "Nvalid_cov", "fraction_modified"):
        numeric = pd.to_numeric(table[column], errors="coerce")
        if numeric.isna().any():
            raise ValueError(
                f"{profile.dataset_id} {stage} sites has non-numeric {column}"
            )
        table[column] = numeric
    for column in ("start", "end", "Nvalid_cov"):
        if not np.allclose(table[column], np.round(table[column]), rtol=0, atol=0):
            raise ValueError(
                f"{profile.dataset_id} {stage} sites has non-integer {column}"
            )
        table[column] = table[column].astype(np.int64)
    if (
        table["start"].lt(0).any()
        or table["end"].le(table["start"]).any()
        or table["Nvalid_cov"].lt(0).any()
        or table["fraction_modified"].lt(0).any()
        or table["fraction_modified"].gt(100).any()
    ):
        raise ValueError(f"{profile.dataset_id} {stage} sites has invalid values")
    strands = table["strand"].astype(str).str.strip()
    if table["strand"].isna().any() or strands.eq("").any():
        raise ValueError(f"{profile.dataset_id} {stage} sites has invalid strand")
    table["strand"] = strands
    columns = [
        *SITE_KEYS,
        "mod_code",
        "Nvalid_cov",
        "fraction_modified",
        *(
            column
            for column in ("transcript_id", "gene_name", "region", "metagene_pos")
            if column in table
        ),
    ]
    return table[columns].copy()


def _summary(stage: str, merged: pd.DataFrame) -> dict[str, object]:
    shared = merged[merged["membership"].eq("shared")]
    n_a = int(merged["membership"].isin(["dataset_a_only", "shared"]).sum())
    n_b = int(merged["membership"].isin(["dataset_b_only", "shared"]).sum())
    n_shared = len(shared)
    n_union = len(merged)
    row: dict[str, object] = {
        "stage": stage,
        "n_sites_a": n_a,
        "n_sites_b": n_b,
        "n_shared_sites": n_shared,
        "n_union_sites": n_union,
        "n_sites_a_only": int(merged["membership"].eq("dataset_a_only").sum()),
        "n_sites_b_only": int(merged["membership"].eq("dataset_b_only").sum()),
        "jaccard": n_shared / n_union if n_union else math.nan,
        "containment_in_a": n_shared / n_a if n_a else math.nan,
        "containment_in_b": n_shared / n_b if n_b else math.nan,
        "shared_fraction_pearson": _correlation(
            shared["fraction_modified_a"], shared["fraction_modified_b"], "pearson"
        ),
        "shared_fraction_spearman": _correlation(
            shared["fraction_modified_a"], shared["fraction_modified_b"], "spearman"
        ),
        "shared_fraction_lins_ccc": _concordance(
            shared["fraction_modified_a"], shared["fraction_modified_b"]
        ),
        "shared_fraction_median_absolute_difference": (
            float(shared["fraction_absolute_difference"].median())
            if n_shared
            else math.nan
        ),
        "shared_coverage_pearson": _correlation(
            shared["Nvalid_cov_a"], shared["Nvalid_cov_b"], "pearson"
        ),
        "shared_coverage_spearman": _correlation(
            shared["Nvalid_cov_a"], shared["Nvalid_cov_b"], "spearman"
        ),
        "shared_coverage_lins_ccc": _concordance(
            shared["Nvalid_cov_a"], shared["Nvalid_cov_b"]
        ),
        "shared_coverage_median_absolute_difference": (
            float(shared["coverage_absolute_difference"].median())
            if n_shared
            else math.nan
        ),
    }
    if {"region_a", "region_b"}.issubset(shared):
        comparable = shared["region_a"].notna() & shared["region_b"].notna()
        row["n_shared_region_comparable"] = int(comparable.sum())
        row["shared_region_agreement_rate"] = (
            float(
                shared.loc[comparable, "region_a"]
                .astype(str)
                .eq(shared.loc[comparable, "region_b"].astype(str))
                .mean()
            )
            if comparable.any()
            else math.nan
        )
    else:
        row["n_shared_region_comparable"] = 0
        row["shared_region_agreement_rate"] = math.nan
    return row


def compare_site_tables(
    first: DatasetProfile,
    second: DatasetProfile,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    memberships = []
    summaries = []
    for stage, first_path, second_path in (
        ("filtered", first.filtered_sites, second.filtered_sites),
        ("metagene", first.metagene_sites, second.metagene_sites),
    ):
        table_a = _load_sites(first_path, first, stage)
        table_b = _load_sites(second_path, second, stage)
        merged = table_a.merge(
            table_b,
            on=SITE_KEYS,
            how="outer",
            suffixes=("_a", "_b"),
            indicator=True,
            validate="one_to_one",
        )
        merged["membership"] = merged.pop("_merge").map(
            {
                "left_only": "dataset_a_only",
                "right_only": "dataset_b_only",
                "both": "shared",
            }
        )
        merged["fraction_delta_b_minus_a"] = (
            merged["fraction_modified_b"] - merged["fraction_modified_a"]
        )
        merged["fraction_absolute_difference"] = merged[
            "fraction_delta_b_minus_a"
        ].abs()
        merged["coverage_delta_b_minus_a"] = (
            merged["Nvalid_cov_b"] - merged["Nvalid_cov_a"]
        )
        merged["coverage_absolute_difference"] = merged[
            "coverage_delta_b_minus_a"
        ].abs()
        merged.insert(0, "stage", stage)
        memberships.append(merged)
        summaries.append(_summary(stage, merged))
    return pd.concat(memberships, ignore_index=True), pd.DataFrame(summaries)
