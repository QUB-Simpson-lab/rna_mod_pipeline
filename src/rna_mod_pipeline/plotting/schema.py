from __future__ import annotations

import numpy as np
import pandas as pd

from ..evidence_quality import add_evidence_quality


REQUIRED_COLUMNS = {
    "database",
    "RBP",
    "mh_or",
    "log2_mh_or",
    "case_overlaps",
    "control_overlaps",
    "estimable",
    "boundary_estimate",
    "is_modification_related",
}


def boolean_series(series: pd.Series) -> pd.Series:
    parsed = series.astype(str).str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )
    if parsed.isna().any():
        raise ValueError(f"Invalid booleans in {series.name}")
    return parsed.astype(bool)


def prepare_results(
    frame: pd.DataFrame,
    fdr_column: str,
    fdr_threshold: float,
) -> pd.DataFrame:
    missing = (REQUIRED_COLUMNS | {fdr_column}) - set(frame.columns)
    if missing:
        raise ValueError(f"Plot input is missing columns: {sorted(missing)}")
    data = frame.copy()
    data["RBP"] = data["RBP"].astype(str)
    data["mh_or"] = pd.to_numeric(data["mh_or"], errors="coerce")
    data["log2_mh_or"] = pd.to_numeric(data["log2_mh_or"], errors="coerce")
    data[fdr_column] = pd.to_numeric(data[fdr_column], errors="coerce")
    data["estimable"] = boolean_series(data["estimable"])
    data["boundary_estimate"] = boolean_series(data["boundary_estimate"])
    data["is_modification_related"] = boolean_series(
        data["is_modification_related"]
    )
    data = add_evidence_quality(
        data,
        odds_column="mh_or",
        fdr_column=fdr_column,
    )
    data["plot_significant"] = (
        data["inference_eligible"]
        & data[fdr_column].notna()
        & data[fdr_column].lt(fdr_threshold)
    )
    data["plot_direction"] = np.select(
        [
            data["plot_significant"] & data["mh_or"].gt(1),
            data["plot_significant"] & data["mh_or"].lt(1),
            data["boundary_estimate"],
            ~data["estimable"],
        ],
        ["enriched", "depleted", "boundary", "not_estimable"],
        default="not_significant",
    )
    return data
