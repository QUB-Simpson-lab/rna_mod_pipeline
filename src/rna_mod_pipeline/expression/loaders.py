from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..io import require_file


@dataclass(frozen=True)
class ExpressionProfile:
    table: pd.DataFrame
    source: str
    cell_line: str
    value_column: str
    threshold: float
    absence_interpretation: str


def load_nanopore(
    path: str | Path,
    *,
    value_column: str = "cpm_0h",
    count_column: str = "count_0h",
    cell_line: str = "MDA-MB-231",
) -> ExpressionProfile:
    frame = pd.read_csv(require_file(path, "Nanopore expression table"), sep="\t")
    required = {"gene", value_column, count_column}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Nanopore expression table is missing: {sorted(missing)}")
    table = frame[["gene", value_column, count_column]].copy()
    table["gene"] = table["gene"].astype(str).str.strip()
    table[value_column] = pd.to_numeric(table[value_column], errors="coerce")
    table[count_column] = pd.to_numeric(table[count_column], errors="coerce")
    table = (
        table.groupby("gene", as_index=False)
        .agg({value_column: "sum", count_column: "sum"})
    )
    table["expression_value"] = table[value_column]
    table["detection_metric"] = table[count_column]
    table["expression_status"] = np.where(
        table[count_column] > 0, "detected", "not_detected_low_depth_uninformative"
    )
    return ExpressionProfile(
        table=table,
        source="nanopore_drs",
        cell_line=cell_line,
        value_column=value_column,
        threshold=0.0,
        absence_interpretation="A zero count is not evidence that the RBP is unexpressed.",
    )


def load_depmap(
    path: str | Path,
    *,
    model_id: str = "ACH-000768",
    threshold: float = 1.0,
) -> ExpressionProfile:
    frame = pd.read_csv(require_file(path, "DepMap expression table"))
    required = {"ModelID"}
    if not required.issubset(frame.columns):
        raise ValueError("DepMap table must contain ModelID")
    selected = frame[frame["ModelID"].astype(str).eq(model_id)]
    if "IsDefaultEntryForModel" in selected:
        defaults = selected[
            selected["IsDefaultEntryForModel"].astype(str).str.lower().isin(
                {"yes", "true", "1"}
            )
        ]
        if len(defaults):
            selected = defaults
    if len(selected) != 1:
        raise ValueError(f"Expected one DepMap row for {model_id}; found {len(selected)}")

    metadata = {
        "",
        "Unnamed: 0",
        "SequencingID",
        "ModelConditionID",
        "ModelID",
        "IsDefaultEntryForMC",
        "IsDefaultEntryForModel",
    }
    row = selected.iloc[0]
    records = []
    for column in frame.columns:
        if column in metadata:
            continue
        match = re.match(r"^(.+?)\s*\(\d+\)$", str(column))
        symbol = (match.group(1) if match else str(column)).strip()
        value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
        if pd.notna(value):
            records.append((symbol, max(float(2**value - 1), 0.0)))
    table = pd.DataFrame(records, columns=["gene", "expression_value"])
    if table["gene"].duplicated().any():
        table = table.groupby("gene", as_index=False)["expression_value"].max()
    table["detection_metric"] = table["expression_value"]
    table["expression_status"] = np.where(
        table["expression_value"] >= threshold, "expressed", "below_tpm_threshold"
    )
    return ExpressionProfile(
        table=table,
        source="depmap_illumina",
        cell_line=model_id,
        value_column="TPM",
        threshold=threshold,
        absence_interpretation=f"Below {threshold:g} TPM is classified as not expressed.",
    )
