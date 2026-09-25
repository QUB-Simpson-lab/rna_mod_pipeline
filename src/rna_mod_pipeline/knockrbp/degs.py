from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from ..config import canonical_rbp
from ..io import require_directory, require_file


DUPLICATE_AUDIT_COLUMNS = [
    "knocked_down_rbp",
    "dataset_id",
    "gene",
    "n_rows",
    "action",
]


def _truth(value: object) -> bool:
    return str(value).strip().lower() in {"yes", "true", "1"}


def _dataset_from_name(path: Path) -> tuple[str, str]:
    match = re.match(r"(.+)_DataSet_(\d+)_degs\.json$", path.name)
    if not match:
        raise ValueError(f"Unexpected KnockRBP filename: {path.name}")
    return canonical_rbp(match.group(1)).upper(), f"DataSet_{match.group(2)}"


def _filtered_rows(
    path: Path,
    *,
    has_pvalues: bool,
    log2fc_cutoff: float,
    padj_cutoff: float,
) -> list[dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    source = payload.get("data") if isinstance(payload, dict) else payload
    if not isinstance(source, list):
        raise ValueError(f"KnockRBP JSON has no row list: {path}")
    rows = []
    for record in source:
        gene = str(record.get("Genesymbol", "")).strip()
        try:
            log2fc = float(record.get("log2FC"))
        except (TypeError, ValueError):
            continue
        padj = pd.to_numeric(pd.Series([record.get("padj")]), errors="coerce").iloc[0]
        if not gene or abs(log2fc) <= log2fc_cutoff:
            continue
        if has_pvalues and (pd.isna(padj) or float(padj) >= padj_cutoff):
            continue
        rows.append(
            {
                "gene": gene,
                "log2fc": log2fc,
                "padj": float(padj) if pd.notna(padj) else np.nan,
                "direction": "up" if log2fc > 0 else "down",
            }
        )
    return rows


def _resolve_duplicates(
    rows: list[dict[str, object]],
    rbp: str,
    dataset_id: str,
) -> tuple[dict[str, dict[str, object]], list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["gene"]), []).append(row)
    genes = {}
    audit = []
    for gene, candidates in sorted(grouped.items()):
        directions = {str(row["direction"]) for row in candidates}
        if len(directions) > 1:
            audit.append(
                {
                    "knocked_down_rbp": rbp,
                    "dataset_id": dataset_id,
                    "gene": gene,
                    "n_rows": len(candidates),
                    "action": "excluded_conflicting_directions",
                }
            )
            continue
        selected = min(
            candidates,
            key=lambda row: (
                float(row["padj"]) if pd.notna(row["padj"]) else float("inf"),
                -abs(float(row["log2fc"])),
            ),
        )
        genes[gene] = selected
        if len(candidates) > 1:
            audit.append(
                {
                    "knocked_down_rbp": rbp,
                    "dataset_id": dataset_id,
                    "gene": gene,
                    "n_rows": len(candidates),
                    "action": "collapsed_same_direction",
                }
            )
    return genes, audit


def load_and_resolve_degs(
    data_dir: str | Path,
    metadata_path: str | Path,
    *,
    allowed_cell_lines: Iterable[str] = ("MDA-MB-231", "MDA-MB-231-LM2"),
    dataset_ids: Iterable[str] | None = None,
    log2fc_cutoff: float = 0.5,
    padj_cutoff: float = 0.05,
) -> tuple[dict[str, dict[str, object]], pd.DataFrame]:
    directory = require_directory(data_dir, "KnockRBP data directory")
    metadata = pd.read_csv(require_file(metadata_path, "KnockRBP metadata"), sep="\t")
    if not {"dataset_id", "rbp", "cell_line", "has_pvalues"}.issubset(metadata):
        raise ValueError("KnockRBP metadata lacks required columns")
    metadata = metadata.set_index("dataset_id", verify_integrity=True)
    allowed = set(allowed_cell_lines)
    requested = set(dataset_ids or [])
    selected: dict[str, dict[str, object]] = {}
    audit_rows = []
    for path in sorted(directory.glob("*_degs.json")):
        rbp, dataset_id = _dataset_from_name(path)
        if dataset_id not in metadata.index:
            raise ValueError(f"Metadata missing for {dataset_id}")
        info = metadata.loc[dataset_id]
        metadata_rbp = canonical_rbp(info["rbp"]).upper()
        if rbp != metadata_rbp:
            raise ValueError(
                f"RBP mismatch for {dataset_id}: filename={rbp}, metadata={metadata_rbp}"
            )
        if requested:
            if dataset_id not in requested:
                continue
        elif str(info["cell_line"]).strip() not in allowed:
            continue
        has_pvalues = _truth(info["has_pvalues"])
        raw = _filtered_rows(
            path,
            has_pvalues=has_pvalues,
            log2fc_cutoff=log2fc_cutoff,
            padj_cutoff=padj_cutoff,
        )
        genes, audit = _resolve_duplicates(raw, rbp, dataset_id)
        audit_rows.extend(audit)
        if dataset_id in selected:
            raise ValueError(f"Duplicate selected dataset ID: {dataset_id}")
        selected[dataset_id] = {
            "rbp": rbp,
            "dataset_id": dataset_id,
            "cell_line": str(info["cell_line"]),
            "has_pvalues": has_pvalues,
            "assay": str(info.get("assay", "")),
            "perturbation": str(info.get("perturbation", "")),
            "treated_n": str(info.get("treated_n", "")),
            "control_n": str(info.get("control_n", "")),
            "genes": genes,
            "n_filtered_rows_before_resolution": len(raw),
            "n_resolved_degs": len(genes),
            "source_file": str(path.resolve()),
        }
    if requested - set(selected):
        missing = sorted(requested - set(selected))
        raise ValueError(f"Requested dataset files were not loaded: {missing}")
    return selected, pd.DataFrame(audit_rows, columns=DUPLICATE_AUDIT_COLUMNS)
