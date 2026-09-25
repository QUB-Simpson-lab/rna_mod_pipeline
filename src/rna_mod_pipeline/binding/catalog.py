from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..config import canonical_rbp
from ..io import require_file


EXPECTED_CATALOG_COUNTS = {
    "ornament": 133,
    "encori": 281,
    "postar3": 216,
}


@dataclass(frozen=True)
class CatalogEntry:
    version: str
    database: str
    rbp_name: str
    canonical_rbp: str
    source_names: tuple[str, ...]


def load_catalog(
    path: str | Path,
    database: str,
    validate_count: bool = True,
) -> list[CatalogEntry]:
    """Load one versioned, validated database panel."""
    if database not in EXPECTED_CATALOG_COUNTS:
        raise ValueError(f"Unsupported binding database: {database}")

    frame = pd.read_csv(require_file(path, "RBP catalog"), sep="\t", dtype=str)
    required = {
        "catalog_version",
        "database",
        "rbp_name",
        "canonical_rbp",
        "source_names",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"RBP catalog is missing columns: {sorted(missing)}")

    selected = frame.loc[frame["database"].eq(database)].copy()
    if selected.empty:
        raise ValueError(f"RBP catalog has no rows for {database}")
    if selected["catalog_version"].nunique() != 1:
        raise ValueError(f"{database} must have exactly one catalog version")
    if selected["rbp_name"].duplicated().any():
        raise ValueError(f"{database} catalog contains duplicate RBP names")
    if selected["canonical_rbp"].duplicated().any():
        duplicates = sorted(
            selected.loc[
                selected["canonical_rbp"].duplicated(keep=False),
                "canonical_rbp",
            ].unique()
        )
        raise ValueError(
            f"{database} catalog contains duplicate canonical RBPs: {duplicates}"
        )

    entries: list[CatalogEntry] = []
    for row in selected.itertuples(index=False):
        rbp_name = str(row.rbp_name).strip()
        canonical = str(row.canonical_rbp).strip().upper()
        expected_canonical = canonical_rbp(rbp_name).upper()
        if canonical != expected_canonical:
            raise ValueError(
                f"{database}/{rbp_name} records canonical_rbp={canonical}; "
                f"expected {expected_canonical}"
            )
        source_field = "" if pd.isna(row.source_names) else str(row.source_names)
        sources = tuple(
            source.strip() for source in source_field.split(";") if source.strip()
        )
        if not rbp_name or not sources:
            raise ValueError(f"{database} catalog contains an empty RBP/source entry")
        entries.append(
            CatalogEntry(
                version=str(row.catalog_version),
                database=database,
                rbp_name=rbp_name,
                canonical_rbp=canonical,
                source_names=sources,
            )
        )

    expected = EXPECTED_CATALOG_COUNTS[database]
    if validate_count and len(entries) != expected:
        raise ValueError(
            f"{database} catalog has {len(entries)} RBPs; expected {expected}"
        )
    return entries
