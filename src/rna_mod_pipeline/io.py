from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd

from .config import BEDMETHYL_COLUMNS, FILTERED_COLUMNS


def require_file(path: str | Path, label: str = "input") -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{label} not found: {resolved}")
    return resolved


def require_directory(path: str | Path, label: str = "directory") -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_dir():
        raise FileNotFoundError(f"{label} not found: {resolved}")
    return resolved


def open_text(path: str | Path):
    resolved = Path(path)
    return gzip.open(resolved, "rt") if resolved.suffix == ".gz" else open(resolved, "rt")


def bedmethyl_chunks(
    path: str | Path,
    chunk_size: int = 2_000_000,
) -> Iterable[pd.DataFrame]:
    return pd.read_csv(
        require_file(path, "bedMethyl file"),
        sep="\t",
        names=BEDMETHYL_COLUMNS,
        header=None,
        chunksize=chunk_size,
        low_memory=False,
    )


def iter_bedmethyl(
    path: str | Path,
) -> Iterator[tuple[str, int, str, str, int, float]]:
    with open_text(require_file(path, "bedMethyl file")) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"Malformed bedMethyl row {line_number}")
            yield (
                fields[0],
                int(fields[1]),
                fields[5],
                fields[3],
                int(fields[9]),
                float(fields[10]),
            )


def load_sites(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(require_file(path, "site table"), sep="\t")
    required = {"chrom", "start"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Site table is missing columns: {sorted(missing)}")
    return frame


def write_filtered_sites(frame: pd.DataFrame, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame.loc[:, FILTERED_COLUMNS].to_csv(destination, sep="\t", index=False)


def write_tsv(frame: pd.DataFrame, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(destination, sep="\t", index=False)
