from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..config import modification as modification_config
from ..io import open_text, require_file


@dataclass(frozen=True)
class BedMethylInspection:
    sampled_rows: int
    complete_file: bool
    modification_codes: tuple[str, ...]
    chromosomes: tuple[str, ...]
    minimum_coverage: int
    maximum_coverage: int
    minimum_fraction: float
    maximum_fraction: float

    def summary(self) -> str:
        scope = "all" if self.complete_file else "first"
        return (
            f"Valid bedMethyl sample: {scope} {self.sampled_rows:,} rows; "
            f"code(s) {', '.join(self.modification_codes)}; "
            f"{len(self.chromosomes):,} chromosome(s); coverage "
            f"{self.minimum_coverage:,}–{self.maximum_coverage:,}; "
            f"modification {self.minimum_fraction:g}–{self.maximum_fraction:g}%."
        )


def inspect_bedmethyl(
    path: str | Path,
    modification: str,
    *,
    maximum_rows: int = 50_000,
) -> BedMethylInspection:
    if maximum_rows < 1:
        raise ValueError("maximum_rows must be positive")
    source = require_file(path, "bedMethyl file")
    expected = {
        value.casefold()
        for value in modification_config(modification).expected_mod_codes
    }
    codes: set[str] = set()
    chromosomes: set[str] = set()
    coverages: list[int] = []
    fractions: list[float] = []
    complete = True
    with open_text(source) as handle:
        for line_number, line in enumerate(handle, 1):
            if line_number > maximum_rows:
                complete = False
                break
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 18:
                raise ValueError(
                    f"Row {line_number} has {len(fields)} columns; expected 18."
                )
            try:
                start, end = int(fields[1]), int(fields[2])
                coverage = int(fields[9])
                fraction = float(fields[10])
            except ValueError as exc:
                raise ValueError(
                    f"Row {line_number} contains an invalid numeric value."
                ) from exc
            if start < 0 or end <= start:
                raise ValueError(
                    f"Row {line_number} has invalid coordinates {start}–{end}."
                )
            if fields[5] not in {"+", "-", "."}:
                raise ValueError(
                    f"Row {line_number} has invalid strand {fields[5]!r}."
                )
            if coverage < 0 or not 0 <= fraction <= 100:
                raise ValueError(
                    f"Row {line_number} has invalid coverage or fraction."
                )
            chromosomes.add(fields[0])
            codes.add(fields[3].strip().casefold())
            coverages.append(coverage)
            fractions.append(fraction)
    if not coverages:
        raise ValueError("The bedMethyl file is empty.")
    unexpected = codes - expected
    if unexpected:
        raise ValueError(
            f"Observed modification code(s) {sorted(unexpected)} do not match "
            f"{modification}; expected one of {sorted(expected)}."
        )
    return BedMethylInspection(
        sampled_rows=len(coverages),
        complete_file=complete,
        modification_codes=tuple(sorted(codes)),
        chromosomes=tuple(sorted(chromosomes)),
        minimum_coverage=min(coverages),
        maximum_coverage=max(coverages),
        minimum_fraction=min(fractions),
        maximum_fraction=max(fractions),
    )
