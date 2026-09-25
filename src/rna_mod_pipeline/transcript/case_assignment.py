from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from rna_mod_pipeline.config import MODIFICATIONS

from .transcripts import PrimaryTranscriptIndex, TranscriptModel


def validate_case_table(
    path: Path,
    modification: str,
    coverage_min: int,
    case_min_fraction: float,
    chromosomes: Optional[set[str]] = None,
) -> tuple[pd.DataFrame, set[tuple[str, int, str]], Counter]:
    sites = pd.read_csv(path, sep="\t", low_memory=False)
    required = {
        "chrom",
        "start",
        "end",
        "strand",
        "mod_code",
        "Nvalid_cov",
        "fraction_modified",
        "transcript_id",
        "gene_name",
        "region",
    }
    missing = required - set(sites.columns)
    if missing:
        raise ValueError(f"High-called site table is missing: {sorted(missing)}")

    flow = Counter(input_high_called_rows=len(sites))
    numeric = ["start", "end", "Nvalid_cov", "fraction_modified"]
    for column in numeric:
        sites[column] = pd.to_numeric(sites[column], errors="coerce")
    values = sites[numeric].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("High-called table has invalid numeric values")
    noninteger = (
        sites["start"].ne(np.floor(sites["start"]))
        | sites["end"].ne(np.floor(sites["end"]))
        | sites["Nvalid_cov"].ne(np.floor(sites["Nvalid_cov"]))
    )
    if noninteger.any():
        raise ValueError("High-called coordinates and coverage must be integers")

    sites["start"] = sites["start"].astype(np.int64)
    sites["end"] = sites["end"].astype(np.int64)
    sites["Nvalid_cov"] = sites["Nvalid_cov"].astype(np.int64)
    sites["chrom"] = sites["chrom"].astype(str).str.strip()
    sites["strand"] = sites["strand"].astype(str).str.strip()
    sites["mod_code"] = sites["mod_code"].astype(str).str.strip().str.lower()
    invalid = (
        sites["start"].lt(0)
        | sites["end"].ne(sites["start"] + 1)
        | ~sites["strand"].isin({"+", "-"})
        | sites["chrom"].eq("")
        | ~sites["fraction_modified"].between(0, 100)
    )
    if invalid.any():
        raise ValueError(f"{int(invalid.sum())} high-called rows are invalid")

    base_key = "m6a" if modification == "m6a_rep2" else modification
    expected_codes = MODIFICATIONS[base_key].expected_mod_codes
    correct_code = sites["mod_code"].isin(expected_codes)
    flow["unexpected_mod_code_rows"] = int((~correct_code).sum())
    if not correct_code.all():
        observed = sorted(sites.loc[~correct_code, "mod_code"].unique())
        raise ValueError(
            "High-called site table contains modification codes outside the "
            f"selected analysis: {observed}; expected {sorted(expected_codes)}"
        )
    eligible = (
        sites["Nvalid_cov"].ge(coverage_min)
        & sites["fraction_modified"].ge(case_min_fraction)
    )
    flow["below_case_threshold_rows"] = int((~eligible).sum())
    sites = sites[eligible].copy()
    if chromosomes:
        sites = sites[sites["chrom"].isin(chromosomes)].copy()
    flow["after_threshold_and_chromosome_filters"] = len(sites)

    all_high = set(
        zip(
            sites["chrom"].astype(str),
            sites["start"].astype(int),
            sites["strand"].astype(str),
        )
    )
    flow["high_coordinates_excluded_from_controls"] = len(all_high)
    mapped = (
        sites["transcript_id"].notna()
        & sites["gene_name"].notna()
        & sites["region"].isin({"5UTR", "CDS", "3UTR"})
    )
    flow["originally_transcript_region_mapped"] = int(mapped.sum())
    flow["originally_unmapped_or_noncanonical_region"] = int((~mapped).sum())

    key = ["chrom", "start", "strand"]
    duplicated = sites[sites.duplicated(key, keep=False)]
    if not duplicated.empty:
        check = [
            "end",
            "mod_code",
            "Nvalid_cov",
            "fraction_modified",
            "transcript_id",
            "gene_name",
            "region",
        ]
        conflicts = sum(
            any(group[column].nunique(dropna=False) > 1 for column in check)
            for _, group in duplicated.groupby(key, sort=False)
        )
        if conflicts:
            raise ValueError(
                f"{conflicts} duplicated high-call coordinates conflict"
            )
        flow["exact_duplicate_high_called_rows_removed"] = int(
            sites.duplicated(key).sum()
        )
        sites = sites.drop_duplicates(key, keep="first")
    sites = sites.sort_values(key, kind="mergesort").reset_index(drop=True)
    flow["deduplicated_threshold_high_called_rows"] = len(sites)
    return sites, all_high, flow


def assign_cases(
    sites: pd.DataFrame,
    index: PrimaryTranscriptIndex,
    models: dict[str, TranscriptModel],
    max_cases: int | None,
    seed: int,
) -> tuple[pd.DataFrame, Counter]:
    rows: list[dict] = []
    flow = Counter()
    for row in sites.itertuples(index=False):
        assigned = index.query(str(row.chrom), str(row.strand), int(row.start))
        if assigned is None:
            flow["primary_assignment_missing"] += 1
            continue
        transcript_id, region = assigned
        model = models[transcript_id]
        original_transcript = (
            str(row.transcript_id) if pd.notna(row.transcript_id) else ""
        )
        original_gene = str(row.gene_name) if pd.notna(row.gene_name) else ""
        original_region = str(row.region) if pd.notna(row.region) else ""
        phase1_mapped = (
            bool(original_transcript)
            and bool(original_gene)
            and original_region in {"5UTR", "CDS", "3UTR"}
        )
        flow["primary_assignment_available"] += 1
        if phase1_mapped and (transcript_id, region) == (
            original_transcript,
            original_region,
        ):
            flow["agrees_with_phase1_assignment"] += 1
        elif phase1_mapped:
            flow["differs_from_phase1_assignment"] += 1
        else:
            flow["phase1_unmapped_newly_assigned"] += 1
        record = row._asdict()
        record.update(
            {
                "phase1_transcript_id": original_transcript,
                "phase1_gene_name": original_gene,
                "phase1_region": original_region,
                "transcript_id": transcript_id,
                "gene_id": model.gene_id,
                "gene_name": model.gene_name,
                "region": region,
            }
        )
        rows.append(record)

    assigned = pd.DataFrame(rows)
    if assigned.empty:
        raise RuntimeError("No high-called sites received a primary transcript")
    if max_cases is not None and len(assigned) > max_cases:
        before = len(assigned)
        assigned = assigned.sample(n=max_cases, random_state=seed, replace=False)
        flow["omitted_by_max_cases"] = before - len(assigned)
    assigned = assigned.sort_values(
        ["chrom", "start", "strand"], kind="mergesort"
    ).reset_index(drop=True)
    flow["retained_assigned_cases"] = len(assigned)
    return assigned, flow
