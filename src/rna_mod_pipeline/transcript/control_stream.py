from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from rna_mod_pipeline.config import MODIFICATIONS
from rna_mod_pipeline.io import open_text

from .transcripts import PrimaryTranscriptIndex


@dataclass
class OpportunityBuffers:
    chroms: list[str]
    positions: list[int]
    strands: list[str]
    coverages: list[int]
    fractions: list[float]
    outcomes: list[bool]
    stratum_codes: list[int]
    case_site_ids: list[int]
    case_records: list[dict]
    selected_expected: dict[tuple[str, int, str], tuple[int, float, str]]


def seed_case_observations(assigned, code_by_key) -> OpportunityBuffers:
    buffers = OpportunityBuffers([], [], [], [], [], [], [], [], [], {})
    for case_id, row in enumerate(assigned.itertuples(index=False)):
        key = (str(row.transcript_id), str(row.region))
        buffers.chroms.append(str(row.chrom))
        buffers.positions.append(int(row.start))
        buffers.strands.append(str(row.strand))
        buffers.coverages.append(int(row.Nvalid_cov))
        buffers.fractions.append(float(row.fraction_modified))
        buffers.outcomes.append(True)
        buffers.stratum_codes.append(code_by_key[key])
        buffers.case_site_ids.append(case_id)
        coordinate = (str(row.chrom), int(row.start), str(row.strand))
        buffers.selected_expected[coordinate] = (
            int(row.Nvalid_cov),
            float(row.fraction_modified),
            str(row.mod_code).lower(),
        )
        buffers.case_records.append(
            {
                "case_site_id": case_id,
                "chrom": row.chrom,
                "start": int(row.start),
                "end": int(row.end),
                "strand": row.strand,
                "Nvalid_cov": int(row.Nvalid_cov),
                "fraction_modified": float(row.fraction_modified),
                "transcript_id": row.transcript_id,
                "gene_id": row.gene_id,
                "gene_name": row.gene_name,
                "region": row.region,
                "phase1_transcript_id": row.phase1_transcript_id,
                "phase1_gene_name": row.phase1_gene_name,
                "phase1_region": row.phase1_region,
            }
        )
    return buffers


def stream_low_fraction_controls(
    bedmethyl_path: Path,
    buffers: OpportunityBuffers,
    all_expected: dict[tuple[str, int, str], tuple[int, float, str]],
    all_high: set[tuple[str, int, str]],
    modification: str,
    coverage_min: int,
    case_min_fraction: float,
    background_max_fraction: float,
    chromosomes: set[str] | None,
    max_bed_rows: int | None,
    transcript_index: PrimaryTranscriptIndex,
    target_strata: set[tuple[str, str]],
    code_by_key: dict[tuple[str, str], int],
) -> tuple[Counter, bool]:
    flow = Counter()
    verified_selected: set[tuple[str, int, str]] = set()
    verified_all: set[tuple[str, int, str]] = set()
    high_only_examples: list[tuple[str, int, str]] = []
    base_key = "m6a" if modification == "m6a_rep2" else modification
    expected_codes = MODIFICATIONS[base_key].expected_mod_codes
    last_start: dict[tuple[str, str], int] = {}
    last_metadata: dict[tuple[str, str], tuple[int, int, float, str]] = {}
    partial = False

    with open_text(bedmethyl_path) as handle:
        for line_number, line in enumerate(handle, 1):
            if max_bed_rows is not None and line_number > max_bed_rows:
                partial = True
                flow["stopped_at_max_bed_rows"] = 1
                break
            flow["bedmethyl_rows_read"] += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= 10:
                flow["malformed_bedmethyl_rows"] += 1
                continue
            chrom = fields[0]
            if chromosomes and chrom not in chromosomes:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
                coverage, fraction = int(fields[9]), float(fields[10])
            except ValueError:
                flow["malformed_numeric_bedmethyl_rows"] += 1
                continue
            strand, mod_code = fields[5], fields[3].strip().lower()
            if (
                start < 0
                or end != start + 1
                or coverage < 0
                or not math.isfinite(fraction)
                or not 0 <= fraction <= 100
                or strand not in {"+", "-"}
            ):
                flow["invalid_bedmethyl_rows"] += 1
                continue
            if mod_code not in expected_codes:
                flow["unexpected_mod_code_rows_in_bedmethyl"] += 1
                raise ValueError(
                    f"Unexpected modification code {mod_code!r} in bedMethyl "
                    f"row {line_number}; expected {sorted(expected_codes)}"
                )
            sort_key = (chrom, strand)
            previous = last_start.get(sort_key)
            metadata = (end, coverage, fraction, mod_code)
            if previous is not None:
                if start < previous:
                    raise ValueError(
                        f"bedMethyl is unsorted within {chrom}/{strand}"
                    )
                if start == previous:
                    if last_metadata.get(sort_key) != metadata:
                        raise ValueError(
                            f"Conflicting bedMethyl duplicate at "
                            f"{chrom}:{start}:{strand}"
                        )
                    flow[
                        "exact_duplicate_expected_code_coordinate_rows_skipped"
                    ] += 1
                    continue
            last_start[sort_key] = start
            last_metadata[sort_key] = metadata
            coordinate = (chrom, start, strand)

            for expected, verified, metric in (
                (
                    all_expected,
                    verified_all,
                    "all_case_bedmethyl_metadata_mismatch",
                ),
                (
                    buffers.selected_expected,
                    verified_selected,
                    "case_bedmethyl_metadata_mismatch",
                ),
            ):
                if coordinate not in expected:
                    continue
                expected_cov, expected_fraction, expected_code = expected[
                    coordinate
                ]
                if (
                    coverage == expected_cov
                    and math.isclose(
                        fraction, expected_fraction, rel_tol=0, abs_tol=1e-9
                    )
                    and mod_code == expected_code
                ):
                    verified.add(coordinate)
                else:
                    flow[metric] += 1

            if coverage >= coverage_min and fraction >= case_min_fraction:
                flow["bedmethyl_rows_meeting_case_threshold"] += 1
                if coordinate not in all_high:
                    flow[
                        "bedmethyl_high_coordinates_missing_from_site_table"
                    ] += 1
                    if len(high_only_examples) < 5:
                        high_only_examples.append(coordinate)

            if coverage < coverage_min:
                continue
            flow["coverage_eligible_rows"] += 1
            if fraction >= background_max_fraction:
                continue
            flow["low_fraction_eligible_rows"] += 1
            if coordinate in all_high:
                flow["excluded_high_called_coordinates"] += 1
                continue
            assignment = transcript_index.query(chrom, strand, start)
            if assignment is None:
                flow["low_called_outside_complete_transcripts"] += 1
                continue
            if assignment not in target_strata:
                flow["low_called_outside_case_strata"] += 1
                continue
            buffers.chroms.append(chrom)
            buffers.positions.append(start)
            buffers.strands.append(strand)
            buffers.coverages.append(coverage)
            buffers.fractions.append(fraction)
            buffers.outcomes.append(False)
            buffers.stratum_codes.append(code_by_key[assignment])
            flow["retained_low_called_controls"] += 1

    missing_selected = set(buffers.selected_expected) - verified_selected
    missing_all = set(all_expected) - verified_all
    flow.update(
        selected_case_coordinates_expected=len(buffers.selected_expected),
        selected_case_coordinates_verified=len(verified_selected),
        selected_case_coordinates_missing_or_mismatched=len(missing_selected),
        all_site_table_high_coordinates_expected=len(all_expected),
        all_site_table_high_coordinates_verified=len(verified_all),
        all_site_table_high_coordinates_missing_or_mismatched=len(missing_all),
    )
    if not partial and (missing_selected or missing_all or high_only_examples):
        raise ValueError(
            "The site table and bedMethyl are not bidirectionally complete: "
            f"selected missing={len(missing_selected)}, "
            f"all missing={len(missing_all)}, "
            "bedMethyl-only high calls="
            f"{flow['bedmethyl_high_coordinates_missing_from_site_table']}"
        )
    return flow, partial
