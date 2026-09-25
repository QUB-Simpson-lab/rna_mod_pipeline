from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from ..io import load_sites, require_file, write_tsv
from ..schemas import validate_modification_codes


@dataclass(frozen=True)
class TranscriptRegions:
    transcript_id: str
    gene_name: str
    gene_id: str
    chrom: str
    strand: str
    utr5: tuple[tuple[int, int], ...]
    cds: tuple[tuple[int, int], ...]
    utr3: tuple[tuple[int, int], ...]

    @property
    def total_length(self) -> int:
        return sum(end - start for group in (self.utr5, self.cds, self.utr3) for start, end in group)


@dataclass(frozen=True)
class IndexedRegions:
    starts: np.ndarray
    ends: np.ndarray
    maximum_end: np.ndarray
    transcript_ids: tuple[str, ...]
    region_types: tuple[str, ...]


def parse_gtf_attributes(value: str) -> dict[str, str]:
    return dict(re.findall(r'(\w+)\s+"([^"]*)"', value))


def build_transcript_regions(
    gtf_path: str | Path,
) -> dict[str, TranscriptRegions]:
    """Load complete GENCODE Basic protein-coding transcript regions."""
    features: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    information: dict[str, dict[str, str]] = {}
    basic: set[str] = set()
    protein_coding: set[str] = set()

    with open(require_file(gtf_path, "GENCODE GTF"), "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue

            if feature == "transcript":
                if attrs.get("transcript_type") == "protein_coding":
                    protein_coding.add(transcript_id)
                if 'tag "basic"' in line:
                    basic.add(transcript_id)
                information[transcript_id] = {
                    "chrom": fields[0],
                    "strand": fields[6],
                    "gene_name": attrs.get("gene_name", ""),
                    "gene_id": attrs.get("gene_id", ""),
                }

            if feature in {"UTR", "CDS", "five_prime_UTR", "three_prime_UTR"}:
                features[transcript_id][feature].append(
                    (int(fields[3]) - 1, int(fields[4]))
                )

    regions: dict[str, TranscriptRegions] = {}
    for transcript_id in sorted(protein_coding & basic):
        info = information.get(transcript_id)
        if info is None:
            continue
        transcript_features = features.get(transcript_id, {})
        cds = transcript_features.get("CDS", [])

        if (
            transcript_features.get("five_prime_UTR")
            and transcript_features.get("three_prime_UTR")
        ):
            utr5 = transcript_features["five_prime_UTR"]
            utr3 = transcript_features["three_prime_UTR"]
        elif transcript_features.get("UTR") and cds:
            cds_min = min(start for start, _ in cds)
            cds_max = max(end for _, end in cds)
            utr5, utr3 = [], []
            for start, end in transcript_features["UTR"]:
                if info["strand"] == "+":
                    if end <= cds_min:
                        utr5.append((start, end))
                    elif start >= cds_max:
                        utr3.append((start, end))
                else:
                    if start >= cds_max:
                        utr5.append((start, end))
                    elif end <= cds_min:
                        utr3.append((start, end))
        else:
            utr5, utr3 = [], []

        if not utr5 or not cds or not utr3:
            continue
        regions[transcript_id] = TranscriptRegions(
            transcript_id=transcript_id,
            gene_name=info["gene_name"],
            gene_id=info["gene_id"],
            chrom=info["chrom"],
            strand=info["strand"],
            utr5=tuple(sorted(utr5)),
            cds=tuple(sorted(cds)),
            utr3=tuple(sorted(utr3)),
        )
    return regions


def build_region_index(
    regions: dict[str, TranscriptRegions],
) -> dict[tuple[str, str], IndexedRegions]:
    raw: dict[tuple[str, str], list[tuple[int, int, str, str]]] = defaultdict(list)
    for transcript_id, transcript in regions.items():
        key = (transcript.chrom, transcript.strand)
        for region_type, intervals in (
            ("5UTR", transcript.utr5),
            ("CDS", transcript.cds),
            ("3UTR", transcript.utr3),
        ):
            raw[key].extend(
                (start, end, transcript_id, region_type)
                for start, end in intervals
            )

    index: dict[tuple[str, str], IndexedRegions] = {}
    for key, values in raw.items():
        values.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
        starts = np.asarray([item[0] for item in values], dtype=np.int64)
        ends = np.asarray([item[1] for item in values], dtype=np.int64)
        index[key] = IndexedRegions(
            starts=starts,
            ends=ends,
            maximum_end=np.maximum.accumulate(ends),
            transcript_ids=tuple(item[2] for item in values),
            region_types=tuple(item[3] for item in values),
        )
    return index


def _scaled_position(
    position: int,
    transcript: TranscriptRegions,
    region_type: str,
) -> float:
    intervals = {
        "5UTR": transcript.utr5,
        "CDS": transcript.cds,
        "3UTR": transcript.utr3,
    }[region_type]
    ordered = sorted(intervals, reverse=transcript.strand == "-")
    region_length = sum(end - start for start, end in ordered)
    cumulative = 0
    for start, end in ordered:
        if start <= position < end:
            cumulative += (
                end - 1 - position
                if transcript.strand == "-"
                else position - start
            )
            break
        cumulative += end - start
    offset = {"5UTR": 0.0, "CDS": 1.0, "3UTR": 2.0}[region_type]
    return offset + cumulative / region_length


def find_region(
    chrom: str,
    position: int,
    strand: str,
    regions: dict[str, TranscriptRegions],
    index: dict[tuple[str, str], IndexedRegions],
) -> tuple[str, str, str, float] | None:
    indexed = index.get((chrom, strand))
    if indexed is None:
        return None

    right = int(np.searchsorted(indexed.starts, position, side="right")) - 1
    candidates: set[tuple[str, str]] = set()
    cursor = right
    while cursor >= 0 and indexed.maximum_end[cursor] > position:
        if indexed.starts[cursor] <= position < indexed.ends[cursor]:
            candidates.add(
                (indexed.transcript_ids[cursor], indexed.region_types[cursor])
            )
        cursor -= 1
    if not candidates:
        return None

    transcript_id, region_type = min(
        candidates,
        key=lambda item: (-regions[item[0]].total_length, item[0], item[1]),
    )
    transcript = regions[transcript_id]
    return (
        transcript_id,
        transcript.gene_name,
        region_type,
        _scaled_position(position, transcript, region_type),
    )


def annotate_metagene(
    sites_path: str | Path,
    gtf_path: str | Path,
    output_path: str | Path,
    expected_modification: str | None = None,
) -> pd.DataFrame:
    sites = load_sites(sites_path)
    if expected_modification is not None:
        validate_modification_codes(
            sites,
            expected_modification,
            "Phase 1 site table",
        )
    required = {"chrom", "start", "strand"}
    missing = required.difference(sites.columns)
    if missing:
        raise ValueError(f"Site table is missing columns: {sorted(missing)}")

    regions = build_transcript_regions(gtf_path)
    index = build_region_index(regions)
    assignments = [
        find_region(
            str(row.chrom),
            int(row.start),
            str(row.strand),
            regions,
            index,
        )
        for row in sites.itertuples(index=False)
    ]

    annotated = sites.copy()
    annotated["transcript_id"] = [
        assignment[0] if assignment else None for assignment in assignments
    ]
    annotated["gene_name"] = [
        assignment[1] if assignment else None for assignment in assignments
    ]
    annotated["region"] = [
        assignment[2] if assignment else "intergenic/intronic"
        for assignment in assignments
    ]
    annotated["metagene_pos"] = [
        assignment[3] if assignment else np.nan for assignment in assignments
    ]
    write_tsv(annotated, output_path)
    return annotated
