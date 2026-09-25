from __future__ import annotations

import heapq
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from rna_mod_pipeline.io import open_text

from .config import REGIONS


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return {
        match.group(1): match.group(2)
        for match in re.finditer(r'(\w+)\s+"([^"]*)"', text)
    }


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    intervals: dict[str, list[tuple[int, int]]]

    def __post_init__(self) -> None:
        for region in REGIONS:
            self.intervals[region] = sorted(self.intervals.get(region, []))
        self.total_length = sum(
            end - start
            for region in REGIONS
            for start, end in self.intervals[region]
        )


def load_complete_basic_transcripts(
    gtf_path: Path,
    chromosomes: Optional[set[str]] = None,
) -> dict[str, TranscriptModel]:
    transcript_info: dict[str, dict[str, str]] = {}
    features: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    protein_coding: set[str] = set()
    basic: set[str] = set()

    with open_text(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {
                "transcript",
                "UTR",
                "CDS",
                "five_prime_UTR",
                "three_prime_UTR",
            }:
                continue
            chrom = fields[0]
            if chromosomes and chrom not in chromosomes:
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue
            if fields[2] == "transcript":
                transcript_info[transcript_id] = {
                    "chrom": chrom,
                    "strand": fields[6],
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", ""),
                }
                if attrs.get("transcript_type") == "protein_coding":
                    protein_coding.add(transcript_id)
                if 'tag "basic"' in line:
                    basic.add(transcript_id)
                continue
            try:
                start, end = int(fields[3]) - 1, int(fields[4])
            except ValueError:
                continue
            if end > start:
                features[transcript_id][fields[2]].append((start, end))

    models: dict[str, TranscriptModel] = {}
    for transcript_id in sorted(protein_coding & basic):
        info = transcript_info.get(transcript_id)
        if not info or info["strand"] not in {"+", "-"}:
            continue
        transcript_features = features.get(transcript_id, {})
        cds = sorted(transcript_features.get("CDS", []))
        if (
            transcript_features.get("five_prime_UTR")
            and transcript_features.get("three_prime_UTR")
        ):
            utr5 = sorted(transcript_features["five_prime_UTR"])
            utr3 = sorted(transcript_features["three_prime_UTR"])
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
        models[transcript_id] = TranscriptModel(
            transcript_id=transcript_id,
            gene_id=info["gene_id"],
            gene_name=info["gene_name"],
            chrom=info["chrom"],
            strand=info["strand"],
            intervals={"5UTR": utr5, "CDS": cds, "3UTR": utr3},
        )
    return models


class PrimaryTranscriptIndex:
    """Assign each genomic base to one longest complete Basic transcript."""

    def __init__(self, models: dict[str, TranscriptModel]):
        grouped: dict[tuple[str, str], list[tuple]] = defaultdict(list)
        transcript_rank = {
            transcript_id: rank
            for rank, transcript_id in enumerate(sorted(models))
        }
        region_rank = {
            region: rank for rank, region in enumerate(sorted(REGIONS))
        }
        for transcript_id, model in models.items():
            for region in REGIONS:
                for start, end in model.intervals[region]:
                    grouped[(model.chrom, model.strand)].append(
                        (
                            int(start),
                            int(end),
                            transcript_id,
                            region,
                            int(model.total_length),
                            transcript_rank[transcript_id],
                            region_rank[region],
                        )
                    )

        self.data: dict[tuple[str, str], tuple[np.ndarray, ...]] = {}
        for key in sorted(grouped):
            records: dict[int, tuple] = {}
            events: list[tuple[int, int, int]] = []
            for interval_id, record in enumerate(grouped[key]):
                start, end, transcript_id, region, total, tx_rank, reg_rank = record
                records[interval_id] = (
                    -total,
                    tx_rank,
                    reg_rank,
                    start,
                    end,
                    interval_id,
                    transcript_id,
                    region,
                )
                events.extend(((start, 1, interval_id), (end, -1, interval_id)))
            events.sort(key=lambda value: (value[0], value[1]))
            active: set[int] = set()
            heap: list[tuple] = []
            starts: list[int] = []
            ends: list[int] = []
            transcripts: list[str] = []
            regions: list[str] = []
            previous: int | None = None
            cursor = 0
            while cursor < len(events):
                position = events[cursor][0]
                while heap and heap[0][5] not in active:
                    heapq.heappop(heap)
                if previous is not None and previous < position and heap:
                    winner = heap[0]
                    transcript_id, region = winner[6], winner[7]
                    if (
                        starts
                        and ends[-1] == previous
                        and transcripts[-1] == transcript_id
                        and regions[-1] == region
                    ):
                        ends[-1] = position
                    else:
                        starts.append(previous)
                        ends.append(position)
                        transcripts.append(transcript_id)
                        regions.append(region)
                same_position = []
                while cursor < len(events) and events[cursor][0] == position:
                    same_position.append(events[cursor])
                    cursor += 1
                for _, action, interval_id in same_position:
                    if action == -1:
                        active.discard(interval_id)
                for _, action, interval_id in same_position:
                    if action == 1:
                        active.add(interval_id)
                        heapq.heappush(heap, records[interval_id])
                previous = position
            self.data[key] = (
                np.asarray(starts, dtype=np.int64),
                np.asarray(ends, dtype=np.int64),
                np.asarray(transcripts, dtype=object),
                np.asarray(regions, dtype=object),
            )

    def query(
        self, chrom: str, strand: str, position: int
    ) -> tuple[str, str] | None:
        indexed = self.data.get((chrom, strand))
        if indexed is None:
            return None
        starts, ends, transcripts, regions = indexed
        index = int(np.searchsorted(starts, position, side="right") - 1)
        if index >= 0 and position < ends[index]:
            return str(transcripts[index]), str(regions[index])
        return None
