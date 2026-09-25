from __future__ import annotations

from pathlib import Path

import pandas as pd
from pyfaidx import Fasta

from ..io import load_sites, require_file, write_tsv
from ..schemas import validate_modification_codes


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def is_drach(sequence: str) -> bool:
    motif = sequence.upper()
    return (
        len(motif) == 5
        and motif[0] in "AGT"
        and motif[1] in "AG"
        and motif[2] == "A"
        and motif[3] == "C"
        and motif[4] in "ACT"
    )


def annotate_drach(
    sites_path: str | Path,
    fasta_path: str | Path,
    output_path: str | Path,
    maximum_non_a_percent: float = 5.0,
    expected_modification: str | None = None,
) -> pd.DataFrame:
    """Add strand-corrected five-mer and DRACH annotations without filtering."""
    sites = load_sites(sites_path)
    if expected_modification is not None:
        validate_modification_codes(
            sites,
            expected_modification,
            "m6A site table",
        )
    required = {"chrom", "start", "strand"}
    missing = required.difference(sites.columns)
    if missing:
        raise ValueError(f"Site table is missing columns: {sorted(missing)}")

    genome = Fasta(str(require_file(fasta_path, "reference FASTA")))
    chromosomes = set(genome.keys())
    motifs: list[str] = []
    flags: list[object] = []
    statuses: list[str] = []
    extracted = 0
    non_a = 0

    for row in sites.itertuples(index=False):
        chrom = str(row.chrom)
        start = int(row.start)
        if (
            chrom not in chromosomes
            or start < 2
            or start + 3 > len(genome[chrom])
        ):
            motifs.append("N/A")
            flags.append(pd.NA)
            statuses.append("unavailable_reference_context")
            continue

        sequence = str(genome[chrom][start - 2 : start + 3]).upper()
        if row.strand == "-":
            sequence = reverse_complement(sequence)
        extracted += 1
        non_a += sequence[2] != "A"
        motifs.append(sequence)
        flags.append(is_drach(sequence))
        statuses.append("assessable")

    if extracted:
        non_a_percent = 100 * non_a / extracted
        if non_a_percent > maximum_non_a_percent:
            raise ValueError(
                f"{non_a_percent:.1f}% of extracted sites have a non-A centre; "
                "check coordinates and reference build"
            )

    annotated = sites.copy()
    annotated["motif_5mer"] = motifs
    annotated["motif_status"] = statuses
    annotated["is_DRACH"] = pd.array(flags, dtype="boolean")
    write_tsv(annotated, output_path)
    return annotated
