from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from ..gui.datasets import DatasetProfile, DatasetRegistry
from ..gui.resources import ResourceProfile, ResourceProfileStore


VERSION = "synthetic-acceptance-v1"
CHROMOSOME = "chrSYNTHETIC"
DATASET_NAME = "SYNTHETIC_m6a"
COUNTS = {"ornament": 133, "encori": 281, "postar3": 216}
SENTINELS = (
    "DEMO_ENRICHED", "DEMO_DEPLETED", "DEMO_NULL", "DEMO_ZERO",
    "DEMO_INFINITY", "DEMO_ALL_BOUND",
)


def _text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8", newline="\n")


def _table(path: Path, header: list[str], rows: list[list]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n",
                            quoting=csv.QUOTE_NONE, quotechar=None)
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def rbp_names(database: str) -> list[str]:
    return [*SENTINELS, *(f"DEMO_N{i:03}" for i in range(COUNTS[database] - 6))]


def bound_counts(name: str, gene: int) -> tuple[int, int]:
    high = 2 if gene < 8 else 3 if gene < 24 else 4
    if name == "DEMO_ENRICHED":
        return high, 1
    if name == "DEMO_DEPLETED":
        return 1, high
    if name == "DEMO_ZERO":
        return 0, 2
    if name == "DEMO_INFINITY":
        return 2, 0
    if name == "DEMO_ALL_BOUND":
        return 4, 4
    return (3, 1) if gene < 16 else (1, 3)


def _reference(root: Path) -> list[dict]:
    sequence = list("C" * 51500)
    records, gtf = [], []
    complement = str.maketrans("ACGT", "TGCA")
    for gene in range(32):
        beginning = 100 + gene * 1600
        strand = "+" if gene < 16 else "-"
        attrs = (
            f'gene_id "DEMO_G{gene:03}"; transcript_id "DEMO_T{gene:03}"; '
            f'gene_name "DEMO_GENE_{gene:03}"; gene_type "protein_coding"; '
            'transcript_type "protein_coding"; tag "basic";'
        )
        gtf.append([CHROMOSOME, "SYNTHETIC", "transcript", beginning + 1,
                    beginning + 1200, ".", strand, ".", attrs])
        regions = ["five_prime_UTR", "CDS", "three_prime_UTR"]
        if strand == "-":
            regions.reverse()
        for block, feature in enumerate(regions):
            left = beginning + 400 * block
            gtf.append([CHROMOSOME, "SYNTHETIC", feature, left + 1,
                        left + 400, ".", strand, "0" if feature == "CDS" else ".", attrs])
            for item in range(8):
                position = left + 25 + item * 40
                motif = "GGACT" if item % 4 < 3 else "TTATT"
                genomic = motif if strand == "+" else motif.translate(complement)[::-1]
                sequence[position - 2:position + 3] = genomic
                records.append({"position": position, "strand": strand,
                                "gene": gene, "item": item, "coverage": 40})
        position = beginning + 1190
        genomic = "GGACT" if strand == "+" else "AGTCC"
        sequence[position - 2:position + 3] = genomic
        records.append({"position": position, "strand": strand,
                        "gene": gene, "item": 0, "coverage": 10})
    fasta = root / "resources/synthetic.fa"
    header = f">{CHROMOSOME}\n"
    _text(fasta, header + "\n".join("".join(sequence[i:i + 60])
                                   for i in range(0, len(sequence), 60)) + "\n")
    _text(Path(f"{fasta}.fai"), f"{CHROMOSOME}\t{len(sequence)}\t{len(header)}\t60\t61\n")
    _table(root / "resources/synthetic.gtf", [], gtf)
    return records


def _bedmethyl(root: Path, records: list[dict]) -> None:
    rows = []
    for row in records:
        position, coverage, item = row["position"], row["coverage"], row["item"]
        fraction = 50 if coverage == 10 else (25 if item < 4 else 0)
        modified = coverage * fraction // 100
        rows.append([CHROMOSOME, position, position + 1, "a", coverage,
                     row["strand"], position, position + 1, "255,0,0", coverage,
                     fraction, modified, coverage - modified, 0, 0, 0, 0, 0])
    _table(root / "inputs/synthetic_m6a.bedmethyl", [], rows)
    _table(root / "inputs/synthetic_m6a_identical_rep.bedmethyl", [], rows)


def _binding(root: Path, records: list[dict]) -> None:
    catalogue, postar = [], []
    eligible = [row for row in records if row["coverage"] == 40]
    for database in COUNTS:
        for name in rbp_names(database):
            catalogue.append([VERSION, database, name, name, name])
            rows = []
            for row in eligible:
                case_count, control_count = bound_counts(name, row["gene"])
                item = row["item"]
                if not ((item < 4 and item < case_count)
                        or (item >= 4 and item - 4 < control_count)):
                    continue
                position = row["position"]
                rows.append([CHROMOSOME, position, position + 1,
                             name, 1, row["strand"]])
                if database == "postar3":
                    postar.append([CHROMOSOME, position, position + 1, "DEMO_PEAK",
                                   row["strand"], name, "SYNTHETIC_CLIP",
                                   "SYNTHETIC_NOT_A_CELL_LINE", "DEMO_NOT_AN_ACCESSION", 1])
            if database != "postar3":
                _table(root / f"resources/{database}/{name}.bed", [], rows)
    _table(root / "resources/synthetic_postar3.txt", [], postar)
    _table(root / "resources/synthetic_rbp_catalog.tsv",
           ["catalog_version", "database", "rbp_name", "canonical_rbp", "source_names"], catalogue)


def prepare_demo(workspace: str | Path) -> Path:
    root = Path(workspace).expanduser().resolve()
    code_root = Path(__file__).resolve().parents[3]
    if root == code_root or root.is_relative_to(code_root):
        raise ValueError("Choose a demo workspace outside the source-code folder")
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise FileExistsError(f"Demo workspace must be new or empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    records = _reference(root)
    _bedmethyl(root, records)
    _binding(root, records)
    profile = ResourceProfile(
        reference_fasta=root / "resources/synthetic.fa",
        annotation_gtf=root / "resources/synthetic.gtf",
        ornament_dir=root / "resources/ornament",
        encori_dir=root / "resources/encori",
        postar3=root / "resources/synthetic_postar3.txt",
        rbp_catalog=root / "resources/synthetic_rbp_catalog.tsv",
        genome_build="SYNTHETIC_NOT_HG38", annotation_release=VERSION,
    )
    ResourceProfileStore(root).save(profile)
    registry = DatasetRegistry(root)
    for name, filename in ((DATASET_NAME, "synthetic_m6a.bedmethyl"),
                           ("SYNTHETIC_identical_rep", "synthetic_m6a_identical_rep.bedmethyl")):
        registry.add(DatasetProfile(
            name=name, modification="m6a", bedmethyl=root / "inputs" / filename,
            fasta=profile.reference_fasta, gtf=profile.annotation_gtf,
            cell_line="SYNTHETIC_NOT_A_CELL_LINE", sample="SOFTWARE_TEST_ONLY",
        ))
    registry.save()
    files = sorted(path for folder in ("inputs", "resources")
                   for path in (root / folder).rglob("*") if path.is_file())
    manifest = {"fixture_version": VERSION, "synthetic_only": True,
                "expected_counts": {"raw_rows": 800, "cases": 384, "controls": 384,
                                    "genes": 32, "strata": 96, "drach_cases": 288},
                "files": {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                          for path in files}}
    _text(root / "SYNTHETIC_FIXTURE.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    _text(root / "READ_ME_SYNTHETIC_ONLY.txt",
          "Artificial SOFTWARE ACCEPTANCE TEST. No human genome, real cell line, "
          "RBP binding observation or biological replicate is represented.\n"
          "All DEMO_* identifiers and interactions are invented test fixtures.\n"
          "Do not use these outputs in the research findings or manuscript.\n")
    return root
