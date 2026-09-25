#!/usr/bin/env python3
"""Create or verify an exact-byte inventory of external pipeline resources."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


REPOSITORY = Path(__file__).resolve().parents[1]
DEFAULT_SNAPSHOT = REPOSITORY / "config/resource_snapshot_20260908.json"
GROUPS = ("reference", "ornament", "encori", "postar3", "knockrbp", "expression")
SOURCES = {
    "reference": {
        "build": "hg38/GRCh38; GENCODE v44 comprehensive CHR annotation",
        "urls": [
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
        ],
    },
    "ornament": {"urls": ["https://rnabiology.ircm.qc.ca/oRNAment/downloads/"]},
    "encori": {"urls": ["https://rnasysu.com/encori/download.php"]},
    "postar3": {"urls": ["https://postar.ncrnalab.org/"]},
    "knockrbp": {"urls": []},
    "expression": {"urls": []},
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    before = path.stat()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(chunk)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"File changed while hashing: {path}")
    return value.hexdigest()


def table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def relative_path(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or "\\" in value or ":" in value:
        raise ValueError(f"Expected a portable relative path: {value}")
    if not path.parts or path.parts == (".",):
        raise ValueError("Empty resource path")
    return path


def inventory_inputs(workspace: Path, catalog: Path) -> list[dict]:
    records = []

    def add(group: str, path: Path, **metadata) -> None:
        if not path.is_file():
            raise FileNotFoundError(path)
        records.append({"group": group, "path": path.relative_to(workspace).as_posix(),
                        **metadata})

    for name in ("hg38.fa", "hg38.fa.fai", "gencode.v44.annotation.gtf"):
        add("reference", workspace / "data" / name)
    rows = table(catalog)
    for group, directory in (("ornament", "HS"), ("encori", "ENCORI")):
        sources = defaultdict(set)
        for row in rows:
            if row["database"] == group:
                for source in row["source_names"].split(";"):
                    sources[source.strip()].add(row["canonical_rbp"])
        for source, names in sorted(sources.items()):
            root = workspace / "data" / directory
            candidates = [root / f"{source}.bed.gz", root / f"{source}.bed"]
            selected = next((path for path in candidates if path.is_file()), None)
            if selected is None:
                raise FileNotFoundError(f"No BED resource for {group}/{source}")
            add(group, selected, source_name=source, canonical_rbps=sorted(names))
    add("postar3", workspace / "data/human.txt")
    root = workspace / "data/knockrbp"
    add("knockrbp", root / "dataset_metadata.tsv")
    datasets = table(root / "dataset_metadata.tsv")
    if len({row["dataset_id"] for row in datasets}) != len(datasets):
        raise ValueError("Duplicate KnockRBP dataset IDs")
    expected_jsons = set()
    for row in datasets:
        name = f"{row['rbp']}_{row['dataset_id']}_degs.json"
        expected_jsons.add(name)
        add("knockrbp", root / name, dataset_id=row["dataset_id"], rbp=row["rbp"],
            cell_line=row["cell_line"],
            selected_by_default=row["cell_line"] in ("MDA-MB-231", "MDA-MB-231-LM2"))
    extras = {path.name for path in root.glob("*_degs.json")} - expected_jsons
    if extras:
        raise ValueError(f"KnockRBP files absent from metadata: {sorted(extras)}")
    for name in ("Kate_231_0h_vs_6h_gene_counts_normalised.tsv",
                 "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"):
        add("expression", workspace / "data" / name)
    return sorted(records, key=lambda row: (row["group"], row["path"]))


def create_snapshot(args) -> int:
    if args.output is None:
        raise ValueError("--create requires --output with a new snapshot filename")
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite snapshot: {args.output}")
    if args.groups or args.maps:
        raise ValueError("--groups and --map are verification options only")
    rows = inventory_inputs(args.workspace, args.catalog)
    last_group = None
    for row in rows:
        if row["group"] != last_group:
            last_group = row["group"]
            print(f"Hashing {last_group} ...", flush=True)
        path = args.workspace / row["path"]
        row["bytes"] = path.stat().st_size
        row["sha256"] = digest(path)
    snapshot = {
        "schema_version": 1,
        "measured_at_utc": datetime.now(timezone.utc).isoformat(),
        "download_date": None,
        "scope": "Exact local external-resource files; excludes raw reads, bedMethyl and patient data",
        "provenance_note": "URLs copied from project RESOURCE_SETUP.md; download dates and unrecorded database releases unavailable",
        "catalog": {"path": "config/rbp_catalog.tsv", "sha256": digest(args.catalog),
                    "versions": sorted({row['catalog_version'] for row in table(args.catalog)})},
        "sources": SOURCES,
        "group_counts": dict(sorted(Counter(row["group"] for row in rows).items())),
        "total_bytes": sum(row["bytes"] for row in rows),
        "files": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(snapshot, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"Created {args.output}: {len(rows)} files; {snapshot['total_bytes']:,} bytes")
    return 0


def mappings(values: list[str]) -> list[tuple[PurePosixPath, Path]]:
    result = []
    for value in values:
        prefix, separator, target = value.partition("=")
        if not separator or not target:
            raise ValueError("--map expects RESOURCE_PATH=LOCAL_PATH")
        result.append((relative_path(prefix), Path(target).expanduser().resolve()))
    if len({prefix for prefix, _ in result}) != len(result):
        raise ValueError("Duplicate --map resource prefix")
    return sorted(result, key=lambda item: len(item[0].parts), reverse=True)


def resolve(path: str, workspace: Path, maps) -> Path:
    relative = relative_path(path)
    for prefix, destination in maps:
        if relative == prefix or prefix in relative.parents:
            return destination.joinpath(*relative.relative_to(prefix).parts)
    return workspace.joinpath(*relative.parts)


def verify_snapshot(args) -> int:
    snapshot = json.loads(args.snapshot.read_text(encoding="utf-8"))
    if snapshot.get("schema_version") != 1:
        raise ValueError("Unsupported resource snapshot schema")
    rows = snapshot["files"]
    if len({row["path"] for row in rows}) != len(rows):
        raise ValueError("Duplicate resource paths in snapshot")
    if any(row["group"] not in GROUPS for row in rows):
        raise ValueError("Unknown resource group in snapshot")
    if Counter(row["group"] for row in rows) != Counter(snapshot["group_counts"]):
        raise ValueError("Snapshot group counts disagree with file records")
    if sum(row["bytes"] for row in rows) != snapshot["total_bytes"]:
        raise ValueError("Snapshot total size disagrees with file records")
    groups = set(args.groups or GROUPS)
    maps = mappings(args.maps)
    selected_paths = [relative_path(row["path"]) for row in rows if row["group"] in groups]
    for prefix, _ in maps:
        if not any(path == prefix or prefix in path.parents for path in selected_paths):
            raise ValueError(f"--map does not match a selected resource: {prefix}")
    failures = 0
    if digest(args.catalog) != snapshot["catalog"]["sha256"]:
        print("FAIL: RBP catalogue differs from snapshot")
        failures += 1
    checked = 0
    for group in sorted(groups):
        selected = [row for row in rows if row["group"] == group]
        if not selected:
            raise ValueError(f"No snapshot files in requested group: {group}")
        print(f"Verifying {group}: {len(selected)} files ...", flush=True)
        for row in selected:
            path = resolve(row["path"], args.workspace, maps)
            try:
                if path.stat().st_size != row["bytes"]:
                    raise ValueError("byte size differs")
                if digest(path) != row["sha256"]:
                    raise ValueError("SHA-256 differs")
                checked += 1
            except (OSError, ValueError) as exc:
                print(f"FAIL: {row['path']}: {exc}")
                failures += 1
    print(f"Verified {checked} files; {failures} failure(s). "
          "Only selected inventory files were checked.")
    return 1 if failures else 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True,
                        help="Workspace containing data/; resources are read only")
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--catalog", type=Path, default=REPOSITORY / "config/rbp_catalog.tsv")
    parser.add_argument("--groups", nargs="+", choices=GROUPS,
                        help="Verify selected groups only; omitted means every group")
    parser.add_argument("--map", dest="maps", action="append", default=[],
                        help="Map inventory path or directory to a local path, e.g. data/HS=D:/RBP/HS")
    parser.add_argument("--create", action="store_true", help="Create a NEW complete inventory")
    parser.add_argument("--output", type=Path, help="New JSON path for --create")
    args = parser.parse_args(argv)
    args.workspace = args.workspace.expanduser().resolve()
    try:
        if args.output and not args.create:
            raise ValueError("--output requires --create")
        return create_snapshot(args) if args.create else verify_snapshot(args)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"Resource inventory error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
