#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.config import MODIFICATIONS
from rna_mod_pipeline.cli import resolve_path
from rna_mod_pipeline.transcript.config import ANALYSES, DATABASE_LABELS
from rna_mod_pipeline.transcript.validation import validate_run
from rna_mod_pipeline.transcript.workflow import RunOptions, run_analysis


RESOURCE_DEFAULTS = {
    "ornament": "data/HS",
    "encori": "data/ENCORI",
    "postar3": "data/human.txt",
}
DEFAULT_RBP_CATALOG = (
    Path(__file__).resolve().parents[1] / "config" / "rbp_catalog.tsv"
)


def _items(values: list[str] | None) -> frozenset[str] | None:
    if not values:
        return None
    result = {
        item.strip()
        for value in values
        for item in value.split(",")
        if item.strip()
    }
    return frozenset(result) or None


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description=(
            "Test RBP overlap after exact transcript and transcript-region "
            "stratification. This is not the strict matched-site workflow."
        )
    )
    command.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[2]),
    )
    command.add_argument(
        "--modification",
        choices=tuple(MODIFICATIONS),
        help="Required unless --validate-run is used.",
    )
    command.add_argument("--sites")
    command.add_argument("--bedmethyl")
    command.add_argument("--gtf", default="data/gencode.v44.annotation.gtf")
    command.add_argument(
        "--rbp-catalog",
        default=DEFAULT_RBP_CATALOG,
    )
    command.add_argument("--ornament-dir")
    command.add_argument("--encori-dir")
    command.add_argument("--postar3")
    command.add_argument(
        "--databases",
        nargs="+",
        choices=tuple(DATABASE_LABELS),
        default=list(DATABASE_LABELS),
    )
    command.add_argument("--output-dir")
    command.add_argument("--coverage-min", type=int, default=20)
    command.add_argument("--case-min-fraction", type=float, default=20.0)
    command.add_argument("--background-max-fraction", type=float, default=20.0)
    command.add_argument("--window", type=int, default=10)
    command.add_argument(
        "--strand-mode", choices=("ignore", "same"), default="ignore"
    )
    command.add_argument("--min-encori-support", type=int, default=1)
    command.add_argument("--fdr", type=float, default=0.05)
    command.add_argument("--top-n", type=int, default=30)
    command.add_argument("--dpi", type=int, default=180)
    command.add_argument("--seed", type=int, default=20260728)
    command.add_argument("--rbps", nargs="+")
    command.add_argument("--chromosomes", nargs="+")
    command.add_argument("--postar-cell-types", nargs="+")
    command.add_argument("--postar-methods", nargs="+")
    command.add_argument("--maximum-rbps-per-database", type=int)
    command.add_argument("--max-cases", type=int)
    command.add_argument("--max-bed-rows", type=int)
    command.add_argument("--max-postar-rows", type=int)
    command.add_argument("--no-plots", action="store_true")
    command.add_argument("--overwrite", action="store_true")
    command.add_argument(
        "--validate-run",
        help="Validate an existing transcript-region output and exit.",
    )
    command.add_argument(
        "--skip-input-hashes",
        action="store_true",
        help="Skip input checksum validation with --validate-run.",
    )
    return command


def main() -> None:
    args = parser().parse_args()
    root = Path(args.project_root).expanduser().resolve()
    if args.skip_input_hashes and not args.validate_run:
        raise ValueError("--skip-input-hashes requires --validate-run")
    if args.validate_run:
        errors = validate_run(
            resolve_path(root, args.validate_run),
            root,
            check_input_hashes=not args.skip_input_hashes,
        )
        if errors:
            raise SystemExit("Validation failed:\n- " + "\n- ".join(errors))
        print("Transcript-region output passed validation.")
        return
    if not args.modification:
        raise SystemExit("--modification is required")
    explicit_resources = {
        "ornament": args.ornament_dir,
        "encori": args.encori_dir,
        "postar3": args.postar3,
    }
    unused_resources = [
        database
        for database, value in explicit_resources.items()
        if value and database not in args.databases
    ]
    if unused_resources:
        raise ValueError(
            "Explicit resource paths were supplied for databases not selected "
            f"by --databases: {', '.join(unused_resources)}"
        )
    resource_paths = {
        database: (
            resolve_path(root, explicit_resources[database])
            if explicit_resources[database]
            else resolve_path(root, RESOURCE_DEFAULTS[database])
        )
        for database in args.databases
    }
    analysis = ANALYSES[args.modification]
    sites = args.sites or analysis.sites
    bedmethyl = args.bedmethyl or analysis.bedmethyl
    output = args.output_dir or (
        root
        / "refactored_outputs"
        / args.modification
        / "transcript_region"
    )
    options = RunOptions(
        project_root=root,
        analysis_id=args.modification,
        modification=analysis.modification,
        sites=resolve_path(root, sites),
        bedmethyl=resolve_path(root, bedmethyl),
        gtf=resolve_path(root, args.gtf),
        rbp_catalog=resolve_path(root, args.rbp_catalog),
        output_dir=resolve_path(root, output),
        databases=tuple(args.databases),
        ornament_dir=resource_paths.get("ornament"),
        encori_dir=resource_paths.get("encori"),
        postar3=resource_paths.get("postar3"),
        coverage_min=args.coverage_min,
        case_min_fraction=args.case_min_fraction,
        background_max_fraction=args.background_max_fraction,
        window=args.window,
        strand_mode=args.strand_mode,
        minimum_encori_support=args.min_encori_support,
        fdr_threshold=args.fdr,
        top_n=args.top_n,
        dpi=args.dpi,
        seed=args.seed,
        rbps=_items(args.rbps),
        chromosomes=_items(args.chromosomes),
        postar_cell_types=_items(args.postar_cell_types),
        postar_methods=_items(args.postar_methods),
        maximum_rbps_per_database=args.maximum_rbps_per_database,
        max_cases=args.max_cases,
        max_bed_rows=args.max_bed_rows,
        max_postar_rows=args.max_postar_rows,
        make_plots=not args.no_plots,
        overwrite=args.overwrite,
    )
    destination = run_analysis(options)
    print(f"Transcript-region analysis written to {destination}")


if __name__ == "__main__":
    main()
