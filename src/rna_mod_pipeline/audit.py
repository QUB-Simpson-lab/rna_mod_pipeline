from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .binding.catalog import EXPECTED_CATALOG_COUNTS, load_catalog
from .binding.sources import binding_source_files
from .config import MODIFICATIONS, canonical_rbp


@dataclass(frozen=True)
class AuditCheck:
    check: str
    passed: bool
    details: str


def _check(label: str, operation) -> AuditCheck:
    try:
        details = operation()
        return AuditCheck(label, True, str(details or "OK"))
    except Exception as exc:
        return AuditCheck(label, False, f"{type(exc).__name__}: {exc}")


def _require_files(root: Path, paths: list[str]) -> str:
    missing = [path for path in paths if not (root / path).is_file()]
    if missing:
        raise FileNotFoundError(f"missing {missing}")
    return f"{len(paths)} required files present"


def _catalog_check(code_root: Path, database: str) -> str:
    entries = load_catalog(code_root / "config" / "rbp_catalog.tsv", database)
    expected = EXPECTED_CATALOG_COUNTS[database]
    if len(entries) != expected:
        raise ValueError(f"{len(entries)} entries; expected {expected}")
    return f"{len(entries)} unique canonical RBPs"


def _resource_check(root: Path, code_root: Path, database: str) -> str:
    entries = load_catalog(code_root / "config" / "rbp_catalog.tsv", database)
    directory = root / "data" / ("HS" if database == "ornament" else "ENCORI")
    paths = binding_source_files(database, directory, entries)
    if any(path.stat().st_size == 0 for path in paths):
        raise ValueError("one or more binding files are empty")
    return f"{len(paths)} source files resolve for {len(entries)} RBPs"


def _phase1_check(root: Path, key: str) -> str:
    config = MODIFICATIONS[key]
    filtered = root / config.filtered_sites
    metagene = root / config.metagene_sites
    _require_files(root, [config.filtered_sites, config.metagene_sites])
    sites = pd.read_csv(filtered, sep="\t")
    required = {
        "chrom",
        "start",
        "strand",
        "mod_code",
        "Nvalid_cov",
        "fraction_modified",
    }
    if not required.issubset(sites):
        raise ValueError(f"filtered table missing {sorted(required - set(sites))}")
    if sites.duplicated(["chrom", "start", "strand"]).any():
        raise ValueError("duplicate genomic site keys")
    if (pd.to_numeric(sites["Nvalid_cov"]) < 20).any():
        raise ValueError("coverage below 20 in default filtered output")
    if (pd.to_numeric(sites["fraction_modified"]) < 20).any():
        raise ValueError("fraction below 20 in default filtered output")
    observed_codes = set(sites["mod_code"].astype(str).str.lower())
    expected_codes = MODIFICATIONS["m6a" if key == "m6a_rep2" else key].expected_mod_codes
    if not observed_codes <= expected_codes:
        raise ValueError(f"unexpected modification codes {sorted(observed_codes)}")
    mapped = pd.read_csv(metagene, sep="\t")
    if len(mapped) != len(sites) or not {
        "transcript_id",
        "gene_name",
        "region",
        "metagene_pos",
    }.issubset(mapped):
        raise ValueError("metagene table is incomplete or changes the site count")
    if config.drach_sites:
        drach = pd.read_csv(root / config.drach_sites, sep="\t")
        if len(drach) != len(sites) or not {"motif_5mer", "is_DRACH"}.issubset(drach):
            raise ValueError("DRACH table is incomplete or changes the site count")
    return f"{len(sites):,} filtered sites; Phase 1 tables align"


def _legacy_loose_paths(root: Path, key: str) -> dict[str, Path]:
    folder = MODIFICATIONS[key].folder
    return {
        "ornament": root / folder / "phase_2_results" / "rbp_enrichment_results.tsv",
        "encori": root / folder / "phase_2b_encori" / "encori_enrichment_results.tsv",
        "postar3": root / folder / "phase_2c_postar3" / "postar3_enrichment_results.tsv",
    }


def _legacy_loose_check(root: Path, code_root: Path, key: str) -> str:
    summaries = []
    for database, path in _legacy_loose_paths(root, key).items():
        table = pd.read_csv(path, sep="\t")
        required = {"RBP", "odds_ratio", "p_value", "fdr"}
        if not required.issubset(table):
            raise ValueError(f"{database} missing {sorted(required - set(table))}")
        observed = {canonical_rbp(value).upper() for value in table["RBP"]}
        expected = {
            entry.canonical_rbp.upper()
            for entry in load_catalog(
                code_root / "config" / "rbp_catalog.tsv",
                database,
            )
        }
        if observed != expected or len(table) != len(expected):
            raise ValueError(
                f"{database} result/catalog mismatch "
                f"(rows={len(table)}, expected={len(expected)})"
            )
        summaries.append(f"{database}={len(table)}")
    return ", ".join(summaries)


def _legacy_transcript_check(root: Path, key: str) -> str:
    legacy_key = "pseudouridine" if key == "pseu" else key
    run = root / "transcript_region_stratified_overlap" / "results" / legacy_key
    table = pd.read_csv(run / "stratified_enrichment_results.tsv", sep="\t")
    required = {
        "database",
        "RBP",
        "mh_or",
        "gene_cluster_pvalue",
        "fdr_within_database",
        "case_overlaps",
        "control_overlaps",
    }
    if not required.issubset(table):
        raise ValueError(f"missing {sorted(required - set(table))}")
    if table.duplicated(["database", "RBP"]).any():
        raise ValueError("duplicate database/RBP results")
    expected = {"oRNAment": 133, "ENCORI": 281, "POSTAR3": 216}
    observed = table.groupby("database")["RBP"].nunique().to_dict()
    if observed != expected:
        raise ValueError(f"RBP panel counts {observed}; expected {expected}")
    for database in ("ornament", "encori", "postar3"):
        for suffix in (".npz", "_rows.tsv"):
            path = run / f"case_binding_matrix_{database}{suffix}"
            if not path.is_file():
                raise FileNotFoundError(path)
    return f"{len(table)} RBP rows with complete database panels"


def _knockrbp_check(root: Path) -> str:
    directory = root / "data" / "knockrbp"
    metadata = pd.read_csv(directory / "dataset_metadata.tsv", sep="\t")
    required = {"dataset_id", "rbp", "cell_line", "has_pvalues"}
    if not required.issubset(metadata):
        raise ValueError(f"metadata missing {sorted(required - set(metadata))}")
    if metadata["dataset_id"].duplicated().any():
        raise ValueError("duplicate dataset_id values")
    indexed = metadata.set_index("dataset_id")
    seen = []
    pattern = re.compile(r"(.+)_DataSet_(\d+)_degs\.json$")
    for path in sorted(directory.glob("*_degs.json")):
        match = pattern.fullmatch(path.name)
        if not match:
            raise ValueError(f"unexpected filename {path.name}")
        dataset_id = f"DataSet_{match.group(2)}"
        if dataset_id not in indexed.index:
            raise ValueError(f"{dataset_id} lacks metadata")
        file_rbp = canonical_rbp(match.group(1)).upper()
        metadata_rbp = canonical_rbp(indexed.loc[dataset_id, "rbp"]).upper()
        if file_rbp != metadata_rbp:
            raise ValueError(f"{dataset_id} RBP mismatch")
        seen.append(dataset_id)
    if set(indexed.index) != set(seen):
        raise ValueError("metadata and DEG file dataset IDs differ")
    per_rbp = metadata.groupby(metadata["rbp"].map(lambda x: canonical_rbp(x).upper())).size()
    repeated = per_rbp[per_rbp > 1].to_dict()
    return f"{len(seen)} datasets; repeated-RBP datasets retained: {repeated}"


def _top_level_function_lengths(path: Path) -> list[tuple[str, int]]:
    tree = ast.parse(path.read_text(), filename=str(path))
    functions = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.append((node.name, node.end_lineno - node.lineno + 1))
        elif isinstance(node, ast.ClassDef):
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    functions.append(
                        (
                            f"{node.name}.{child.name}",
                            child.end_lineno - child.lineno + 1,
                        )
                    )
    return functions


def _code_structure_check(code_root: Path) -> str:
    expected = {
        "phase1_filter.py",
        "phase1_drach.py",
        "phase1_metagene.py",
        "loose_overlap.py",
        "transcript_region_overlap.py",
        "plot_enrichment.py",
        "cross_database.py",
        "integrate_expression.py",
        "knockrbp_validation.py",
        "string_analysis.py",
        "compare_modifications.py",
        "compare_datasets.py",
        "pipeline_gui.py",
        "validate_handoff.py",
    }
    scripts = {path.name for path in (code_root / "scripts").glob("*.py")}
    if expected - scripts:
        raise ValueError(f"missing command wrappers {sorted(expected - scripts)}")
    too_long = {
        path.name: len(path.read_text().splitlines())
        for path in (code_root / "scripts").glob("*.py")
        if len(path.read_text().splitlines()) > 300
    }
    if too_long:
        raise ValueError(f"command wrappers exceed 300 lines: {too_long}")

    package_root = code_root / "src" / "rna_mod_pipeline"
    package_files = sorted(package_root.rglob("*.py"))
    oversized_modules = {
        path.relative_to(code_root).as_posix(): len(path.read_text().splitlines())
        for path in package_files
        if len(path.read_text().splitlines()) > 500
    }
    if oversized_modules:
        raise ValueError(f"package modules exceed 500 lines: {oversized_modules}")

    oversized_functions = {}
    production_files = [*package_files, *sorted((code_root / "scripts").glob("*.py"))]
    transcript_root = package_root / "transcript"
    function_count = 0
    for path in production_files:
        threshold = 180 if transcript_root in path.parents else 250
        for name, length in _top_level_function_lengths(path):
            function_count += 1
            if length > threshold:
                key = f"{path.relative_to(code_root).as_posix()}:{name}"
                oversized_functions[key] = f"{length} lines (limit {threshold})"
    if oversized_functions:
        raise ValueError(
            "production functions exceed structural limits: "
            f"{oversized_functions}"
        )
    return (
        f"{len(expected)} runnable wrappers; {len(package_files)} package modules "
        f"and {function_count} production functions within structural limits"
    )


def run_handoff_audit(
    project_root: str | Path,
    *,
    include_legacy_results: bool = True,
) -> list[AuditCheck]:
    root = Path(project_root).expanduser().resolve()
    code_root = Path(__file__).resolve().parents[2]
    checks = [
        _check(
            "Core project inputs",
            lambda: _require_files(
                root,
                [
                    "data/hg38.fa",
                    "data/gencode.v44.annotation.gtf",
                    "data/human.txt",
                    "data/Kate_231_0h_vs_6h_gene_counts_normalised.tsv",
                    "data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                    "data/knockrbp/dataset_metadata.tsv",
                    *(definition.bedmethyl for definition in MODIFICATIONS.values()),
                ],
            ),
        ),
        _check("Command structure", lambda: _code_structure_check(code_root)),
        *[
            _check(
                f"{database} catalogue",
                lambda database=database: _catalog_check(code_root, database),
            )
            for database in EXPECTED_CATALOG_COUNTS
        ],
        *[
            _check(
                f"{database} binding resources",
                lambda database=database: _resource_check(root, code_root, database),
            )
            for database in ("ornament", "encori")
        ],
        _check(
            "POSTAR3 resource",
            lambda: (
                f"{(root / 'data' / 'human.txt').stat().st_size:,} bytes"
                if (root / "data" / "human.txt").stat().st_size > 0
                else (_ for _ in ()).throw(ValueError("POSTAR3 file is empty"))
            ),
        ),
        _check("KnockRBP dataset mapping", lambda: _knockrbp_check(root)),
    ]
    if include_legacy_results:
        for key in MODIFICATIONS:
            checks.extend(
                [
                    _check(
                        f"{key} Phase 1 baseline",
                        lambda key=key: _phase1_check(root, key),
                    ),
                    _check(
                        f"{key} loose-result baseline",
                        lambda key=key: _legacy_loose_check(root, code_root, key),
                    ),
                    _check(
                        f"{key} transcript-result baseline",
                        lambda key=key: _legacy_transcript_check(root, key),
                    ),
                ]
            )
    return checks
