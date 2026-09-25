from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from scipy.special import betainc

from .generate import COUNTS, DATASET_NAME, SENTINELS, VERSION, rbp_names


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _number(observed: str, expected: float, label: str) -> None:
    value = float(observed) if observed else math.nan
    same = math.isnan(value) if math.isnan(expected) else math.isclose(
        value, expected, rel_tol=2e-8, abs_tol=1e-300)
    _assert(same, f"{label}: expected {expected!r}, observed {value!r}")


def fisher_reference(a: int, c: int, size: int = 384) -> float:
    bound = a + c
    denominator = math.comb(size * 2, size)
    probabilities = {
        x: math.comb(bound, x) * math.comb(size * 2 - bound, size - x) / denominator
        for x in range(max(0, bound - size), min(size, bound) + 1)
    }
    observed = probabilities[a]
    return min(1.0, math.fsum(value for value in probabilities.values()
                             if value <= observed * (1 + 1e-12)))


def _bh(values: dict[str, float], *, declared: bool) -> dict[str, float]:
    finite = {name: value for name, value in values.items() if math.isfinite(value)}
    size = len(values) if declared else len(finite)
    ordered = sorted(finite, key=lambda name: finite[name])
    result = {name: math.nan for name in values}
    for i, name in enumerate(ordered):
        result[name] = min(1.0, *(finite[other] * size / (j + 1)
                                for j, other in enumerate(ordered) if j >= i))
    return result


def analytical_reference() -> dict[str, dict[str, float]]:
    cluster_se = math.sqrt(8 / 279)
    statistic = math.log(9) / cluster_se
    cluster_p = float(betainc(31 / 2, 0.5, 31 / (31 + statistic**2)))
    specifications = {
        "DEMO_ENRICHED": (288, 96, 9.0, cluster_p),
        "DEMO_DEPLETED": (96, 288, 1 / 9, cluster_p),
        "DEMO_NULL": (192, 192, 1.0, 1.0),
        "DEMO_ZERO": (0, 192, 0.0, math.nan),
        "DEMO_INFINITY": (192, 0, math.inf, math.nan),
        "DEMO_ALL_BOUND": (384, 384, math.nan, math.nan),
    }
    return {name: {"cases": a, "controls": c, "or": odds,
                   "fisher_p": fisher_reference(a, c), "cluster_p": p,
                   "cluster_se": cluster_se if name in SENTINELS[:2] else math.nan}
            for name, (a, c, odds, p) in specifications.items()}


def _check_phase1(root: Path, output: Path) -> None:
    filtered = _rows(output / "phase1/filtered_m6A.tsv")
    drach = _rows(output / "phase1/filtered_m6A_drach.tsv")
    metagene = _rows(output / "phase1/filtered_m6A_metagene.tsv")
    with (root / "inputs/synthetic_m6a.bedmethyl").open() as handle:
        raw_count = sum(1 for _ in handle)
    _assert(raw_count == 800, f"Raw bedMethyl rows: {raw_count}, expected 800")
    _assert(len(filtered) == len(drach) == len(metagene) == 384,
            "All three Phase 1 tables must contain 384 cases; DRACH must not filter cases")
    _assert(sum(row["is_DRACH"].lower() == "true" for row in drach) == 288,
            "Expected 288 DRACH and 96 non-DRACH cases")
    _assert(Counter(row["region"] for row in metagene) == {"5UTR": 128, "CDS": 128, "3UTR": 128},
            "Expected 128 cases in each transcript region")
    _assert(Counter(row["strand"] for row in metagene) == {"+": 192, "-": 192},
            "Expected 192 cases on each strand")
    _assert(len({row["transcript_id"] for row in metagene}) == 32, "Expected 32 complete transcripts")
    _assert(all(row["motif_5mer"] in {"GGACT", "TTATT"} for row in drach),
            "Strand-corrected reference motifs differ from the synthetic design")


def _check_loose(output: Path, references: dict) -> list[dict]:
    filenames = {"ornament": "rbp_enrichment_results.tsv", "encori": "encori_enrichment_results.tsv",
                 "postar3": "postar3_enrichment_results.tsv"}
    summary = []
    for database, filename in filenames.items():
        table = _rows(output / "loose" / database / filename)
        _assert(len(table) == COUNTS[database], f"{database}: incomplete loose panel")
        _assert({row["RBP"] for row in table} == set(rbp_names(database)), f"{database}: wrong RBP identities")
        pvalues = {row["RBP"]: references.get(row["RBP"], references["DEMO_NULL"])["fisher_p"]
                   for row in table}
        adjusted = _bh(pvalues, declared=False)
        for row in table:
            name = row["RBP"]
            reference = references.get(name, references["DEMO_NULL"])
            for column, key in (("mod_overlaps", "cases"), ("bg_overlaps", "controls"),
                                ("odds_ratio", "or"), ("p_value", "fisher_p")):
                _number(row[column], reference[key], f"loose/{database}/{name}/{column}")
            _number(row["fdr"], adjusted[name], f"loose/{database}/{name}/FDR")
            _assert(row["partial_input"].lower() == "false", "Unexpected partial-input flag")
            if name in SENTINELS[:2]:
                _assert(row["inference_eligible"].lower() == "true", f"{name}: failed loose quality policy")
                summary.append({"design": "loose", "database": database, "RBP": name,
                                "odds_ratio": float(row["odds_ratio"]), "fdr": float(row["fdr"])})
    return summary


def _check_transcript(output: Path, references: dict) -> list[dict]:
    folder = output / "transcript_region"
    strata = _rows(folder / "transcript_region_strata.tsv")
    _assert(len(strata) == 96, "Expected 96 transcript-region strata")
    _assert(all(row["n_case_sites"] == row["n_control_sites"] == "4" for row in strata),
            "Every synthetic stratum must contain four cases and four controls")
    summary = []
    labels = {"oRNAment": "ornament", "ENCORI": "encori", "POSTAR3": "postar3"}
    for filename, multiplier in (("stratified_enrichment_results.tsv", 1),
                                 ("region_specific_results.tsv", 3)):
        table = _rows(folder / filename)
        _assert(len(table) == 630 * multiplier, f"Unexpected row count in {filename}")
        keys = {(row["database"], row["RBP"], row["region_filter"]) for row in table}
        _assert(len(keys) == len(table), f"Duplicate tests in {filename}")
        regions = ("all",) if multiplier == 1 else ("5UTR", "CDS", "3UTR")
        expected_keys = {(label, name, region) for label, database in labels.items()
                         for name in rbp_names(database) for region in regions}
        _assert(keys == expected_keys, f"Incomplete RBP-by-region panel in {filename}")
        for label, database in labels.items():
            _assert({row["RBP"] for row in table if row["database"] == label} == set(rbp_names(database)),
                    f"{database}: wrong transcript RBP panel")
        for row in table:
            name, database = row["RBP"], labels[row["database"]]
            reference = references.get(name, references["DEMO_NULL"])
            for column, expected in (("n_case_sites", 384 / multiplier),
                                     ("n_control_sites", 384 / multiplier),
                                     ("n_strata_total", 96 / multiplier), ("n_genes", 32),
                                     ("case_overlaps", reference["cases"] / multiplier),
                                     ("control_overlaps", reference["controls"] / multiplier),
                                     ("mh_or", reference["or"]),
                                     ("gene_cluster_pvalue", reference["cluster_p"])):
                _number(row[column], expected, f"transcript/{database}/{name}/{column}")
            fdr_column = "fdr_within_database" if multiplier == 1 else "fdr_within_database_region"
            expected_fdr = min(1.0, reference["cluster_p"] * COUNTS[database] / 2)
            if math.isnan(reference["cluster_p"]):
                expected_fdr = math.nan
            _number(row[fdr_column], expected_fdr, f"transcript/{database}/{name}/FDR")
            if name in SENTINELS[:2]:
                _number(row["gene_cluster_log_or_se"], reference["cluster_se"], "cluster SE")
                _assert(row["direction_stable_after_each_gene_removed"].lower() == "true", "Direction is not stable")
                _assert(row["inference_eligible"].lower() == "true", "Sentinel failed transcript quality policy")
                _assert(float(row[fdr_column]) < 0.05, "Expected significant finite sentinel")
                if multiplier == 1:
                    summary.append({"design": "transcript-region", "database": database,
                                    "RBP": name, "odds_ratio": float(row["mh_or"]),
                                    "fdr": float(row[fdr_column])})
            _assert(row["partial_input"].lower() == "false", "Unexpected partial-input flag")
    return summary


def _check_crossdb(output: Path) -> None:
    for design in ("loose", "transcript-region"):
        path = output / "cross_database" / f"{design}_all/triple_direction_concordant.tsv"
        rows = _rows(path)
        names = {row["RBP"] for row in rows}
        _assert(len(rows) == 2 and names == set(SENTINELS[:2]),
                f"{design}: expected only finite enriched/depleted sentinels as primary triple consensus; got {names}")


def _check_plots(output: Path) -> int:
    from PIL import Image

    required = [output / "phase1/plots" / name for name in (
        "01_coverage_distribution.png", "02_fraction_modified_distribution.png",
        "03_chromosome_distribution.png", "04_threshold_sensitivity.png",
        "04a_joint_threshold_heatmap.png", "05_filtered_chromosome_distribution.png",
        "06_filtered_fraction_distribution.png", "07_coverage_vs_fraction.png",
        "08_drach_pie.png", "09_top_motifs.png", "12_metagene_profile.png", "14_region_distribution.png")]
    for database in COUNTS:
        required.extend(output / "loose" / database / "plots" / f"{database}_{suffix}.png"
                        for suffix in ("all_rbps_overview", "top_enriched", "top_depleted", "machinery_heatmap"))
        required.extend(output / "transcript_region" / f"plots_{database}" / f"{database}_{suffix}.png"
                        for suffix in ("all_rbps_overview", "top_enriched", "top_depleted", "machinery_heatmap"))
        required.extend(output / "transcript_region/regional_top_plots" / f"{database}_{region}_top_{direction}.png"
                        for region in ("5utr", "cds", "3utr") for direction in ("enriched", "depleted"))
    for design in ("loose", "transcript-region"):
        required.extend(output / "cross_database" / f"{design}_all/plots" / f"{name}.png"
                        for name in ("validation_summary", "pairwise_scatter", "pairwise_shared_top",
                                     "pairwise_agreement", "database_membership_counts", "triple_validation"))
    for path in required:
        _assert(path.is_file(), f"Required plot not found: {path}")
    images = list(output.rglob("*.png"))
    for path in images:
        with Image.open(path) as image:
            _assert(image.width > 100 and image.height > 100, f"Invalid image dimensions: {path}")
            image.verify()
    return len(images)


def check_demo(workspace: str | Path, *, plots: bool = True) -> dict:
    root = Path(workspace).expanduser().resolve()
    manifest = json.loads((root / "SYNTHETIC_FIXTURE.json").read_text(encoding="utf-8"))
    _assert(manifest.get("fixture_version") == VERSION and manifest.get("synthetic_only") is True,
            "Unrecognised synthetic fixture")
    destination = root / "acceptance_check.json"
    destination.write_text(json.dumps({"passed": False, "status": "checks_not_completed",
                                      "synthetic_only": True}) + "\n", encoding="utf-8")
    for relative, digest in manifest["files"].items():
        path = root / relative
        _assert(path.resolve().is_relative_to(root), "Unsafe fixture path")
        _assert(hashlib.sha256(path.read_bytes()).hexdigest() == digest, f"Input changed: {relative}")
    output = root / "refactored_outputs/datasets" / DATASET_NAME
    references = analytical_reference()
    _check_phase1(root, output)
    summary = _check_loose(output, references) + _check_transcript(output, references)
    _check_crossdb(output)
    plot_count = _check_plots(output) if plots else 0
    result = {"passed": True, "synthetic_only": True, "fixture_version": VERSION,
              "plots_checked": plot_count, "numeric_results": summary,
              "scope": "Phase 1, all three loose panels, transcript-region panels, cross-database consensus"}
    destination.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result
