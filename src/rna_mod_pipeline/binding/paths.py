from __future__ import annotations

import json
from pathlib import Path


ENRICHMENT_FILENAMES = {
    "ornament": "rbp_enrichment_results.tsv",
    "encori": "encori_enrichment_results.tsv",
    "postar3": "postar3_enrichment_results.tsv",
}


def loose_defaults(
    root: Path,
    analysis: str,
    database: str,
) -> tuple[Path, Path, Path]:
    sources = {
        "ornament": root / "data" / "HS",
        "encori": root / "data" / "ENCORI",
        "postar3": root / "data" / "human.txt",
    }
    output = root / "refactored_outputs" / analysis / "loose" / database
    return sources[database], output, output / "plots"


def default_phase1_sites(root: Path, analysis: str, short_name: str) -> Path:
    return (
        root
        / "refactored_outputs"
        / analysis
        / "phase1"
        / f"filtered_{short_name}_metagene.tsv"
    )


def loose_output_paths(
    output_dir: Path,
    short_name: str,
    database: str,
) -> tuple[Path, Path, Path]:
    enrichment = ENRICHMENT_FILENAMES[database]
    annotation = {
        "ornament": f"filtered_{short_name}_rbp_annotated.tsv",
        "encori": f"filtered_{short_name}_encori_annotated.tsv",
        "postar3": f"filtered_{short_name}_postar3_annotated.tsv",
    }[database]
    return (
        output_dir / enrichment,
        output_dir / annotation,
        output_dir / "run_manifest.json",
    )


def validate_loose_output_ownership(
    output_dir: Path,
    modification: str,
    database: str,
) -> None:
    foreign = [
        output_dir / filename
        for key, filename in ENRICHMENT_FILENAMES.items()
        if key != database and (output_dir / filename).exists()
    ]
    if foreign:
        raise ValueError(
            "Loose output directory already contains another database run: "
            + ", ".join(str(path) for path in foreign)
        )
    manifest_path = output_dir / "run_manifest.json"
    if not manifest_path.is_file():
        return
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(
            f"Cannot validate existing loose-output manifest: {manifest_path}"
        ) from exc
    parameters = manifest.get("parameters", {})
    observed = (
        manifest.get("workflow"),
        parameters.get("modification"),
        parameters.get("database"),
    )
    expected = ("loose_overlap", modification, database)
    if observed != expected:
        raise ValueError(
            "Loose output directory belongs to a different run "
            f"{observed!r}, not {expected!r}: {output_dir}"
        )
