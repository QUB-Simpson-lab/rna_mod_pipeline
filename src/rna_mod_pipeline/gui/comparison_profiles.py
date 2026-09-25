from __future__ import annotations

from pathlib import Path

from .datasets import DatasetProfile, _write_json


def _transcript_enrichment(output_root: Path) -> dict[str, str]:
    table = output_root / "transcript_region" / "stratified_enrichment_results.tsv"
    if not table.is_file():
        return {}
    import pandas as pd

    database_values = pd.read_csv(
        table,
        sep="\t",
        usecols=["database"],
    )["database"]
    supported = {"ornament", "encori", "postar3"}
    return {
        database: str(table)
        for raw in database_values.dropna().astype(str).unique()
        if (database := raw.strip().casefold()) in supported
    }


def comparison_profile_payload(profile: DatasetProfile) -> dict[str, object]:
    if profile.output_root is None or profile.fasta is None or profile.gtf is None:
        raise ValueError("Comparison profiles require resolved dataset paths")
    short = {"m6a": "m6A", "m5c": "m5C", "pseu": "pseU"}[profile.modification]
    phase = profile.output_root / "phase1"
    names = {
        "ornament": "rbp_enrichment_results.tsv",
        "encori": "encori_enrichment_results.tsv",
        "postar3": "postar3_enrichment_results.tsv",
    }
    loose = {
        database: str(path)
        for database, filename in names.items()
        if (path := profile.output_root / "loose" / database / filename).is_file()
    }
    enrichment: dict[str, dict[str, str]] = {}
    if loose:
        enrichment["loose"] = loose
    transcript = _transcript_enrichment(profile.output_root)
    if transcript:
        enrichment["transcript-region"] = transcript
    metadata = dict(profile.sample_metadata or {})
    reference = {
        "genome_build": metadata.get("genome_build", "hg38"),
        "annotation_release": metadata.get("annotation_release", "gencode.v44"),
    }
    for key in ("fasta_sha256", "gtf_sha256"):
        if metadata.get(key):
            reference[key] = metadata[key]
    return {
        "schema_version": 1,
        "dataset_id": profile.name,
        "analysis_id": profile.modification,
        "modification": profile.modification,
        "reference": reference,
        "tables": {
            "filtered_sites": str(phase / f"filtered_{short}.tsv"),
            "metagene_sites": str(phase / f"filtered_{short}_metagene.tsv"),
            "enrichment": enrichment,
        },
    }


def comparison_profile_errors(profile: DatasetProfile) -> tuple[str, ...]:
    payload = comparison_profile_payload(profile)
    tables = payload["tables"]
    assert isinstance(tables, dict)
    errors = []
    for key in ("filtered_sites", "metagene_sites"):
        path = Path(str(tables[key])).expanduser().resolve()
        if not path.is_file():
            errors.append(f"Required {key} table not found: {path}")
    return tuple(errors)


def write_comparison_profile(profile: DatasetProfile) -> Path:
    payload = comparison_profile_payload(profile)
    assert profile.output_root is not None
    return _write_json(
        profile.output_root / "dataset_comparison_profile.json",
        payload,
    )


def prepare_comparison_profiles(
    first: DatasetProfile,
    second: DatasetProfile,
    project_root: str | Path,
) -> tuple[Path, Path]:
    left = first.resolved(project_root)
    right = second.resolved(project_root)
    errors = [*comparison_profile_errors(left), *comparison_profile_errors(right)]
    if errors:
        raise ValueError("\n".join(errors))
    return write_comparison_profile(left), write_comparison_profile(right)
