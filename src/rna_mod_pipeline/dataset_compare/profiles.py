from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ..config import MODIFICATIONS
from ..crossdb.io import DATABASES
from ..io import require_file


BIOLOGICAL_MODIFICATION = {
    "m6a": "m6a",
    "m6a_rep2": "m6a",
    "m5c": "m5c",
    "pseu": "pseu",
}
DESIGNS = ("loose", "transcript-region")
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")


@dataclass(frozen=True)
class ReferenceIdentity:
    genome_build: str
    annotation_release: str
    fasta_sha256: str | None = None
    gtf_sha256: str | None = None


@dataclass(frozen=True)
class DatasetProfile:
    profile_path: Path
    dataset_id: str
    analysis_id: str
    modification: str
    reference: ReferenceIdentity
    filtered_sites: Path
    metagene_sites: Path
    enrichment: dict[tuple[str, str], Path]

    @property
    def input_paths(self) -> tuple[Path, ...]:
        paths = [
            self.profile_path,
            self.filtered_sites,
            self.metagene_sites,
            *self.enrichment.values(),
        ]
        return tuple(dict.fromkeys(path.resolve() for path in paths))


def dataset_profile_document(
    *,
    dataset_id: str,
    analysis_id: str,
    genome_build: str,
    annotation_release: str,
    filtered_sites: str,
    metagene_sites: str,
    enrichment: dict[str, dict[str, str]] | None = None,
    fasta_sha256: str | None = None,
    gtf_sha256: str | None = None,
) -> dict:
    """Return the JSON-serializable DatasetProfile document used by the GUI/CLI."""
    normalized_analysis = analysis_id.strip().lower()
    if normalized_analysis not in BIOLOGICAL_MODIFICATION:
        raise ValueError(
            f"analysis_id must be one of {sorted(BIOLOGICAL_MODIFICATION)}"
        )
    for label, value in (
        ("dataset_id", dataset_id),
        ("genome_build", genome_build),
        ("annotation_release", annotation_release),
        ("filtered_sites", filtered_sites),
        ("metagene_sites", metagene_sites),
    ):
        if not str(value).strip():
            raise ValueError(f"{label} must be non-empty")
    normalized_enrichment: dict[str, dict[str, str]] = {}
    for design, databases in (enrichment or {}).items():
        normalized_design = str(design).strip().lower()
        if normalized_design not in DESIGNS or not isinstance(databases, dict):
            raise ValueError(f"Invalid enrichment design: {design!r}")
        normalized_databases: dict[str, str] = {}
        for database, path in databases.items():
            normalized_database = str(database).strip().lower()
            if normalized_database not in DATABASES or not str(path).strip():
                raise ValueError(
                    f"Invalid {normalized_design} enrichment entry: {database!r}"
                )
            if normalized_database in normalized_databases:
                raise ValueError(
                    f"Duplicate normalized database key: {normalized_database}"
                )
            normalized_databases[normalized_database] = str(path)
        normalized_enrichment[normalized_design] = normalized_databases
    reference = {
        "genome_build": str(genome_build),
        "annotation_release": str(annotation_release),
    }
    for key, value in (
        ("fasta_sha256", fasta_sha256),
        ("gtf_sha256", gtf_sha256),
    ):
        if value:
            if not _SHA256.fullmatch(str(value).strip()):
                raise ValueError(f"{key} must be a 64-character SHA-256 digest")
            reference[key] = str(value).strip().lower()
    return {
        "schema_version": 1,
        "dataset_id": str(dataset_id),
        "analysis_id": normalized_analysis,
        "modification": BIOLOGICAL_MODIFICATION[normalized_analysis],
        "reference": reference,
        "tables": {
            "filtered_sites": str(filtered_sites),
            "metagene_sites": str(metagene_sites),
            "enrichment": normalized_enrichment,
        },
    }


def enrichment_pair_inventory(
    first: DatasetProfile,
    second: DatasetProfile,
) -> dict[str, tuple[tuple[str, str], ...]]:
    first_pairs = set(first.enrichment)
    second_pairs = set(second.enrichment)
    return {
        "dataset_a": tuple(sorted(first_pairs)),
        "dataset_b": tuple(sorted(second_pairs)),
        "common": tuple(sorted(first_pairs & second_pairs)),
        "dataset_a_only": tuple(sorted(first_pairs - second_pairs)),
        "dataset_b_only": tuple(sorted(second_pairs - first_pairs)),
    }


def format_pairs(pairs: tuple[tuple[str, str], ...]) -> str:
    return ";".join(f"{design}/{database}" for design, database in pairs)


def _text(mapping: dict, key: str, label: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} requires a non-empty {key!r} string")
    return value.strip()


def _checksum(mapping: dict, key: str, label: str) -> str | None:
    value = mapping.get(key)
    if value in (None, ""):
        return None
    if not isinstance(value, str) or not _SHA256.fullmatch(value.strip()):
        raise ValueError(f"{label}.{key} must be a 64-character SHA-256 digest")
    return value.strip().lower()


def _table_path(base: Path, value: object, label: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty path string")
    path = Path(value).expanduser()
    resolved = path.resolve() if path.is_absolute() else (base / path).resolve()
    return require_file(resolved, label)


def _enrichment_paths(
    base: Path,
    tables: dict,
    label: str,
) -> dict[tuple[str, str], Path]:
    raw = tables.get("enrichment", {})
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise ValueError(f"{label}.tables.enrichment must be an object")
    unknown_designs = set(raw) - set(DESIGNS)
    if unknown_designs:
        raise ValueError(
            f"{label} has unknown enrichment designs: {sorted(unknown_designs)}"
        )
    paths: dict[tuple[str, str], Path] = {}
    for design, database_paths in raw.items():
        if not isinstance(database_paths, dict):
            raise ValueError(f"{label} enrichment {design!r} must be an object")
        normalized = {}
        for key, value in database_paths.items():
            normalized_key = str(key).strip().lower()
            if normalized_key in normalized:
                raise ValueError(
                    f"{label} {design} repeats normalized database "
                    f"{normalized_key!r}"
                )
            normalized[normalized_key] = value
        unknown = set(normalized) - set(DATABASES)
        if unknown:
            raise ValueError(
                f"{label} {design} has unknown databases: {sorted(unknown)}"
            )
        for database, value in normalized.items():
            paths[(design, database)] = _table_path(
                base,
                value,
                f"{label} {design}/{database} enrichment table",
            )
    return paths


def load_profile(path: str | Path) -> DatasetProfile:
    profile_path = require_file(path, "dataset profile").resolve()
    try:
        raw = json.loads(profile_path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid dataset profile JSON: {profile_path}") from exc
    if not isinstance(raw, dict):
        raise ValueError("Dataset profile must contain a JSON object")
    if raw.get("schema_version") != 1:
        raise ValueError("Dataset profile schema_version must equal 1")

    label = f"dataset profile {profile_path}"
    dataset_id = _text(raw, "dataset_id", label)
    analysis_id = _text(raw, "analysis_id", label).lower()
    if analysis_id not in MODIFICATIONS:
        raise ValueError(
            f"{label} analysis_id must be one of {sorted(MODIFICATIONS)}"
        )
    modification = _text(raw, "modification", label).lower()
    biological = BIOLOGICAL_MODIFICATION[analysis_id]
    if modification != biological:
        raise ValueError(
            f"{label} modification {modification!r} is incompatible with "
            f"analysis_id {analysis_id!r} ({biological!r})"
        )

    reference_raw = raw.get("reference")
    if not isinstance(reference_raw, dict):
        raise ValueError(f"{label} requires a reference object")
    reference = ReferenceIdentity(
        genome_build=_text(reference_raw, "genome_build", f"{label}.reference"),
        annotation_release=_text(
            reference_raw, "annotation_release", f"{label}.reference"
        ),
        fasta_sha256=_checksum(reference_raw, "fasta_sha256", f"{label}.reference"),
        gtf_sha256=_checksum(reference_raw, "gtf_sha256", f"{label}.reference"),
    )

    tables = raw.get("tables")
    if not isinstance(tables, dict):
        raise ValueError(f"{label} requires a tables object")
    base = profile_path.parent
    return DatasetProfile(
        profile_path=profile_path,
        dataset_id=dataset_id,
        analysis_id=analysis_id,
        modification=modification,
        reference=reference,
        filtered_sites=_table_path(
            base, tables.get("filtered_sites"), f"{label} filtered-sites table"
        ),
        metagene_sites=_table_path(
            base, tables.get("metagene_sites"), f"{label} metagene-sites table"
        ),
        enrichment=_enrichment_paths(base, tables, label),
    )


def _digest_status(first: str | None, second: str | None) -> str:
    if first is None and second is None:
        return "not_provided"
    if first is None:
        return "unverified_missing_in_dataset_a"
    if second is None:
        return "unverified_missing_in_dataset_b"
    if first != second:
        return "mismatch"
    return "verified_equal"


def validate_compatibility(
    first: DatasetProfile,
    second: DatasetProfile,
) -> pd.DataFrame:
    if first.dataset_id == second.dataset_id:
        raise ValueError("Dataset profiles must use distinct dataset_id values")
    if first.modification != second.modification:
        raise ValueError(
            "Dataset biological modifications differ: "
            f"{first.modification!r} versus {second.modification!r}"
        )
    if first.reference.genome_build.casefold() != second.reference.genome_build.casefold():
        raise ValueError(
            "Dataset reference genome builds differ: "
            f"{first.reference.genome_build!r} versus "
            f"{second.reference.genome_build!r}"
        )
    if (
        first.reference.annotation_release.casefold()
        != second.reference.annotation_release.casefold()
    ):
        raise ValueError(
            "Dataset annotation releases differ: "
            f"{first.reference.annotation_release!r} versus "
            f"{second.reference.annotation_release!r}"
        )
    inventory = enrichment_pair_inventory(first, second)
    fasta_status = _digest_status(
        first.reference.fasta_sha256, second.reference.fasta_sha256
    )
    gtf_status = _digest_status(first.reference.gtf_sha256, second.reference.gtf_sha256)
    if "mismatch" in (fasta_status, gtf_status):
        raise ValueError(
            "Dataset reference checksum mismatch: "
            f"FASTA={fasta_status}, GTF={gtf_status}"
        )
    return pd.DataFrame(
        [
            {
                "comparison_type": "same_modification_dataset_robustness",
                "dataset_a": first.dataset_id,
                "dataset_b": second.dataset_id,
                "analysis_id_a": first.analysis_id,
                "analysis_id_b": second.analysis_id,
                "modification": first.modification,
                "genome_build": first.reference.genome_build,
                "annotation_release": first.reference.annotation_release,
                "fasta_identity_status": fasta_status,
                "gtf_identity_status": gtf_status,
                "n_enrichment_pairs_a": len(inventory["dataset_a"]),
                "n_enrichment_pairs_b": len(inventory["dataset_b"]),
                "n_enrichment_pairs_compared": len(inventory["common"]),
                "enrichment_pairs_a": format_pairs(inventory["dataset_a"]),
                "enrichment_pairs_b": format_pairs(inventory["dataset_b"]),
                "enrichment_pairs_compared": format_pairs(inventory["common"]),
                "enrichment_pairs_a_only": format_pairs(
                    inventory["dataset_a_only"]
                ),
                "enrichment_pairs_b_only": format_pairs(
                    inventory["dataset_b_only"]
                ),
                "independent_clip_replication": False,
                "interpretation": (
                    "Robustness to dataset/callset variation; shared RBP resources "
                    "are not independent CLIP replication."
                ),
            }
        ]
    )
