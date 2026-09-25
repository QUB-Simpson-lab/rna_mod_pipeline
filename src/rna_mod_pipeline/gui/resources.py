from __future__ import annotations

import json
import os
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping

from .datasets import _write_json
from .profile_paths import PATH_FORMAT, load_path, serialize_path, validate_path_format


RESOURCE_PROFILE_SCHEMA_VERSION = 1
DEFAULT_RESOURCE_PROFILE = Path("refactored_outputs") / "resource_profile.json"

RESOURCE_KINDS = {
    "reference_fasta": "file",
    "annotation_gtf": "file",
    "ornament_dir": "directory",
    "encori_dir": "directory",
    "postar3": "file",
    "nanopore_expression": "file",
    "depmap_expression": "file",
    "knockrbp_dir": "directory",
    "knockrbp_metadata": "file",
    "rbp_catalog": "file",
}

RESOURCE_LABELS = {
    "reference_fasta": "Reference FASTA",
    "annotation_gtf": "GENCODE GTF",
    "ornament_dir": "oRNAment directory",
    "encori_dir": "ENCORI directory",
    "postar3": "POSTAR3 human table",
    "nanopore_expression": "Nanopore expression table",
    "depmap_expression": "DepMap expression table",
    "knockrbp_dir": "KnockRBP dataset directory",
    "knockrbp_metadata": "KnockRBP metadata table",
    "rbp_catalog": "RBP catalogue",
}


def fasta_index_checks(fasta: str | Path) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Inspect pyfaidx's adjacent index without creating or changing it."""
    reference = Path(fasta).expanduser().resolve()
    index = Path(f"{reference}.fai")
    errors: list[str] = []
    warnings: list[str] = []
    if not reference.is_file():
        return (f"Reference FASTA is not an existing file: {reference}",), ()
    if index.exists():
        if (
            not index.is_file()
            or not os.access(index, os.R_OK)
            or index.stat().st_size == 0
        ):
            errors.append(f"FASTA index is not a readable, non-empty file: {index}")
        elif index.stat().st_mtime < reference.stat().st_mtime:
            warnings.append(
                f"FASTA index is older than the reference and pyfaidx may rebuild it: {index}"
            )
            if not os.access(reference.parent, os.W_OK):
                errors.append(
                    "The FASTA directory is read-only, so the stale index cannot be rebuilt; "
                    f"provide an up-to-date adjacent index: {index}"
                )
    elif os.access(reference.parent, os.W_OK):
        warnings.append(
            "No adjacent FASTA index was found. pyfaidx will create it when DRACH "
            f"analysis starts: {index}"
        )
    else:
        errors.append(
            "No adjacent FASTA index was found and the reference directory is "
            f"not writable. Provide {index} or use a writable reference copy."
        )
    return tuple(errors), tuple(warnings)


def _path(value: object) -> Path | None:
    if value in (None, ""):
        return None
    return Path(str(value)).expanduser()


def _resolve(path: Path | None, root: Path) -> Path | None:
    if path is None:
        return None
    return (path if path.is_absolute() else root / path).resolve()


def _portable(path: Path | None, root: Path) -> str | None:
    if path is None:
        return None
    absolute = _resolve(path, root)
    assert absolute is not None
    return serialize_path(absolute, root)


@dataclass(frozen=True)
class ResourceProfile:
    reference_fasta: Path | None = None
    annotation_gtf: Path | None = None
    ornament_dir: Path | None = None
    encori_dir: Path | None = None
    postar3: Path | None = None
    nanopore_expression: Path | None = None
    depmap_expression: Path | None = None
    knockrbp_dir: Path | None = None
    knockrbp_metadata: Path | None = None
    rbp_catalog: Path | None = None
    genome_build: str = ""
    annotation_release: str = ""

    def __post_init__(self) -> None:
        for key in RESOURCE_KINDS:
            object.__setattr__(self, key, _path(getattr(self, key)))
        object.__setattr__(self, "genome_build", str(self.genome_build).strip())
        object.__setattr__(
            self, "annotation_release", str(self.annotation_release).strip()
        )

    def resolved(self, project_root: str | Path) -> ResourceProfile:
        root = Path(project_root).expanduser().resolve()
        return replace(
            self,
            **{key: _resolve(getattr(self, key), root) for key in RESOURCE_KINDS},
        )

    def path(self, key: str) -> Path | None:
        if key not in RESOURCE_KINDS:
            raise KeyError(f"Unknown resource: {key}")
        return getattr(self, key)

    def require(self, key: str) -> Path:
        value = self.path(key)
        if value is None:
            raise ValueError(
                f"Configure {RESOURCE_LABELS[key]} in Resources before running this workflow"
            )
        return value

    def validation_errors(self, required: tuple[str, ...] = ()) -> tuple[str, ...]:
        unknown = set(required).difference(RESOURCE_KINDS)
        if unknown:
            raise KeyError(f"Unknown required resources: {sorted(unknown)}")
        errors: list[str] = []
        for key, kind in RESOURCE_KINDS.items():
            value = self.path(key)
            if value is None:
                if key in required:
                    errors.append(f"{RESOURCE_LABELS[key]} is not configured")
                continue
            exists = value.is_dir() if kind == "directory" else value.is_file()
            if not exists:
                expected = "directory" if kind == "directory" else "file"
                errors.append(
                    f"{RESOURCE_LABELS[key]} is not an existing {expected}: {value}"
                )
        if self.reference_fasta is not None and self.reference_fasta.is_file():
            index_errors, _warnings = fasta_index_checks(self.reference_fasta)
            errors.extend(index_errors)
        return tuple(errors)

    def validation_warnings(self) -> tuple[str, ...]:
        if self.reference_fasta is None or not self.reference_fasta.is_file():
            return ()
        _errors, warnings = fasta_index_checks(self.reference_fasta)
        return warnings

    def to_dict(self, project_root: str | Path) -> dict[str, object]:
        root = Path(project_root).expanduser().resolve()
        return {
            "schema_version": RESOURCE_PROFILE_SCHEMA_VERSION,
            "path_format": PATH_FORMAT,
            "resources": {
                key: _portable(self.path(key), root) for key in RESOURCE_KINDS
            },
            "reference": {
                "genome_build": self.genome_build,
                "annotation_release": self.annotation_release,
            },
        }

    @classmethod
    def from_dict(
        cls,
        payload: Mapping[str, object],
        project_root: str | Path,
    ) -> ResourceProfile:
        if payload.get("schema_version") != RESOURCE_PROFILE_SCHEMA_VERSION:
            raise ValueError("Unsupported resource-profile schema version")
        resources = payload.get("resources", {})
        reference = payload.get("reference", {})
        if not isinstance(resources, Mapping):
            raise ValueError("Resource profile must contain a resources object")
        if not isinstance(reference, Mapping):
            raise ValueError("Resource profile reference must be an object")
        path_format = payload.get("path_format")
        validate_path_format(path_format)
        profile = cls(
            **{
                key: load_path(resources.get(key), project_root, path_format)
                for key in RESOURCE_KINDS
            },
            genome_build=str(reference.get("genome_build", "")),
            annotation_release=str(reference.get("annotation_release", "")),
        )
        return profile.resolved(project_root)


def discover_project_resources(
    project_root: str | Path,
    code_root: str | Path | None = None,
) -> ResourceProfile:
    root = Path(project_root).expanduser().resolve()
    code = Path(
        code_root or Path(__file__).resolve().parents[3]
    ).expanduser().resolve()
    candidates = {
        "reference_fasta": root / "data/hg38.fa",
        "annotation_gtf": root / "data/gencode.v44.annotation.gtf",
        "ornament_dir": root / "data/HS",
        "encori_dir": root / "data/ENCORI",
        "postar3": root / "data/human.txt",
        "nanopore_expression": root
        / "data/Kate_231_0h_vs_6h_gene_counts_normalised.tsv",
        "depmap_expression": root
        / "data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
        "knockrbp_dir": root / "data/knockrbp",
        "knockrbp_metadata": root / "data/knockrbp/dataset_metadata.tsv",
        "rbp_catalog": code / "config/rbp_catalog.tsv",
    }
    values = {
        key: path
        for key, path in candidates.items()
        if (
            path.is_dir()
            if RESOURCE_KINDS[key] == "directory"
            else path.is_file()
        )
    }
    # A discovered profile is a convenience for this checkout, not a dependency.
    return ResourceProfile(
        **values,
        genome_build="hg38" if "reference_fasta" in values else "",
        annotation_release=(
            "GENCODE v44" if "annotation_gtf" in values else ""
        ),
    )


class ResourceProfileStore:
    def __init__(
        self,
        project_root: str | Path,
        path: str | Path | None = None,
        code_root: str | Path | None = None,
    ) -> None:
        self.project_root = Path(project_root).expanduser().resolve()
        self.code_root = Path(
            code_root or Path(__file__).resolve().parents[3]
        ).expanduser().resolve()
        candidate = Path(path).expanduser() if path else DEFAULT_RESOURCE_PROFILE
        self.path = _resolve(candidate, self.project_root)
        assert self.path is not None

    def load(self) -> ResourceProfile:
        if not self.path.is_file():
            return discover_project_resources(self.project_root, self.code_root)
        try:
            payload = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Cannot read resource profile {self.path}: {exc}") from exc
        if not isinstance(payload, Mapping):
            raise ValueError("Resource profile must be a JSON object")
        return ResourceProfile.from_dict(payload, self.project_root)

    def save(self, profile: ResourceProfile) -> Path:
        resolved = profile.resolved(self.project_root)
        errors = resolved.validation_errors()
        if errors:
            raise ValueError("\n".join(errors))
        return _write_json(self.path, resolved.to_dict(self.project_root))
