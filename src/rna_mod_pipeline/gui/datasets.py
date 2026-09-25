from __future__ import annotations
import json
import os
import re
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping

from .profile_paths import PATH_FORMAT, load_path, serialize_path, validate_path_format

MODIFICATIONS = frozenset({"m6a", "m5c", "pseu"})
REGISTRY_SCHEMA_VERSION = 1
DEFAULT_REGISTRY = Path("refactored_outputs") / "dataset_registry.json"
def _path(value: str | Path | None) -> Path | None:
    return Path(value).expanduser() if value is not None else None
def _absolute(path: Path, root: Path) -> Path:
    return (path if path.is_absolute() else root / path).resolve()
def _portable(path: Path | None, root: Path | None) -> str | None:
    if path is None:
        return None
    if root is None:
        return serialize_path(path, None)
    absolute = _absolute(path, root)
    return serialize_path(absolute, root)
def _slug(name: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._-")
    return value or "dataset"
def _write_json(path: Path, payload: Mapping[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.",
        suffix=".tmp", delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path
@dataclass(frozen=True)
class DatasetProfile:
    name: str
    modification: str
    bedmethyl: Path
    fasta: Path | None = None
    gtf: Path | None = None
    sample_metadata: Mapping[str, str] | None = None
    output_root: Path | None = None
    sample: str | None = None
    cell_line: str | None = None
    replicate: str | None = None

    def __post_init__(self) -> None:
        name = str(self.name).strip()
        modification = str(self.modification).strip().lower()
        if not name or "\n" in name or "\r" in name:
            raise ValueError("Dataset name must be non-empty and single-line")
        if modification not in MODIFICATIONS:
            raise ValueError(
                "Dataset modification must be one of m6a, m5c, or pseu"
            )
        metadata = {
            str(key): str(value)
            for key, value in dict(self.sample_metadata or {}).items()
            if value is not None
        }
        for key, value in (
            ("sample", self.sample),
            ("cell_line", self.cell_line),
            ("replicate", self.replicate),
        ):
            if value is not None:
                metadata[key] = str(value)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "modification", modification)
        object.__setattr__(self, "bedmethyl", _path(self.bedmethyl))
        object.__setattr__(self, "fasta", _path(self.fasta))
        object.__setattr__(self, "gtf", _path(self.gtf))
        object.__setattr__(self, "output_root", _path(self.output_root))
        object.__setattr__(self, "sample_metadata", metadata)
        object.__setattr__(self, "sample", metadata.get("sample"))
        object.__setattr__(self, "cell_line", metadata.get("cell_line"))
        object.__setattr__(self, "replicate", metadata.get("replicate"))
    def resolved(self, project_root: str | Path) -> DatasetProfile:
        root = Path(project_root).expanduser().resolve()
        profile = replace(
            self,
            bedmethyl=_absolute(self.bedmethyl, root),
            fasta=_absolute(self.fasta or Path("data/hg38.fa"), root),
            gtf=_absolute(
                self.gtf or Path("data/gencode.v44.annotation.gtf"), root
            ),
            output_root=_absolute(
                self.output_root
                or Path("refactored_outputs") / "datasets" / _slug(self.name),
                root,
            ),
        )
        inputs = (profile.bedmethyl, profile.fasta, profile.gtf)
        if any(path == profile.output_root or path.is_relative_to(profile.output_root)
               for path in inputs):
            raise ValueError("Dataset output root cannot contain an input file")
        return profile
    def to_dict(self, project_root: str | Path | None = None) -> dict[str, object]:
        root = (
            Path(project_root).expanduser().resolve()
            if project_root is not None
            else None
        )
        return {
            "path_format": PATH_FORMAT,
            "name": self.name,
            "modification": self.modification,
            "bedmethyl": _portable(self.bedmethyl, root),
            "fasta": _portable(self.fasta, root),
            "gtf": _portable(self.gtf, root),
            "sample_metadata": dict(self.sample_metadata or {}),
            "output_root": _portable(self.output_root, root),
        }
    @classmethod
    def from_dict(
        cls,
        record: Mapping[str, object],
        project_root: str | Path | None = None,
    ) -> DatasetProfile:
        required = {"name", "modification", "bedmethyl"}
        missing = required.difference(record)
        if missing:
            raise ValueError(f"Dataset profile is missing: {sorted(missing)}")
        metadata = record.get("sample_metadata")
        if metadata is not None and not isinstance(metadata, Mapping):
            raise ValueError("sample_metadata must be a JSON object")
        path_format = record.get("path_format")
        validate_path_format(path_format)
        paths = {
            key: load_path(record.get(key), project_root, path_format)
            for key in ("bedmethyl", "fasta", "gtf", "output_root")
        }
        if paths["bedmethyl"] is None:
            raise ValueError("Dataset bedmethyl path must be non-empty")
        return cls(
            name=str(record["name"]),
            modification=str(record["modification"]),
            bedmethyl=paths["bedmethyl"],
            fasta=paths["fasta"],
            gtf=paths["gtf"],
            sample_metadata=metadata,
            output_root=paths["output_root"],
        )
class DatasetRegistry:
    def __init__(
        self,
        project_root: str | Path,
        profiles: tuple[DatasetProfile, ...] = (),
        registry_path: str | Path | None = None,
    ) -> None:
        self.project_root = Path(project_root).expanduser().resolve()
        self.path = (
            Path(registry_path).expanduser().resolve()
            if registry_path
            else self.project_root / DEFAULT_REGISTRY
        )
        self._profiles: dict[str, DatasetProfile] = {}
        for profile in profiles:
            self.add(profile)
    @classmethod
    def load(
        cls,
        project_root: str | Path,
        registry_path: str | Path | None = None,
    ) -> DatasetRegistry:
        registry = cls(project_root, registry_path=registry_path)
        if not registry.path.is_file():
            return registry
        payload = json.loads(registry.path.read_text(encoding="utf-8"))
        if payload.get("schema_version") != REGISTRY_SCHEMA_VERSION:
            raise ValueError("Unsupported dataset-registry schema version")
        records = payload.get("datasets")
        if not isinstance(records, list):
            raise ValueError("Dataset registry must contain a datasets list")
        for record in records:
            if not isinstance(record, Mapping):
                raise ValueError("Every dataset record must be a JSON object")
            registry.add(DatasetProfile.from_dict(record, registry.project_root))
        return registry
    def list(self) -> tuple[DatasetProfile, ...]:
        return tuple(self._profiles[key] for key in sorted(self._profiles))
    def get(self, name: str) -> DatasetProfile:
        try:
            return self._profiles[name.casefold()]
        except KeyError as exc:
            raise KeyError(f"Unknown dataset: {name}") from exc
    def add(self, profile: DatasetProfile) -> DatasetProfile:
        resolved = profile.resolved(self.project_root)
        key = resolved.name.casefold()
        if key in self._profiles:
            raise ValueError(f"Dataset name already exists: {resolved.name}")
        for existing in self._profiles.values():
            first, second = resolved.output_root, existing.output_root
            if first == second or first.is_relative_to(second) or second.is_relative_to(first):
                raise ValueError(
                    "Dataset output roots must be separate and non-nested"
                )
        self._profiles[key] = resolved
        return resolved
    def replace(self, original_name: str, profile: DatasetProfile) -> DatasetProfile:
        original_key = original_name.casefold()
        if original_key not in self._profiles:
            raise KeyError(f"Unknown dataset: {original_name}")
        resolved = profile.resolved(self.project_root)
        replacement_key = resolved.name.casefold()
        if replacement_key != original_key and replacement_key in self._profiles:
            raise ValueError(f"Dataset name already exists: {resolved.name}")
        for key, existing in self._profiles.items():
            if key == original_key:
                continue
            first, second = resolved.output_root, existing.output_root
            if first == second or first.is_relative_to(second) or second.is_relative_to(first):
                raise ValueError(
                    "Dataset output roots must be separate and non-nested"
                )
        updated = dict(self._profiles)
        del updated[original_key]
        updated[replacement_key] = resolved
        self._profiles = updated
        return resolved
    def remove(self, name: str) -> DatasetProfile:
        try:
            return self._profiles.pop(name.casefold())
        except KeyError as exc:
            raise KeyError(f"Unknown dataset: {name}") from exc
    def save(self, path: str | Path | None = None) -> Path:
        destination = (
            Path(path).expanduser().resolve() if path is not None else self.path
        )
        payload = {
            "schema_version": REGISTRY_SCHEMA_VERSION,
            "datasets": [
                profile.to_dict(self.project_root) for profile in self.list()
            ],
        }
        _write_json(destination, payload)
        self.path = destination
        return destination
