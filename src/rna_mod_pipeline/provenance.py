from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy
import pandas
import scipy

from . import __version__


def sha256_file(path: str | Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def portable_path(path: str | Path, project_root: str | Path) -> str:
    resolved = Path(path).expanduser().resolve()
    root = Path(project_root).expanduser().resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return resolved.as_posix()


def file_record(path: str | Path, project_root: str | Path) -> dict[str, Any]:
    resolved = Path(path).expanduser().resolve()
    return {
        "path": portable_path(resolved, project_root),
        "bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def _package_version(distribution: str) -> str:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "not_installed"


def _git_record(root: Path) -> tuple[str, str]:
    try:
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout
        return commit, str(bool(dirty)).lower()
    except (OSError, subprocess.SubprocessError):
        return "not_available", "not_available"


def software_record(project_root: str | Path | None = None) -> dict[str, str]:
    record = {
        "rna_mod_pipeline": __version__,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
        "matplotlib": _package_version("matplotlib"),
        "networkx": _package_version("networkx"),
        "pyfaidx": _package_version("pyfaidx"),
    }
    if project_root is not None:
        root = Path(project_root).expanduser().resolve()
        commit, dirty = _git_record(root)
        record["git_commit"] = commit
        record["git_worktree_dirty"] = dirty
        code_root = Path(__file__).resolve().parents[2]
        code_commit, code_dirty = _git_record(code_root)
        record["code_git_commit"] = code_commit
        record["code_git_worktree_dirty"] = code_dirty
    return record


def constraint_records(project_root: str | Path) -> list[dict[str, Any]]:
    """Record available handoff constraint files without requiring a checkout."""
    code_root = Path(__file__).resolve().parents[2]
    names = (
        "requirements-lock.txt",
        "requirements-gui-lock.txt",
        "requirements-build-lock.txt",
    )
    return [
        file_record(path, project_root)
        for name in names
        if (path := code_root / name).is_file()
    ]


def make_manifest(
    workflow: str,
    project_root: str | Path,
    inputs: Iterable[str | Path],
    parameters: dict[str, Any],
    outputs: Iterable[str | Path] = (),
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "workflow": workflow,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": ".",
        "parameters": parameters,
        "inputs": [file_record(path, project_root) for path in inputs],
        "outputs": [file_record(path, project_root) for path in outputs],
        "software": software_record(project_root),
        "dependency_constraints": constraint_records(project_root),
    }


def remap_manifest_output_paths(
    manifest: dict[str, Any],
    output_paths: Iterable[str | Path],
    staging_root: str | Path,
    destination_root: str | Path,
    project_root: str | Path,
) -> None:
    """Replace staged output paths with their final published locations."""
    paths = [Path(path).expanduser().resolve() for path in output_paths]
    records = manifest.get("outputs", [])
    if len(records) != len(paths):
        raise ValueError("Manifest output records do not match staged output paths")
    stage = Path(staging_root).expanduser().resolve()
    destination = Path(destination_root).expanduser().resolve()
    for record, path in zip(records, paths):
        try:
            relative = path.relative_to(stage)
        except ValueError as exc:
            raise ValueError(f"Staged output is outside staging root: {path}") from exc
        record["path"] = portable_path(destination / relative, project_root)


def remap_manifest_output_files(
    manifest: dict[str, Any],
    staged_paths: Iterable[str | Path],
    final_paths: Iterable[str | Path],
    project_root: str | Path,
) -> None:
    """Replace staged output paths with explicitly paired final file paths."""
    staged = [Path(path).expanduser().resolve() for path in staged_paths]
    final = [Path(path).expanduser().resolve() for path in final_paths]
    records = manifest.get("outputs", [])
    if len(records) != len(staged) or len(staged) != len(final):
        raise ValueError("Manifest output records do not match staged/final files")
    for record, staged_path, final_path in zip(records, staged, final):
        recorded = resolve_record_path(record, project_root).resolve()
        if recorded != staged_path:
            raise ValueError(
                f"Manifest output order does not match staged file: {staged_path}"
            )
        record["path"] = portable_path(final_path, project_root)


def write_manifest(manifest: dict[str, Any], path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary.replace(destination)


def resolve_record_path(record: dict[str, Any], project_root: str | Path) -> Path:
    path = Path(record["path"])
    if path.is_absolute():
        return path
    return Path(project_root).expanduser().resolve() / path


def validate_file_record(
    record: dict[str, Any],
    project_root: str | Path,
) -> list[str]:
    path = resolve_record_path(record, project_root)
    errors: list[str] = []
    if not path.is_file():
        return [f"missing file: {path}"]
    if path.stat().st_size != int(record["bytes"]):
        errors.append(f"size mismatch: {path}")
    if sha256_file(path) != record["sha256"]:
        errors.append(f"checksum mismatch: {path}")
    return errors


def validate_manifest(path: str | Path, project_root: str | Path) -> list[str]:
    manifest = json.loads(Path(path).read_text())
    errors: list[str] = []
    for record in [
        *manifest.get("inputs", []),
        *manifest.get("outputs", []),
        *manifest.get("dependency_constraints", []),
    ]:
        errors.extend(validate_file_record(record, project_root))
    return errors
