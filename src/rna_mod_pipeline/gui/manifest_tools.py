from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class OutputRecord:
    path: Path
    bytes: int | None
    sha256: str | None


@dataclass(frozen=True)
class ManifestView:
    path: Path
    workflow: str
    created_utc: str
    parameters: dict[str, Any]
    outputs: tuple[OutputRecord, ...]
    raw: dict[str, Any]


def locate_manifest(selection: str | Path) -> Path:
    path = Path(selection).expanduser().resolve()
    return path / "run_manifest.json" if path.is_dir() else path


def load_manifest(selection: str | Path, project_root: str | Path) -> ManifestView:
    path = locate_manifest(selection)
    raw = json.loads(path.read_text())
    root = Path(project_root).expanduser().resolve()
    outputs = []
    for record in raw.get("outputs", []):
        value = Path(str(record.get("path", ""))).expanduser()
        resolved = value.resolve() if value.is_absolute() else (root / value).resolve()
        outputs.append(
            OutputRecord(
                resolved,
                int(record["bytes"]) if "bytes" in record else None,
                str(record["sha256"]) if record.get("sha256") else None,
            )
        )
    return ManifestView(
        path=path,
        workflow=str(raw.get("workflow", "unknown")),
        created_utc=str(raw.get("created_utc", "")),
        parameters=dict(raw.get("parameters", {})),
        outputs=tuple(outputs),
        raw=raw,
    )


def fast_validate(view: ManifestView) -> list[str]:
    errors = []
    if not view.workflow or view.workflow == "unknown":
        errors.append("manifest has no workflow identity")
    if not view.outputs:
        errors.append("manifest contains no output records")
    for record in view.outputs:
        if not record.path.is_file():
            errors.append(f"missing output: {record.path}")
        elif record.bytes is not None and record.path.stat().st_size != record.bytes:
            errors.append(f"size mismatch: {record.path}")
    return errors


def validate_completed_manifest(
    manifest_path: str | Path,
    project_root: str | Path,
) -> tuple[ManifestView | None, list[str]]:
    path = Path(manifest_path).expanduser().resolve()
    if not path.is_file():
        return None, [f"expected run manifest is missing: {path}"]
    try:
        view = load_manifest(path, project_root)
    except Exception as exc:
        return None, [f"cannot read completed run manifest: {exc}"]
    return view, fast_validate(view)


def text_preview(path: str | Path, limit: int = 200_000) -> str:
    source = Path(path)
    with source.open("rb") as handle:
        content = handle.read(limit + 1)
    truncated = len(content) > limit
    text = content[:limit].decode("utf-8", errors="replace")
    return text + ("\n\n[Preview truncated]" if truncated else "")
