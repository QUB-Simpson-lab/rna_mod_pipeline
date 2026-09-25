from __future__ import annotations

import os
from pathlib import Path, PurePath, PurePosixPath, PureWindowsPath


PATH_FORMAT = "posix-relative"


def validate_path_format(value: object) -> None:
    if value not in (None, PATH_FORMAT):
        raise ValueError(f"Unsupported saved path format: {value}")


def serialize_path(path: PurePath | None, root: PurePath | None) -> str | None:
    if path is None:
        return None
    if root is not None:
        try:
            return path.relative_to(root).as_posix()
        except ValueError:
            return str(path)
    return str(path) if path.is_absolute() else path.as_posix()


def load_path(
    value: object,
    project_root: str | Path | None = None,
    path_format: object = None,
) -> Path | None:
    validate_path_format(path_format)
    if value in (None, ""):
        return None
    text = str(value)
    native = Path(text).expanduser()
    root = Path(project_root).expanduser().resolve() if project_root else None
    if os.name == "nt":
        if PurePosixPath(text).is_absolute() and not native.is_absolute():
            raise ValueError(f"Reselect this absolute path on this computer: {text}")
        if path_format == PATH_FORMAT and "\\" in text and not native.is_absolute():
            raise ValueError(f"Windows cannot represent this literal backslash path: {text}")
        return native
    if native.is_absolute():
        return native
    if root is not None and (root / native).exists():
        return native
    windows = PureWindowsPath(text)
    if windows.drive or text.startswith("\\"):
        raise ValueError(
            f"Reselect this Windows drive or network path on this computer: {text}. "
            "Use a new workspace, or move the saved profile JSON aside before "
            "opening the GUI to configure local paths."
        )
    if "\\" not in text or path_format == PATH_FORMAT:
        return native
    migrated = Path(windows.as_posix())
    if root is not None and (root / migrated).exists():
        return migrated
    raise ValueError(
        f"Cannot safely interpret the legacy backslash path: {text}. "
        "It may be a Windows relative path or a literal filename. "
        "Reselect and save this path on its original computer, or use a new "
        "workspace and configure local paths."
    )
