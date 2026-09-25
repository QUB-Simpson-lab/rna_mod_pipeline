"""Optional graphical launcher for the command-line pipeline.

Importing this package never requires PySide6. Call ``launch_gui`` to start the
desktop application.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


def launch_gui(
    project_root: str | Path,
    argv: Sequence[str] | None = None,
) -> int:
    """Import the Qt application only when the desktop interface is launched."""
    from .app import launch_gui as _launch_gui

    return _launch_gui(project_root, argv)


__all__ = ["launch_gui"]
