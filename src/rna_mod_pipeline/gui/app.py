from __future__ import annotations

import sys
from pathlib import Path
from typing import Sequence


INSTALL_MESSAGE = (
    "The optional graphical launcher requires PySide6.\n"
    "Install it in the pipeline virtual environment with:\n"
    "  python -m pip install PySide6\n"
    "The command-line pipeline does not require PySide6."
)


def launch_gui(
    project_root: str | Path,
    argv: Sequence[str] | None = None,
) -> int:
    try:
        from PySide6.QtWidgets import QApplication, QMessageBox
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.startswith("PySide6"):
            print(INSTALL_MESSAGE, file=sys.stderr)
            return 2
        raise

    from .main_window import MainWindow

    application = QApplication.instance() or QApplication(list(argv or sys.argv))
    application.setApplicationName("RNA Modification Pipeline")
    application.setOrganizationName("Queen's University Belfast")
    try:
        window = MainWindow(Path(project_root).expanduser().resolve())
    except Exception as exc:
        QMessageBox.critical(
            None,
            "Unable to start launcher",
            str(exc),
        )
        return 1
    window.show()
    return application.exec()
