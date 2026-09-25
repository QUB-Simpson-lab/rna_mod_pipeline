from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFileDialog, QHBoxLayout, QLineEdit, QPushButton, QWidget


class PathEdit(QWidget):
    changed = Signal(str)

    def __init__(self, mode: str = "file", parent: QWidget | None = None):
        super().__init__(parent)
        self.mode = mode
        self.edit = QLineEdit()
        self.button = QPushButton("Browse…")
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.edit, 1)
        layout.addWidget(self.button)
        self.edit.textChanged.connect(self.changed.emit)
        self.button.clicked.connect(self._browse)

    def text(self) -> str:
        return self.edit.text().strip()

    def setText(self, value: str) -> None:
        self.edit.setText(value)

    def _browse(self) -> None:
        current = self.text()
        start = str(Path(current).expanduser().parent) if current else ""
        if self.mode in {"directory", "dir", "output"}:
            selected = QFileDialog.getExistingDirectory(self, "Select directory", start)
        else:
            selected, _ = QFileDialog.getOpenFileName(self, "Select file", start)
        if selected:
            self.setText(selected)
