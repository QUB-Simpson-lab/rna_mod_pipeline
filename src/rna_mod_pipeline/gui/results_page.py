from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QHBoxLayout,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .manifest_tools import ManifestView, fast_validate, load_manifest, text_preview
from .image_preview import ImagePreview
from .widgets import PathEdit


class ResultsPage(QWidget):
    def __init__(self, project_root: Path, parent=None):
        super().__init__(parent)
        self.project_root = project_root
        self.view: ManifestView | None = None
        self.path = PathEdit("directory")
        self.load_button = QPushButton("Load manifest")
        self.open_button = QPushButton("Open output folder")
        self.validate_button = QPushButton("Fast validation")
        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Output", "Size", "Status"])
        self.tree.header().setStretchLastSection(True)
        self.text_preview = QPlainTextEdit()
        self.text_preview.setReadOnly(True)
        self.image = ImagePreview()
        self.preview_splitter = QSplitter(Qt.Orientation.Vertical)
        self.preview_splitter.addWidget(self.text_preview)
        self.preview_splitter.addWidget(self.image)
        self.text_preview.hide()

        controls = QHBoxLayout()
        controls.addWidget(self.path, 1)
        controls.addWidget(self.load_button)
        controls.addWidget(self.validate_button)
        controls.addWidget(self.open_button)
        left = QSplitter(Qt.Orientation.Vertical)
        left.addWidget(self.summary)
        left.addWidget(self.tree)
        body = QSplitter()
        body.addWidget(left)
        body.addWidget(self.preview_splitter)
        body.setSizes([600, 700])
        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(body)

        self.load_button.clicked.connect(self._load_selected)
        self.open_button.clicked.connect(self._open_folder)
        self.validate_button.clicked.connect(self._validate)
        self.tree.currentItemChanged.connect(self._preview)

    def load_output(self, value: str | Path) -> None:
        self.path.setText(str(value))
        self._load_selected()

    def _load_selected(self) -> None:
        self.view = None
        self.tree.clear()
        self.summary.clear()
        self._clear_preview()
        try:
            self.view = load_manifest(self.path.text(), self.project_root)
        except Exception as exc:
            QMessageBox.critical(self, "Cannot load manifest", str(exc))
            return
        summary = {
            "workflow": self.view.workflow,
            "created_utc": self.view.created_utc,
            "manifest": str(self.view.path),
            "parameters": self.view.parameters,
        }
        self.summary.setPlainText(json.dumps(summary, indent=2, default=str))
        for record in self.view.outputs:
            exists = record.path.is_file()
            size = record.path.stat().st_size if exists else record.bytes
            item = QTreeWidgetItem(
                [
                    record.path.name,
                    f"{size:,}" if size is not None else "",
                    "available" if exists else "missing",
                ]
            )
            item.setData(0, Qt.ItemDataRole.UserRole, str(record.path))
            self.tree.addTopLevelItem(item)
        self.tree.resizeColumnToContents(0)

    def _validate(self) -> None:
        if self.view is None:
            self._load_selected()
        if self.view is None:
            return
        errors = fast_validate(self.view)
        if errors:
            QMessageBox.warning(
                self, "Validation problems", "\n".join(f"• {item}" for item in errors)
            )
        else:
            QMessageBox.information(
                self,
                "Fast validation passed",
                "All manifest outputs exist and their recorded sizes match.",
            )

    def _open_folder(self) -> None:
        if self.view is None:
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.view.path.parent)))

    def _clear_preview(self) -> None:
        self.text_preview.clear()
        self.text_preview.hide()
        self.image.clear()
        self.image.show()

    def _preview(self, current: QTreeWidgetItem | None, _previous) -> None:
        self._clear_preview()
        if current is None:
            return
        path = Path(str(current.data(0, Qt.ItemDataRole.UserRole)))
        is_image = path.suffix.lower() in {".png", ".jpg", ".jpeg", ".bmp"}
        if is_image and path.is_file():
            self.image.set_image(path)
            return
        self.image.hide()
        self.text_preview.show()
        if not path.is_file():
            self.text_preview.setPlainText("Output file is missing.")
        elif path.suffix.lower() in {".tsv", ".csv", ".txt", ".md", ".json"}:
            try:
                self.text_preview.setPlainText(text_preview(path))
            except Exception as exc:
                self.text_preview.setPlainText(f"Preview failed: {exc}")
        else:
            self.text_preview.setPlainText(
                f"{path}\n\nBinary or unsupported preview format."
            )
