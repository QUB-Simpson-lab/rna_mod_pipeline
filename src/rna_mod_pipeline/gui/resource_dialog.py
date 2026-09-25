from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QVBoxLayout,
)

from .resources import (
    RESOURCE_KINDS,
    RESOURCE_LABELS,
    ResourceProfile,
)
from .widgets import PathEdit


class ResourceProfileDialog(QDialog):
    def __init__(
        self,
        project_root: Path,
        profile: ResourceProfile,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.project_root = project_root
        self.setWindowTitle("Pipeline resources")
        self.setMinimumWidth(760)
        self.paths: dict[str, PathEdit] = {}

        introduction = QLabel(
            "Configure shared references and RBP resources once. Paths inside "
            "the project are saved relatively so the project can be moved. "
            "Leave resources you do not use blank."
        )
        introduction.setWordWrap(True)
        form = QFormLayout()
        for key, kind in RESOURCE_KINDS.items():
            widget = PathEdit(kind)
            value = profile.path(key)
            if value is not None:
                widget.setText(str(value))
            self.paths[key] = widget
            form.addRow(RESOURCE_LABELS[key], widget)
        self.genome_build = QLineEdit(profile.genome_build)
        self.annotation_release = QLineEdit(profile.annotation_release)
        form.addRow("Genome build label", self.genome_build)
        form.addRow("Annotation release label", self.annotation_release)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addWidget(introduction)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def profile(self) -> ResourceProfile:
        return ResourceProfile(
            **{
                key: widget.text() or None
                for key, widget in self.paths.items()
            },
            genome_build=self.genome_build.text(),
            annotation_release=self.annotation_release.text(),
        ).resolved(self.project_root)

    def accept(self) -> None:
        profile = self.profile()
        errors = profile.validation_errors()
        if errors:
            QMessageBox.warning(
                self,
                "Invalid resources",
                "\n".join(errors),
            )
            return
        warnings = profile.validation_warnings()
        if warnings:
            QMessageBox.information(
                self,
                "Resource note",
                "\n".join(warnings),
            )
        super().accept()
