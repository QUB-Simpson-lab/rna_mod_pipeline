from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from .bedmethyl_inspection import inspect_bedmethyl
from .datasets import DatasetProfile
from .resources import ResourceProfile
from .widgets import PathEdit


def resolve_dialog_path(project_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else project_root / path).resolve()


class AddDatasetDialog(QDialog):
    def __init__(
        self,
        project_root: Path,
        parent=None,
        resources: ResourceProfile | None = None,
        initial: DatasetProfile | None = None,
        clone: bool = False,
    ):
        super().__init__(parent)
        self.project_root = project_root
        self.setWindowTitle(
            "Clone dataset" if clone else "Edit dataset" if initial else "Add dataset"
        )
        self.setMinimumWidth(620)

        self.name = QLineEdit()
        self.modification = QComboBox()
        self.modification.addItem("m6A", "m6a")
        self.modification.addItem("m5C", "m5c")
        self.modification.addItem("Pseudouridine (Ψ)", "pseu")
        self.bedmethyl = PathEdit("file")
        self.fasta = PathEdit("file")
        self.gtf = PathEdit("file")
        self.genome_build = QLineEdit(
            resources.genome_build if resources else "hg38"
        )
        self.annotation_release = QLineEdit(
            resources.annotation_release if resources else "GENCODE v44"
        )
        self.sample = QLineEdit()
        self.cell_line = QLineEdit()
        self.replicate = QLineEdit()
        self.output_root = PathEdit("directory")
        self.inspect_button = QPushButton("Inspect bedMethyl")
        self.inspection = QLabel(
            "The file will be checked before the dataset is registered."
        )
        self.inspection.setWordWrap(True)
        if resources is None:
            self.fasta.setText(str(project_root / "data" / "hg38.fa"))
            self.gtf.setText(str(project_root / "data" / "gencode.v44.annotation.gtf"))
        else:
            if resources.reference_fasta:
                self.fasta.setText(str(resources.reference_fasta))
            if resources.annotation_gtf:
                self.gtf.setText(str(resources.annotation_gtf))
        if initial is not None:
            metadata = dict(initial.sample_metadata or {})
            self.name.setText(initial.name + " copy" if clone else initial.name)
            index = self.modification.findData(initial.modification)
            self.modification.setCurrentIndex(max(0, index))
            self.bedmethyl.setText(str(initial.bedmethyl))
            if initial.fasta:
                self.fasta.setText(str(initial.fasta))
            if initial.gtf:
                self.gtf.setText(str(initial.gtf))
            self.genome_build.setText(metadata.get("genome_build", ""))
            self.annotation_release.setText(metadata.get("annotation_release", ""))
            self.sample.setText(metadata.get("sample", ""))
            self.cell_line.setText(metadata.get("cell_line", ""))
            self.replicate.setText(metadata.get("replicate", ""))
            if initial.output_root and not clone:
                self.output_root.setText(str(initial.output_root))

        form = QFormLayout()
        form.addRow("Dataset name *", self.name)
        form.addRow("Modification *", self.modification)
        form.addRow("bedMethyl *", self.bedmethyl)
        form.addRow("Reference FASTA *", self.fasta)
        form.addRow("GENCODE GTF *", self.gtf)
        form.addRow("Genome build", self.genome_build)
        form.addRow("Annotation release", self.annotation_release)
        form.addRow("Sample", self.sample)
        form.addRow("Cell line", self.cell_line)
        form.addRow("Replicate", self.replicate)
        form.addRow("Output root", self.output_root)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self.inspect_button)
        layout.addWidget(self.inspection)
        layout.addWidget(buttons)
        self.inspect_button.clicked.connect(self._inspect)

    def accept(self) -> None:
        if not self.name.text().strip():
            QMessageBox.warning(self, "Missing name", "Enter a dataset name.")
            return
        path = resolve_dialog_path(self.project_root, self.bedmethyl.text())
        if not path.is_file():
            QMessageBox.warning(
                self, "Invalid bedMethyl", "Select an existing bedMethyl file."
            )
            return
        for label, raw in (("Reference FASTA", self.fasta.text()), ("GTF", self.gtf.text())):
            if not raw or not resolve_dialog_path(self.project_root, raw).is_file():
                QMessageBox.warning(
                    self, f"Invalid {label}", f"Select an existing file for {label}."
                )
                return
        if not self._inspect():
            return
        super().accept()

    def _inspect(self) -> bool:
        try:
            result = inspect_bedmethyl(
                resolve_dialog_path(self.project_root, self.bedmethyl.text()),
                self.modification.currentData(),
            )
        except Exception as exc:
            self.inspection.setText(f"Inspection failed: {exc}")
            QMessageBox.warning(self, "Invalid bedMethyl", str(exc))
            return False
        self.inspection.setText(result.summary())
        return True

    def values(self) -> dict[str, str | None]:
        optional = lambda widget: widget.text() or None
        return {
            "name": self.name.text().strip(),
            "modification": self.modification.currentData(),
            "bedmethyl": self.bedmethyl.text(),
            "fasta": optional(self.fasta),
            "gtf": optional(self.gtf),
            "sample_metadata": {
                key: value
                for key, value in {
                    "sample": self.sample.text().strip(),
                    "cell_line": self.cell_line.text().strip(),
                    "replicate": self.replicate.text().strip(),
                    "genome_build": self.genome_build.text().strip(),
                    "annotation_release": self.annotation_release.text().strip(),
                }.items()
                if value
            },
            "output_root": optional(self.output_root),
        }
