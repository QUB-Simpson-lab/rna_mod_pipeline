from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QListWidget,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .adapters import (
    add_dataset,
    dataset_details,
    dataset_label,
    dataset_profiles,
    remove_dataset,
    replace_dataset,
)
from .dataset_dialog import AddDatasetDialog
from .resource_dialog import ResourceProfileDialog
from .resources import ResourceProfile, ResourceProfileStore, discover_project_resources


class DatasetPage(QWidget):
    datasets_changed = Signal()
    resources_changed = Signal(object)

    def __init__(
        self,
        project_root: Path,
        registry: Any,
        parent=None,
        resource_store: ResourceProfileStore | None = None,
        resources: ResourceProfile | None = None,
    ):
        super().__init__(parent)
        self.project_root = project_root
        self.registry = registry
        self.resource_store = resource_store or ResourceProfileStore(project_root)
        self.resources = resources or discover_project_resources(project_root)
        self.list_widget = QListWidget()
        self.details = QTableWidget(0, 2)
        self.details.setHorizontalHeaderLabels(["Field", "Value"])
        self.details.horizontalHeader().setStretchLastSection(True)
        self.add_button = QPushButton("Add dataset…")
        self.edit_button = QPushButton("Edit…")
        self.clone_button = QPushButton("Clone…")
        self.remove_button = QPushButton("Remove")
        self.resources_button = QPushButton("Resources…")
        self.check_resources_button = QPushButton("Check resources…")
        self.reload_button = QPushButton("Reload")

        heading = QLabel(f"Project: {project_root}")
        heading.setWordWrap(True)
        introduction = QLabel(
            "Add one bedMethyl callset at a time. The historical four-file "
            "project layout is not required; every dataset receives its own "
            "protected output workspace."
        )
        introduction.setWordWrap(True)
        buttons = QHBoxLayout()
        buttons.addWidget(self.add_button)
        buttons.addWidget(self.edit_button)
        buttons.addWidget(self.clone_button)
        buttons.addWidget(self.remove_button)
        buttons.addWidget(self.resources_button)
        buttons.addWidget(self.check_resources_button)
        buttons.addWidget(self.reload_button)
        buttons.addStretch(1)
        layout = QVBoxLayout(self)
        layout.addWidget(heading)
        layout.addWidget(introduction)
        layout.addLayout(buttons)
        layout.addWidget(self.list_widget, 2)
        layout.addWidget(self.details, 1)

        self.add_button.clicked.connect(self._add)
        self.edit_button.clicked.connect(self._edit)
        self.clone_button.clicked.connect(self._clone)
        self.remove_button.clicked.connect(self._remove)
        self.resources_button.clicked.connect(self._resources)
        self.check_resources_button.clicked.connect(self._check_resources)
        self.reload_button.clicked.connect(self.reload)
        self.list_widget.currentRowChanged.connect(self._show_details)
        self.reload()

    def profiles(self) -> list[Any]:
        return dataset_profiles(self.registry)

    def reload(self) -> None:
        selected = self.list_widget.currentRow()
        self.list_widget.clear()
        for profile in self.profiles():
            self.list_widget.addItem(dataset_label(profile))
        if self.list_widget.count():
            self.list_widget.setCurrentRow(max(0, min(selected, self.list_widget.count() - 1)))
        else:
            self.details.setRowCount(0)
        self.datasets_changed.emit()

    def _show_details(self, row: int) -> None:
        profiles = self.profiles()
        if row < 0 or row >= len(profiles):
            self.details.setRowCount(0)
            return
        values = dataset_details(profiles[row])
        self.details.setRowCount(len(values))
        for index, (name, value) in enumerate(values.items()):
            self.details.setItem(index, 0, QTableWidgetItem(name))
            self.details.setItem(index, 1, QTableWidgetItem(value))

    def _add(self) -> None:
        dialog = AddDatasetDialog(
            self.project_root,
            self,
            resources=self.resources,
        )
        if not dialog.exec():
            return
        try:
            add_dataset(self.registry, dialog.values())
        except Exception as exc:
            QMessageBox.critical(self, "Could not add dataset", str(exc))
            return
        self.reload()
        self.list_widget.setCurrentRow(self.list_widget.count() - 1)

    def _selected(self) -> Any | None:
        row = self.list_widget.currentRow()
        profiles = self.profiles()
        return profiles[row] if 0 <= row < len(profiles) else None

    def _edit(self) -> None:
        profile = self._selected()
        if profile is None:
            QMessageBox.information(self, "No dataset", "Select a dataset first.")
            return
        dialog = AddDatasetDialog(
            self.project_root,
            self,
            resources=self.resources,
            initial=profile,
        )
        if not dialog.exec():
            return
        try:
            replace_dataset(self.registry, profile.name, dialog.values())
        except Exception as exc:
            QMessageBox.critical(self, "Could not edit dataset", str(exc))
            return
        self.reload()

    def _clone(self) -> None:
        profile = self._selected()
        if profile is None:
            QMessageBox.information(self, "No dataset", "Select a dataset first.")
            return
        dialog = AddDatasetDialog(
            self.project_root,
            self,
            resources=self.resources,
            initial=profile,
            clone=True,
        )
        if not dialog.exec():
            return
        try:
            add_dataset(self.registry, dialog.values())
        except Exception as exc:
            QMessageBox.critical(self, "Could not clone dataset", str(exc))
            return
        self.reload()

    def _remove(self) -> None:
        profile = self._selected()
        if profile is None:
            QMessageBox.information(self, "No dataset", "Select a dataset first.")
            return
        answer = QMessageBox.question(
            self,
            "Remove dataset?",
            f"Remove {profile.name!r} from the registry?\n\n"
            "No bedMethyl files or analysis outputs will be deleted.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        try:
            remove_dataset(self.registry, profile.name)
        except Exception as exc:
            QMessageBox.critical(self, "Could not remove dataset", str(exc))
            return
        self.reload()

    def _resources(self) -> None:
        dialog = ResourceProfileDialog(
            self.project_root,
            self.resources,
            self,
        )
        if not dialog.exec():
            return
        profile = dialog.profile()
        try:
            self.resource_store.save(profile)
        except Exception as exc:
            QMessageBox.critical(self, "Could not save resources", str(exc))
            return
        self.resources = profile
        self.resources_changed.emit(profile)

    def _check_resources(self) -> None:
        from ..doctor import run_resource_doctor

        checks = run_resource_doctor(self.project_root)
        failures = [check for check in checks if check.failed]
        warnings = [check for check in checks if check.status == "WARN"]
        details = "\n".join(
            f"[{check.status}] {check.check}: {check.details}" for check in checks
        )
        details += (
            "\n\nCommand-line equivalent:\n"
            f"rna-mod-doctor --project-root {self.project_root}"
        )
        message = QMessageBox(self)
        message.setDetailedText(details)
        if failures:
            message.setIcon(QMessageBox.Icon.Warning)
            message.setWindowTitle("Resource check failed")
            message.setText(f"{len(failures)} resource check(s) failed.")
        elif warnings:
            message.setIcon(QMessageBox.Icon.Information)
            message.setWindowTitle("Resources pass with warnings")
            message.setText(
                f"Resources passed with {len(warnings)} warning(s)."
            )
        else:
            message.setIcon(QMessageBox.Icon.Information)
            message.setWindowTitle("Resources passed")
            message.setText("All configured resource checks passed.")
        message.exec()
