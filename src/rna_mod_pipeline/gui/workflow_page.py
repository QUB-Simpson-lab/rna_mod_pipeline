from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from .adapters import (
    BuildPlan,
    build_plan,
    dataset_label,
    dataset_profiles,
    workflow_specs,
)
from .form_controls import DynamicForm
from .preflight_adapter import preflight_errors
from .comparison_profiles import prepare_comparison_profiles
from .resources import ResourceProfile
from .state import format_argv
from .workflow_registry import prerequisites_for


class WorkflowPage(QWidget):
    run_requested = Signal(object)
    cancel_requested = Signal()

    def __init__(
        self,
        project_root: Path,
        registry: Any,
        parent=None,
        resources: ResourceProfile | None = None,
    ):
        super().__init__(parent)
        self.project_root = project_root
        self.registry = registry
        self.resources = resources
        self.specs = workflow_specs()
        self.current_plan: BuildPlan | None = None

        self.dataset = QComboBox()
        self.second_dataset = QComboBox()
        self.additional_datasets = QListWidget()
        self.additional_datasets.setSelectionMode(
            QAbstractItemView.SelectionMode.MultiSelection
        )
        self.additional_datasets.setMaximumHeight(130)
        self.workflow = QComboBox()
        for spec in self.specs:
            self.workflow.addItem(spec.label, spec.key)
        self.description = QLabel()
        self.description.setWordWrap(True)
        self.form = DynamicForm()
        self.preview = QPlainTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setMaximumBlockCount(500)
        self.validation = QLabel()
        self.validation.setWordWrap(True)
        self.validate_button = QPushButton("Validate")
        self.run_button = QPushButton("Run")
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.status = QLabel("Status: idle")
        self.progress = QProgressBar()
        self.progress.setTextVisible(False)
        self.progress.setVisible(False)
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(10_000)

        selectors = QFormLayout()
        selectors.addRow("Dataset", self.dataset)
        selectors.addRow("Second dataset", self.second_dataset)
        self.second_dataset_label = selectors.labelForField(self.second_dataset)
        selectors.addRow("Additional datasets", self.additional_datasets)
        self.additional_datasets_label = selectors.labelForField(
            self.additional_datasets
        )
        selectors.addRow("Workflow", self.workflow)
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.addLayout(selectors)
        left_layout.addWidget(self.description)
        left_layout.addWidget(self.form)
        left_layout.addWidget(QLabel("Exact command"))
        left_layout.addWidget(self.preview)
        left_layout.addWidget(self.validation)
        buttons = QHBoxLayout()
        buttons.addWidget(self.validate_button)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        buttons.addStretch(1)
        left_layout.addLayout(buttons)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.addWidget(self.status)
        right_layout.addWidget(self.progress)
        right_layout.addWidget(self.log)
        splitter = QSplitter()
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([650, 650])
        layout = QVBoxLayout(self)
        layout.addWidget(splitter)

        self.workflow.currentIndexChanged.connect(self._build_form)
        self.dataset.currentIndexChanged.connect(self._update_preview)
        self.second_dataset.currentIndexChanged.connect(self._update_preview)
        self.additional_datasets.itemSelectionChanged.connect(self._update_preview)
        self.validate_button.clicked.connect(self._validate)
        self.run_button.clicked.connect(self._run)
        self.cancel_button.clicked.connect(self.cancel_requested.emit)
        self.refresh_datasets()
        self._build_form()

    def refresh_datasets(self) -> None:
        selected = self.dataset.currentText()
        second_selected = self.second_dataset.currentText()
        additional_selected = {
            item.text() for item in self.additional_datasets.selectedItems()
        }
        self.dataset.clear()
        self.second_dataset.clear()
        self.additional_datasets.clear()
        for profile in dataset_profiles(self.registry):
            label = dataset_label(profile)
            self.dataset.addItem(label, profile)
            self.second_dataset.addItem(label, profile)
            item = QListWidgetItem(label)
            item.setData(Qt.ItemDataRole.UserRole, profile)
            self.additional_datasets.addItem(item)
            if label in additional_selected:
                item.setSelected(True)
        index = self.dataset.findText(selected)
        self.dataset.setCurrentIndex(index if index >= 0 else 0)
        second_index = self.second_dataset.findText(second_selected)
        if second_index < 0 and self.second_dataset.count() > 1:
            second_index = 1
        self.second_dataset.setCurrentIndex(max(0, second_index))
        self._update_preview()

    def set_running(self, running: bool) -> None:
        self.run_button.setEnabled(not running)
        self.validate_button.setEnabled(not running)
        self.cancel_button.setEnabled(running)
        self.dataset.setEnabled(not running)
        self.second_dataset.setEnabled(not running)
        self.additional_datasets.setEnabled(not running)
        self.workflow.setEnabled(not running)
        self.form.setEnabled(not running)
        self.progress.setVisible(running)
        self.progress.setRange(0, 0 if running else 1)

    def set_status(self, status: str) -> None:
        self.status.setText(f"Status: {status}")

    def append_log(self, message: str) -> None:
        self.log.appendPlainText(message)

    def _build_form(self) -> None:
        if not self.specs or self.workflow.currentIndex() < 0:
            return
        spec = self.specs[self.workflow.currentIndex()]
        comparison = spec.key in {
            "compare_datasets",
            "dataset_comparison",
            "compare_modifications",
        }
        pairwise = spec.key in {"compare_datasets", "dataset_comparison"}
        multiple = spec.key == "compare_modifications"
        self.second_dataset.setVisible(pairwise)
        self.second_dataset_label.setVisible(pairwise)
        self.additional_datasets.setVisible(multiple)
        self.additional_datasets_label.setVisible(multiple)
        self.form.populate(spec.fields, self._update_preview)
        self._update_preview()

    def _selection(self) -> tuple[str, Any, dict[str, Any]]:
        if self.workflow.currentIndex() < 0:
            raise ValueError("No workflow is selected")
        dataset = self.dataset.currentData()
        if dataset is None:
            raise ValueError("Add and select a dataset first")
        values = self.form.values()
        workflow = self.workflow.currentData()
        if workflow in {"compare_datasets", "dataset_comparison"}:
            second = self.second_dataset.currentData()
            if second is None:
                raise ValueError("Select a second dataset")
            if second is dataset:
                raise ValueError("Choose two different datasets")
            values["second_dataset"] = second
        elif workflow == "compare_modifications":
            comparisons = [
                item.data(Qt.ItemDataRole.UserRole)
                for item in self.additional_datasets.selectedItems()
            ]
            comparisons = [item for item in comparisons if item is not dataset]
            if not comparisons:
                raise ValueError("Select at least one additional modification dataset")
            values["comparison_profiles"] = comparisons
        return workflow, dataset, values

    def set_resources(self, resources: ResourceProfile) -> None:
        self.resources = resources
        self._update_preview()

    def _update_description(
        self,
        key: str,
        dataset: Any,
        values: dict[str, Any],
    ) -> None:
        spec = next(item for item in self.specs if item.key == key)
        required = prerequisites_for(key, dataset.modification, values)
        labels = {
            item.key: item.label for item in self.specs
        }
        suffix = ""
        if required:
            suffix = "\nPrerequisites: " + ", ".join(
                labels.get(item, item) for item in required
            )
        self.description.setText(spec.description + suffix)

    def _update_preview(self, *_args) -> None:
        try:
            key, dataset, values = self._selection()
            self._update_description(key, dataset, values)
            self.current_plan = build_plan(
                key,
                dataset,
                values,
                project_root=self.project_root,
                resource_profile=self.resources,
            )
            self.preview.setPlainText(format_argv(self.current_plan.argv))
            self.validation.setText("")
        except Exception as exc:
            self.current_plan = None
            self.preview.setPlainText(f"Command unavailable: {exc}")

    def _validate_plan(self) -> tuple[BuildPlan, list[str]]:
        key, dataset, values = self._selection()
        missing = self.form.missing_required()
        if missing:
            raise ValueError("Required values are missing: " + ", ".join(missing))
        plan = build_plan(
            key,
            dataset,
            values,
            project_root=self.project_root,
            resource_profile=self.resources,
        )
        errors, warnings = preflight_errors(key, dataset, values, plan)
        if errors:
            raise ValueError("\n".join(errors))
        return plan, [*plan.warnings, *warnings]

    def _validate(self) -> None:
        try:
            plan, warnings = self._validate_plan()
        except Exception as exc:
            self.validation.setText(f"Validation failed: {exc}")
            return
        message = "Validation passed."
        if warnings:
            message += "\nWarnings:\n- " + "\n- ".join(warnings)
        self.validation.setText(message)
        self.preview.setPlainText(format_argv(plan.argv))

    def _run(self) -> None:
        try:
            plan, warnings = self._validate_plan()
        except Exception as exc:
            QMessageBox.critical(self, "Cannot start workflow", str(exc))
            return
        if warnings:
            message = "Preflight warnings:\n\n" + "\n".join(f"• {item}" for item in warnings)
            answer = QMessageBox.question(
                self,
                "Run with warnings?",
                message,
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return
        try:
            key, dataset, values = self._selection()
            if key in {"compare_datasets", "dataset_comparison"}:
                prepare_comparison_profiles(
                    dataset,
                    values["second_dataset"],
                    self.project_root,
                )
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Cannot prepare workflow",
                str(exc),
            )
            return
        self.log.clear()
        self.run_requested.emit(plan)
