from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QStandardPaths
from PySide6.QtWidgets import QMainWindow, QMessageBox, QTabWidget

from .adapters import BuildPlan, load_dataset_registry
from .dataset_page import DatasetPage
from .manifest_tools import validate_completed_manifest
from .results_page import ResultsPage
from .resources import ResourceProfileStore
from .runner import ProcessRunner
from .workflow_page import WorkflowPage


class MainWindow(QMainWindow):
    def __init__(self, project_root: Path):
        super().__init__()
        self.project_root = project_root
        self.registry = load_dataset_registry(project_root)
        self.resource_store = ResourceProfileStore(project_root)
        self.resource_profile = self.resource_store.load()
        self.current_plan: BuildPlan | None = None
        self.runner = ProcessRunner(self)

        self.datasets = DatasetPage(
            project_root,
            self.registry,
            resource_store=self.resource_store,
            resources=self.resource_profile,
        )
        self.workflows = WorkflowPage(
            project_root,
            self.registry,
            resources=self.resource_profile,
        )
        self.results = ResultsPage(project_root)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.datasets, "Projects and datasets")
        self.tabs.addTab(self.workflows, "Run workflow")
        self.tabs.addTab(self.results, "Results")
        self.setCentralWidget(self.tabs)
        self.setWindowTitle("RNA Modification Pipeline")
        self.resize(1320, 820)

        self.datasets.datasets_changed.connect(self.workflows.refresh_datasets)
        self.datasets.resources_changed.connect(self.workflows.set_resources)
        self.workflows.run_requested.connect(self._run)
        self.workflows.cancel_requested.connect(self.runner.cancel)
        self.runner.log_received.connect(self.workflows.append_log)
        self.runner.state_changed.connect(self._state_changed)
        self.runner.completed.connect(self._completed)

    def _run(self, plan: BuildPlan) -> None:
        self.current_plan = plan
        log_root = Path(
            QStandardPaths.writableLocation(
                QStandardPaths.StandardLocation.AppLocalDataLocation
            )
        ) / "runs"
        try:
            self.runner.start(
                plan.argv,
                self.project_root,
                log_root,
                plan.output_dir,
            )
        except Exception as exc:
            QMessageBox.critical(self, "Could not start workflow", str(exc))

    def _state_changed(self, state: str) -> None:
        running = state in {"starting", "running", "cancelling"}
        self.workflows.set_running(running)
        self.workflows.set_status(state)

    def _completed(self, summary: dict[str, object]) -> None:
        state = str(summary.get("state", "unknown"))
        log_path = summary.get("log_path")
        self.workflows.append_log(f"Job finished with status: {state}")
        if log_path:
            self.workflows.append_log(f"Run log: {log_path}")
        if state != "succeeded" or not self.current_plan:
            return
        manifest = self.current_plan.manifest_path
        if manifest is None:
            return
        _view, errors = validate_completed_manifest(manifest, self.project_root)
        if errors:
            self.workflows.set_status("validation failed")
            self.workflows.append_log(
                "Post-run validation failed:\n- " + "\n- ".join(errors)
            )
            return
        self.workflows.set_status("succeeded and validated")
        self.results.load_output(manifest)
        self.tabs.setCurrentWidget(self.results)

    def closeEvent(self, event) -> None:
        if not self.runner.busy:
            event.accept()
            return
        answer = QMessageBox.question(
            self,
            "Job still running",
            "Cancel the running job and close the launcher?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if answer == QMessageBox.StandardButton.Yes:
            self.runner.cancel()
        event.ignore()
