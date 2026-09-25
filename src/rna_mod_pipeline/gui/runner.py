from __future__ import annotations

import os
import signal
from pathlib import Path
from typing import Sequence

from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, QTimer, Signal

from .state import JobJournal, RunState, format_argv


class ProcessRunner(QObject):
    log_received = Signal(str)
    state_changed = Signal(str)
    completed = Signal(object)

    def __init__(self, parent: QObject | None = None):
        super().__init__(parent)
        self.process = QProcess(self)
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.SeparateChannels)
        self.process.readyReadStandardOutput.connect(self._read_stdout)
        self.process.readyReadStandardError.connect(self._read_stderr)
        self.process.started.connect(self._started)
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(self._error)
        self._force_stop_timer = QTimer(self)
        self._force_stop_timer.setSingleShot(True)
        self._force_stop_timer.timeout.connect(self._force_stop)
        self.journal: JobJournal | None = None
        self._cancel_requested = False
        self._finalized = False

    @property
    def busy(self) -> bool:
        return self.process.state() != QProcess.ProcessState.NotRunning

    def start(
        self,
        argv: Sequence[str],
        cwd: str | Path,
        log_root: str | Path,
        output_dir: str | Path | None = None,
    ) -> None:
        if self.busy:
            raise RuntimeError("A pipeline job is already running")
        if not argv:
            raise ValueError("Cannot run an empty command")
        self._force_stop_timer.stop()
        self.journal = JobJournal.create(log_root, argv, cwd, output_dir)
        self.journal.transition(RunState.STARTING)
        self._cancel_requested = False
        self._finalized = False
        self.state_changed.emit(RunState.STARTING.value)
        self.log_received.emit(f"$ {format_argv(argv)}")
        self.process.setWorkingDirectory(str(Path(cwd).expanduser().resolve()))
        self.process.setProcessEnvironment(QProcessEnvironment.systemEnvironment())
        self.process.setProgram(str(argv[0]))
        self.process.setArguments([str(value) for value in argv[1:]])
        self.process.start()

    def cancel(self) -> None:
        if not self.busy or self._cancel_requested:
            return
        self._cancel_requested = True
        if self.journal and self.journal.state == RunState.RUNNING:
            self.journal.transition(RunState.CANCELLING)
        self.state_changed.emit(RunState.CANCELLING.value)
        self.log_received.emit("Cancellation requested; waiting for graceful cleanup…")
        pid = int(self.process.processId())
        if os.name == "posix" and pid:
            try:
                os.kill(pid, signal.SIGINT)
            except OSError:
                self.process.terminate()
        else:
            self.process.terminate()
        self._force_stop_timer.start(5000)

    def _started(self) -> None:
        if self.journal and self.journal.state == RunState.STARTING:
            self.journal.transition(RunState.RUNNING)
            if self._cancel_requested:
                self.journal.transition(RunState.CANCELLING)
        self.state_changed.emit(
            RunState.CANCELLING.value
            if self._cancel_requested
            else RunState.RUNNING.value
        )

    def _emit_chunk(self, stream: str, content: bytes) -> None:
        text = content.decode("utf-8", errors="replace")
        if self.journal:
            self.journal.append(stream, text)
        for line in text.splitlines():
            self.log_received.emit(f"[{stream}] {line}")

    def _read_stdout(self) -> None:
        self._emit_chunk("stdout", bytes(self.process.readAllStandardOutput()))

    def _read_stderr(self) -> None:
        self._emit_chunk("stderr", bytes(self.process.readAllStandardError()))

    def _finished(self, exit_code: int, _status) -> None:
        self._read_stdout()
        self._read_stderr()
        state = (
            RunState.INTERRUPTED
            if self._cancel_requested
            else RunState.SUCCEEDED
            if exit_code == 0
            else RunState.FAILED
        )
        self._finalize(state, exit_code)

    def _error(self, error: QProcess.ProcessError) -> None:
        if error == QProcess.ProcessError.FailedToStart:
            self.log_received.emit(f"Process failed to start: {self.process.errorString()}")
            self._finalize(RunState.FAILED, None)

    def _force_stop(self) -> None:
        if self.busy:
            self.log_received.emit("Grace period expired; force-stopping the process.")
            self.process.kill()

    def _finalize(self, state: RunState, return_code: int | None) -> None:
        if self._finalized:
            return
        self._finalized = True
        self._force_stop_timer.stop()
        if self.journal:
            self.journal.finish(state, return_code)
            summary = self.journal.summary()
        else:
            summary = {"state": state.value, "return_code": return_code}
        self.state_changed.emit(state.value)
        self.completed.emit(summary)
