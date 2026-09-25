from __future__ import annotations

import json
import os
import shlex
import subprocess
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Sequence


class RunState(str, Enum):
    IDLE = "idle"
    STARTING = "starting"
    RUNNING = "running"
    CANCELLING = "cancelling"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    INTERRUPTED = "interrupted"


_TRANSITIONS = {
    RunState.IDLE: {RunState.STARTING},
    RunState.STARTING: {RunState.RUNNING, RunState.FAILED, RunState.INTERRUPTED},
    RunState.RUNNING: {
        RunState.CANCELLING,
        RunState.SUCCEEDED,
        RunState.FAILED,
        RunState.INTERRUPTED,
    },
    RunState.CANCELLING: {RunState.INTERRUPTED, RunState.FAILED},
    RunState.SUCCEEDED: {RunState.STARTING},
    RunState.FAILED: {RunState.STARTING},
    RunState.INTERRUPTED: {RunState.STARTING},
}


def format_argv(argv: Sequence[str], platform: str | None = None) -> str:
    values = [str(value) for value in argv]
    target = (platform or os.name).casefold()
    if target in {"nt", "windows", "win32"}:
        return subprocess.list2cmdline(values)
    return shlex.join(values)


@dataclass
class JobJournal:
    directory: Path
    argv: list[str]
    cwd: str
    output_dir: str | None = None
    job_id: str = field(default_factory=lambda: uuid.uuid4().hex)
    state: RunState = RunState.IDLE
    started_utc: str | None = None
    finished_utc: str | None = None
    return_code: int | None = None

    @classmethod
    def create(
        cls,
        log_root: str | Path,
        argv: Sequence[str],
        cwd: str | Path,
        output_dir: str | Path | None = None,
    ) -> "JobJournal":
        job_id = uuid.uuid4().hex
        directory = Path(log_root).expanduser().resolve() / job_id
        directory.mkdir(parents=True, exist_ok=False)
        journal = cls(
            directory=directory,
            argv=[str(value) for value in argv],
            cwd=str(Path(cwd).expanduser().resolve()),
            output_dir=(
                str(Path(output_dir).expanduser().resolve())
                if output_dir is not None
                else None
            ),
            job_id=job_id,
        )
        journal._write_metadata()
        (directory / "command.txt").write_text(
            format_argv(journal.argv) + "\n",
            encoding="utf-8",
        )
        return journal

    @property
    def log_path(self) -> Path:
        return self.directory / "run.log"

    def transition(self, state: RunState) -> None:
        if state not in _TRANSITIONS[self.state]:
            raise ValueError(f"Invalid run-state transition: {self.state} -> {state}")
        self.state = state
        now = datetime.now(timezone.utc).isoformat()
        if state == RunState.STARTING:
            self.started_utc = now
        if state in {RunState.SUCCEEDED, RunState.FAILED, RunState.INTERRUPTED}:
            self.finished_utc = now
        self._write_metadata()

    def append(self, stream: str, text: str) -> None:
        if not text:
            return
        timestamp = datetime.now(timezone.utc).isoformat()
        with self.log_path.open("a", encoding="utf-8") as handle:
            for line in text.splitlines():
                handle.write(f"{timestamp}\t{stream}\t{line}\n")

    def finish(self, state: RunState, return_code: int | None) -> None:
        self.return_code = return_code
        self.transition(state)

    def summary(self) -> dict[str, object]:
        return {
            "job_id": self.job_id,
            "state": self.state.value,
            "argv": self.argv,
            "cwd": self.cwd,
            "output_dir": self.output_dir,
            "started_utc": self.started_utc,
            "finished_utc": self.finished_utc,
            "return_code": self.return_code,
            "log_path": str(self.log_path),
        }

    def _write_metadata(self) -> None:
        path = self.directory / "job.json"
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(
            json.dumps(self.summary(), indent=2) + "\n",
            encoding="utf-8",
        )
        temporary.replace(path)
