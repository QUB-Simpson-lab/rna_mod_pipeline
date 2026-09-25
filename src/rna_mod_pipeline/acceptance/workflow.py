from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

from ..gui.command_builder import build_command
from ..gui.datasets import DatasetRegistry
from ..gui.resources import ResourceProfileStore
from .generate import DATASET_NAME, VERSION


def demo_commands(workspace: Path, *, plots: bool = True) -> list[tuple[str, list[str]]]:
    resources = ResourceProfileStore(workspace).load()
    dataset = DatasetRegistry.load(workspace).get(DATASET_NAME)
    steps = [("phase1_filter", {}), ("phase1_drach", {}), ("phase1_metagene", {})]
    steps += [("loose_overlap", {"database": db}) for db in ("ornament", "encori", "postar3")]
    steps += [("transcript_region_overlap", {}), ("transcript_validate", {})]
    steps += [("cross_database", {"design": design}) for design in ("loose", "transcript-region")]
    commands = []
    for workflow, settings in steps:
        if not plots and workflow != "transcript_validate":
            settings["no_plots" if workflow == "transcript_region_overlap" else "skip_plots"] = True
        command = build_command(dataset, workflow, settings, project_root=workspace,
                                resource_profile=resources)
        label = workflow + "_" + str(settings.get("database", settings.get("design", "")))
        commands.append((label.rstrip("_"), command))
    return commands


def run_demo(workspace: str | Path, *, plots: bool = True) -> Path:
    root = Path(workspace).expanduser().resolve()
    fixture = json.loads((root / "SYNTHETIC_FIXTURE.json").read_text(encoding="utf-8"))
    if fixture.get("fixture_version") != VERSION or not fixture.get("synthetic_only"):
        raise ValueError("Workspace is not a recognised synthetic acceptance fixture")
    output = DatasetRegistry.load(root).get(DATASET_NAME).output_root
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Demo results already exist; use a new demo workspace for an automatic run")
    logs = root / "acceptance_logs"
    logs.mkdir(exist_ok=True)
    commands = demo_commands(root, plots=plots)
    (logs / "commands.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")
    environment = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONDONTWRITEBYTECODE="1",
                       MPLCONFIGDIR=str(root / ".matplotlib"))
    for index, (label, command) in enumerate(commands, 1):
        rendered = subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)
        print(f"[{index}/{len(commands)}] SYNTHETIC: {label}\n{rendered}", flush=True)
        with (logs / f"{index:02}_{label}.log").open("w", encoding="utf-8") as handle:
            result = subprocess.run(command, cwd=root, env=environment, stdout=handle,
                                    stderr=subprocess.STDOUT, check=False)
        if result.returncode:
            raise RuntimeError(f"{label} failed with exit {result.returncode}; see {handle.name}")
    (logs / "run_settings.json").write_text(
        json.dumps({"plots": plots, "python": sys.version, "platform": sys.platform}, indent=2) + "\n",
        encoding="utf-8")
    return output
