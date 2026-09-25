#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.gui import launch_gui


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description="Launch the optional RNA-modification pipeline desktop interface"
    )
    command.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[2]),
        help="Workspace directory used for data, datasets, and outputs.",
    )
    return command


def main() -> None:
    args = parser().parse_args()
    raise SystemExit(launch_gui(args.project_root))


if __name__ == "__main__":
    main()
