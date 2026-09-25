#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.acceptance.generate import prepare_demo
from rna_mod_pipeline.acceptance.workflow import run_demo


def main() -> None:
    parser = argparse.ArgumentParser(description="Create artificial software-test inputs; no biological data.")
    parser.add_argument("--workspace", required=True, help="New/empty folder outside the source checkout.")
    parser.add_argument("--run", action="store_true", help="Also execute the core CLI acceptance workflow.")
    parser.add_argument("--skip-plots", action="store_true", help="Faster numerical test; does not verify PNG generation.")
    args = parser.parse_args()
    if args.skip_plots and not args.run:
        parser.error("--skip-plots requires --run")
    root = prepare_demo(args.workspace)
    print(f"Synthetic inputs and GUI profiles prepared: {root}", flush=True)
    if args.run:
        run_demo(root, plots=not args.skip_plots)
    print("SOFTWARE TEST ONLY: do not interpret DEMO_* data as research findings.")


if __name__ == "__main__":
    main()
