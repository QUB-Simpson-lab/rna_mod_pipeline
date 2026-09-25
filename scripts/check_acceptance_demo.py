#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.acceptance.checks import check_demo


def main() -> None:
    parser = argparse.ArgumentParser(description="Check synthetic acceptance results against analytical expectations.")
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--skip-plots", action="store_true", help="Check numbers only, without requiring PNGs.")
    args = parser.parse_args()
    report = check_demo(args.workspace, plots=not args.skip_plots)
    print(f"PASS: synthetic numerical checks; {report['plots_checked']} PNGs verified.")
    print(f"Report: {Path(args.workspace).resolve() / 'acceptance_check.json'}")
    print("This verifies software behaviour, not biological validity or native GUI interaction.")


if __name__ == "__main__":
    main()
