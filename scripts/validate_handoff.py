#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rna_mod_pipeline.audit import run_handoff_audit


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(
        description="Validate the historical project and current regression baselines."
    )
    command.add_argument(
        "--project-root",
        default=str(Path(__file__).resolve().parents[2]),
    )
    command.add_argument(
        "--skip-legacy-results",
        action="store_true",
        help="Check code and inputs without requiring the existing result folders.",
    )
    return command


def main() -> None:
    args = parser().parse_args()
    checks = run_handoff_audit(
        args.project_root,
        include_legacy_results=not args.skip_legacy_results,
    )
    width = max(len(check.check) for check in checks)
    for check in checks:
        status = "PASS" if check.passed else "FAIL"
        print(f"[{status}] {check.check:<{width}}  {check.details}")
    failures = [check for check in checks if not check.passed]
    print(f"\n{len(checks) - len(failures)}/{len(checks)} checks passed.")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
