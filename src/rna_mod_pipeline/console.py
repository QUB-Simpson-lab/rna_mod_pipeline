from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from .gui.resources import RESOURCE_KINDS


def _default_project_root() -> Path:
    checkout = Path(__file__).resolve().parents[2]
    if (checkout / "pyproject.toml").is_file():
        return checkout.parent
    current = Path.cwd().resolve()
    if (current / "refactored_code").is_dir():
        return current
    if current.name == "refactored_code":
        return current.parent
    return current


def gui_main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Launch the RNA-modification pipeline desktop interface."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=_default_project_root(),
        help="Workspace directory used for data, datasets, and outputs.",
    )
    args = parser.parse_args(argv)
    from .gui import launch_gui

    return launch_gui(args.project_root)


def audit_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate the historical project and regression baselines."
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=_default_project_root(),
    )
    parser.add_argument(
        "--skip-legacy-results",
        action="store_true",
        help=(
            "Check code and historical project inputs without requiring existing "
            "result folders; this is not a portable installation-only check."
        ),
    )
    return parser


def audit_main(argv: Sequence[str] | None = None) -> int:
    parser = audit_parser()
    args = parser.parse_args(argv)
    from .audit import run_handoff_audit

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
    return int(bool(failures))


def doctor_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Check portable references and analysis resources without requiring "
            "the historical project datasets or results."
        )
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        default=_default_project_root(),
        help="Directory used for datasets, saved resources, and outputs.",
    )
    parser.add_argument(
        "--resource-profile",
        type=Path,
        help="Optional resource-profile JSON; defaults to the GUI profile.",
    )
    parser.add_argument(
        "--require",
        action="append",
        choices=tuple(RESOURCE_KINDS),
        default=[],
        help="Require an additional named resource; may be repeated.",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Require every supported shared resource.",
    )
    return parser


def doctor_main(argv: Sequence[str] | None = None) -> int:
    parser = doctor_parser()
    args = parser.parse_args(argv)
    from .doctor import CORE_RESOURCES, run_resource_doctor

    required = tuple(RESOURCE_KINDS) if args.full else (*CORE_RESOURCES, *args.require)
    checks = run_resource_doctor(
        args.project_root,
        profile_path=args.resource_profile,
        required_resources=required,
    )
    width = max(len(check.check) for check in checks)
    for check in checks:
        print(f"[{check.status}] {check.check:<{width}}  {check.details}")
    failures = [check for check in checks if check.failed]
    passed = sum(check.status == "PASS" for check in checks)
    skipped = sum(check.status == "SKIP" for check in checks)
    warnings = sum(check.status == "WARN" for check in checks)
    print(
        f"\n{passed} passed; {warnings} warning(s); {skipped} skipped; "
        f"{len(failures)} failed."
    )
    return int(bool(failures))
