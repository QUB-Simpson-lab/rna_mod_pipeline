from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from .comparison_profiles import (
    comparison_profile_errors,
    comparison_profile_payload,
)
from .datasets import DatasetProfile
from .resources import fasta_index_checks


SHELL_TOKENS = frozenset({"|", "||", "&&", ";", ">", ">>", "<", "2>", "&"})
INPUT_FLAGS = {
    "phase1_filter.py": ("--input",),
    "phase1_drach.py": ("--input", "--fasta"),
    "phase1_metagene.py": ("--input", "--gtf"),
    "loose_overlap.py": (
        "--sites", "--bedmethyl", "--binding-source", "--catalog",
    ),
    "transcript_region_overlap.py": (
        "--sites", "--bedmethyl", "--gtf", "--rbp-catalog",
        "--ornament-dir", "--encori-dir", "--postar3", "--validate-run",
    ),
    "plot_enrichment.py": ("--results", "--regional-results"),
    "cross_database.py": ("--ornament", "--encori", "--postar3"),
    "integrate_expression.py": (
        "--expression", "--enrichment", "--metagene", "--gene-lengths",
    ),
    "knockrbp_validation.py": (
        "--knockrbp-dir", "--metadata", "--encori-annotated",
        "--postar3-annotated", "--encori-enrichment",
        "--postar3-enrichment", "--transcript-run-dir", "--cross-database",
    ),
    "string_analysis.py": ("--cross-database",),
    "compare_datasets.py": ("--left-profile", "--right-profile"),
}


@dataclass(frozen=True)
class PreflightResult:
    errors: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return not self.errors


def _values(argv: Sequence[str], flag: str) -> tuple[str, ...]:
    return tuple(
        argv[index + 1]
        for index, value in enumerate(argv[:-1])
        if value == flag and not argv[index + 1].startswith("--")
    )


def _value(argv: Sequence[str], flag: str) -> str | None:
    values = _values(argv, flag)
    return values[-1] if values else None


def _resolve(raw: str, root: Path) -> Path:
    path = Path(raw).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def _group_paths(
    argv: Sequence[str], flag: str, width: int, path_index: int, root: Path
) -> tuple[Path, ...]:
    paths = []
    for index, value in enumerate(argv):
        if value != flag:
            continue
        group = argv[index + 1:index + 1 + width]
        if len(group) == width:
            paths.append(_resolve(group[path_index], root))
    return tuple(paths)


def _profile_table_errors(profile_path: Path) -> list[str]:
    try:
        payload = json.loads(profile_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"Cannot read comparison profile {profile_path}: {exc}"]
    tables = payload.get("tables")
    if not isinstance(tables, Mapping):
        return [f"Comparison profile has no tables object: {profile_path}"]
    errors = []
    for key in ("filtered_sites", "metagene_sites"):
        raw = tables.get(key)
        if not isinstance(raw, str) or not raw:
            errors.append(f"{profile_path} does not declare {key}")
            continue
        path = _resolve(raw, profile_path.parent)
        if not path.is_file():
            errors.append(f"Required {key} table not found: {path}")
    return errors


def _comparison_errors(profile_paths: tuple[Path, ...]) -> list[str]:
    errors = [
        message
        for profile in profile_paths
        if profile.is_file()
        for message in _profile_table_errors(profile)
    ]
    if errors or len(profile_paths) != 2 or not all(p.is_file() for p in profile_paths):
        return errors
    try:
        from rna_mod_pipeline.dataset_compare.profiles import (
            load_profile,
            validate_compatibility,
        )

        validate_compatibility(
            load_profile(profile_paths[0]), load_profile(profile_paths[1])
        )
    except (OSError, ValueError) as exc:
        errors.append(f"Dataset comparison is incompatible: {exc}")
    return errors


def _output_errors(
    argv: Sequence[str], root: Path, inputs: tuple[Path, ...]
) -> list[str]:
    overwrite = "--overwrite" in argv
    errors = []
    outputs = []
    for flag in ("--output", "--manifest"):
        if raw := _value(argv, flag):
            path = _resolve(raw, root)
            outputs.append(path)
            if path.exists() and not overwrite:
                errors.append(f"Output already exists; enable overwrite or change it: {path}")
    if raw := _value(argv, "--output-dir"):
        path = _resolve(raw, root)
        outputs.append(path)
        if path.exists() and not path.is_dir():
            errors.append(f"Output directory is an existing file: {path}")
        elif path.is_dir() and any(path.iterdir()) and not overwrite:
            errors.append(f"Output directory is not empty: {path}")
    for output in outputs:
        for source in inputs:
            if output == source:
                errors.append(f"Input and output resolve to the same path: {output}")
            elif source.is_dir() and output.is_relative_to(source):
                errors.append(f"Output is inside an input directory: {output}")
            elif output.is_dir() and source.is_relative_to(output):
                errors.append(f"Output directory contains an input: {source}")
    return errors


def preflight_command(
    command: Sequence[str],
    deferred_inputs: Sequence[str | Path] = (),
) -> PreflightResult:
    if isinstance(command, (str, bytes)):
        return PreflightResult(("Command must be an argv sequence, not shell text.",))
    argv = tuple(str(value) for value in command)
    errors: list[str] = []
    warnings: list[str] = []
    if len(argv) < 2:
        return PreflightResult(("Command must contain an executable and script.",))
    dangerous = sorted(SHELL_TOKENS.intersection(argv))
    if dangerous:
        errors.append(f"Shell control tokens are not allowed in argv: {dangerous}")
    executable = Path(argv[0]).expanduser()
    if (executable.is_absolute() or executable.parent != Path(".")):
        if not executable.is_file():
            errors.append(f"Python executable not found: {executable}")
    elif shutil.which(argv[0]) is None:
        errors.append(f"Executable is not available: {argv[0]}")
    script = Path(argv[1]).expanduser()
    if not script.is_file():
        errors.append(f"Workflow script not found: {script}")
    root = _resolve(_value(argv, "--project-root") or ".", Path.cwd())
    deferred = {_resolve(str(path), root) for path in deferred_inputs}
    input_paths = tuple(
        _resolve(raw, root)
        for flag in INPUT_FLAGS.get(script.name, ())
        for raw in _values(argv, flag)
    )
    if script.name == "compare_modifications.py":
        input_paths += _group_paths(argv, "--input", 3, 2, root)
        input_paths += _group_paths(argv, "--metagene", 2, 1, root)
    missing_allowed = "--allow-missing" in argv
    database_flags = {"--ornament", "--encori", "--postar3"}
    for flag in INPUT_FLAGS.get(script.name, ()):
        for raw in _values(argv, flag):
            path = _resolve(raw, root)
            if (
                not path.exists()
                and path not in deferred
                and not (missing_allowed and flag in database_flags)
            ):
                errors.append(f"Required input not found ({flag}): {path}")
    if script.name == "compare_modifications.py" and not missing_allowed:
        for path in input_paths:
            if not path.exists():
                errors.append(f"Required comparison input not found: {path}")
    if script.name == "string_analysis.py" and "--offline" in argv:
        cache = _value(argv, "--cache-dir")
        if not cache or not _resolve(cache, root).is_dir():
            errors.append("Offline STRING analysis requires an existing --cache-dir")
    if script.name == "phase1_drach.py":
        if raw_fasta := _value(argv, "--fasta"):
            fasta = _resolve(raw_fasta, root)
            if fasta.is_file():
                index_errors, index_warnings = fasta_index_checks(fasta)
                errors.extend(index_errors)
                warnings.extend(index_warnings)
    if script.name == "compare_datasets.py":
        profiles = tuple(
            _resolve(raw, root)
            for flag in ("--left-profile", "--right-profile")
            for raw in _values(argv, flag)
        )
        if not profiles or not all(path in deferred for path in profiles):
            errors.extend(_comparison_errors(profiles))
    errors.extend(_output_errors(argv, root, input_paths))
    if "--overwrite" in argv:
        warnings.append("Overwrite is enabled; CLI ownership checks still apply.")
    return PreflightResult(tuple(dict.fromkeys(errors)), tuple(warnings))


def validate_request(
    workflow: str,
    dataset: DatasetProfile,
    options: Mapping[str, object],
    command: Sequence[str],
) -> PreflightResult:
    warnings: list[str] = []
    deferred_inputs: tuple[Path, ...] = ()
    if workflow == "compare_datasets":
        second = options.get("second_dataset")
        if not isinstance(second, DatasetProfile):
            return PreflightResult(("Select a valid second dataset.",))
        root = _resolve(_value(tuple(map(str, command)), "--project-root") or ".", Path.cwd())
        first_resolved = dataset.resolved(root)
        second_resolved = second.resolved(root)
        if first_resolved.name.casefold() == second_resolved.name.casefold():
            return PreflightResult(("Compared datasets must be distinct.",))
        if first_resolved.modification != second_resolved.modification:
            return PreflightResult(
                ("Compared datasets must use the same modification.",)
            )
        first_reference = comparison_profile_payload(first_resolved)["reference"]
        second_reference = comparison_profile_payload(second_resolved)["reference"]
        assert isinstance(first_reference, Mapping)
        assert isinstance(second_reference, Mapping)
        for key, label in (
            ("genome_build", "genome build"),
            ("annotation_release", "annotation release"),
        ):
            first_value = str(first_reference.get(key, "")).strip().casefold()
            second_value = str(second_reference.get(key, "")).strip().casefold()
            if first_value != second_value:
                return PreflightResult(
                    (f"Compared datasets use different {label} values.",)
                )
        for key, label in (("fasta_sha256", "FASTA"), ("gtf_sha256", "GTF")):
            first_digest = first_reference.get(key)
            second_digest = second_reference.get(key)
            if first_digest and second_digest and first_digest != second_digest:
                return PreflightResult(
                    (f"Compared datasets have different {label} checksums.",)
                )
        for label, first_path, second_path in (
            ("FASTA", first_resolved.fasta, second_resolved.fasta),
            ("GTF", first_resolved.gtf, second_resolved.gtf),
        ):
            if first_path != second_path:
                warnings.append(
                    f"Compared datasets use different {label} paths. "
                    "Reference identity is only verified when profile checksums "
                    "are supplied."
                )
        source_errors = [
            *comparison_profile_errors(first_resolved),
            *comparison_profile_errors(second_resolved),
        ]
        if source_errors:
            return PreflightResult(tuple(dict.fromkeys(source_errors)), tuple(warnings))
        deferred_inputs = (
            first_resolved.output_root / "dataset_comparison_profile.json",
            second_resolved.output_root / "dataset_comparison_profile.json",
        )
    result = preflight_command(command, deferred_inputs=deferred_inputs)
    return PreflightResult(
        result.errors,
        tuple(dict.fromkeys([*warnings, *result.warnings])),
    )
