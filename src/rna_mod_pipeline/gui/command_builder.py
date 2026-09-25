from __future__ import annotations
import sys
from dataclasses import replace
from pathlib import Path
from typing import Mapping, Sequence
from .command_downstream import build_downstream
from .datasets import DatasetProfile
from .resources import ResourceProfile
from .workflow_registry import FORWARDED_OPTIONS, WorkflowRegistry
SHORT_NAMES = {"m6a": "m6A", "m5c": "m5C", "pseu": "pseU"}
DATABASES = ("ornament", "encori", "postar3")
ENRICHMENT_NAMES = {
    "ornament": "rbp_enrichment_results.tsv",
    "encori": "encori_enrichment_results.tsv",
    "postar3": "postar3_enrichment_results.tsv",
}
def _paths(profile: DatasetProfile) -> dict[str, Path]:
    assert profile.output_root and profile.fasta and profile.gtf
    short = SHORT_NAMES[profile.modification]
    phase = profile.output_root / "phase1"
    filtered = phase / f"filtered_{short}.tsv"
    drach = phase / f"filtered_{short}_drach.tsv"
    return {"phase": phase, "filtered": filtered, "drach": drach,
            "metagene": phase / f"filtered_{short}_metagene.tsv",
            "phase_plots": phase / "plots",
            "transcript": profile.output_root / "transcript_region"}
def _value(options: Mapping[str, object], key: str, default: Path | str) -> str:
    return str(options.get(key) or default)
def _resource(
    options: Mapping[str, object],
    option_key: str,
    profile: ResourceProfile | None,
    resource_key: str,
    legacy_default: Path,
) -> Path | str:
    if options.get(option_key):
        return str(options[option_key])
    if profile is None:
        return legacy_default
    return profile.require(resource_key)
def _add(command: list[str], flag: str, value: object) -> None:
    if value is None or value is False:
        return
    if isinstance(value, str) and not value.strip():
        return
    if isinstance(value, (list, tuple, set, frozenset)) and not value:
        return
    command.append(flag)
    if value is True:
        return
    if isinstance(value, (list, tuple, set, frozenset)):
        command.extend(str(item) for item in value)
    else:
        command.append(str(value))
def _forward(
    command: list[str],
    workflow: str,
    options: Mapping[str, object],
    consumed: set[str],
) -> None:
    for key in FORWARDED_OPTIONS[workflow]:
        if key in consumed or key not in options:
            continue
        value = options[key]
        if key in {"cell_types", "methods"} and isinstance(value, Sequence) and not isinstance(value, str):
            value = ",".join(str(item) for item in value)
        _add(command, "--" + key.replace("_", "-"), value)
def _enrichment(profile: DatasetProfile, design: str, database: str, context: str) -> Path:
    assert profile.output_root
    if design == "loose":
        return profile.output_root / "loose" / database / ENRICHMENT_NAMES[database]
    table = "stratified_enrichment_results.tsv" if context == "all" else "region_specific_results.tsv"
    return profile.output_root / "transcript_region" / table
def build_command(
    profile: DatasetProfile | str,
    workflow: str | DatasetProfile,
    options: Mapping[str, object] | None = None,
    refactored_root: str | Path | None = None,
    project_root: str | Path | None = None,
    resource_profile: ResourceProfile | None = None,
) -> list[str]:
    settings = dict(options or {})
    if isinstance(profile, str):
        profile, workflow = workflow, profile
    if not isinstance(profile, DatasetProfile) or not isinstance(workflow, str):
        raise TypeError("build_command requires a DatasetProfile and workflow key")
    code_root = Path(refactored_root or Path(__file__).resolve().parents[3]).expanduser().resolve()
    data_root = Path(project_root or code_root.parent).expanduser().resolve()
    resources = (
        resource_profile.resolved(data_root)
        if resource_profile is not None
        else None
    )
    replacements: dict[str, Path] = {}
    if profile.fasta is None and resources and resources.reference_fasta:
        replacements["fasta"] = resources.reference_fasta
    if profile.gtf is None and resources and resources.annotation_gtf:
        replacements["gtf"] = resources.annotation_gtf
    profile = replace(profile, **replacements) if replacements else profile
    if resources is not None:
        if workflow == "phase1_drach" and profile.fasta is None:
            resources.require("reference_fasta")
        if workflow in {"phase1_metagene", "transcript_region_overlap"} and profile.gtf is None:
            resources.require("annotation_gtf")
    dataset = profile.resolved(data_root)
    spec = WorkflowRegistry().get(workflow)
    if dataset.modification not in spec.modifications:
        raise ValueError(f"{workflow} does not support {dataset.modification}")
    command = [
        str(settings.pop("python_executable", sys.executable)),
        str(code_root / "scripts" / spec.script),
        "--project-root", str(data_root),
    ]
    consumed = {"python_executable"}
    p = _paths(dataset)
    phase_plots = _value(settings, "plot_dir", p["phase_plots"])

    if workflow == "phase1_filter":
        command += ["--modification", dataset.modification, "--input", str(dataset.bedmethyl),
                    "--output", _value(settings, "output", p["filtered"]),
                    "--plot-dir", phase_plots, "--manifest",
                    _value(settings, "manifest", p["phase"] / "phase1_filter_manifest.json")]
        consumed |= {"output", "plot_dir", "manifest"}
    elif workflow == "phase1_drach":
        command += ["--modification", "m6a", "--input", _value(settings, "input", p["filtered"]),
                    "--fasta", str(dataset.fasta), "--output", _value(settings, "output", p["drach"]),
                    "--plot-dir", phase_plots, "--manifest",
                    _value(settings, "manifest", p["phase"] / "phase1_drach_manifest.json")]
        consumed |= {"input", "output", "plot_dir", "manifest"}
    elif workflow == "phase1_metagene":
        source = p["drach"] if dataset.modification == "m6a" else p["filtered"]
        command += ["--modification", dataset.modification, "--input", _value(settings, "input", source),
                    "--gtf", str(dataset.gtf), "--output", _value(settings, "output", p["metagene"]),
                    "--plot-dir", phase_plots, "--manifest",
                    _value(settings, "manifest", p["phase"] / "phase1_metagene_manifest.json")]
        consumed |= {"input", "output", "plot_dir", "manifest"}
    elif workflow == "loose_overlap":
        database = str(settings.get("database", ""))
        if database not in DATABASES:
            raise ValueError("loose_overlap requires a database")
        resource_key = {
            "ornament": "ornament_dir",
            "encori": "encori_dir",
            "postar3": "postar3",
        }[database]
        legacy_source = {
            "ornament": data_root / "data/HS",
            "encori": data_root / "data/ENCORI",
            "postar3": data_root / "data/human.txt",
        }[database]
        binding_source = _resource(
            settings,
            "binding_source",
            resources,
            resource_key,
            legacy_source,
        )
        output = dataset.output_root / "loose" / database
        command += ["--modification", dataset.modification, "--database", database,
                    "--sites", _value(settings, "sites", p["metagene"]),
                    "--bedmethyl", str(dataset.bedmethyl), "--binding-source",
                    str(binding_source), "--catalog",
                    _value(
                        settings,
                        "catalog",
                        resources.rbp_catalog
                        if resources and resources.rbp_catalog
                        else code_root / "config/rbp_catalog.tsv",
                    ),
                    "--output-dir", _value(settings, "output_dir", output), "--plot-dir",
                    _value(settings, "plot_dir", output / "plots")]
        consumed |= {"database", "sites", "binding_source", "catalog", "output_dir", "plot_dir"}
    elif workflow == "transcript_region_overlap":
        databases = tuple(settings.get("databases", DATABASES))
        command += ["--modification", dataset.modification, "--sites",
                    _value(settings, "sites", p["metagene"]), "--bedmethyl", str(dataset.bedmethyl),
                    "--gtf", str(dataset.gtf), "--rbp-catalog",
                    _value(
                        settings,
                        "rbp_catalog",
                        resources.rbp_catalog
                        if resources and resources.rbp_catalog
                        else code_root / "config/rbp_catalog.tsv",
                    ),
                    "--databases", *databases, "--output-dir",
                    _value(settings, "output_dir", p["transcript"])]
        source_options = {
            "ornament": (
                "--ornament-dir", "ornament_dir", data_root / "data/HS"
            ),
            "encori": (
                "--encori-dir", "encori_dir", data_root / "data/ENCORI"
            ),
            "postar3": (
                "--postar3", "postar3", data_root / "data/human.txt"
            ),
        }
        for database in databases:
            if database not in source_options:
                raise ValueError(f"Unknown binding database: {database}")
            flag, resource_key, legacy_source = source_options[database]
            option_key = flag[2:].replace("-", "_")
            source = _resource(
                settings,
                option_key,
                resources,
                resource_key,
                legacy_source,
            )
            _add(command, flag, source)
        consumed |= {"databases", "sites", "rbp_catalog", "output_dir", "ornament_dir", "encori_dir", "postar3"}
    elif workflow == "transcript_validate":
        command += ["--validate-run", _value(settings, "validate_run", p["transcript"])]
        consumed.add("validate_run")
    elif workflow == "plot_enrichment":
        command += ["--modification", dataset.modification, "--results",
                    _value(settings, "results", p["transcript"] / "stratified_enrichment_results.tsv"),
                    "--regional-results", _value(settings, "regional_results", p["transcript"] / "region_specific_results.tsv"),
                    "--output-dir", _value(settings, "output_dir", p["transcript"] / "replotted")]
        consumed |= {"results", "regional_results", "output_dir"}
    else:
        command, consumed = build_downstream(
            command,
            dataset,
            workflow,
            settings,
            code_root,
            data_root,
            p,
            resources,
        )
    _forward(command, workflow, settings, consumed)
    return command
