from __future__ import annotations

from pathlib import Path
from typing import Mapping

from .command_dataset_compare import build_dataset_comparison
from .datasets import DatasetProfile
from .resources import ResourceProfile


DATABASES = ("ornament", "encori", "postar3")
SHORT_NAMES = {"m6a": "m6A", "m5c": "m5C", "pseu": "pseU"}
ENRICHMENT_NAMES = {
    "ornament": "rbp_enrichment_results.tsv",
    "encori": "encori_enrichment_results.tsv",
    "postar3": "postar3_enrichment_results.tsv",
}


def _value(options: Mapping[str, object], key: str, default: Path | str) -> str:
    return str(options.get(key) or default)


def _add(command: list[str], flag: str, value: object) -> None:
    if value is None or value is False:
        return
    if isinstance(value, str) and not value.strip():
        return
    if isinstance(value, (list, tuple, set, frozenset)) and not value:
        return
    command.append(flag)
    if value is not True:
        if isinstance(value, (list, tuple, set, frozenset)):
            command.extend(str(item) for item in value)
        else:
            command.append(str(value))


def _resource(
    settings: Mapping[str, object],
    option_key: str,
    resources: ResourceProfile | None,
    resource_key: str,
    legacy_default: Path,
) -> Path | str:
    if settings.get(option_key):
        return str(settings[option_key])
    if resources is None:
        return legacy_default
    return resources.require(resource_key)


def _enrichment(
    profile: DatasetProfile,
    design: str,
    database: str,
    context: str,
) -> Path:
    assert profile.output_root
    if design == "loose":
        return (
            profile.output_root / "loose" / database / ENRICHMENT_NAMES[database]
        )
    table = (
        "stratified_enrichment_results.tsv"
        if context == "all"
        else "region_specific_results.tsv"
    )
    return profile.output_root / "transcript_region" / table


def build_downstream(
    command: list[str],
    dataset: DatasetProfile,
    workflow: str,
    settings: Mapping[str, object],
    code_root: Path,
    project_root: Path,
    paths: Mapping[str, Path],
    resources: ResourceProfile | None = None,
) -> tuple[list[str], set[str]]:
    del code_root
    assert dataset.output_root
    consumed: set[str] = set()
    design = str(settings.get("design", "loose"))
    context = str(settings.get("context", "all"))
    if design == "loose" and context != "all":
        raise ValueError(
            "Loose-design downstream workflows support only all context; "
            "use transcript-region design for 5UTR, CDS, or 3UTR."
        )
    databases = tuple(settings.get("databases", DATABASES))
    if workflow == "cross_database":
        output = dataset.output_root / "cross_database" / f"{design}_{context}"
        command += [
            "--modification", dataset.modification, "--design", design,
            "--context", context, "--databases", *databases, "--output-dir",
            _value(settings, "output_dir", output),
        ]
        for database in databases:
            _add(
                command,
                f"--{database}",
                settings.get(database)
                or _enrichment(dataset, design, database, context),
            )
        consumed |= {"design", "context", "databases", "output_dir", *DATABASES}
    elif workflow == "integrate_expression":
        source = str(settings.get("source", ""))
        if source not in {"nanopore", "depmap"}:
            raise ValueError("Expression source must be nanopore or depmap")
        database = settings.get("database")
        if design == "loose" and database not in DATABASES:
            raise ValueError("Loose expression integration requires a database")
        expression = _resource(
            settings,
            "expression",
            resources,
            "nanopore_expression" if source == "nanopore" else "depmap_expression",
            project_root
            / "data"
            / (
                "Kate_231_0h_vs_6h_gene_counts_normalised.tsv"
                if source == "nanopore"
                else "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
            ),
        )
        enrichment = _enrichment(
            dataset, design, str(database or "encori"), context
        )
        output = dataset.output_root / "expression" / (
            f"{source}_{design}_{database or 'all'}_{context}"
        )
        command += [
            "--modification", dataset.modification, "--source", source,
            "--design", design, "--context", context, "--expression",
            str(expression), "--enrichment",
            _value(settings, "enrichment", enrichment), "--metagene",
            _value(settings, "metagene", paths["metagene"]), "--output-dir",
            _value(settings, "output_dir", output),
        ]
        if database:
            _add(command, "--database", database)
        consumed |= {
            "source", "design", "context", "database", "expression",
            "enrichment", "metagene", "output_dir",
        }
    elif workflow == "knockrbp_validation":
        output = dataset.output_root / "knockrbp" / f"{design}_{context}"
        command += [
            "--modification", dataset.modification, "--design", design,
            "--context", context, "--knockrbp-dir",
            str(
                _resource(
                    settings,
                    "knockrbp_dir",
                    resources,
                    "knockrbp_dir",
                    project_root / "data/knockrbp",
                )
            ),
            "--metadata",
            str(
                _resource(
                    settings,
                    "metadata",
                    resources,
                    "knockrbp_metadata",
                    project_root / "data/knockrbp/dataset_metadata.tsv",
                )
            ),
            "--output-dir", _value(settings, "output_dir", output),
        ]
        if design == "loose":
            short = SHORT_NAMES[dataset.modification]
            command += [
                "--encori-annotated",
                _value(
                    settings,
                    "encori_annotated",
                    dataset.output_root
                    / f"loose/encori/filtered_{short}_encori_annotated.tsv",
                ),
                "--postar3-annotated",
                _value(
                    settings,
                    "postar3_annotated",
                    dataset.output_root
                    / f"loose/postar3/filtered_{short}_postar3_annotated.tsv",
                ),
            ]
        else:
            command += [
                "--transcript-run-dir",
                _value(settings, "transcript_run_dir", paths["transcript"]),
            ]
        if settings.get("cross_database"):
            _add(command, "--cross-database", settings["cross_database"])
        elif settings.get("include_regulatory_network"):
            _add(
                command,
                "--cross-database",
                dataset.output_root
                / "cross_database"
                / f"{design}_{context}"
                / "cross_database_results.tsv",
            )
        consumed |= {
            "design", "context", "knockrbp_dir", "metadata", "output_dir",
            "encori_annotated", "postar3_annotated", "transcript_run_dir",
            "cross_database", "include_regulatory_network",
        }
    elif workflow == "string_analysis":
        cross = (
            dataset.output_root
            / "cross_database"
            / f"{design}_{context}"
            / "cross_database_results.tsv"
        )
        output = dataset.output_root / "string" / f"{design}_{context}"
        command += [
            "--modification", dataset.modification, "--design", design,
            "--context", context, "--cross-database",
            _value(settings, "cross_database", cross), "--cache-dir",
            _value(settings, "cache_dir", dataset.output_root / "string_cache"),
            "--output-dir", _value(settings, "output_dir", output),
        ]
        consumed |= {
            "design", "context", "cross_database", "cache_dir", "output_dir"
        }
    elif workflow == "compare_modifications":
        profiles = [dataset, *(settings.get("comparison_profiles", ()) or ())]
        resolved = {item.modification: item.resolved(project_root) for item in profiles}
        if len(resolved) < 2 or len(resolved) != len(profiles):
            raise ValueError("Comparison requires distinct modifications")
        command += [
            "--modifications", *resolved, "--databases", *databases,
            "--design", design, "--context", context, "--output-dir",
            _value(
                settings,
                "output_dir",
                dataset.output_root / "cross_modification",
            ),
        ]
        for modification, item in resolved.items():
            for database in databases:
                command += [
                    "--input", modification, database,
                    str(_enrichment(item, design, database, context)),
                ]
            if not settings.get("skip_phase1_context"):
                command += [
                    "--metagene",
                    modification,
                    str(
                        item.output_root
                        / "phase1"
                        / f"filtered_{SHORT_NAMES[modification]}_metagene.tsv"
                    ),
                ]
        consumed |= {
            "comparison_profiles", "databases", "design", "context", "output_dir"
        }
    elif workflow == "compare_datasets":
        return build_dataset_comparison(command, dataset, settings, project_root)
    elif workflow != "validate_handoff":
        raise KeyError(f"Unsupported workflow: {workflow}")
    return command, consumed
