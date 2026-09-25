from __future__ import annotations

import inspect
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping


@dataclass(frozen=True)
class FieldSpec:
    key: str
    label: str
    kind: str = "text"
    required: bool = False
    choices: tuple[str, ...] = ()
    default: Any = None
    help_text: str = ""
    advanced: bool = False


@dataclass(frozen=True)
class WorkflowSpec:
    key: str
    label: str
    description: str = ""
    fields: tuple[FieldSpec, ...] = ()


@dataclass(frozen=True)
class BuildPlan:
    argv: tuple[str, ...]
    output_dir: Path | None = None
    manifest_path: Path | None = None
    warnings: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)


def _get(value: Any, *names: str, default: Any = None) -> Any:
    for name in names:
        if isinstance(value, Mapping) and name in value:
            return value[name]
        if hasattr(value, name):
            return getattr(value, name)
    return default


def load_dataset_registry(project_root: str | Path) -> Any:
    from .datasets import DatasetRegistry

    root = Path(project_root).expanduser().resolve()
    loader = getattr(DatasetRegistry, "load", None)
    if callable(loader):
        return loader(root)
    return DatasetRegistry(root)


def dataset_profiles(registry: Any) -> list[Any]:
    for name in ("list", "all", "profiles", "datasets"):
        value = getattr(registry, name, None)
        if callable(value):
            return list(value())
        if value is not None:
            return list(value.values()) if isinstance(value, Mapping) else list(value)
    try:
        return list(registry)
    except TypeError as exc:
        raise TypeError("DatasetRegistry does not expose an iterable dataset list") from exc


def dataset_label(profile: Any) -> str:
    name = str(_get(profile, "name", "dataset_name", "id", default="Unnamed dataset"))
    modification = _get(profile, "modification", "mod_type", default="")
    return f"{name} ({modification})" if modification else name


def dataset_details(profile: Any) -> dict[str, str]:
    keys = (
        ("name", ("name", "dataset_name", "id")),
        ("modification", ("modification", "mod_type")),
        ("bedMethyl", ("bedmethyl", "bedmethyl_path", "input_path")),
        ("FASTA", ("fasta", "fasta_path")),
        ("GTF", ("gtf", "gtf_path")),
        ("sample metadata", ("sample_metadata", "metadata")),
        ("output root", ("output_root", "output_dir")),
    )
    return {
        label: str(value)
        for label, names in keys
        if (value := _get(profile, *names)) not in (None, "")
    }


def add_dataset(registry: Any, values: Mapping[str, Any]) -> Any:
    from .datasets import DatasetProfile

    aliases = {
        "bedmethyl_path": values.get("bedmethyl"),
        "fasta_path": values.get("fasta"),
        "gtf_path": values.get("gtf"),
        "metadata": values.get("sample_metadata"),
        **values,
    }
    signature = inspect.signature(DatasetProfile)
    accepts_kwargs = any(
        item.kind == inspect.Parameter.VAR_KEYWORD
        for item in signature.parameters.values()
    )
    kwargs = (
        dict(values)
        if accepts_kwargs
        else {
            name: aliases[name]
            for name in signature.parameters
            if name in aliases and aliases[name] not in (None, "")
        }
    )
    profile = DatasetProfile(**kwargs)
    method = getattr(registry, "add", None) or getattr(registry, "add_dataset", None)
    if not callable(method):
        raise TypeError("DatasetRegistry must provide add(profile) or add_dataset(profile)")
    added = method(profile)
    saver = getattr(registry, "save", None)
    if callable(saver):
        saver()
    return added if added is not None else profile


def replace_dataset(
    registry: Any,
    original_name: str,
    values: Mapping[str, Any],
) -> Any:
    from .datasets import DatasetProfile

    profile = DatasetProfile(**dict(values))
    method = getattr(registry, "replace", None)
    if not callable(method):
        raise TypeError("DatasetRegistry must provide replace(name, profile)")
    updated = method(original_name, profile)
    saver = getattr(registry, "save", None)
    if callable(saver):
        saver()
    return updated


def remove_dataset(registry: Any, name: str) -> Any:
    method = getattr(registry, "remove", None)
    if not callable(method):
        raise TypeError("DatasetRegistry must provide remove(name)")
    removed = method(name)
    saver = getattr(registry, "save", None)
    if callable(saver):
        saver()
    return removed


def _normalise_field(value: Any) -> FieldSpec:
    choices = _get(value, "choices", "options", default=()) or ()
    return FieldSpec(
        key=str(_get(value, "key", "name")),
        label=str(_get(value, "label", "title", "key", "name")),
        kind=str(_get(value, "kind", "type", default="text")),
        required=bool(_get(value, "required", default=False)),
        choices=tuple(str(item) for item in choices),
        default=_get(value, "default"),
        help_text=str(_get(value, "help_text", "help", "description", default="")),
        advanced=bool(_get(value, "advanced", default=False)),
    )


def workflow_specs() -> list[WorkflowSpec]:
    from . import workflow_registry as module
    from .form_specs import fields_for

    source = None
    for name in ("list_workflows", "get_workflows"):
        candidate = getattr(module, name, None)
        if callable(candidate):
            source = candidate()
            break
    if source is None:
        source = getattr(module, "WORKFLOWS", None)
    if source is None:
        registry = getattr(module, "WorkflowRegistry", None)
        if registry is not None:
            instance = registry()
            source = getattr(instance, "all", lambda: instance)()
    if source is None:
        raise RuntimeError("workflow_registry exposes no workflow definitions")
    values: Iterable[Any] = source.values() if isinstance(source, Mapping) else source
    result = []
    for value in values:
        fields = _get(value, "fields", "parameters", default=()) or fields_for(
            str(_get(value, "key", "name"))
        )
        result.append(
            WorkflowSpec(
                key=str(_get(value, "key", "name")),
                label=str(_get(value, "label", "title", "key", "name")),
                description=str(_get(value, "description", "help", default="")),
                fields=tuple(_normalise_field(item) for item in fields),
            )
        )
    if not any(item.key == "compare_datasets" for item in result):
        result.append(
            WorkflowSpec(
                key="compare_datasets",
                label="Compare two datasets",
                description=(
                    "Compare two independently processed datasets as a robustness "
                    "analysis. This does not pool raw replicates."
                ),
                fields=tuple(
                    _normalise_field(item) for item in fields_for("compare_datasets")
                ),
            )
        )
    return result


def build_plan(
    workflow_key: str,
    dataset: Any,
    values: Mapping[str, Any],
    project_root: str | Path | None = None,
    resource_profile: Any | None = None,
) -> BuildPlan:
    from .command_builder import build_command

    kwargs: dict[str, Any] = {"resource_profile": resource_profile}
    if project_root is not None:
        kwargs["project_root"] = project_root
    result = build_command(workflow_key, dataset, dict(values), **kwargs)
    argv = _get(result, "argv", "command", default=result)
    if isinstance(argv, str) or not isinstance(argv, Iterable):
        raise TypeError("build_command must return an argv sequence, not a shell string")
    output = _get(result, "output_dir", "output", default=values.get("output_dir"))
    arguments = [str(item) for item in argv]
    manifest = _get(result, "manifest_path", "manifest")
    if not output:
        for flag in ("--output-dir", "--manifest"):
            if flag in arguments and arguments.index(flag) + 1 < len(arguments):
                candidate = Path(arguments[arguments.index(flag) + 1])
                output = candidate if flag == "--output-dir" else candidate.parent
                break
    if not manifest:
        if "--manifest" in arguments and arguments.index("--manifest") + 1 < len(arguments):
            manifest = arguments[arguments.index("--manifest") + 1]
        elif output:
            manifest = Path(output) / "run_manifest.json"
    warnings = _get(result, "warnings", default=()) or ()
    metadata = _get(result, "metadata", default={}) or {}
    return BuildPlan(
        argv=tuple(str(item) for item in argv),
        output_dir=Path(output).expanduser().resolve() if output else None,
        manifest_path=Path(manifest).expanduser().resolve() if manifest else None,
        warnings=tuple(str(item) for item in warnings),
        metadata=metadata,
    )
