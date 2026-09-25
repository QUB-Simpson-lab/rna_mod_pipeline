from __future__ import annotations

from pathlib import Path
from typing import Mapping

from .datasets import DatasetProfile


def _profiles(
    dataset: DatasetProfile,
    settings: Mapping[str, object],
    project_root: Path,
) -> tuple[DatasetProfile, DatasetProfile]:
    second = settings.get("second_dataset")
    if not isinstance(second, DatasetProfile):
        raise ValueError("compare_datasets requires second_dataset")
    first = dataset.resolved(project_root)
    second = second.resolved(project_root)
    if first.name.casefold() == second.name.casefold():
        raise ValueError("Dataset comparison requires two distinct datasets")
    if first.modification != second.modification:
        raise ValueError("Dataset comparison requires the same modification")
    assert first.output_root and second.output_root
    if (
        first.output_root == second.output_root
        or first.output_root.is_relative_to(second.output_root)
        or second.output_root.is_relative_to(first.output_root)
    ):
        raise ValueError("Compared datasets require separate output roots")
    return first, second


def build_dataset_comparison(
    command: list[str],
    dataset: DatasetProfile,
    settings: Mapping[str, object],
    project_root: Path,
) -> tuple[list[str], set[str]]:
    first, second = _profiles(dataset, settings, project_root)
    assert first.output_root and second.output_root
    left = first.output_root / "dataset_comparison_profile.json"
    right = second.output_root / "dataset_comparison_profile.json"
    default = first.output_root.parent / (
        f"comparison_{first.output_root.name}_vs_{second.output_root.name}"
    )
    output = Path(str(settings.get("output_dir") or default)).expanduser()
    output = (output if output.is_absolute() else project_root / output).resolve()
    if any(
        output == root or output.is_relative_to(root)
        for root in (first.output_root, second.output_root)
    ):
        raise ValueError("Comparison output must be outside both dataset outputs")
    command += [
        "--left-profile", str(left),
        "--right-profile", str(right),
        "--output-dir", str(output),
    ]
    return command, {"second_dataset", "output_dir"}
