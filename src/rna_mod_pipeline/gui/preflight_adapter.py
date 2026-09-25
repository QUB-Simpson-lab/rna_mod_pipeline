from __future__ import annotations

from typing import Any, Mapping

from .adapters import BuildPlan
from .preflight import validate_request


def preflight_errors(
    workflow_key: str,
    dataset: Any,
    values: Mapping[str, Any],
    plan: BuildPlan,
) -> tuple[list[str], list[str]]:
    result = validate_request(
        workflow_key,
        dataset,
        dict(values),
        list(plan.argv),
    )
    return (
        [str(item) for item in result.errors],
        [str(item) for item in result.warnings],
    )
