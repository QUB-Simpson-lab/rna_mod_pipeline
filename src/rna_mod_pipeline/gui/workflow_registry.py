from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping


FORWARDED_OPTIONS = {
    "phase1_filter": ("min_coverage", "min_fraction", "chunk_size", "skip_plots", "overwrite"),
    "phase1_drach": ("skip_plots", "overwrite"),
    "phase1_metagene": ("skip_plots", "overwrite"),
    "loose_overlap": ("window", "background_min_coverage", "background_max_fraction",
        "case_min_coverage", "case_min_fraction", "fdr_threshold", "minimum_clip_experiments",
        "cell_types", "methods", "site_annotation", "top_n", "allow_site_subset",
        "allow_empty_resources", "skip_plots", "overwrite"),
    "transcript_region_overlap": ("coverage_min", "case_min_fraction",
        "background_max_fraction", "window", "strand_mode", "min_encori_support",
        "fdr", "top_n", "dpi", "seed", "rbps", "chromosomes", "postar_cell_types",
        "postar_methods", "maximum_rbps_per_database", "max_cases", "max_bed_rows",
        "max_postar_rows", "no_plots", "overwrite"),
    "transcript_validate": ("skip_input_hashes",),
    "plot_enrichment": ("fdr", "top_n", "dpi", "overwrite"),
    "cross_database": ("fdr", "allow_missing", "skip_plots", "plot_top_n", "plot_dpi", "overwrite"),
    "integrate_expression": ("model_id", "tpm_threshold", "nanopore_value_column",
        "nanopore_count_column", "gene_lengths", "overwrite"),
    "knockrbp_validation": ("cell_lines", "dataset_ids", "encori_enrichment",
        "postar3_enrichment", "log2fc_cutoff", "padj_cutoff",
        "regulatory_plot_top_n", "overwrite"),
    "string_analysis": ("directions", "top_n", "fdr", "offline", "retries", "api_delay", "overwrite"),
    "compare_modifications": ("skip_phase1_context", "allow_missing", "fdr", "overwrite"),
    "compare_datasets": ("fdr", "dpi", "skip_plots", "overwrite"),
    "validate_handoff": ("skip_legacy_results",),
}


@dataclass(frozen=True)
class WorkflowSpec:
    key: str
    label: str
    script: str
    category: str
    requires: tuple[str, ...] = ()
    modifications: tuple[str, ...] = ("m6a", "m5c", "pseu")
    dataset_local: bool = True
    description: str = ""


_WORKFLOWS = (
    WorkflowSpec(
        "phase1_filter",
        "Phase 1 · Filter bedMethyl",
        "phase1_filter.py",
        "Phase 1",
        description="Apply coverage and modification-fraction thresholds.",
    ),
    WorkflowSpec(
        "phase1_drach",
        "Phase 1 · DRACH annotation",
        "phase1_drach.py",
        "Phase 1",
        requires=("phase1_filter",),
        modifications=("m6a",),
        description="Annotate m6A calls without removing non-DRACH sites.",
    ),
    WorkflowSpec(
        "phase1_metagene",
        "Phase 1 · Transcript annotation",
        "phase1_metagene.py",
        "Phase 1",
        requires=("phase1_filter",),
        description="Assign calls to complete transcripts and regions.",
    ),
    WorkflowSpec(
        "loose_overlap",
        "RBP overlap · Loose",
        "loose_overlap.py",
        "RBP overlap",
        requires=("phase1_metagene",),
        description="Run one binding database without transcript matching.",
    ),
    WorkflowSpec(
        "transcript_region_overlap",
        "RBP overlap · Transcript-region",
        "transcript_region_overlap.py",
        "RBP overlap",
        requires=("phase1_metagene",),
        description="Stratify cases and controls by transcript and region.",
    ),
    WorkflowSpec(
        "transcript_validate",
        "Validate transcript-region run",
        "transcript_region_overlap.py",
        "Validation",
        requires=("transcript_region_overlap",),
        description="Validate tables, matrices, checksums, and manifests.",
    ),
    WorkflowSpec(
        "plot_enrichment",
        "Regenerate transcript plots",
        "plot_enrichment.py",
        "Plotting",
        requires=("transcript_region_overlap",),
    ),
    WorkflowSpec(
        "cross_database",
        "Cross-database comparison",
        "cross_database.py",
        "Integration",
        requires=("loose_overlap",),
    ),
    WorkflowSpec(
        "integrate_expression",
        "Expression integration",
        "integrate_expression.py",
        "Integration",
        requires=("phase1_metagene",),
    ),
    WorkflowSpec(
        "knockrbp_validation",
        "KnockRBP validation",
        "knockrbp_validation.py",
        "Validation",
        requires=("loose_overlap",),
    ),
    WorkflowSpec(
        "string_analysis",
        "STRING analysis",
        "string_analysis.py",
        "Functional analysis",
        requires=("cross_database",),
    ),
    WorkflowSpec(
        "compare_modifications",
        "Cross-modification comparison",
        "compare_modifications.py",
        "Integration",
        requires=("cross_database",),
        dataset_local=False,
    ),
    WorkflowSpec(
        "compare_datasets",
        "Compare two datasets",
        "compare_datasets.py",
        "Validation",
        requires=("phase1_metagene",),
        description="Same-modification robustness comparison without pooling.",
    ),
    WorkflowSpec(
        "validate_handoff",
        "Validate historical project",
        "validate_handoff.py",
        "Validation",
        dataset_local=False,
    ),
)


class WorkflowRegistry:
    def __init__(self, workflows: Iterable[WorkflowSpec] = _WORKFLOWS) -> None:
        items = tuple(workflows)
        self._workflows = {workflow.key: workflow for workflow in items}
        if len(self._workflows) != len(items):
            raise ValueError("Workflow keys must be unique")

    def list(self) -> tuple[WorkflowSpec, ...]:
        return tuple(self._workflows.values())

    def get(self, key: str) -> WorkflowSpec:
        try:
            return self._workflows[key]
        except KeyError as exc:
            raise KeyError(f"Unknown workflow: {key}") from exc

    def keys(self) -> tuple[str, ...]:
        return tuple(self._workflows)


_DEFAULT_REGISTRY = WorkflowRegistry()


def list_workflows() -> tuple[WorkflowSpec, ...]:
    return _DEFAULT_REGISTRY.list()


def prerequisites_for(
    workflow: str,
    modification: str | None = None,
    options: Mapping[str, object] | None = None,
) -> tuple[str, ...]:
    """Return the prerequisites implied by the selected workflow options."""
    settings = options or {}
    if workflow == "phase1_metagene":
        return ("phase1_drach",) if modification == "m6a" else ("phase1_filter",)
    if workflow in {
        "cross_database",
        "integrate_expression",
        "knockrbp_validation",
    }:
        if workflow == "knockrbp_validation" and settings.get(
            "include_regulatory_network"
        ):
            return ("cross_database",)
        design = str(settings.get("design", "loose"))
        return (
            ("transcript_region_overlap",)
            if design == "transcript-region"
            else ("loose_overlap",)
        )
    return _DEFAULT_REGISTRY.get(workflow).requires
