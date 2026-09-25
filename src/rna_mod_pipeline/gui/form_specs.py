from __future__ import annotations


def _field(
    key: str,
    label: str,
    kind: str = "text",
    default=None,
    choices=(),
    required: bool = False,
    help_text: str = "",
    advanced: bool = False,
) -> dict[str, object]:
    return {
        "key": key,
        "label": label,
        "kind": kind,
        "default": default,
        "choices": choices,
        "required": required,
        "help_text": help_text,
        "advanced": advanced,
    }


COMMON = {
    "overwrite": _field(
        "overwrite",
        "Overwrite matching completed run",
        "bool",
        False,
        advanced=True,
    ),
}


WORKFLOW_FIELDS = {
    "phase1_filter": (
        _field("min_coverage", "Minimum coverage", "integer", 20),
        _field("min_fraction", "Minimum modification (%)", "float", 20),
        _field("skip_plots", "Skip plots", "bool", False),
        _field(
            "chunk_size",
            "Rows per input chunk",
            "integer",
            2_000_000,
            advanced=True,
            help_text="Memory/performance setting; it does not change the filter rule.",
        ),
        COMMON["overwrite"],
    ),
    "phase1_drach": (
        _field("skip_plots", "Skip plots", "bool", False),
        COMMON["overwrite"],
    ),
    "phase1_metagene": (
        _field("skip_plots", "Skip plots", "bool", False),
        COMMON["overwrite"],
    ),
    "loose_overlap": (
        _field(
            "database",
            "Binding database",
            choices=("ornament", "encori", "postar3"),
            required=True,
        ),
        _field("window", "Overlap window (nt)", "integer", 10),
        _field("case_min_coverage", "Case minimum coverage", "integer", 20),
        _field("case_min_fraction", "Case minimum modification (%)", "float", 20),
        _field("background_min_coverage", "Control minimum coverage", "integer", 20),
        _field("background_max_fraction", "Control maximum modification (%)", "float", 20),
        _field("minimum_clip_experiments", "Minimum ENCORI experiments", "integer", 1),
        _field("fdr_threshold", "FDR threshold", "float", 0.05, advanced=True),
        _field(
            "cell_types",
            "POSTAR3 cell types (comma-separated)",
            "list",
            "",
            advanced=True,
        ),
        _field(
            "methods",
            "POSTAR3 methods (comma-separated)",
            "list",
            "",
            advanced=True,
        ),
        _field(
            "site_annotation",
            "Site annotation",
            default="all",
            choices=("all", "significant", "none"),
        ),
        _field("skip_plots", "Skip plots", "bool", False),
        _field("top_n", "Top RBPs per plot", "integer", 30, advanced=True),
        _field(
            "allow_site_subset",
            "Allow an intentional case-site subset",
            "bool",
            False,
            advanced=True,
            help_text=(
                "Required for DRACH-only or another intentional subset; the run is "
                "labelled partial-input and excluded from primary inference."
            ),
        ),
        COMMON["overwrite"],
    ),
    "transcript_region_overlap": (
        _field(
            "databases",
            "Databases (comma-separated)",
            "list",
            "ornament, encori, postar3",
        ),
        _field("coverage_min", "Minimum coverage", "integer", 20),
        _field("case_min_fraction", "Case minimum modification (%)", "float", 20),
        _field("background_max_fraction", "Control maximum modification (%)", "float", 20),
        _field("window", "Overlap window (nt)", "integer", 10),
        _field("strand_mode", "Strand mode", default="ignore", choices=("ignore", "same")),
        _field("min_encori_support", "Minimum ENCORI support", "integer", 1, advanced=True),
        _field("fdr", "FDR threshold", "float", 0.05, advanced=True),
        _field("top_n", "Top RBPs per plot", "integer", 30, advanced=True),
        _field("dpi", "Plot DPI", "integer", 180, advanced=True),
        _field("seed", "Reproducibility seed", "integer", 20260728, advanced=True),
        _field("rbps", "Limit to RBPs (comma-separated)", "list", "", advanced=True),
        _field(
            "chromosomes",
            "Limit to chromosomes (comma-separated)",
            "list",
            "",
            advanced=True,
        ),
        _field(
            "postar_cell_types",
            "POSTAR3 cell types (comma-separated)",
            "list",
            "",
            advanced=True,
        ),
        _field(
            "postar_methods",
            "POSTAR3 methods (comma-separated)",
            "list",
            "",
            advanced=True,
        ),
        _field("no_plots", "Skip plots", "bool", False),
        COMMON["overwrite"],
    ),
    "transcript_validate": (
        _field("skip_input_hashes", "Skip input checksums", "bool", False),
    ),
    "plot_enrichment": (
        _field("fdr", "FDR threshold", "float", 0.05),
        _field("top_n", "Top RBPs per panel", "integer", 30),
        _field("dpi", "Plot DPI", "integer", 180),
        COMMON["overwrite"],
    ),
    "cross_database": (
        _field("design", "Design", default="loose", choices=("loose", "transcript-region")),
        _field("context", "Context", default="all", choices=("all", "5UTR", "CDS", "3UTR")),
        _field("databases", "Databases (comma-separated)", "list", "ornament, encori, postar3"),
        _field("fdr", "FDR threshold", "float", 0.05),
        _field("allow_missing", "Allow missing database", "bool", False),
        _field("skip_plots", "Skip plots", "bool", False),
        _field("plot_top_n", "Top RBPs in comparison plots", "integer", 20, advanced=True),
        _field("plot_dpi", "Comparison plot DPI", "integer", 180, advanced=True),
        COMMON["overwrite"],
    ),
    "integrate_expression": (
        _field("source", "Expression source", choices=("nanopore", "depmap"), required=True),
        _field("design", "Design", default="loose", choices=("loose", "transcript-region")),
        _field("context", "Context", default="all", choices=("all", "5UTR", "CDS", "3UTR")),
        _field(
            "database",
            "RBP database (required for loose design)",
            choices=("", "ornament", "encori", "postar3"),
        ),
        _field("model_id", "DepMap model ID", default="ACH-000768", advanced=True),
        _field("tpm_threshold", "DepMap TPM threshold", "float", 1.0, advanced=True),
        _field(
            "nanopore_value_column",
            "Nanopore normalised-value column",
            default="cpm_0h",
            advanced=True,
        ),
        _field(
            "nanopore_count_column",
            "Nanopore count column",
            default="count_0h",
            advanced=True,
        ),
        _field("gene_lengths", "Optional gene-length table", "file", "", advanced=True),
        COMMON["overwrite"],
    ),
    "knockrbp_validation": (
        _field("design", "Design", default="loose", choices=("loose", "transcript-region")),
        _field("context", "Context", default="all", choices=("all", "5UTR", "CDS", "3UTR")),
        _field("log2fc_cutoff", "|log2FC| cutoff", "float", 0.5),
        _field("padj_cutoff", "Adjusted-p cutoff", "float", 0.05),
        _field(
            "include_regulatory_network",
            "Include the cross-database RBP regulatory network",
            "bool",
            True,
            help_text=(
                "Uses the expected completed cross-database result for this "
                "design and context. Uncheck to run target/DEG overlap only."
            ),
        ),
        _field(
            "cell_lines",
            "Included cell lines (comma-separated)",
            "list",
            "MDA-MB-231, MDA-MB-231-LM2",
            advanced=True,
        ),
        _field("dataset_ids", "Limit to dataset IDs", "list", "", advanced=True),
        _field(
            "regulatory_plot_top_n",
            "Regulatory heatmap RBPs",
            "integer",
            40,
            advanced=True,
        ),
        _field(
            "cross_database",
            "Optional cross-database result table",
            "file",
            "",
            advanced=True,
            help_text="Adds quality-aware modification-association labels to the RBP network.",
        ),
        COMMON["overwrite"],
    ),
    "string_analysis": (
        _field("design", "Design", default="loose", choices=("loose", "transcript-region")),
        _field("context", "Context", default="all", choices=("all", "5UTR", "CDS", "3UTR")),
        _field("directions", "Directions (comma-separated)", "list", "enriched, depleted"),
        _field("top_n", "Top RBPs", "integer", 10),
        _field("fdr", "FDR threshold", "float", 0.05),
        _field("offline", "Offline/cache-only", "bool", False),
        _field("retries", "STRING request attempts", "integer", 3, advanced=True),
        _field("api_delay", "Delay between STRING requests (seconds)", "float", 1.0, advanced=True),
        COMMON["overwrite"],
    ),
    "compare_modifications": (
        _field("databases", "Databases (comma-separated)", "list", "ornament, encori, postar3"),
        _field("design", "Design", default="loose", choices=("loose", "transcript-region")),
        _field("context", "Context", default="all", choices=("all", "5UTR", "CDS", "3UTR")),
        _field("fdr", "FDR threshold", "float", 0.05),
        _field(
            "skip_phase1_context",
            "Skip transcript/site context summaries",
            "bool",
            False,
            advanced=True,
        ),
        _field("allow_missing", "Allow missing database inputs", "bool", False, advanced=True),
        COMMON["overwrite"],
    ),
    "compare_datasets": (
        _field("fdr", "FDR threshold", "float", 0.05),
        _field("dpi", "Plot DPI", "integer", 180),
        _field("skip_plots", "Skip plots", "bool", False),
        COMMON["overwrite"],
    ),
    "validate_handoff": (
        _field("skip_legacy_results", "Skip legacy-result checks", "bool", False),
    ),
}


def fields_for(workflow: str) -> tuple[dict[str, object], ...]:
    return tuple(WORKFLOW_FIELDS.get(workflow, ()))
