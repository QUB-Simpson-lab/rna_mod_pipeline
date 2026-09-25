from __future__ import annotations

from collections.abc import Mapping

from ..evidence_quality import evidence_quality_manifest
from .loose import LooseOverlapResult


def loose_manifest_parameters(
    result: LooseOverlapResult,
    *,
    modification: str,
    database: str,
    catalog_version: str,
    options: Mapping[str, object],
) -> dict[str, object]:
    """Build the scientific and filtering metadata for one loose run."""
    return {
        "analysis": "loose_all_context",
        "modification": modification,
        "database": database,
        "catalog_version": catalog_version,
        "case_sites": result.case_count,
        "qualifying_raw_case_sites": result.qualifying_raw_case_count,
        "omitted_qualifying_case_sites": result.omitted_qualifying_case_count,
        "allow_site_subset": options["allow_site_subset"],
        "background_sites": result.background_count,
        "rbps_tested": len(result.enrichment),
        "window_bp": options["window"],
        "strand_policy": "ignored",
        "background_min_coverage": options["background_min_coverage"],
        "background_max_fraction_percent": options["background_max_fraction"],
        "case_min_coverage": options["case_min_coverage"],
        "case_min_fraction_percent": options["case_min_fraction"],
        "fdr_threshold": options["fdr_threshold"],
        "minimum_clip_experiments": options["minimum_clip_experiments"],
        "cell_types": options["cell_types"],
        "methods": options["methods"],
        "site_annotation": options["site_annotation"],
        "allow_empty_resources": options["allow_empty_resources"],
        "ornament_coordinate_policy": (
            "historical loose interpretation; start==end rows retained"
            if database == "ornament"
            else "not_applicable"
        ),
        "empty_resources": int(
            result.enrichment["resource_status"].eq("resource_empty").sum()
        ),
        "site_annotation_significance_scope": (
            "quality-qualified significant RBPs; raw statistical calls "
            "remain in the enrichment table"
        ),
        "evidence_quality": evidence_quality_manifest(
            result.enrichment,
            fdr_threshold=float(options["fdr_threshold"]),
        ),
    }
