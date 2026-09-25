from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from ..config import KNOWN_RELATED, modification as modification_config
from ..evidence_quality import add_significance_calls
from ..intervals import IntervalIndex, group_positions
from ..io import iter_bedmethyl, load_sites
from ..stats import bh_adjust, fisher_test
from .sources import BindingTrack


@dataclass(frozen=True)
class LooseOverlapResult:
    enrichment: pd.DataFrame
    annotated_sites: pd.DataFrame
    case_count: int
    background_count: int
    qualifying_raw_case_count: int
    omitted_qualifying_case_count: int
    case_subset_allowed: bool


def _overlap_mask(
    positions: np.ndarray,
    chromosome_groups: dict[str, np.ndarray],
    intervals: IntervalIndex,
    window: int,
) -> np.ndarray:
    mask = np.zeros(len(positions), dtype=bool)
    for chrom, indexes in chromosome_groups.items():
        indexed = intervals.get(chrom)
        if indexed is None or len(indexed.starts) == 0:
            continue
        selected = positions[indexes]
        right = np.searchsorted(
            indexed.starts,
            selected + window + 1,
            side="left",
        )
        candidates = right > 0
        if candidates.any():
            furthest_end = indexed.maximum_end[right[candidates] - 1]
            mask[indexes[candidates]] = (
                furthest_end > selected[candidates] - window
            )
    return mask


def _background_positions(
    bedmethyl_path: str | Path,
    sites: pd.DataFrame,
    minimum_coverage: int,
    maximum_fraction: float,
    case_minimum_coverage: int,
    case_minimum_fraction: float,
    expected_mod_codes: frozenset[str],
    allow_site_subset: bool,
) -> tuple[dict[str, np.ndarray], int, int]:
    case_keys = set(zip(sites["chrom"].astype(str), sites["start"].astype(int)))
    grouped: dict[str, list[int]] = defaultdict(list)
    qualifying_case_keys: set[tuple[str, int]] = set()
    observed_case_keys: set[tuple[str, int]] = set()
    expected_codes = {code.lower() for code in expected_mod_codes}
    for chrom, start, _strand, mod_code, coverage, fraction in iter_bedmethyl(
        bedmethyl_path
    ):
        normalized_code = str(mod_code).strip().lower()
        if normalized_code not in expected_codes:
            raise ValueError(
                "bedMethyl modification code does not match the selected "
                f"analysis: observed {normalized_code!r}, expected one of "
                f"{sorted(expected_codes)}"
            )
        key = (chrom, start)
        if key in case_keys:
            observed_case_keys.add(key)
        if (
            coverage >= case_minimum_coverage
            and fraction >= case_minimum_fraction
        ):
            qualifying_case_keys.add(key)
        if (
            coverage >= minimum_coverage
            and fraction < maximum_fraction
            and key not in case_keys
        ):
            grouped[chrom].append(start)
    missing_cases = case_keys.difference(observed_case_keys)
    if missing_cases:
        preview = sorted(missing_cases)[:5]
        raise ValueError(
            f"{len(missing_cases)} supplied case sites are absent from the raw "
            f"bedMethyl file; examples: {preview}"
        )
    omitted = qualifying_case_keys.difference(case_keys)
    if omitted and not allow_site_subset:
        raise ValueError(
            f"The site table omits {len(omitted)} raw calls meeting the case "
            "thresholds. Use the complete Phase 1 table, or pass "
            "--allow-site-subset for an intentional subset such as DRACH-only."
        )
    return {
        chrom: np.asarray(positions, dtype=np.int64)
        for chrom, positions in grouped.items()
    }, len(qualifying_case_keys), len(omitted)


def _overlap_count(
    grouped_positions: dict[str, np.ndarray],
    intervals: IntervalIndex,
    window: int,
) -> int:
    count = 0
    for chrom, positions in grouped_positions.items():
        indexed = intervals.get(chrom)
        if indexed is None or len(indexed.starts) == 0:
            continue
        right = np.searchsorted(
            indexed.starts,
            positions + window + 1,
            side="left",
        )
        candidates = right > 0
        if candidates.any():
            furthest_end = indexed.maximum_end[right[candidates] - 1]
            count += int(
                np.sum(furthest_end > positions[candidates] - window)
            )
    return count


def _source_columns(database: str, track: BindingTrack) -> dict[str, object]:
    if database == "ornament":
        return {
            "n_source_tracks": len(track.source_names),
            "source_tracks": ",".join(track.source_names),
        }
    if database == "encori":
        return {"data_source": "ENCORI_CLIP-seq"}
    if database == "postar3":
        return {
            "data_source": "POSTAR3_CLIP-seq",
            "n_experiments": track.metadata.get("n_experiments", 0),
            "n_cell_types": track.metadata.get("n_cell_types", 0),
            "methods": track.metadata.get("methods", ""),
        }
    raise ValueError(f"Unsupported binding database: {database}")


def _estimate_reason(
    odds_ratio: float,
    case_bound: int,
    case_count: int,
    control_bound: int,
    background_count: int,
    resource_status: str,
) -> str:
    if resource_status == "resource_empty":
        return "resource_empty"
    if np.isnan(odds_ratio):
        if case_bound == 0 and control_bound == 0:
            return "no_case_or_control_overlap"
        if case_bound == case_count and control_bound == background_count:
            return "all_case_and_control_sites_bound"
        return "degenerate_contingency_table"
    if odds_ratio == 0:
        return "zero_case_overlap" if case_bound == 0 else "all_control_sites_bound"
    if np.isposinf(odds_ratio):
        return "zero_control_overlap" if control_bound == 0 else "all_case_sites_bound"
    return ""


def run_loose_overlap(
    sites_path: str | Path,
    bedmethyl_path: str | Path,
    tracks: Iterable[BindingTrack],
    modification: str,
    database: str,
    window: int = 10,
    background_min_coverage: int = 20,
    background_max_fraction: float = 20.0,
    fdr_threshold: float = 0.05,
    site_annotation: str = "all",
    case_min_coverage: int = 20,
    case_min_fraction: float = 20.0,
    allow_site_subset: bool = False,
) -> LooseOverlapResult:
    """Compare all modification sites with the loose same-base background."""
    if modification == "m6a_rep2":
        machinery_key = "m6a"
    else:
        machinery_key = modification
    if machinery_key not in KNOWN_RELATED:
        raise ValueError(f"Unsupported modification: {modification}")
    if database not in {"ornament", "encori", "postar3"}:
        raise ValueError(f"Unsupported binding database: {database}")
    if window < 0:
        raise ValueError("window must be non-negative")
    if background_min_coverage < 0:
        raise ValueError("background_min_coverage must be non-negative")
    if not 0 <= background_max_fraction <= 100:
        raise ValueError("background_max_fraction must lie in [0, 100]")
    if not 0 < fdr_threshold <= 1:
        raise ValueError("fdr_threshold must lie in (0, 1]")
    if case_min_coverage < 0:
        raise ValueError("case_min_coverage must be non-negative")
    if not 0 <= case_min_fraction <= 100:
        raise ValueError("case_min_fraction must lie in [0, 100]")
    if background_max_fraction > case_min_fraction:
        raise ValueError(
            "background_max_fraction cannot exceed case_min_fraction; "
            "case and background definitions must not overlap"
        )
    if site_annotation not in {"all", "significant", "none"}:
        raise ValueError("site_annotation must be all, significant, or none")

    sites = load_sites(sites_path)
    required_site_columns = {
        "chrom",
        "start",
        "strand",
        "mod_code",
        "Nvalid_cov",
        "fraction_modified",
    }
    missing_columns = required_site_columns.difference(sites.columns)
    if missing_columns:
        raise ValueError(
            f"Site table is missing provenance columns: {sorted(missing_columns)}"
        )
    if sites.duplicated(["chrom", "start"]).any():
        raise ValueError("Site table contains duplicate genomic case positions")
    expected_codes = modification_config(modification).expected_mod_codes
    site_codes = {
        str(value).strip().lower()
        for value in sites["mod_code"].dropna().unique()
    }
    if not site_codes or not site_codes.issubset(
        {code.lower() for code in expected_codes}
    ):
        raise ValueError(
            "Site-table modification codes do not match the selected "
            f"analysis: observed {sorted(site_codes)}"
        )
    case_coverage = pd.to_numeric(sites["Nvalid_cov"], errors="coerce")
    case_fraction = pd.to_numeric(sites["fraction_modified"], errors="coerce")
    if case_coverage.isna().any() or case_fraction.isna().any():
        raise ValueError("Case coverage and modification fraction must be numeric")
    if case_coverage.lt(case_min_coverage).any():
        raise ValueError("Site table contains cases below case_min_coverage")
    if case_fraction.lt(case_min_fraction).any():
        raise ValueError("Site table contains cases below case_min_fraction")

    background, qualifying_raw_count, omitted_case_count = _background_positions(
        bedmethyl_path,
        sites,
        minimum_coverage=background_min_coverage,
        maximum_fraction=background_max_fraction,
        case_minimum_coverage=case_min_coverage,
        case_minimum_fraction=case_min_fraction,
        expected_mod_codes=expected_codes,
        allow_site_subset=allow_site_subset,
    )
    case_chromosomes = sites["chrom"].astype(str).to_numpy()
    case_positions = sites["start"].astype(np.int64).to_numpy()
    case_groups = group_positions(case_chromosomes)
    del case_chromosomes
    case_count = len(sites)
    background_count = sum(len(positions) for positions in background.values())
    if case_count == 0 or background_count == 0:
        raise ValueError("Cases and background must both contain at least one site")

    rows: list[dict[str, object]] = []
    case_masks: dict[str, np.ndarray] = {}
    known = KNOWN_RELATED[machinery_key]
    for track in tracks:
        case_mask = _overlap_mask(
            case_positions,
            case_groups,
            track.intervals,
            window,
        )
        control_bound = _overlap_count(
            background,
            track.intervals,
            window,
        )
        case_bound = int(case_mask.sum())
        tested = fisher_test(
            case_bound,
            case_count,
            control_bound,
            background_count,
        )
        odds_ratio, p_value = tested.odds_ratio, tested.p_value

        case_masks[track.rbp_name] = case_mask
        row = {
            "RBP": track.rbp_name,
            "modification": modification,
            "database": database,
            "analysis_design": "loose_all_context",
            "partial_input": bool(allow_site_subset and omitted_case_count > 0),
            "input_subset_status": (
                "intentional_case_subset"
                if allow_site_subset and omitted_case_count > 0
                else "complete_qualifying_case_set"
            ),
            "total_binding_sites": track.interval_count,
            "mod_overlaps": case_bound,
            "mod_overlap_pct": round(100 * case_bound / case_count, 2),
            "bg_overlaps": control_bound,
            "bg_overlap_pct": round(100 * control_bound / background_count, 2),
            "odds_ratio": odds_ratio,
            "p_value": p_value,
            "is_known_related": track.canonical_rbp in known,
            "resource_status": track.metadata.get("resource_status", "loaded"),
            "estimate_status": (
                "non_estimable"
                if np.isnan(odds_ratio)
                else (
                    "boundary"
                    if odds_ratio == 0 or np.isinf(odds_ratio)
                    else "estimable"
                )
            ),
            "estimate_reason": _estimate_reason(
                odds_ratio, case_bound, case_count, control_bound,
                background_count, track.metadata.get("resource_status", "loaded"),
            ),
        }
        row.update(_source_columns(database, track))
        rows.append(row)
        del track

    enrichment = pd.DataFrame(rows)
    if enrichment.empty:
        raise ValueError("No RBP binding tracks were supplied")
    enrichment["fdr"] = bh_adjust(enrichment["p_value"])
    enrichment = add_significance_calls(
        enrichment,
        odds_column="odds_ratio",
        fdr_column="fdr",
        fdr_threshold=fdr_threshold,
    )
    enrichment["direction"] = enrichment["statistical_direction"]
    enrichment = enrichment.sort_values(
        ["p_value", "RBP"],
        kind="mergesort",
    ).reset_index(drop=True)

    annotated = sites.copy()
    annotated["analysis_modification"] = modification
    annotated["analysis_database"] = database
    annotated["analysis_design"] = "loose_all_context"
    significant = set(
        enrichment.loc[enrichment["significant"], "RBP"]
    )
    if site_annotation == "all":
        prefix = {
            "ornament": "ornament",
            "encori": "encori",
            "postar3": "postar3",
        }[database]
        binary = pd.DataFrame(
            {
                f"{prefix}_{rbp_name}": case_masks[rbp_name].astype(np.int8)
                for rbp_name in enrichment["RBP"]
            },
            index=annotated.index,
        )
        annotated = pd.concat([annotated, binary], axis=1)

    if site_annotation != "none":
        lists: list[str] = []
        counts: list[int] = []
        ordered_names = sorted(significant)
        for site_index in range(case_count):
            bound = [
                name for name in ordered_names if case_masks[name][site_index]
            ]
            lists.append(",".join(bound))
            counts.append(len(bound))
        annotated["overlapping_RBPs"] = lists
        annotated["n_overlapping_RBPs"] = counts

    return LooseOverlapResult(
        enrichment=enrichment,
        annotated_sites=annotated,
        case_count=case_count,
        background_count=background_count,
        qualifying_raw_case_count=qualifying_raw_count,
        omitted_qualifying_case_count=omitted_case_count,
        case_subset_allowed=allow_site_subset,
    )
