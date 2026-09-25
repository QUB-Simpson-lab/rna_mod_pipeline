from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from .config import modification as modification_config


ENRICHMENT_REQUIRED = {
    "RBP",
    "mod_overlaps",
    "bg_overlaps",
    "odds_ratio",
    "p_value",
    "fdr",
}

MODIFICATION_LABELS = {
    "m6a": "m6a",
    "m6a_rep2": "m6a_rep2",
    "m5c": "m5c",
    "pseu": "pseu",
    "pseudouridine": "pseu",
    "ψ": "pseu",
}

DESIGN_LABELS = {
    "loose": "loose",
    "loose_all_context": "loose",
    "transcript": "transcript-region",
    "transcript-region": "transcript-region",
    "transcript_region": "transcript-region",
}


def require_columns(
    frame: pd.DataFrame,
    columns: Iterable[str],
    table_name: str,
) -> None:
    missing = set(columns).difference(frame.columns)
    if missing:
        raise ValueError(f"{table_name} is missing columns: {sorted(missing)}")


def validate_unique(frame: pd.DataFrame, columns: Iterable[str], table_name: str) -> None:
    keys = list(columns)
    duplicate = frame.duplicated(keys, keep=False)
    if duplicate.any():
        values = frame.loc[duplicate, keys].head(5).to_dict("records")
        raise ValueError(f"{table_name} has duplicate keys {keys}: {values}")


def validate_enrichment(frame: pd.DataFrame, table_name: str = "enrichment table") -> None:
    require_columns(frame, ENRICHMENT_REQUIRED, table_name)
    validate_unique(frame, ["RBP"], table_name)
    for column in ("mod_overlaps", "bg_overlaps"):
        numeric = pd.to_numeric(frame[column], errors="coerce")
        if numeric.isna().any() or (numeric < 0).any():
            raise ValueError(f"{table_name} has invalid {column} values")
    for column in ("p_value", "fdr"):
        numeric = pd.to_numeric(frame[column], errors="coerce")
        finite = numeric.dropna()
        if ((finite < 0) | (finite > 1)).any():
            raise ValueError(f"{table_name} has invalid {column} values")


def _normalized_values(
    frame: pd.DataFrame,
    column: str,
    aliases: dict[str, str],
    table_name: str,
) -> set[str]:
    values = {
        aliases.get(str(value).strip().lower(), str(value).strip().lower())
        for value in frame[column].dropna()
    }
    if not values:
        raise ValueError(f"{table_name} has an empty {column} label")
    return values


def validate_result_identity(
    frame: pd.DataFrame,
    *,
    expected_modification: str | None = None,
    expected_database: str | None = None,
    expected_design: str | None = None,
    context: str = "all",
    table_name: str = "result table",
) -> dict[str, str]:
    """Reject embedded labels or schemas that conflict with requested analysis."""
    status = {
        "modification": "unverified_no_embedded_label",
        "database": "unverified_no_embedded_label",
        "design": "unverified_no_embedded_label",
        "context": "unverified_no_embedded_label",
    }
    if expected_modification is not None:
        expected = MODIFICATION_LABELS.get(
            expected_modification.lower(), expected_modification.lower()
        )
        biological_expected = "m6a" if expected == "m6a_rep2" else expected
        if "analysis_id" in frame:
            observed_analysis = _normalized_values(
                frame, "analysis_id", MODIFICATION_LABELS, table_name
            )
            if observed_analysis != {expected}:
                raise ValueError(
                    f"{table_name} analysis_id {sorted(observed_analysis)} does "
                    f"not match requested modification {expected!r}"
                )
            status["modification"] = "verified_analysis_id"
        modification_column = next(
            (
                column
                for column in ("modification", "analysis_modification")
                if column in frame
            ),
            None,
        )
        if modification_column is not None:
            observed = _normalized_values(
                frame, modification_column, MODIFICATION_LABELS, table_name
            )
            expected_label = (
                biological_expected if "analysis_id" in frame else expected
            )
            if observed != {expected_label}:
                raise ValueError(
                    f"{table_name} modification label {sorted(observed)} does "
                    f"not match requested {expected_label!r}"
                )
            if "analysis_id" not in frame:
                status["modification"] = "verified_embedded_label"

    database_column = next(
        (
            column
            for column in ("database", "analysis_database")
            if column in frame
        ),
        None,
    )
    if expected_database is not None and database_column is not None:
        expected = expected_database.lower()
        observed = {
            str(value).strip().lower()
            for value in frame[database_column].dropna()
        }
        if observed != {expected}:
            raise ValueError(
                f"{table_name} database label {sorted(observed)} does not "
                f"match requested {expected!r}"
            )
        status["database"] = "verified_embedded_label"

    if expected_design is not None:
        expected = DESIGN_LABELS.get(expected_design.lower(), expected_design.lower())
        design_column = next(
            (column for column in ("analysis_design", "design") if column in frame),
            None,
        )
        if design_column is not None:
            observed = _normalized_values(
                frame, design_column, DESIGN_LABELS, table_name
            )
            if observed != {expected}:
                raise ValueError(
                    f"{table_name} design label {sorted(observed)} does not "
                    f"match requested {expected!r}"
                )
            status["design"] = "verified_embedded_label"
        else:
            inferred = (
                "transcript-region"
                if {"mh_or", "region_filter"}.issubset(frame.columns)
                else "loose"
                if "odds_ratio" in frame.columns
                else None
            )
            if inferred is not None and inferred != expected:
                raise ValueError(
                    f"{table_name} schema implies {inferred!r}, not requested "
                    f"{expected!r}"
                )
            if inferred is not None:
                status["design"] = "verified_schema"

    requested_context = context.upper()
    if "context_region" in frame:
        observed = {
            str(value).strip().upper()
            for value in frame["context_region"].dropna()
        }
        if observed != {requested_context}:
            raise ValueError(
                f"{table_name} context label {sorted(observed)} does not "
                f"match requested {context!r}"
            )
        status["context"] = "verified_context_label"
    elif "region_filter" in frame:
        observed = {
            str(value).strip().upper()
            for value in frame["region_filter"].dropna()
        }
        if requested_context not in observed:
            raise ValueError(
                f"{table_name} has no requested context {context!r}; "
                f"available: {sorted(observed)}"
            )
        status["context"] = "verified_region_filter"
    elif requested_context == "ALL":
        status["context"] = "implicit_all_context"
    else:
        raise ValueError(
            f"{table_name} has no region_filter column for context {context!r}"
        )
    return status


def validate_modification_codes(
    frame: pd.DataFrame,
    expected_modification: str,
    table_name: str,
) -> str:
    if "mod_code" not in frame:
        return "unverified_no_modification_code"
    observed = {
        str(value).strip().lower()
        for value in frame["mod_code"].dropna().unique()
    }
    expected = {
        value.lower()
        for value in modification_config(expected_modification).expected_mod_codes
    }
    if not observed or not observed.issubset(expected):
        raise ValueError(
            f"{table_name} modification codes {sorted(observed)} do not match "
            f"requested {expected_modification!r} ({sorted(expected)})"
        )
    if expected_modification in {"m6a", "m6a_rep2"}:
        return "verified_m6a_code_family_replicate_identity_unavailable"
    return "verified_modification_code"
