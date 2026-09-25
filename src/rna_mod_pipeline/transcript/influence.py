from __future__ import annotations

import math

import numpy as np

from .opportunities import OpportunityDesign


def leave_one_gene_out_diagnostics(
    design: OpportunityDesign,
    numerator: float,
    denominator: float,
    log_or: float,
    estimable: bool,
    direction: str,
    gene_numerator: np.ndarray,
    gene_denominator: np.ndarray,
    gene_influence: np.ndarray,
    analysis_id: str,
    modification: str,
    database: str,
    rbp: str,
    region_filter: str | None,
) -> tuple[float, bool, int, int, list[dict]]:
    maximum_influence = (
        float(np.nanmax(np.abs(gene_influence)))
        if np.isfinite(gene_influence).any()
        else math.nan
    )
    stable = False
    reversals = 0
    boundaries = 0
    rows: list[dict] = []
    if not estimable:
        return maximum_influence, stable, reversals, boundaries, rows

    loo_numerator = numerator - gene_numerator
    loo_denominator = denominator - gene_denominator
    valid = (loo_numerator > 0) & (loo_denominator > 0)
    loo_or = np.full(design.n_genes, np.nan)
    loo_or[valid] = loo_numerator[valid] / loo_denominator[valid]
    loo_direction = np.where(
        loo_or > 1,
        "enriched",
        np.where(loo_or < 1, "depleted", "no_effect"),
    )
    reversed_mask = valid & (loo_direction != direction)
    reversals = int(reversed_mask.sum())
    boundaries = int((~valid).sum())
    stable = bool(valid.all() and not reversed_mask.any())
    delta = np.abs(np.log(loo_or) - log_or)
    order = np.argsort(
        np.where(np.isfinite(delta), delta, -1), kind="mergesort"
    )[::-1][:10]
    for gene_index in order:
        if not np.isfinite(delta[gene_index]):
            continue
        rows.append(
            {
                "analysis_id": analysis_id,
                "modification": modification,
                "database": database,
                "RBP": rbp,
                "region_filter": region_filter or "all",
                "gene_id": str(design.genes[gene_index]),
                "gene_name": str(design.gene_names[gene_index]),
                "mh_numerator_contribution": gene_numerator[gene_index],
                "mh_denominator_contribution": gene_denominator[gene_index],
                "cluster_influence": gene_influence[gene_index],
                "leave_one_gene_out_or": loo_or[gene_index],
                "absolute_leave_one_out_delta_log_or": delta[gene_index],
                "direction_reversal": bool(reversed_mask[gene_index]),
            }
        )
    return maximum_influence, stable, reversals, boundaries, rows
