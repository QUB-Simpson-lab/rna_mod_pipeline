from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from ..config import expression_gene_symbol
from ..evidence_quality import add_evidence_quality
from ..io import require_file
from ..stats import bh_adjust
from .loaders import ExpressionProfile


def annotate_rbp_expression(
    enrichment_path: str | Path,
    profile: ExpressionProfile,
    *,
    database: str | None = None,
    context: str = "all",
) -> pd.DataFrame:
    enrichment = pd.read_csv(require_file(enrichment_path, "RBP enrichment table"), sep="\t")
    if "RBP" not in enrichment:
        raise ValueError("RBP enrichment table must contain RBP")
    enrichment = enrichment.copy()
    database_column = next(
        (
            column
            for column in ("database", "analysis_database")
            if column in enrichment
        ),
        None,
    )
    if database and database_column:
        enrichment = enrichment[
            enrichment[database_column].astype(str).str.lower() == database.lower()
        ].copy()
    if context != "all":
        if "region_filter" not in enrichment:
            raise ValueError(
                "Regional expression integration requires region_filter"
            )
        enrichment = enrichment[
            enrichment["region_filter"].astype(str).str.upper()
            == context.upper()
        ].copy()
    if enrichment.empty:
        raise ValueError("No RBP enrichment rows remain in the selected scope")
    enrichment = add_evidence_quality(enrichment)
    enrichment["gene_symbol"] = enrichment["RBP"].map(expression_gene_symbol)
    enrichment["expression_mapping_status"] = np.where(
        enrichment["gene_symbol"].isna(), "explicitly_unmappable", "mapped"
    )
    expression = profile.table.rename(columns={"gene": "gene_symbol"})
    result = enrichment.merge(expression, on="gene_symbol", how="left", validate="many_to_one")
    result["expression_source"] = profile.source
    result["expression_cell_line"] = profile.cell_line
    result["expression_threshold"] = profile.threshold
    result["expression_status"] = result["expression_status"].fillna(
        "gene_not_available_in_expression_source"
    )
    result["absence_interpretation"] = profile.absence_interpretation
    return result


def _correlation(x: pd.Series, y: pd.Series, method: str) -> tuple[int, float, float]:
    valid = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(valid) < 3 or valid["x"].nunique() < 2 or valid["y"].nunique() < 2:
        return len(valid), math.nan, math.nan
    result = (
        stats.spearmanr(valid["x"], valid["y"])
        if method == "spearman"
        else stats.pearsonr(valid["x"], valid["y"])
    )
    return len(valid), float(result.statistic), float(result.pvalue)


def modification_expression_association(
    metagene_path: str | Path,
    profile: ExpressionProfile,
    *,
    gene_lengths_path: str | Path | None = None,
    gene_column: str = "gene_name",
    context: str = "all",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sites = pd.read_csv(require_file(metagene_path, "metagene table"), sep="\t")
    if gene_column not in sites:
        raise ValueError(f"Metagene table must contain {gene_column}")
    if context != "all":
        if "region" not in sites:
            raise ValueError("Regional expression integration requires region")
        sites = sites[
            sites["region"].astype(str).str.upper() == context.upper()
        ].copy()
    genes = sites[gene_column].astype("string").str.strip()
    valid = genes.notna() & ~genes.str.lower().isin({"", "nan", "none", "<na>"})
    counts = (
        genes[valid].value_counts().rename_axis("gene").rename("modification_site_count")
        .reset_index()
    )
    gene_table = profile.table.merge(counts, on="gene", how="left")
    gene_table["modification_site_count"] = (
        gene_table["modification_site_count"].fillna(0).astype(int)
    )
    metrics = ["modification_site_count"]
    if gene_lengths_path:
        lengths = pd.read_csv(require_file(gene_lengths_path, "gene-length table"), sep="\t")
        if not {"gene", "length_nt"}.issubset(lengths):
            raise ValueError("Gene-length table must contain gene and length_nt")
        lengths["length_nt"] = pd.to_numeric(lengths["length_nt"], errors="coerce")
        lengths = lengths.groupby("gene", as_index=False)["length_nt"].max()
        gene_table = gene_table.merge(lengths, on="gene", how="left")
        gene_table["modification_sites_per_kb"] = (
            gene_table["modification_site_count"]
            / pd.to_numeric(gene_table["length_nt"], errors="coerce")
            * 1000
        )
        metrics.append("modification_sites_per_kb")

    expressed = gene_table["expression_status"].isin({"detected", "expressed"})
    modified = gene_table["modification_site_count"].gt(0)
    subsets = (
        (
            "detected_or_expressed_and_modified",
            "primary_legacy_compatible",
            "Genes with at least one modification site and expression status "
            "detected or expressed.",
            gene_table[modified & expressed],
        ),
        (
            "all_expression_genes",
            "exploratory_zero_inclusive",
            "All genes represented in the expression source; genes without "
            "modification sites have a site count of zero.",
            gene_table,
        ),
        (
            "genes_with_at_least_one_modification",
            "exploratory_modified_only",
            "All genes with at least one modification site, regardless of "
            "expression detection or threshold status.",
            gene_table[modified],
        ),
    )

    rows = []
    for metric in metrics:
        for subset_name, role, definition, subset in subsets:
            for method in ("spearman", "pearson"):
                n, coefficient, p_value = _correlation(
                    np.log2(subset["expression_value"].astype(float) + 1),
                    subset[metric],
                    method,
                )
                rows.append(
                    {
                        "expression_source": profile.source,
                        "cell_line": profile.cell_line,
                        "subset": subset_name,
                        "estimand_role": role,
                        "estimand_definition": definition,
                        "modification_metric": metric,
                        "expression_transform": "log2(value+1)",
                        "method": method,
                        "legacy_parity_expected": (
                            role == "primary_legacy_compatible"
                            and metric == "modification_site_count"
                            and method == "spearman"
                        ),
                        "n_genes": n,
                        "coefficient": coefficient,
                        "p_value": p_value,
                    }
                )
    summary = pd.DataFrame(rows)
    summary["fdr_within_run"] = bh_adjust(summary["p_value"])
    return gene_table, summary
