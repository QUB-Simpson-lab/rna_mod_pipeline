from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .client import CachedStringClient, StringAPIError

SPECIES = 9606
CALLER = "rna_modification_qub_refactored"
CONFIDENCE = {"medium": 400, "high": 700}


def resolve_symbols(client: CachedStringClient, symbols: list[str]) -> pd.DataFrame:
    if not symbols:
        return pd.DataFrame(
            columns=["query_symbol", "string_id", "preferred_name", "resolved"]
        )
    response = client.call(
        "get_string_ids",
        {
            "identifiers": "\r".join(sorted(set(symbols))),
            "species": SPECIES,
            "limit": 1,
            "caller_identity": CALLER,
        },
    )
    if not isinstance(response, list):
        raise StringAPIError("STRING identifier response was not a list")
    by_query = {}
    for record in response:
        query = str(record.get("queryItem", ""))
        if query in by_query:
            raise StringAPIError(f"Duplicate STRING resolution for {query}")
        by_query[query] = record
    return pd.DataFrame(
        [
            {
                "query_symbol": symbol,
                "string_id": str(by_query.get(symbol, {}).get("stringId", "")),
                "preferred_name": str(by_query.get(symbol, {}).get("preferredName", "")),
                "resolved": bool(by_query.get(symbol, {}).get("stringId")),
            }
            for symbol in sorted(set(symbols))
        ]
    )


def _interactions(response: object) -> pd.DataFrame:
    columns = [
        "protein_A",
        "protein_B",
        "combined_score",
        "neighborhood_score",
        "fusion_score",
        "cooccurrence_score",
        "coexpression_score",
        "experimental_score",
        "database_score",
        "textmining_score",
    ]
    if not isinstance(response, list):
        raise StringAPIError("STRING network response was not a list")
    return pd.DataFrame(
        [
            {
                "protein_A": edge.get("preferredName_A", ""),
                "protein_B": edge.get("preferredName_B", ""),
                "combined_score": edge.get("score", np.nan),
                "neighborhood_score": edge.get("nscore", np.nan),
                "fusion_score": edge.get("fscore", np.nan),
                "cooccurrence_score": edge.get("pscore", np.nan),
                "coexpression_score": edge.get("ascore", np.nan),
                "experimental_score": edge.get("escore", np.nan),
                "database_score": edge.get("dscore", np.nan),
                "textmining_score": edge.get("tscore", np.nan),
            }
            for edge in response
        ],
        columns=columns,
    )


def _enrichment(response: object, fdr_threshold: float) -> pd.DataFrame:
    columns = [
        "category",
        "term",
        "description",
        "p_value",
        "fdr",
        "gene_count",
        "background_count",
        "genes",
    ]
    if not isinstance(response, list):
        raise StringAPIError("STRING enrichment response was not a list")
    rows = []
    for term in response:
        fdr = pd.to_numeric(pd.Series([term.get("fdr")]), errors="coerce").iloc[0]
        if pd.isna(fdr) or float(fdr) >= fdr_threshold:
            continue
        rows.append(
            {
                "category": term.get("category", ""),
                "term": term.get("term", ""),
                "description": term.get("description", ""),
                "p_value": term.get("p_value", np.nan),
                "fdr": float(fdr),
                "gene_count": term.get("number_of_genes", np.nan),
                "background_count": term.get("number_of_genes_in_background", np.nan),
                "genes": ",".join(term.get("preferredNames", []) or []),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _ppi(response: object, direction: str, confidence: str) -> pd.DataFrame:
    columns = [
        "direction",
        "confidence",
        "n_nodes",
        "n_edges",
        "expected_edges",
        "average_node_degree",
        "local_clustering_coefficient",
        "ppi_enrichment_p_value",
    ]
    record = response[0] if isinstance(response, list) and response else response
    if not isinstance(record, dict):
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(
        [
            {
                "direction": direction,
                "confidence": confidence,
                "n_nodes": record.get("number_of_nodes", np.nan),
                "n_edges": record.get("number_of_edges", np.nan),
                "expected_edges": record.get("expected_number_of_edges", np.nan),
                "average_node_degree": record.get("average_node_degree", np.nan),
                "local_clustering_coefficient": record.get(
                    "local_clustering_coefficient", np.nan
                ),
                "ppi_enrichment_p_value": record.get("p_value", np.nan),
            }
        ]
    )[columns]


def run_string_groups(
    selection: pd.DataFrame,
    client: CachedStringClient,
    output_dir: str | Path,
    *,
    requested_directions: tuple[str, ...] = ("enriched", "depleted"),
    fdr_threshold: float = 0.05,
) -> tuple[list[Path], pd.DataFrame]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    output_files: list[Path] = []
    summaries = []
    for direction in requested_directions:
        group = (
            selection[selection["direction"].eq(direction)].copy()
            if "direction" in selection
            else pd.DataFrame(columns=["RBP"])
        )
        symbols = list(group["RBP"].astype(str))
        resolved = resolve_symbols(client, symbols)
        resolution_path = output / f"string_ids_{direction}.tsv"
        resolved.to_csv(resolution_path, sep="\t", index=False)
        output_files.append(resolution_path)
        identifiers = list(resolved.loc[resolved["resolved"], "string_id"].astype(str))
        if len(identifiers) < 2:
            summaries.append(
                {
                    "direction": direction,
                    "selected_rbps": len(symbols),
                    "resolved_proteins": len(identifiers),
                    "status": "insufficient_resolved_proteins",
                }
            )
            for confidence in CONFIDENCE:
                path = output / f"interactions_{direction}_{confidence}.tsv"
                _interactions([]).to_csv(path, sep="\t", index=False)
                output_files.append(path)
            ppi_path = output / f"ppi_enrichment_{direction}.tsv"
            pd.DataFrame(
                columns=[
                    "direction",
                    "confidence",
                    "n_nodes",
                    "n_edges",
                    "expected_edges",
                    "average_node_degree",
                    "local_clustering_coefficient",
                    "ppi_enrichment_p_value",
                ]
            ).to_csv(ppi_path, sep="\t", index=False)
            output_files.append(ppi_path)
            path = output / f"functional_enrichment_{direction}.tsv"
            _enrichment([], fdr_threshold).to_csv(path, sep="\t", index=False)
            output_files.append(path)
            continue
        common = {
            "identifiers": "\r".join(identifiers),
            "species": SPECIES,
            "caller_identity": CALLER,
        }
        ppi_tables = []
        for confidence, required_score in CONFIDENCE.items():
            parameters = {**common, "required_score": required_score}
            interactions = _interactions(client.call("network", parameters))
            interaction_path = output / f"interactions_{direction}_{confidence}.tsv"
            interactions.to_csv(interaction_path, sep="\t", index=False)
            output_files.append(interaction_path)
            ppi_tables.append(
                _ppi(
                    client.call("ppi_enrichment", parameters),
                    direction,
                    confidence,
                )
            )
        ppi = pd.concat(ppi_tables, ignore_index=True)
        ppi_path = output / f"ppi_enrichment_{direction}.tsv"
        ppi.to_csv(ppi_path, sep="\t", index=False)
        output_files.append(ppi_path)
        enrichment = _enrichment(
            client.call("enrichment", common), fdr_threshold
        )
        enrichment_path = output / f"functional_enrichment_{direction}.tsv"
        enrichment.to_csv(enrichment_path, sep="\t", index=False)
        output_files.append(enrichment_path)
        summaries.append(
            {
                "direction": direction,
                "selected_rbps": len(symbols),
                "resolved_proteins": len(identifiers),
                "status": "analysed",
                "significant_enrichment_terms": len(enrichment),
            }
        )
    return output_files, pd.DataFrame(summaries)
