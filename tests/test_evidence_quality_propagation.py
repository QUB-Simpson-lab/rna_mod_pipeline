from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from rna_mod_pipeline.binding.loose import run_loose_overlap
from rna_mod_pipeline.binding.manifest import loose_manifest_parameters
from rna_mod_pipeline.binding.plots import (
    _machinery_labels,
    _top_effect_rows,
    plot_loose_results,
)
from rna_mod_pipeline.binding.sources import BindingTrack
from rna_mod_pipeline.crossdb.analysis import (
    compare_databases,
    recompute_consensus_fields,
)
from rna_mod_pipeline.crossdb.plot_data import pairwise_plot_source
from rna_mod_pipeline.crossmod.analysis import compare_modifications
from rna_mod_pipeline.crossmod.plots import create_crossmod_plots
from rna_mod_pipeline.dataset_compare.enrichment import _compare_pair
from rna_mod_pipeline.dataset_compare.plots import plot_enrichment_concordance
from rna_mod_pipeline.evidence_quality import (
    add_evidence_quality,
    add_significance_calls,
    evidence_quality_manifest,
)
from rna_mod_pipeline.expression.analysis import annotate_rbp_expression
from rna_mod_pipeline.expression.loaders import ExpressionProfile
from rna_mod_pipeline.intervals import build_interval_index
from rna_mod_pipeline.knockrbp.analysis import regulatory_network
from rna_mod_pipeline.knockrbp.plots import regulatory_plot_source
from rna_mod_pipeline.stringdb.selection import (
    consensus_selection_audit,
    select_consensus_rbps,
)


def database_table(database: str, *, low_quality: bool = True) -> pd.DataFrame:
    rows = [
        {
            "RBP": "GOOD",
            "odds_ratio": 2.0,
            "fdr": 0.01,
            "p_value": 0.001,
            "database": database,
            "source_file": f"{database}.tsv",
            "partial_input": False,
            "estimable": True,
            "boundary_estimate": False,
            "n_gene_clusters_for_robust_variance": 30,
            "low_gene_cluster_count": False,
            "direction_stable_after_each_gene_removed": True,
            "leave_one_gene_out_direction_reversals": 0,
            "leave_one_gene_out_boundaries": 0,
            "maximum_absolute_gene_influence": 0.1,
        },
        {
            "RBP": "LOW",
            "odds_ratio": 3.0,
            "fdr": 0.005,
            "p_value": 0.0005,
            "database": database,
            "source_file": f"{database}.tsv",
            "partial_input": False,
            "estimable": True,
            "boundary_estimate": False,
            "n_gene_clusters_for_robust_variance": 10 if low_quality else 30,
            "low_gene_cluster_count": low_quality,
            "direction_stable_after_each_gene_removed": True,
            "leave_one_gene_out_direction_reversals": 0,
            "leave_one_gene_out_boundaries": 0,
            "maximum_absolute_gene_influence": 0.2,
        },
    ]
    return add_evidence_quality(pd.DataFrame(rows))


def cross_database_table() -> pd.DataFrame:
    combined, _ = compare_databases(
        {
            database: database_table(database)
            for database in ("ornament", "encori", "postar3")
        }
    )
    return combined


class CrossDatabaseQualityTests(unittest.TestCase):
    def test_primary_and_raw_calls_are_both_retained(self) -> None:
        combined = cross_database_table().set_index("RBP")
        self.assertTrue(bool(combined.loc["GOOD", "triple_significant"]))
        self.assertTrue(
            bool(combined.loc["GOOD", "triple_statistically_significant"])
        )
        self.assertFalse(bool(combined.loc["LOW", "triple_significant"]))
        self.assertTrue(
            bool(combined.loc["LOW", "triple_statistically_significant"])
        )
        self.assertEqual(
            combined.loc["LOW", "consensus_evidence_quality"],
            "sensitivity_only_consensus_quality_excluded",
        )
        self.assertIn(
            "low_gene_cluster_count",
            combined.loc["LOW", "consensus_quality_reasons"],
        )

    def test_pairwise_primary_correlation_flag_excludes_warning_rows(self) -> None:
        combined, pairwise = compare_databases(
            {
                database: database_table(database)
                for database in ("ornament", "encori", "postar3")
            }
        )
        source = pairwise_plot_source(combined, pairwise, 0.05)
        pair = source[source["pair"].eq("ornament_vs_encori")].set_index("RBP")
        self.assertTrue(
            bool(pair.loc["GOOD", "eligible_for_primary_correlation"])
        )
        self.assertFalse(
            bool(pair.loc["LOW", "eligible_for_primary_correlation"])
        )


class StringSelectionQualityTests(unittest.TestCase):
    def test_quality_excluded_candidates_remain_in_audit_but_not_primary_selection(self) -> None:
        cross = cross_database_table()
        audit = consensus_selection_audit(
            cross, top_n=10, directions=("enriched",)
        ).set_index("RBP")
        self.assertEqual(audit.loc["GOOD", "selection_status"], "selected_primary")
        self.assertEqual(
            audit.loc["LOW", "selection_status"],
            "sensitivity_only_quality_excluded",
        )
        selected = select_consensus_rbps(
            cross, top_n=10, directions=("enriched",)
        )
        self.assertEqual(list(selected["RBP"]), ["GOOD"])

    def test_legacy_table_without_diagnostics_remains_selectable(self) -> None:
        legacy = pd.DataFrame(
            {
                "RBP": ["LEGACY"],
                "encori_odds_ratio": [2.0],
                "encori_fdr": [0.01],
                "encori_direction": ["Enriched"],
                "postar3_odds_ratio": [2.5],
                "postar3_fdr": [0.02],
                "postar3_direction": ["Enriched"],
            }
        )
        selected = select_consensus_rbps(
            legacy, top_n=10, directions=("enriched",)
        )
        self.assertEqual(list(selected["RBP"]), ["LEGACY"])
        self.assertEqual(
            selected.loc[0, "encori_evidence_quality"],
            "primary_eligible_core_only",
        )


class CrossModificationQualityTests(unittest.TestCase):
    def test_primary_patterns_exclude_warning_rows_but_raw_patterns_retain_them(self) -> None:
        tables = {
            ("m6a", "encori"): database_table("encori", low_quality=True),
            ("m5c", "encori"): database_table("encori", low_quality=False),
        }
        matrix, patterns, _, overlaps = compare_modifications(
            tables,
            modifications=("m6a", "m5c"),
            databases=("encori",),
        )
        low_m6a = matrix[
            matrix["RBP"].eq("LOW") & matrix["modification"].eq("m6a")
        ].iloc[0]
        self.assertEqual(low_m6a["statistical_state"], "significant_enriched")
        self.assertEqual(low_m6a["state"], "quality_excluded")
        low_pattern = patterns[patterns["RBP"].eq("LOW")].iloc[0]
        self.assertEqual(low_pattern["n_modifications_significant"], 1)
        self.assertEqual(
            low_pattern["n_modifications_statistically_significant"], 2
        )
        overlap = overlaps.iloc[0]
        self.assertEqual(overlap["n_shared_significant"], 1)
        self.assertEqual(overlap["n_shared_statistically_significant"], 2)

    def test_primary_plots_use_eligible_metrics_and_retain_raw_sources(self) -> None:
        tables = {
            ("m6a", "encori"): database_table("encori", low_quality=True),
            ("m5c", "encori"): database_table("encori", low_quality=False),
        }
        matrix, patterns, correlations, _ = compare_modifications(
            tables,
            modifications=("m6a", "m5c"),
            databases=("encori",),
        )
        with tempfile.TemporaryDirectory() as temporary:
            paths = create_crossmod_plots(
                matrix, patterns, correlations, temporary
            )
            names = {path.name for path in paths}
            self.assertIn("crossmod_effect_correlations.png", names)
            self.assertIn(
                "crossmod_effect_correlations_raw_sensitivity.png", names
            )
            correlation_source = pd.read_csv(
                Path(temporary) / "crossmod_plot_correlation_source.tsv",
                sep="\t",
            )
            self.assertEqual(
                set(correlation_source["primary_plot_metric"]),
                {"inference_eligible_coefficient"},
            )
            effect_source = pd.read_csv(
                Path(temporary) / "crossmod_plot_effect_source.tsv", sep="\t"
            )
            low_m6a = effect_source[
                effect_source["RBP"].eq("LOW")
                & effect_source["modification"].eq("m6a")
            ].iloc[0]
            self.assertFalse(bool(low_m6a["inference_eligible"]))
            self.assertEqual(low_m6a["state"], "quality_excluded")


class OtherDownstreamPropagationTests(unittest.TestCase):
    def test_regulatory_plot_source_is_bounded_and_deterministic(self) -> None:
        table = pd.DataFrame(
            {
                "affected_rbp": [f"RBP{i:02d}" for i in range(45)],
                "knocked_down_rbp": "KDRBP",
                "dataset_id": "dataset",
                "cell_line": "MDA-MB-231",
                "log2fc": [float(i) / 10 for i in range(45)],
            }
        )
        source = regulatory_plot_source(table, top_n=40)
        self.assertEqual(source["affected_rbp"].nunique(), 40)
        self.assertIn("RBP44", set(source["affected_rbp"]))
        self.assertNotIn("RBP00", set(source["affected_rbp"]))
        self.assertEqual(
            source.loc[source["regulatory_plot_rank"].eq(1), "affected_rbp"].iloc[0],
            "RBP44",
        )

    def test_knockrbp_edges_retain_consensus_quality(self) -> None:
        cross = cross_database_table()
        deg_data = {
            "dataset": {
                "rbp": "KDRBP",
                "cell_line": "MDA-MB-231",
                "genes": {
                    "GOOD": {"log2fc": 1.0, "padj": 0.01, "direction": "up"},
                    "LOW": {"log2fc": -1.0, "padj": 0.02, "direction": "down"},
                },
            }
        }
        network = regulatory_network(deg_data, cross).set_index("affected_rbp")
        self.assertEqual(
            network.loc["GOOD", "consensus_evidence_quality"],
            "primary_eligible_consensus",
        )
        self.assertEqual(
            network.loc["LOW", "consensus_evidence_quality"],
            "sensitivity_only_consensus_quality_excluded",
        )
        self.assertFalse(bool(network.loc["LOW", "encori_inference_eligible"]))

    def test_knockrbp_recomputes_stale_aggregate_consensus(self) -> None:
        cross = cross_database_table()
        low = cross["RBP"].eq("LOW")
        cross.loc[low, "triple_significant"] = True
        cross.loc[low, "n_databases_significant"] = 3
        cross.loc[low, "consensus_evidence_quality"] = (
            "primary_eligible_consensus"
        )
        deg_data = {
            "dataset": {
                "rbp": "KDRBP",
                "cell_line": "MDA-MB-231",
                "genes": {
                    "LOW": {"log2fc": 1.0, "padj": 0.01, "direction": "up"}
                },
            }
        }
        observed = regulatory_network(deg_data, cross).iloc[0]
        self.assertFalse(bool(observed["triple_significant"]))
        self.assertEqual(int(observed["n_databases_significant"]), 0)
        self.assertEqual(
            observed["consensus_evidence_quality"],
            "sensitivity_only_consensus_quality_excluded",
        )

    def test_knockrbp_recompute_preserves_missing_database_as_not_covered(self) -> None:
        cross, _ = compare_databases(
            {
                database: database_table(database)
                for database in ("encori", "postar3")
            },
            allow_missing=True,
        )
        good = cross["RBP"].eq("GOOD")
        cross.loc[good, "n_databases_significant"] = 0
        cross.loc[good, "consensus_direction"] = "Conflicting"
        cross.loc[good, "consensus_evidence_quality"] = "stale"

        recomputed = recompute_consensus_fields(cross).set_index("RBP")
        self.assertEqual(int(recomputed.loc["GOOD", "n_databases_present"]), 2)
        self.assertEqual(
            recomputed.loc["GOOD", "ornament_primary_direction"],
            "Unavailable",
        )
        self.assertEqual(
            recomputed.loc["GOOD", "ornament_evidence_quality"],
            "not_covered",
        )

        deg_data = {
            "dataset": {
                "rbp": "KDRBP",
                "cell_line": "MDA-MB-231",
                "genes": {
                    "GOOD": {
                        "log2fc": 1.0,
                        "padj": 0.01,
                        "direction": "up",
                    }
                },
            }
        }
        observed = regulatory_network(deg_data, cross).iloc[0]
        self.assertEqual(int(observed["n_databases_present"]), 2)
        self.assertFalse(bool(observed["ornament_present"]))
        self.assertEqual(int(observed["n_databases_significant"]), 2)
        self.assertEqual(observed["consensus_direction"], "Enriched")
        self.assertEqual(
            observed["consensus_evidence_quality"],
            "primary_eligible_consensus",
        )

    def test_expression_annotation_propagates_quality_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "enrichment.tsv"
            table = database_table("encori")
            table.to_csv(path, sep="\t", index=False)
            expression = pd.DataFrame(
                {
                    "gene": ["GOOD", "LOW"],
                    "expression_value": [5.0, 3.0],
                    "detection_metric": [5.0, 3.0],
                    "expression_status": ["expressed", "expressed"],
                }
            )
            profile = ExpressionProfile(
                table=expression,
                source="synthetic",
                cell_line="MDA-MB-231",
                value_column="TPM",
                threshold=1.0,
                absence_interpretation="synthetic",
            )
            annotated = annotate_rbp_expression(path, profile, database="encori")
            observed = annotated.set_index("RBP")
            self.assertTrue(bool(observed.loc["GOOD", "inference_eligible"]))
            self.assertFalse(bool(observed.loc["LOW", "inference_eligible"]))

    def test_dataset_comparison_keeps_raw_and_primary_agreement_separate(self) -> None:
        first = database_table("encori", low_quality=False)
        second = database_table("encori", low_quality=True)
        for table in (first, second):
            for label in ("modification", "database", "design", "context"):
                table[f"input_{label}_status"] = "verified"
        detail, summary = _compare_pair(
            first, second, "transcript-region", "encori", 0.05
        )
        low = detail[detail["RBP"].eq("LOW")].iloc[0]
        self.assertTrue(bool(low["both_significant"]))
        self.assertFalse(bool(low["both_primary_significant"]))
        self.assertEqual(summary["n_significant_both"], 2)
        self.assertEqual(summary["n_primary_significant_both"], 1)
        self.assertTrue(
            bool(
                detail.loc[
                    detail["RBP"].eq("GOOD"),
                    "eligible_for_primary_effect_comparison",
                ].iloc[0]
            )
        )
        self.assertFalse(
            bool(low["eligible_for_primary_effect_comparison"])
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "dataset_concordance.png"
            plot_enrichment_concordance(
                detail,
                path,
                "dataset_a",
                "dataset_b",
                100,
            )
            self.assertTrue(path.is_file())


class LoosePrimarySignificanceTests(unittest.TestCase):
    def test_partial_run_retains_raw_call_but_excludes_site_annotation_and_top_plot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case_positions = list(range(1_000, 1_020))
            omitted_positions = list(range(1_020, 1_040))
            background_positions = list(range(2_000, 2_100))
            sites = pd.DataFrame(
                {
                    "chrom": "chr1",
                    "start": case_positions,
                    "strand": "+",
                    "mod_code": "a",
                    "Nvalid_cov": 20,
                    "fraction_modified": 30.0,
                }
            )
            sites_path = root / "sites.tsv"
            sites.to_csv(sites_path, sep="\t", index=False)
            bedmethyl_path = root / "calls.bedmethyl"
            rows = []
            for position in [*case_positions, *omitted_positions]:
                rows.append(
                    [
                        "chr1",
                        position,
                        position + 1,
                        "a",
                        0,
                        "+",
                        position,
                        position + 1,
                        "0",
                        20,
                        30.0,
                    ]
                )
            for position in background_positions:
                rows.append(
                    [
                        "chr1",
                        position,
                        position + 1,
                        "a",
                        0,
                        "+",
                        position,
                        position + 1,
                        "0",
                        20,
                        0.0,
                    ]
                )
            pd.DataFrame(rows).to_csv(
                bedmethyl_path, sep="\t", index=False, header=False
            )
            bound = [*case_positions[:15], *background_positions[:10]]
            track = BindingTrack(
                rbp_name="TESTRBP",
                canonical_rbp="TESTRBP",
                intervals=build_interval_index(
                    {"chr1": [(position, position + 1) for position in bound]}
                ),
                source_names=("TESTRBP",),
            )
            result = run_loose_overlap(
                sites_path,
                bedmethyl_path,
                [track],
                modification="m6a",
                database="ornament",
                window=0,
                site_annotation="significant",
                allow_site_subset=True,
            )
            manifest_parameters = loose_manifest_parameters(
                result,
                modification="m6a",
                database="ornament",
                catalog_version="test_catalog",
                options={
                    "allow_site_subset": True,
                    "window": 0,
                    "background_min_coverage": 20,
                    "background_max_fraction": 20.0,
                    "case_min_coverage": 20,
                    "case_min_fraction": 20.0,
                    "fdr_threshold": 0.05,
                    "minimum_clip_experiments": 1,
                    "cell_types": None,
                    "methods": None,
                    "site_annotation": "significant",
                    "allow_empty_resources": False,
                },
            )
            plot_dir = root / "plots"
            plot_loose_results(
                result.enrichment,
                "m6A",
                "ornament",
                plot_dir,
            )
            self.assertEqual(len(list(plot_dir.glob("*.png"))), 4)
        row = result.enrichment.iloc[0]
        self.assertTrue(bool(row["statistically_significant"]))
        self.assertFalse(bool(row["significant"]))
        self.assertEqual(row["statistical_direction"], "Enriched")
        self.assertEqual(row["primary_direction"], "Quality excluded")
        self.assertEqual(manifest_parameters["case_sites"], 20)
        self.assertEqual(manifest_parameters["omitted_qualifying_case_sites"], 20)
        self.assertEqual(
            manifest_parameters["evidence_quality"][
                "n_raw_significant_quality_excluded"
            ],
            1,
        )
        self.assertTrue(result.annotated_sites["overlapping_RBPs"].eq("").all())
        self.assertTrue(result.annotated_sites["n_overlapping_RBPs"].eq(0).all())
        top = _top_effect_rows(result.enrichment, enriched=True, top_n=30)
        self.assertTrue(top.empty)
        self.assertEqual(_machinery_labels(result.enrichment), ["TESTRBP"])

    def test_quality_manifest_separates_raw_and_primary_counts(self) -> None:
        calls = add_significance_calls(
            pd.DataFrame(
                [
                    {"RBP": "GOOD", "odds_ratio": 2.0, "fdr": 0.01},
                    {
                        "RBP": "PART",
                        "odds_ratio": 3.0,
                        "fdr": 0.01,
                        "partial_input": True,
                    },
                ]
            )
        )
        summary = evidence_quality_manifest(calls)
        self.assertEqual(summary["n_raw_statistically_significant"], 2)
        self.assertEqual(summary["n_primary_significant"], 1)
        self.assertEqual(summary["n_raw_significant_quality_excluded"], 1)


if __name__ == "__main__":
    unittest.main()
