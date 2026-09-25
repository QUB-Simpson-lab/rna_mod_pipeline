from __future__ import annotations

import math
import unittest

import pandas as pd

from rna_mod_pipeline.evidence_quality import (
    add_evidence_quality,
    add_prefixed_evidence_quality,
    classify_evidence_quality,
)


def robust_record(**updates: object) -> dict[str, object]:
    record: dict[str, object] = {
        "odds_ratio": 2.0,
        "fdr": 0.01,
        "partial_input": False,
        "estimable": True,
        "boundary_estimate": False,
        "n_gene_clusters_for_robust_variance": 25,
        "low_gene_cluster_count": False,
        "direction_stable_after_each_gene_removed": True,
        "leave_one_gene_out_direction_reversals": 0,
        "leave_one_gene_out_boundaries": 0,
        "maximum_absolute_gene_influence": 0.25,
    }
    record.update(updates)
    return record


class EvidenceQualityClassificationTests(unittest.TestCase):
    def test_complete_robust_estimate_is_primary_eligible(self) -> None:
        quality = classify_evidence_quality(robust_record())
        self.assertTrue(quality["estimate_eligible"])
        self.assertTrue(quality["inference_eligible"])
        self.assertEqual(quality["evidence_quality"], "primary_eligible_robust")
        self.assertEqual(quality["input_completeness"], "complete")
        self.assertEqual(quality["gene_cluster_quality"], "adequate")
        self.assertEqual(quality["leave_one_gene_out_quality"], "stable")

    def test_partial_input_is_retained_but_not_inference_eligible(self) -> None:
        quality = classify_evidence_quality(robust_record(partial_input=True))
        self.assertFalse(quality["inference_eligible"])
        self.assertIn("partial_input", quality["inference_exclusion_reasons"])
        self.assertEqual(quality["input_completeness"], "partial")

    def test_zero_and_infinite_odds_ratios_are_boundary_estimates(self) -> None:
        for odds in (0.0, math.inf):
            with self.subTest(odds=odds):
                quality = classify_evidence_quality(
                    robust_record(odds_ratio=odds, estimable=False)
                )
                self.assertFalse(quality["inference_eligible"])
                self.assertEqual(quality["estimate_quality"], "boundary")
                self.assertIn(
                    "boundary_estimate", quality["inference_exclusion_reasons"]
                )

    def test_low_clusters_and_leave_one_gene_out_failures_are_sensitivity_only(self) -> None:
        low = classify_evidence_quality(
            robust_record(
                n_gene_clusters_for_robust_variance=19,
                low_gene_cluster_count=False,
            )
        )
        self.assertTrue(low["estimate_eligible"])
        self.assertFalse(low["inference_eligible"])
        self.assertEqual(
            low["evidence_quality"], "sensitivity_only_robustness_warning"
        )
        self.assertIn("low_gene_cluster_count", low["evidence_quality_reasons"])

        unstable = classify_evidence_quality(
            robust_record(
                direction_stable_after_each_gene_removed=False,
                leave_one_gene_out_direction_reversals=1,
            )
        )
        self.assertFalse(unstable["inference_eligible"])
        self.assertIn(
            "leave_one_gene_out_direction_reversal",
            unstable["inference_exclusion_reasons"],
        )

    def test_missing_robustness_diagnostics_are_backward_compatible(self) -> None:
        quality = classify_evidence_quality({"odds_ratio": 1.5, "fdr": 0.02})
        self.assertTrue(quality["inference_eligible"])
        self.assertEqual(
            quality["evidence_quality"], "primary_eligible_core_only"
        )
        self.assertEqual(
            quality["input_completeness"], "not_reported_assumed_complete"
        )
        self.assertIn(
            "gene_cluster_diagnostics_unavailable",
            quality["evidence_quality_reasons"],
        )

        nan_diagnostics = classify_evidence_quality(
            {
                "odds_ratio": 1.5,
                "fdr": 0.02,
                "n_gene_clusters_for_robust_variance": math.nan,
                "leave_one_gene_out_direction_reversals": math.nan,
                "leave_one_gene_out_boundaries": math.nan,
            }
        )
        self.assertEqual(
            nan_diagnostics["quality_diagnostics_available"], "none"
        )

    def test_single_gene_influence_has_no_arbitrary_magnitude_exclusion(self) -> None:
        quality = classify_evidence_quality(
            robust_record(maximum_absolute_gene_influence=0.99)
        )
        self.assertTrue(quality["inference_eligible"])
        self.assertEqual(
            quality["single_gene_influence_status"],
            "available_direction_stable",
        )

    def test_finite_fdr_is_required_but_significance_is_not(self) -> None:
        nonsignificant = classify_evidence_quality(robust_record(fdr=0.9))
        self.assertTrue(nonsignificant["inference_eligible"])
        missing = classify_evidence_quality(robust_record(fdr=math.nan))
        self.assertFalse(missing["inference_eligible"])
        self.assertIn("fdr_unavailable", missing["inference_exclusion_reasons"])
        for invalid_value in (-0.01, 1.01):
            with self.subTest(invalid_fdr=invalid_value):
                invalid = classify_evidence_quality(
                    robust_record(fdr=invalid_value)
                )
                self.assertFalse(invalid["inference_eligible"])
                self.assertIn(
                    "invalid_fdr", invalid["inference_exclusion_reasons"]
                )


class EvidenceQualityFrameTests(unittest.TestCase):
    def test_frame_and_prefixed_helpers_preserve_rows_and_add_fields(self) -> None:
        frame = pd.DataFrame(
            [robust_record(RBP="GOOD"), robust_record(RBP="PART", partial_input=True)]
        )
        qualified = add_evidence_quality(frame)
        self.assertEqual(list(qualified["RBP"]), ["GOOD", "PART"])
        self.assertEqual(list(qualified["inference_eligible"]), [True, False])

        prefixed = frame.rename(
            columns={column: f"encori_{column}" for column in frame if column != "RBP"}
        )
        qualified_prefixed = add_prefixed_evidence_quality(prefixed, "encori_")
        self.assertEqual(
            list(qualified_prefixed["encori_inference_eligible"]),
            [True, False],
        )


if __name__ == "__main__":
    unittest.main()
