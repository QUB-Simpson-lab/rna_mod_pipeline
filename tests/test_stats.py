from __future__ import annotations

import math
import unittest

import numpy as np

from rna_mod_pipeline.stats import bh_adjust, direction, fisher_test, safe_log2


class BenjaminiHochbergTests(unittest.TestCase):
    def test_adjustment_preserves_order_ties_and_missing_values(self) -> None:
        observed = bh_adjust([0.01, 0.04, 0.03, 0.002, math.nan])
        np.testing.assert_allclose(observed[:4], [0.02, 0.04, 0.04, 0.008])
        self.assertTrue(math.isnan(observed[4]))

    def test_invalid_probability_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "between 0 and 1"):
            bh_adjust([0.1, 1.01])


class FisherTests(unittest.TestCase):
    def test_table_uses_bound_and_unbound_counts(self) -> None:
        result = fisher_test(2, 10, 3, 10)
        self.assertAlmostEqual(result.odds_ratio, 14 / 24)
        self.assertGreaterEqual(result.p_value, 0)
        self.assertLessEqual(result.p_value, 1)

    def test_exact_zero_and_infinite_boundaries_are_retained(self) -> None:
        self.assertEqual(fisher_test(0, 10, 5, 10).odds_ratio, 0)
        self.assertTrue(math.isinf(fisher_test(5, 10, 0, 10).odds_ratio))

    def test_empty_comparison_returns_unavailable(self) -> None:
        result = fisher_test(0, 0, 0, 10)
        self.assertTrue(math.isnan(result.odds_ratio))
        self.assertTrue(math.isnan(result.p_value))

    def test_inconsistent_counts_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Bound counts"):
            fisher_test(11, 10, 0, 10)


class ClassificationTests(unittest.TestCase):
    def test_fdr_boundary_is_not_significant(self) -> None:
        self.assertEqual(direction(2.0, 0.05), "Not significant")

    def test_finite_significant_effects_are_directional(self) -> None:
        self.assertEqual(direction(2.0, 0.049), "Enriched")
        self.assertEqual(direction(0.5, 0.049), "Depleted")

    def test_log_transform_rejects_zero_and_nonfinite_values(self) -> None:
        self.assertTrue(math.isnan(safe_log2(0)))
        self.assertTrue(math.isnan(safe_log2(math.inf)))
        self.assertEqual(safe_log2(4), 2)


if __name__ == "__main__":
    unittest.main()
