from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from rna_mod_pipeline.binding.loose import _estimate_reason, run_loose_overlap
from rna_mod_pipeline.binding.plots import plot_all_rbps
from rna_mod_pipeline.binding.sources import BindingTrack
from rna_mod_pipeline.intervals import build_interval_index


class LooseDegenerateTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="loose_degenerate_")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        cases, controls = [10, 20, 30, 40], [110, 120, 130, 140]
        sites = pd.DataFrame({
            "chrom": "chr1", "start": cases, "strand": "+", "mod_code": "a",
            "Nvalid_cov": 20, "fraction_modified": 30.0,
        })
        sites_path = self.root / "sites.tsv"
        sites.to_csv(sites_path, sep="\t", index=False)
        raw = [
            ["chr1", pos, pos + 1, "a", 0, "+", pos, pos + 1, "0", 20,
             30.0 if pos in cases else 0.0]
            for pos in cases + controls
        ]
        bedmethyl = self.root / "calls.bedmethyl"
        pd.DataFrame(raw).to_csv(bedmethyl, sep="\t", index=False, header=False)
        positions = {
            "ALL": cases + controls, "NONE": [1_000], "EMPTY": [],
            "SAT_CONTROL": cases[:1] + controls,
            "SAT_CASE": cases + controls[:1],
            "ZERO_CASE": controls[:1], "ZERO_CONTROL": cases[:1],
            "FINITE": cases[:3] + controls[:1],
        }
        tracks = [
            BindingTrack(
                rbp_name=name, canonical_rbp=name, source_names=(name,),
                intervals=build_interval_index({"chr1": [(pos, pos + 1) for pos in bound]}),
                metadata={"resource_status": "resource_empty" if name == "EMPTY" else "loaded"},
            )
            for name, bound in positions.items()
        ]
        self.result = run_loose_overlap(
            sites_path, bedmethyl, tracks, "m6a", "ornament", window=0,
        )
        self.rows = self.result.enrichment.set_index("RBP")

    def test_all_bound_none_bound_and_empty_are_distinguished(self) -> None:
        for name, count, reason in (
            ("ALL", 4, "all_case_and_control_sites_bound"),
            ("NONE", 0, "no_case_or_control_overlap"),
            ("EMPTY", 0, "resource_empty"),
        ):
            with self.subTest(name=name):
                row = self.rows.loc[name]
                self.assertTrue(np.isnan(row.odds_ratio))
                self.assertEqual(row.mod_overlaps, count)
                self.assertEqual(row.bg_overlaps, count)
                self.assertEqual(row.estimate_status, "non_estimable")
                self.assertEqual(row.estimate_reason, reason)
                self.assertEqual(row.p_value, 1.0)
                self.assertEqual(row.fdr, 1.0)

    def test_boundary_reasons_distinguish_saturation_from_zero_overlap(self) -> None:
        expected = {
            "SAT_CONTROL": "all_control_sites_bound",
            "SAT_CASE": "all_case_sites_bound",
            "ZERO_CASE": "zero_case_overlap",
            "ZERO_CONTROL": "zero_control_overlap",
        }
        for name, reason in expected.items():
            with self.subTest(name=name):
                self.assertEqual(self.rows.loc[name, "estimate_reason"], reason)
                self.assertEqual(self.rows.loc[name, "estimate_status"], "boundary")

    def test_numerical_results_and_annotations_are_preserved(self) -> None:
        expected = {
            "SAT_CONTROL": (1, 4, 0.0, 1 / 7, 4 / 7),
            "SAT_CASE": (4, 1, np.inf, 1 / 7, 4 / 7),
            "ZERO_CASE": (0, 1, 0.0, 1.0, 1.0),
            "ZERO_CONTROL": (1, 0, np.inf, 1.0, 1.0),
            "FINITE": (3, 1, 9.0, 17 / 35, 1.0),
        }
        self.assertEqual(self.result.case_count, 4)
        self.assertEqual(self.result.background_count, 4)
        for name, (case, control, odds, p_value, fdr) in expected.items():
            with self.subTest(name=name):
                row = self.rows.loc[name]
                self.assertEqual(row.mod_overlaps, case)
                self.assertEqual(row.bg_overlaps, control)
                self.assertEqual(row.odds_ratio, odds)
                self.assertAlmostEqual(row.p_value, p_value)
                self.assertAlmostEqual(row.fdr, fdr)
        annotated = self.result.annotated_sites
        self.assertTrue(annotated["ornament_ALL"].eq(1).all())
        self.assertTrue(annotated[["ornament_NONE", "ornament_EMPTY"]].eq(0).all().all())

    def footer(self, names, stale_reason=False) -> str:
        frame = self.rows.loc[names].reset_index()
        if stale_reason:
            frame["estimate_reason"] = "no_case_or_control_overlap"
        texts = []
        with patch.object(Figure, "savefig", lambda figure, *args, **kwargs:
                          texts.extend(text.get_text() for text in figure.texts)):
            plot_all_rbps(frame, "m6A", "ornament", self.root / "plots")
        return " ".join(texts)

    def test_footer_explains_each_non_estimable_reason(self) -> None:
        all_bound = self.footer(["ALL"])
        self.assertIn("all case and comparison sites bound; no unbound sites", all_bound)
        self.assertNotIn("zero overlaps", all_bound)
        self.assertIn("zero overlaps in both groups", self.footer(["NONE"]))
        self.assertIn("empty binding resource", self.footer(["EMPTY"]))

    def test_mixed_and_stale_reasons_do_not_claim_zero_overlap(self) -> None:
        for footer in (self.footer(["ALL", "NONE"]), self.footer(["ALL"], stale_reason=True)):
            self.assertNotIn("zero overlaps in both groups", footer)
            self.assertIn("see estimate_reason", footer)

    def test_other_nan_degeneracy_has_an_explicit_fallback(self) -> None:
        self.assertEqual(
            _estimate_reason(np.nan, 2, 4, 1, 4, "loaded"),
            "degenerate_contingency_table",
        )


if __name__ == "__main__":
    unittest.main()
