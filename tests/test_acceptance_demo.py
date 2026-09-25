from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from pyfaidx import Fasta

from rna_mod_pipeline.acceptance.checks import analytical_reference, fisher_reference
from rna_mod_pipeline.acceptance.generate import CHROMOSOME, COUNTS, prepare_demo
from rna_mod_pipeline.acceptance.workflow import demo_commands
from rna_mod_pipeline.binding.catalog import load_catalog
from rna_mod_pipeline.binding.sources import binding_source_files
from rna_mod_pipeline.gui.datasets import DatasetRegistry
from rna_mod_pipeline.gui.resources import ResourceProfileStore
from rna_mod_pipeline.phase1.drach import reverse_complement
from rna_mod_pipeline.transcript.transcripts import load_complete_basic_transcripts


class AcceptanceFixtureTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="rna_demo_test_")
        self.addCleanup(self.temporary.cleanup)
        self.root = prepare_demo(Path(self.temporary.name) / "demo with spaces")

    def test_complete_transcripts_and_strand_correct_reference(self) -> None:
        models = load_complete_basic_transcripts(self.root / "resources/synthetic.gtf")
        self.assertEqual(len(models), 32)
        self.assertEqual({model.total_length for model in models.values()}, {1200})
        self.assertEqual({model.strand for model in models.values()}, {"+", "-"})
        with Fasta(str(self.root / "resources/synthetic.fa")) as genome:
            for gene in (0, 15, 16, 31):
                position = 100 + gene * 1600 + 25
                sequence = str(genome[CHROMOSOME][position - 2:position + 3])
                if gene >= 16:
                    sequence = reverse_complement(sequence)
                self.assertEqual(sequence, "GGACT")

    def test_all_declared_resources_exist_without_empty_resource_override(self) -> None:
        profile = ResourceProfileStore(self.root).load()
        self.assertEqual(profile.validation_errors(), ())
        for database, expected in COUNTS.items():
            entries = load_catalog(profile.rbp_catalog, database)
            self.assertEqual(len(entries), expected)
            source = profile.postar3 if database == "postar3" else self.root / "resources" / database
            self.assertTrue(all(path.stat().st_size > 0
                                for path in binding_source_files(database, source, entries)))

    def test_portable_profiles_build_commands_for_external_workspace(self) -> None:
        registry = DatasetRegistry.load(self.root)
        self.assertEqual(len(registry.list()), 2)
        for profile in registry.list():
            self.assertTrue(profile.bedmethyl.is_relative_to(self.root))
            self.assertTrue(profile.output_root.is_relative_to(self.root))
        commands = demo_commands(self.root, plots=False)
        self.assertEqual(len(commands), 10)
        for _name, command in commands:
            self.assertIn(str(self.root), command)
            self.assertNotIn("--allow-site-subset", command)
            self.assertNotIn("--allow-empty-resources", command)

    def test_inputs_are_deterministic_and_existing_work_is_protected(self) -> None:
        with self.assertRaises(FileExistsError):
            prepare_demo(self.root)
        other = prepare_demo(Path(self.temporary.name) / "second")
        first = json.loads((self.root / "SYNTHETIC_FIXTURE.json").read_text())
        second = json.loads((other / "SYNTHETIC_FIXTURE.json").read_text())
        self.assertEqual(first, second)

    def test_analytical_reference_has_known_exact_test_and_finite_cluster_variance(self) -> None:
        self.assertAlmostEqual(fisher_reference(2, 0, 2), 1 / 3)
        reference = analytical_reference()
        enriched = reference["DEMO_ENRICHED"]
        self.assertAlmostEqual(enriched["cluster_se"], math.sqrt(8 / 279))
        self.assertGreater(enriched["cluster_p"], 0)
        self.assertLess(enriched["cluster_p"], 1e-10)
        self.assertEqual(reference["DEMO_NULL"]["fisher_p"], 1)


if __name__ == "__main__":
    unittest.main()
