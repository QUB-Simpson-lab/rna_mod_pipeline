from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1] / "scripts/verify_resources.py"
SPEC = importlib.util.spec_from_file_location("resource_inventory", SOURCE)
assert SPEC is not None and SPEC.loader is not None
inventory = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inventory)


class ResourceSnapshotTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="resource_inventory_test_")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.catalog = self.root / "catalog.tsv"
        self.catalog.write_text(
            "catalog_version\tdatabase\trbp_name\tcanonical_rbp\tsource_names\n"
            "1\tornament\tSRSF1\tSRSF1\tSRSF1;SF2\n"
            "1\tencori\tYTHDF2\tYTHDF2\tYTHDF2\n",
            encoding="utf-8",
        )
        names = (
            "hg38.fa", "hg38.fa.fai", "gencode.v44.annotation.gtf", "human.txt",
            "HS/SRSF1.bed.gz", "HS/SRSF1.bed", "HS/SF2.bed.gz",
            "ENCORI/YTHDF2.bed", "knockrbp/YTHDF2_DataSet_8_degs.json",
            "Kate_231_0h_vs_6h_gene_counts_normalised.tsv",
            "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
        )
        for name in names:
            path = self.root / "data" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"abc")
        (self.root / "data/knockrbp/dataset_metadata.tsv").write_text(
            "dataset_id\trbp\tcell_line\nDataSet_8\tYTHDF2\tMDA-MB-231-LM2\n",
            encoding="utf-8",
        )
        self.snapshot = self.root / "snapshot.json"
        self.base = ["--workspace", str(self.root), "--catalog", str(self.catalog)]
        self.assertEqual(self.run_main(["--create", "--output", str(self.snapshot)]), 0)
        self.base += ["--snapshot", str(self.snapshot)]

    def run_main(self, args=()) -> int:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return inventory.main(self.base + list(args))

    def test_create_verify_and_aliases(self) -> None:
        self.assertEqual(self.run_main(), 0)
        document = json.loads(self.snapshot.read_text(encoding="utf-8"))
        beds = [row for row in document["files"] if row["group"] == "ornament"]
        self.assertEqual(len(beds), 2)
        self.assertTrue(all(row["path"].endswith(".bed.gz") for row in beds))
        self.assertTrue(all(row["canonical_rbps"] == ["SRSF1"] for row in beds))

    def test_subset_ignores_unrequested_missing(self) -> None:
        (self.root / "data/hg38.fa").unlink()
        self.assertEqual(self.run_main(["--groups", "ornament"]), 0)
        self.assertEqual(self.run_main(), 1)

    def test_same_size_corruption_detected(self) -> None:
        (self.root / "data/human.txt").write_bytes(b"xyz")
        self.assertEqual(self.run_main(), 1)

    def test_size_corruption_detected(self) -> None:
        (self.root / "data/human.txt").write_bytes(b"x")
        self.assertEqual(self.run_main(), 1)

    def test_catalog_change_detected(self) -> None:
        self.catalog.write_text(self.catalog.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        self.assertEqual(self.run_main(), 1)

    def test_remapping_directory(self) -> None:
        destination = self.root / "Elsewhere with spaces"
        (self.root / "data/HS").rename(destination)
        self.assertEqual(self.run_main(["--map", f"data/HS={destination}"]), 0)

    def test_remapping_exact_file(self) -> None:
        destination = self.root / "Another name.fa"
        (self.root / "data/hg38.fa").rename(destination)
        self.assertEqual(self.run_main(["--map", f"data/hg38.fa={destination}"]), 0)

    def test_reject_unused_mapping(self) -> None:
        self.assertEqual(self.run_main(["--map", f"data/typo={self.root}"]), 2)

    def test_reject_duplicate_mapping(self) -> None:
        self.assertEqual(self.run_main([
            "--map", f"data/HS={self.root}", "--map", f"data/HS={self.root}"
        ]), 2)

    def test_refuse_snapshot_overwrite(self) -> None:
        before = self.snapshot.read_bytes()
        self.assertEqual(self.run_main(["--create", "--output", str(self.snapshot)]), 2)
        self.assertEqual(self.snapshot.read_bytes(), before)

    def test_reject_unsafe_snapshot_path(self) -> None:
        document = json.loads(self.snapshot.read_text(encoding="utf-8"))
        document["files"][0]["path"] = "../unsafe"
        self.snapshot.write_text(json.dumps(document), encoding="utf-8")
        self.assertEqual(self.run_main(), 2)

    def test_reject_invalid_totals(self) -> None:
        document = json.loads(self.snapshot.read_text(encoding="utf-8"))
        document["total_bytes"] += 1
        self.snapshot.write_text(json.dumps(document), encoding="utf-8")
        self.assertEqual(self.run_main(), 2)

    def test_reject_unregistered_knockrbp(self) -> None:
        (self.root / "data/knockrbp/X_DataSet_20_degs.json").write_text("{}", encoding="utf-8")
        self.assertEqual(self.run_main([
            "--create", "--output", str(self.root / "other.json")
        ]), 2)


if __name__ == "__main__":
    unittest.main()
