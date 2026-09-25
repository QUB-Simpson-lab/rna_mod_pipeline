from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from rna_mod_pipeline.doctor import run_resource_doctor
from rna_mod_pipeline.gui.resources import ResourceProfile, ResourceProfileStore
from rna_mod_pipeline.provenance import make_manifest, validate_manifest


CODE_ROOT = Path(__file__).resolve().parents[1]


def touch(path: Path, text: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


class ResourceDoctorTests(unittest.TestCase):
    def test_store_discovers_catalog_independently_of_checkout_name(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            profile = ResourceProfileStore(root).load()
            self.assertEqual(
                profile.rbp_catalog,
                (CODE_ROOT / "config/rbp_catalog.tsv").resolve(),
            )

    def test_portable_doctor_uses_profile_without_legacy_datasets(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = touch(root / "resources/reference.fa", ">chr1\nA\n")
            touch(Path(f"{fasta}.fai"), "chr1\t1\t6\t1\t2\n")
            gtf = touch(root / "resources/annotation.gtf")
            store = ResourceProfileStore(root)
            store.save(
                ResourceProfile(
                    reference_fasta=fasta,
                    annotation_gtf=gtf,
                    rbp_catalog=CODE_ROOT / "config/rbp_catalog.tsv",
                )
            )
            checks = run_resource_doctor(root)
            self.assertFalse([check for check in checks if check.failed])
            self.assertTrue(
                any(
                    check.check == "RBP catalogue panels"
                    and "ornament=133" in check.details
                    for check in checks
                )
            )
            self.assertTrue(any(check.status == "SKIP" for check in checks))


class ProvenanceTests(unittest.TestCase):
    def test_manifest_records_constraints_and_extended_versions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = touch(root / "input.tsv")
            output = touch(root / "output.tsv")
            manifest = make_manifest(
                "test",
                root,
                [source],
                {"threshold": 20},
                [output],
            )
            self.assertIn("matplotlib", manifest["software"])
            self.assertIn("networkx", manifest["software"])
            self.assertIn("pyfaidx", manifest["software"])
            self.assertIn("code_git_commit", manifest["software"])
            self.assertEqual(len(manifest["dependency_constraints"]), 3)
            path = root / "manifest.json"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            self.assertEqual(validate_manifest(path, root), [])


if __name__ == "__main__":
    unittest.main()
