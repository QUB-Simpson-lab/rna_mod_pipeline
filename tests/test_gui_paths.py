from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path, PureWindowsPath

from rna_mod_pipeline.gui.datasets import DatasetProfile, DatasetRegistry
from rna_mod_pipeline.gui.profile_paths import PATH_FORMAT, load_path, serialize_path
from rna_mod_pipeline.gui.resources import ResourceProfile


def touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("example\n", encoding="utf-8")
    return path


class SerializedPathTests(unittest.TestCase):
    def test_windows_helpers_stop_if_their_source_directory_is_unavailable(self) -> None:
        source = Path(__file__).resolve().parents[1]
        for name in ("install_windows.bat", "run_windows.bat", "test_demo_windows.bat",
                     "run_demo_windows.bat", "reset_environment_windows.bat"):
            with self.subTest(helper=name):
                content = (source / name).read_bytes()
                self.assertNotIn(b"\n", content.replace(b"\r\n", b""))
                text = content.decode("utf-8").replace("\r\n", "\n")
                guard = text.index('cd /d "%~dp0" || (')
                end = text.index("\n)", guard)
                self.assertIn("exit /b 1", text[guard:end])
                self.assertLess(end, text.index(".venv"))

    def test_windows_relative_serialization_uses_forward_slashes(self) -> None:
        root = PureWindowsPath(r"C:\RNA project")
        path = root / "resources" / "reference.fa"
        self.assertEqual(serialize_path(path, root), "resources/reference.fa")
        self.assertEqual(
            serialize_path(PureWindowsPath(r"inputs\sample.bedmethyl"), None),
            "inputs/sample.bedmethyl",
        )

    def test_external_windows_paths_keep_their_absolute_location(self) -> None:
        root = PureWindowsPath(r"C:\RNA project")
        for path in (
            PureWindowsPath(r"C:\external\ref.fa"),
            PureWindowsPath(r"D:\resources\ref.fa"),
            PureWindowsPath(r"\\server\share\ref.fa"),
        ):
            self.assertEqual(serialize_path(path, root), str(path))

    def test_saved_dataset_and_resources_relocate_together(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary).resolve() / "original workspace"
            moved = Path(temporary).resolve() / "moved workspace"
            bed = touch(first / "inputs/sample.bedmethyl")
            fasta = touch(first / "resources/reference.fa")
            gtf = touch(first / "resources/genes.gtf")
            dataset = DatasetProfile("sample", "m6a", bed, fasta, gtf).resolved(first)
            record = dataset.to_dict(first)
            resources = ResourceProfile(reference_fasta=fasta, annotation_gtf=gtf)
            payload = resources.to_dict(first)
            self.assertEqual(record["path_format"], PATH_FORMAT)
            self.assertEqual(payload["path_format"], PATH_FORMAT)
            first.rename(moved)
            restored = DatasetProfile.from_dict(record, moved).resolved(moved)
            restored_resources = ResourceProfile.from_dict(payload, moved)
            self.assertEqual(restored.bedmethyl, moved / "inputs/sample.bedmethyl")
            self.assertEqual(restored.fasta, moved / "resources/reference.fa")
            self.assertEqual(restored.gtf, moved / "resources/genes.gtf")
            self.assertEqual(restored_resources.reference_fasta, restored.fasta)
            self.assertEqual(restored_resources.annotation_gtf, restored.gtf)
            self.assertEqual(
                restored.output_root, moved / "refactored_outputs/datasets/sample"
            )
            self.assertFalse(restored.output_root.exists())

    def test_unknown_path_format_is_not_silently_interpreted(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unsupported saved path format"):
            load_path("resources/ref.fa", path_format="future-format")


@unittest.skipIf(os.name == "nt", "POSIX interpretation of transferred Windows profiles")
class LegacyWindowsPathTests(unittest.TestCase):
    def test_existing_legacy_dataset_paths_migrate_without_rewriting_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            bed = touch(root / "inputs/calls.bedmethyl")
            fasta = touch(root / "resources/ref.fa")
            gtf = touch(root / "resources/genes.gtf")
            output = root / "results/sample"
            output.mkdir(parents=True)
            registry_path = root / "registry.json"
            payload = {
                "schema_version": 1,
                "datasets": [{
                    "name": "sample", "modification": "m6a",
                    "bedmethyl": r"inputs\calls.bedmethyl",
                    "fasta": r"resources\ref.fa",
                    "gtf": r"resources\genes.gtf",
                    "output_root": r"results\sample",
                }],
            }
            original = json.dumps(payload)
            registry_path.write_text(original, encoding="utf-8")
            registry = DatasetRegistry.load(root, registry_path)
            restored = registry.get("sample")
            self.assertEqual(restored.bedmethyl, bed)
            self.assertEqual(restored.fasta, fasta)
            self.assertEqual(restored.gtf, gtf)
            self.assertEqual(restored.output_root, output)
            self.assertEqual(registry_path.read_text(encoding="utf-8"), original)
            registry.save()
            saved = json.loads(registry_path.read_text(encoding="utf-8"))
            self.assertEqual(saved["datasets"][0]["bedmethyl"], "inputs/calls.bedmethyl")

    def test_legacy_resource_paths_migrate_when_resolved_files_exist(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            fasta = touch(root / "resources/ref.fa")
            payload = {
                "schema_version": 1,
                "resources": {"reference_fasta": r"resources\ref.fa"},
            }
            profile = ResourceProfile.from_dict(payload, root)
            self.assertEqual(profile.reference_fasta, fasta)

    def test_existing_literal_backslash_filename_has_priority(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            literal = touch(root / r"inputs\calls.bedmethyl")
            touch(root / "inputs/calls.bedmethyl")
            value = load_path(r"inputs\calls.bedmethyl", root)
            self.assertEqual(root / value, literal)
            self.assertEqual(load_path(str(literal), root), literal)

    def test_new_posix_profile_preserves_missing_literal_backslash_filename(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            value = r"inputs/calls\sample.bedmethyl"
            touch(root / "inputs/calls/sample.bedmethyl")
            self.assertEqual(load_path(value, root, PATH_FORMAT), Path(value))

    def test_ambiguous_missing_legacy_path_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError, "Cannot safely interpret"):
                load_path(r"results\sample", temporary)

    def test_foreign_drives_and_network_paths_are_never_prefixed_with_root(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            for value in (r"C:\data\ref.fa", "C:/data/ref.fa", r"C:ref.fa",
                          r"\\server\share\ref.fa"):
                with self.subTest(value=value):
                    with self.assertRaisesRegex(ValueError, "Reselect this Windows"):
                        load_path(value, temporary, PATH_FORMAT)


@unittest.skipUnless(os.name == "nt", "Requires native Windows path semantics")
class NativeWindowsPathTests(unittest.TestCase):
    def test_legacy_windows_and_new_relative_paths_are_accepted(self) -> None:
        self.assertEqual(load_path(r"inputs\calls.bedmethyl"), Path("inputs/calls.bedmethyl"))
        self.assertEqual(
            load_path("resources/ref.fa", path_format=PATH_FORMAT),
            Path("resources/ref.fa"),
        )
        self.assertEqual(load_path(r"D:\resources\ref.fa"), Path(r"D:\resources\ref.fa"))

    def test_foreign_posix_absolute_and_literal_backslash_paths_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Reselect this absolute path"):
            load_path("/Users/example/ref.fa", path_format=PATH_FORMAT)
        with self.assertRaisesRegex(ValueError, "literal backslash path"):
            load_path(r"inputs/calls\sample.bedmethyl", path_format=PATH_FORMAT)


if __name__ == "__main__":
    unittest.main()
