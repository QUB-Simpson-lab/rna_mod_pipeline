from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from rna_mod_pipeline.cli import (
    StagedOutputDirectory,
    StagedOutputFiles,
    project_protected_trees,
    validate_output_location,
)


def _manifest(path: Path, workflow: str, outputs: list[Path]) -> None:
    path.write_text(
        json.dumps(
            {
                "workflow": workflow,
                "parameters": {"dataset": "sample"},
                "outputs": [{"path": str(output)} for output in outputs],
            }
        )
    )


class DirectoryPublicationTests(unittest.TestCase):
    def test_overwrite_requires_the_expected_workflow_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "result"
            destination.mkdir()
            (destination / "run_manifest.json").write_text(
                json.dumps({"workflow": "another-workflow"})
            )
            with self.assertRaisesRegex(ValueError, "another-workflow"):
                StagedOutputDirectory(
                    destination,
                    overwrite=True,
                    expected_workflow="example",
                )

    def test_successful_overwrite_replaces_complete_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "result"
            destination.mkdir()
            (destination / "run_manifest.json").write_text(
                json.dumps({"workflow": "example"})
            )
            (destination / "old.txt").write_text("old")

            transaction = StagedOutputDirectory(
                destination,
                overwrite=True,
                expected_workflow="example",
            )
            (transaction.stage / "run_manifest.json").write_text(
                json.dumps({"workflow": "example"})
            )
            (transaction.stage / "new.txt").write_text("new")
            transaction.publish()

            self.assertEqual((destination / "new.txt").read_text(), "new")
            self.assertFalse((destination / "old.txt").exists())

    def test_backup_cleanup_failure_is_nonfatal_after_publication(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "result"
            destination.mkdir()
            (destination / "run_manifest.json").write_text(
                json.dumps({"workflow": "example"})
            )
            (destination / "old.txt").write_text("old")
            transaction = StagedOutputDirectory(
                destination,
                overwrite=True,
                expected_workflow="example",
            )
            (transaction.stage / "run_manifest.json").write_text(
                json.dumps({"workflow": "example"})
            )
            (transaction.stage / "new.txt").write_text("new")
            real_rmtree = shutil.rmtree

            def fail_for_backup(path, *args, **kwargs):
                if ".backup." in Path(path).name:
                    raise OSError("simulated cleanup failure")
                return real_rmtree(path, *args, **kwargs)

            with patch(
                "rna_mod_pipeline.cli.shutil.rmtree",
                side_effect=fail_for_backup,
            ):
                with self.assertWarnsRegex(RuntimeWarning, "publication succeeded"):
                    transaction.publish()

            self.assertEqual((destination / "new.txt").read_text(), "new")
            self.assertTrue(transaction.cleanup_warnings)
            backups = list(root.glob(".result.backup.*"))
            self.assertEqual(len(backups), 1)
            real_rmtree(backups[0])


class FilePublicationTests(unittest.TestCase):
    def _transaction(self, root: Path) -> tuple[StagedOutputFiles, Path, Path]:
        output = root / "result.tsv"
        manifest = root / "run_manifest.json"
        output.write_text("old data")
        _manifest(manifest, "example", [output, manifest])
        transaction = StagedOutputFiles(
            [output, manifest],
            overwrite=True,
            expected_workflow="example",
            ownership_manifest=manifest,
            project_root=root,
            expected_parameters={"dataset": "sample"},
        )
        transaction.path(output).write_text("new data")
        transaction.path(manifest).write_text("new manifest")
        return transaction, output, manifest

    def test_backup_cleanup_failure_is_nonfatal_after_publication(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            transaction, output, manifest = self._transaction(root)
            real_unlink = Path.unlink

            def fail_for_backup(path, *args, **kwargs):
                if ".backup." in path.name:
                    raise OSError("simulated cleanup failure")
                return real_unlink(path, *args, **kwargs)

            with patch.object(Path, "unlink", fail_for_backup):
                with self.assertWarnsRegex(RuntimeWarning, "publication succeeded"):
                    transaction.publish()

            self.assertEqual(output.read_text(), "new data")
            self.assertEqual(manifest.read_text(), "new manifest")
            self.assertTrue(transaction.cleanup_warnings)
            for backup in root.glob(".*.backup.*"):
                real_unlink(backup)

    def test_publication_failure_restores_every_old_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            transaction, output, manifest = self._transaction(root)
            original_replace = transaction._replace

            def fail_for_second_publication(source: Path, destination: Path) -> None:
                if (
                    source.parent.name.startswith(".rna-mod-staging.")
                    and destination == manifest.resolve()
                ):
                    raise OSError("simulated publication failure")
                original_replace(source, destination)

            with patch.object(transaction, "_replace", side_effect=fail_for_second_publication):
                with self.assertRaisesRegex(OSError, "publication failure"):
                    transaction.publish()

            self.assertEqual(output.read_text(), "old data")
            self.assertIn('"workflow": "example"', manifest.read_text())
            self.assertFalse(list(root.glob(".*.backup.*")))

    def test_manifest_owned_stale_file_blocks_partial_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "result.tsv"
            stale = root / "old_plot.png"
            manifest = root / "run_manifest.json"
            output.write_text("old data")
            stale.write_text("old plot")
            _manifest(manifest, "example", [output, stale, manifest])

            with self.assertRaisesRegex(ValueError, "owns files not selected"):
                StagedOutputFiles(
                    [output, manifest],
                    overwrite=True,
                    expected_workflow="example",
                    ownership_manifest=manifest,
                    project_root=root,
                    expected_parameters={"dataset": "sample"},
                )

            self.assertEqual(output.read_text(), "old data")
            self.assertEqual(stale.read_text(), "old plot")


class OutputProtectionTests(unittest.TestCase):
    def test_actual_checkout_is_protected_independently_of_folder_name(self) -> None:
        checkout = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            protected = project_protected_trees(project)
            self.assertIn(checkout, protected)
            with self.assertRaisesRegex(ValueError, "protected tree"):
                validate_output_location(
                    checkout / "generated-output",
                    [],
                    protected_trees=protected,
                )

    def test_output_cannot_contain_an_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "input.tsv"
            source.write_text("input")
            with self.assertRaisesRegex(ValueError, "replace input"):
                validate_output_location(root, [source])

    def test_output_cannot_intersect_a_protected_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            protected = root / "data"
            with self.assertRaisesRegex(ValueError, "protected tree"):
                validate_output_location(protected / "new", [], protected_trees=[protected])


if __name__ == "__main__":
    unittest.main()
