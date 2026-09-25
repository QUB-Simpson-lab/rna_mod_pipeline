from __future__ import annotations

import json
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rna_mod_pipeline.gui.command_builder import build_command
from rna_mod_pipeline.gui.comparison_profiles import prepare_comparison_profiles
from rna_mod_pipeline.gui.datasets import DatasetProfile, DatasetRegistry
from rna_mod_pipeline.gui.dataset_dialog import resolve_dialog_path
from rna_mod_pipeline.gui.preflight import validate_request
from rna_mod_pipeline.gui.resources import (
    ResourceProfile,
    ResourceProfileStore,
    fasta_index_checks,
)
from rna_mod_pipeline.gui.state import format_argv
from rna_mod_pipeline.gui.workflow_registry import prerequisites_for


CODE_ROOT = Path(__file__).resolve().parents[1]


def touch(path: Path, text: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


class ResourceProfileTests(unittest.TestCase):
    def test_round_trip_is_versioned_and_project_paths_are_relative(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = touch(root / "resources/reference.fa", ">chr1\nA\n")
            touch(Path(f"{fasta}.fai"), "chr1\t1\t6\t1\t2\n")
            gtf = touch(root / "resources/annotation.gtf")
            ornament = root / "resources/ornament"
            ornament.mkdir()
            external = touch(root.parent / f"{root.name}_postar3.txt")
            profile = ResourceProfile(
                reference_fasta=fasta,
                annotation_gtf=gtf,
                ornament_dir=ornament,
                postar3=external,
                genome_build="hg38",
                annotation_release="GENCODE v44",
            )
            store = ResourceProfileStore(root)
            path = store.save(profile)
            payload = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema_version"], 1)
            self.assertEqual(
                payload["resources"]["reference_fasta"],
                "resources/reference.fa",
            )
            self.assertTrue(Path(payload["resources"]["postar3"]).is_absolute())
            loaded = store.load()
            self.assertEqual(loaded.reference_fasta, fasta.resolve())
            self.assertEqual(loaded.postar3, external.resolve())
            external.unlink()

    def test_fasta_index_check_is_read_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fasta = touch(Path(temporary) / "reference.fa", ">chr1\nA\n")
            errors, warnings = fasta_index_checks(fasta)
            self.assertFalse(errors)
            self.assertTrue(any("will create" in item for item in warnings))
            self.assertFalse(Path(f"{fasta}.fai").exists())


class DatasetRegistryTests(unittest.TestCase):
    def test_portable_save_replace_and_non_destructive_remove(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bed = touch(root / "inputs/calls.bedmethyl")
            fasta = touch(root / "resources/ref.fa")
            gtf = touch(root / "resources/genes.gtf")
            registry = DatasetRegistry(root)
            added = registry.add(
                DatasetProfile("sample", "m6a", bed, fasta=fasta, gtf=gtf)
            )
            registry.save()
            payload = json.loads(registry.path.read_text(encoding="utf-8"))
            self.assertEqual(
                payload["datasets"][0]["bedmethyl"], "inputs/calls.bedmethyl"
            )
            replacement = registry.replace(
                "sample",
                DatasetProfile(
                    "sample-renamed",
                    "m6a",
                    bed,
                    fasta=fasta,
                    gtf=gtf,
                ),
            )
            marker = touch(replacement.output_root / "keep.txt")
            removed = registry.remove("sample-renamed")
            self.assertEqual(removed.name, "sample-renamed")
            self.assertTrue(marker.is_file())
            self.assertEqual(registry.list(), ())


class CommandBuilderTests(unittest.TestCase):
    def _dataset(self, root: Path) -> DatasetProfile:
        return DatasetProfile(
            "sample",
            "m6a",
            touch(root / "calls.bedmethyl"),
            fasta=touch(root / "ref.fa"),
            gtf=touch(root / "genes.gtf"),
        )

    def test_resource_profile_and_explicit_override_precedence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dataset = self._dataset(root)
            encori = root / "resources/encori"
            encori.mkdir(parents=True)
            resources = ResourceProfile(encori_dir=encori)
            command = build_command(
                dataset,
                "loose_overlap",
                {"database": "encori"},
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=resources,
            )
            self.assertEqual(
                command[command.index("--binding-source") + 1],
                str(encori.resolve()),
            )
            explicit = root / "override"
            command = build_command(
                dataset,
                "loose_overlap",
                {"database": "encori", "binding_source": explicit},
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=ResourceProfile(),
            )
            self.assertEqual(command[command.index("--binding-source") + 1], str(explicit))

    def test_profile_missing_selected_resource_fails_but_legacy_default_remains(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dataset = self._dataset(root)
            with self.assertRaisesRegex(ValueError, "Configure ENCORI"):
                build_command(
                    dataset,
                    "loose_overlap",
                    {"database": "encori"},
                    refactored_root=CODE_ROOT,
                    project_root=root,
                    resource_profile=ResourceProfile(),
                )
            command = build_command(
                dataset,
                "loose_overlap",
                {"database": "encori"},
                refactored_root=CODE_ROOT,
                project_root=root,
            )
            self.assertEqual(
                command[command.index("--binding-source") + 1],
                str(root.resolve() / "data/ENCORI"),
            )

    def test_transcript_builder_requires_only_selected_database_resources(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dataset = self._dataset(root)
            ornament = root / "ornament"
            ornament.mkdir()
            command = build_command(
                dataset,
                "transcript_region_overlap",
                {"databases": ["ornament"]},
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=ResourceProfile(ornament_dir=ornament),
            )
            self.assertIn("--ornament-dir", command)
            self.assertNotIn("--encori-dir", command)
            self.assertNotIn("--postar3", command)

    def test_empty_optional_values_do_not_create_dangling_flags(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dataset = self._dataset(root)
            postar3 = touch(root / "resources/human.txt")
            command = build_command(
                dataset,
                "loose_overlap",
                {
                    "database": "postar3",
                    "cell_types": [],
                    "methods": [],
                },
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=ResourceProfile(postar3=postar3),
            )
            self.assertNotIn("--cell-types", command)
            self.assertNotIn("--methods", command)

    def test_downstream_resources_and_three_modification_comparison(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dataset = self._dataset(root)
            expression = touch(root / "resources/expression.tsv")
            knock_dir = root / "resources/knockrbp"
            knock_dir.mkdir(parents=True)
            metadata = touch(knock_dir / "metadata.tsv")
            resources = ResourceProfile(
                nanopore_expression=expression,
                knockrbp_dir=knock_dir,
                knockrbp_metadata=metadata,
            )
            expression_command = build_command(
                dataset,
                "integrate_expression",
                {"source": "nanopore", "database": "encori"},
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=resources,
            )
            self.assertEqual(
                expression_command[expression_command.index("--expression") + 1],
                str(expression.resolve()),
            )
            knock_command = build_command(
                dataset,
                "knockrbp_validation",
                {},
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=resources,
            )
            self.assertEqual(
                knock_command[knock_command.index("--knockrbp-dir") + 1],
                str(knock_dir.resolve()),
            )
            network_command = build_command(
                dataset,
                "knockrbp_validation",
                {
                    "include_regulatory_network": True,
                    "regulatory_plot_top_n": 25,
                },
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=resources,
            )
            expected_cross = (
                dataset.resolved(root).output_root
                / "cross_database/loose_all/cross_database_results.tsv"
            )
            self.assertEqual(
                network_command[network_command.index("--cross-database") + 1],
                str(expected_cross),
            )
            self.assertEqual(
                network_command[
                    network_command.index("--regulatory-plot-top-n") + 1
                ],
                "25",
            )
            other_profiles = [
                DatasetProfile(
                    modification,
                    modification,
                    dataset.bedmethyl,
                    fasta=dataset.fasta,
                    gtf=dataset.gtf,
                )
                for modification in ("m5c", "pseu")
            ]
            comparison = build_command(
                dataset,
                "compare_modifications",
                {"comparison_profiles": other_profiles},
                refactored_root=CODE_ROOT,
                project_root=root,
            )
            start = comparison.index("--modifications") + 1
            self.assertEqual(comparison[start:start + 3], ["m6a", "m5c", "pseu"])

    def test_reference_profile_populates_unset_dataset_reference(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bed = touch(root / "calls.bedmethyl")
            fasta = touch(root / "portable/reference.fa")
            touch(Path(f"{fasta}.fai"))
            profile = DatasetProfile("sample", "m6a", bed)
            command = build_command(
                profile,
                "phase1_drach",
                {},
                refactored_root=CODE_ROOT,
                project_root=root,
                resource_profile=ResourceProfile(reference_fasta=fasta),
            )
            self.assertEqual(
                command[command.index("--fasta") + 1],
                str(fasta.resolve()),
            )


class PreflightTests(unittest.TestCase):
    def test_dataset_comparison_validation_does_not_write_profiles(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bed = touch(root / "calls.bedmethyl")
            fasta = touch(root / "ref.fa")
            gtf = touch(root / "genes.gtf")
            profiles = []
            for name in ("first", "second"):
                output = root / "outputs" / name
                profile = DatasetProfile(
                    name,
                    "m6a",
                    bed,
                    fasta=fasta,
                    gtf=gtf,
                    output_root=output,
                    sample_metadata={
                        "genome_build": "hg38",
                        "annotation_release": "GENCODE v44",
                    },
                ).resolved(root)
                touch(output / "phase1/filtered_m6A.tsv")
                touch(output / "phase1/filtered_m6A_metagene.tsv")
                profiles.append(profile)
            command = build_command(
                profiles[0],
                "compare_datasets",
                {"second_dataset": profiles[1]},
                refactored_root=CODE_ROOT,
                project_root=root,
            )
            result = validate_request(
                "compare_datasets",
                profiles[0],
                {"second_dataset": profiles[1]},
                command,
            )
            self.assertTrue(result.ok, result.errors)
            targets = [
                item.output_root / "dataset_comparison_profile.json"
                for item in profiles
            ]
            self.assertFalse(any(path.exists() for path in targets))
            prepare_comparison_profiles(profiles[0], profiles[1], root)
            self.assertTrue(all(path.is_file() for path in targets))

    def test_drach_preflight_reports_missing_index_without_creating_it(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bed = touch(root / "calls.bedmethyl")
            fasta = touch(root / "reference.fa", ">chr1\nA\n")
            gtf = touch(root / "genes.gtf")
            profile = DatasetProfile(
                "sample", "m6a", bed, fasta=fasta, gtf=gtf
            ).resolved(root)
            touch(profile.output_root / "phase1/filtered_m6A.tsv")
            command = build_command(
                profile,
                "phase1_drach",
                {},
                refactored_root=CODE_ROOT,
                project_root=root,
            )
            result = validate_request("phase1_drach", profile, {}, command)
            self.assertTrue(result.ok, result.errors)
            self.assertTrue(any("pyfaidx will create" in item for item in result.warnings))
            self.assertFalse(Path(f"{fasta}.fai").exists())


class InterfaceHelpersTests(unittest.TestCase):
    def test_transcript_cli_catalog_default_is_checkout_relative(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            checkout = Path(temporary) / "rna_mod_pipeline"
            (checkout / "scripts").mkdir(parents=True)
            (checkout / "config").mkdir()
            shutil.copy2(
                CODE_ROOT / "scripts/transcript_region_overlap.py",
                checkout / "scripts/transcript_region_overlap.py",
            )
            shutil.copytree(CODE_ROOT / "src", checkout / "src")
            original_path = sys.path[:]
            try:
                namespace = runpy.run_path(
                    str(checkout / "scripts/transcript_region_overlap.py")
                )
            finally:
                sys.path[:] = original_path
            args = namespace["parser"]().parse_args(
                ["--project-root", str(Path(temporary) / "workspace")]
            )
            self.assertEqual(
                Path(args.rbp_catalog),
                (checkout / "config/rbp_catalog.tsv").resolve(),
            )

    def test_dataset_dialog_paths_are_relative_to_project_not_process_cwd(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            expected = root / "inputs/calls.bedmethyl"
            self.assertEqual(
                resolve_dialog_path(root, "inputs/calls.bedmethyl"),
                expected.resolve(),
            )

    def test_dynamic_prerequisites(self) -> None:
        self.assertEqual(
            prerequisites_for("phase1_metagene", "m6a"),
            ("phase1_drach",),
        )
        self.assertEqual(
            prerequisites_for("phase1_metagene", "m5c"),
            ("phase1_filter",),
        )
        self.assertEqual(
            prerequisites_for(
                "knockrbp_validation",
                "m6a",
                {"design": "transcript-region"},
            ),
            ("transcript_region_overlap",),
        )
        self.assertEqual(
            prerequisites_for(
                "knockrbp_validation",
                "m6a",
                {"design": "loose", "include_regulatory_network": True},
            ),
            ("cross_database",),
        )

    def test_loose_downstream_rejects_regional_context(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bed = touch(root / "calls.bedmethyl")
            fasta = touch(root / "ref.fa")
            gtf = touch(root / "genes.gtf")
            dataset = DatasetProfile(
                "sample", "m6a", bed, fasta=fasta, gtf=gtf
            )
            with self.assertRaisesRegex(ValueError, "only all context"):
                build_command(
                    dataset,
                    "cross_database",
                    {"design": "loose", "context": "CDS"},
                    refactored_root=CODE_ROOT,
                    project_root=root,
                )

    def test_command_preview_is_platform_aware(self) -> None:
        argv = ["python", "folder with spaces/script.py", "plain"]
        self.assertEqual(format_argv(argv, "nt"), subprocess.list2cmdline(argv))
        self.assertEqual(
            format_argv(argv, "posix"),
            "python 'folder with spaces/script.py' plain",
        )


if __name__ == "__main__":
    unittest.main()
