from __future__ import annotations

import json
import shutil
import tempfile
import uuid
import warnings
from collections.abc import Iterable
from pathlib import Path


LEGACY_PROJECT_TREES = (
    "data",
    "m6a",
    "m5c",
    "pseudouridine",
    "m6a_rep2",
    "m6a_frac10",
    "m6a_frac30",
    "m5c_frac10",
    "m5c_frac30",
    "pseu_frac10",
    "pseu_frac30",
    "cross_validation",
    "comparison",
    "target_gene_analysis",
    "target_gene_comparison",
    "transcript_region_stratified_overlap",
    "all_rbp_matched_screen",
    "matched_site_validation",
    "matched_site_validation_corrected",
)


def project_protected_trees(project_root: str | Path) -> tuple[Path, ...]:
    root = Path(project_root).expanduser().resolve()
    checkout = Path(__file__).resolve().parents[2]
    candidates = (
        checkout,
        root / "refactored_code",
        *(root / name for name in LEGACY_PROJECT_TREES),
    )
    return tuple(dict.fromkeys(path.resolve() for path in candidates))


def resolve_path(project_root: str | Path, value: str | Path) -> Path:
    root = Path(project_root).expanduser().resolve()
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def validate_output_location(
    output: str | Path,
    inputs: Iterable[str | Path],
    protected_roots: Iterable[str | Path] = (),
    protected_trees: Iterable[str | Path] = (),
) -> Path:
    destination = Path(output).expanduser().resolve()
    exact_protected = {
        Path.home().resolve(),
        Path(destination.anchor).resolve(),
        *(Path(item).expanduser().resolve() for item in protected_roots),
    }
    if destination in exact_protected:
        raise ValueError(f"Refusing to use protected output path: {destination}")
    for value in protected_trees:
        tree = Path(value).expanduser().resolve()
        if (
            destination == tree
            or tree in destination.parents
            or destination in tree.parents
        ):
            raise ValueError(
                f"Output path intersects protected tree {tree}: {destination}"
            )
    for value in inputs:
        source = Path(value).expanduser().resolve()
        if destination == source or destination in source.parents:
            raise ValueError(
                f"Output path {destination} would replace input {source}"
            )
        if source.is_dir() and source in destination.parents:
            raise ValueError(
                f"Output path {destination} is inside input directory {source}"
            )
    return destination


def csv_values(value: str) -> list[str]:
    values = [item.strip() for item in value.split(",") if item.strip()]
    if not values:
        raise ValueError("at least one value is required")
    return values


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise ValueError("value must be greater than zero")
    return parsed


def percentage(value: str) -> float:
    parsed = float(value)
    if parsed < 0 or parsed > 100:
        raise ValueError("percentage must be between 0 and 100")
    return parsed


def validate_workflow_directory(
    path: str | Path,
    expected_workflow: str | None,
) -> None:
    destination = Path(path).expanduser().resolve()
    if not expected_workflow:
        raise ValueError(
            "Refusing broad directory overwrite without an expected workflow"
        )
    manifest_path = destination / "run_manifest.json"
    if not manifest_path.is_file():
        raise ValueError(
            f"Refusing to overwrite an unmarked directory: {destination}"
        )
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot validate output manifest: {manifest_path}") from exc
    observed = manifest.get("workflow")
    if observed != expected_workflow:
        raise ValueError(
            f"Output directory belongs to workflow {observed!r}, "
            f"not {expected_workflow!r}: {destination}"
        )


class StagedOutputDirectory:
    def __init__(
        self,
        destination: str | Path,
        overwrite: bool,
        expected_workflow: str,
        inputs: Iterable[str | Path] = (),
        protected_roots: Iterable[str | Path] = (),
        protected_trees: Iterable[str | Path] = (),
    ):
        self.destination = validate_output_location(
            destination,
            inputs,
            protected_roots,
            protected_trees,
        )
        if self.destination.exists() and not self.destination.is_dir():
            raise FileExistsError(
                f"Output path is an existing file: {self.destination}"
            )
        nonempty = self.destination.exists() and any(self.destination.iterdir())
        if nonempty and not overwrite:
            raise FileExistsError(
                f"Output directory is not empty: {self.destination}. "
                "Choose another directory or pass --overwrite."
            )
        if nonempty and overwrite:
            validate_workflow_directory(self.destination, expected_workflow)
        self.destination.parent.mkdir(parents=True, exist_ok=True)
        self.stage = Path(
            tempfile.mkdtemp(
                prefix=f".{self.destination.name}.staging.",
                dir=self.destination.parent,
            )
        )
        self._published = False
        self._cleanup_warnings: list[str] = []

    @property
    def cleanup_warnings(self) -> tuple[str, ...]:
        """Cleanup problems recorded after a successful publication."""
        return tuple(self._cleanup_warnings)

    def _warn_cleanup(self, path: Path, exc: Exception) -> None:
        message = (
            "Output publication succeeded, but an obsolete backup could not "
            f"be removed: {path} ({type(exc).__name__}: {exc}). The new output "
            "is live; remove the backup manually after checking it."
        )
        self._cleanup_warnings.append(message)
        warnings.warn(message, RuntimeWarning, stacklevel=2)

    def publish(self) -> None:
        backup: Path | None = None
        if self.destination.exists():
            backup = self.destination.with_name(
                f".{self.destination.name}.backup.{uuid.uuid4().hex}"
            )
            self.destination.replace(backup)
        try:
            self.stage.replace(self.destination)
            self._published = True
        except Exception:
            if backup is not None and not self.destination.exists():
                backup.replace(self.destination)
            raise
        if backup is not None:
            try:
                shutil.rmtree(backup)
            except OSError as exc:
                self._warn_cleanup(backup, exc)

    def __enter__(self) -> Path:
        return self.stage

    def __exit__(self, exc_type, exc, traceback) -> None:
        if not self._published and self.stage.exists():
            shutil.rmtree(self.stage)


def _manifest_output_paths(
    manifest_path: Path,
    project_root: Path,
    expected_workflow: str,
    expected_parameters: dict[str, object],
) -> set[Path]:
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot validate output manifest: {manifest_path}") from exc
    observed = manifest.get("workflow")
    if observed != expected_workflow:
        raise ValueError(
            f"Output files belong to workflow {observed!r}, "
            f"not {expected_workflow!r}: {manifest_path}"
        )
    parameters = manifest.get("parameters", {})
    mismatches = {
        key: (parameters.get(key), expected)
        for key, expected in expected_parameters.items()
        if parameters.get(key) != expected
    }
    if mismatches:
        raise ValueError(
            "Output manifest belongs to a different run identity: "
            f"{mismatches}"
        )
    owned: set[Path] = set()
    for record in manifest.get("outputs", []):
        raw_path = record.get("path")
        if not isinstance(raw_path, str) or not raw_path:
            continue
        value = Path(raw_path)
        owned.add(
            value.expanduser().resolve()
            if value.is_absolute()
            else (project_root / value).resolve()
        )
    return owned


class StagedOutputFiles:
    """Publish a known set of files together without replacing their directories."""

    def __init__(
        self,
        destinations: Iterable[str | Path],
        overwrite: bool,
        expected_workflow: str,
        ownership_manifest: str | Path,
        project_root: str | Path,
        expected_parameters: dict[str, object] | None = None,
        inputs: Iterable[str | Path] = (),
        protected_roots: Iterable[str | Path] = (),
        protected_trees: Iterable[str | Path] = (),
    ):
        self.project_root = Path(project_root).expanduser().resolve()
        self.destinations = tuple(
            validate_output_location(
                value,
                inputs,
                protected_roots,
                protected_trees,
            )
            for value in destinations
        )
        if not self.destinations:
            raise ValueError("At least one output file is required")
        if len(set(self.destinations)) != len(self.destinations):
            raise ValueError("Output file paths must be unique")
        for destination in self.destinations:
            if destination.exists() and not destination.is_file():
                raise FileExistsError(
                    f"Output path is not a regular file: {destination}"
                )

        self.ownership_manifest = Path(ownership_manifest).expanduser().resolve()
        if self.ownership_manifest not in self.destinations:
            raise ValueError("The ownership manifest must be one of the output files")
        existing = {path for path in self.destinations if path.exists()}
        if existing and not overwrite:
            joined = "\n  ".join(str(path) for path in sorted(existing))
            raise FileExistsError(
                f"Output files already exist:\n  {joined}\n"
                "Choose other paths or pass --overwrite."
            )
        if existing and overwrite:
            if not self.ownership_manifest.is_file():
                raise ValueError(
                    "Refusing to overwrite output files without their workflow "
                    f"manifest: {self.ownership_manifest}"
                )
            owned = _manifest_output_paths(
                self.ownership_manifest,
                self.project_root,
                expected_workflow,
                expected_parameters or {},
            )
            stale_owned = {
                path
                for path in owned.difference(self.destinations)
                if path.exists()
            }
            if stale_owned:
                joined = "\n  ".join(str(path) for path in sorted(stale_owned))
                raise ValueError(
                    "The existing run owns files not selected by this rerun:\n  "
                    f"{joined}\nUse the same output/plot selection or choose "
                    "fresh output paths; no files were changed."
                )
            unowned = existing.difference(
                {self.ownership_manifest, *owned}
            )
            if unowned:
                joined = "\n  ".join(str(path) for path in sorted(unowned))
                raise ValueError(
                    f"Refusing to overwrite files not owned by the manifest:\n  {joined}"
                )

        self._transaction_id = uuid.uuid4().hex
        self._staging_dirs: dict[Path, Path] = {}
        for parent in {path.parent for path in self.destinations}:
            parent.mkdir(parents=True, exist_ok=True)
            self._staging_dirs[parent] = Path(
                tempfile.mkdtemp(
                    prefix=".rna-mod-staging.",
                    dir=parent,
                )
            )
        self._staged = {
            destination: self._staging_dirs[destination.parent] / destination.name
            for destination in self.destinations
        }
        self._backups: dict[Path, Path] = {}
        self._published = False
        self._recovery_failed = False
        self._cleanup_warnings: list[str] = []

    @property
    def cleanup_warnings(self) -> tuple[str, ...]:
        """Cleanup problems recorded after a successful publication."""
        return tuple(self._cleanup_warnings)

    def path(self, destination: str | Path) -> Path:
        resolved = Path(destination).expanduser().resolve()
        try:
            return self._staged[resolved]
        except KeyError as exc:
            raise KeyError(f"Unregistered output file: {resolved}") from exc

    def directory(self, destination: str | Path) -> Path:
        resolved = Path(destination).expanduser().resolve()
        try:
            return self._staging_dirs[resolved]
        except KeyError as exc:
            raise KeyError(f"No output files are registered in: {resolved}") from exc

    @staticmethod
    def _replace(source: Path, destination: Path) -> None:
        source.replace(destination)

    def _cleanup(self, remove_backups: bool = True) -> None:
        for directory in self._staging_dirs.values():
            if directory.exists():
                shutil.rmtree(directory)
        if remove_backups:
            for backup in self._backups.values():
                if backup.exists():
                    backup.unlink()

    def _cleanup_after_publish(self) -> None:
        paths = [*self._staging_dirs.values(), *self._backups.values()]
        for path in paths:
            try:
                if path.is_dir():
                    shutil.rmtree(path)
                elif path.exists():
                    path.unlink()
            except OSError as exc:
                message = (
                    "Output publication succeeded, but transaction cleanup "
                    f"could not remove {path} ({type(exc).__name__}: {exc}). "
                    "The new output is live; remove this transaction artifact "
                    "manually after checking it."
                )
                self._cleanup_warnings.append(message)
                warnings.warn(message, RuntimeWarning, stacklevel=2)

    def publish(self) -> None:
        missing = [
            destination
            for destination, staged in self._staged.items()
            if not staged.is_file()
        ]
        if missing:
            joined = "\n  ".join(str(path) for path in missing)
            raise FileNotFoundError(f"Staged output files are missing:\n  {joined}")

        published: list[Path] = []
        try:
            for destination in self.destinations:
                if destination.exists():
                    backup = destination.with_name(
                        f".{destination.name}.backup.{self._transaction_id}"
                    )
                    self._replace(destination, backup)
                    self._backups[destination] = backup
            for destination in self.destinations:
                self._replace(self._staged[destination], destination)
                published.append(destination)
        except Exception as publish_error:
            for destination in reversed(published):
                if destination.exists():
                    destination.unlink()
            recovery_errors: list[Exception] = []
            for destination, backup in reversed(tuple(self._backups.items())):
                try:
                    if backup.exists():
                        self._replace(backup, destination)
                except Exception as exc:
                    recovery_errors.append(exc)
            self._recovery_failed = bool(recovery_errors)
            self._cleanup(remove_backups=not self._recovery_failed)
            if recovery_errors:
                raise RuntimeError(
                    "Output publication failed and automatic restoration was "
                    "incomplete; backup files were preserved"
                ) from publish_error
            raise

        self._published = True
        self._cleanup_after_publish()

    def __enter__(self) -> "StagedOutputFiles":
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        if not self._published:
            self._cleanup(remove_backups=not self._recovery_failed)
