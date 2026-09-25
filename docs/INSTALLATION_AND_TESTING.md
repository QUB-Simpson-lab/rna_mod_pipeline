# Installation, constraints, and testing

## Supported environment

The supported Python range is 3.10 through 3.12 on current Windows, macOS, and
Linux runners. Python 3.11 is recommended. Python 3.13 and later are rejected
until the scientific and GUI dependency set is tested there.

`requirements-lock.txt` constrains the complete observed core runtime
dependency closure, and `requirements-gui-lock.txt` constrains the optional
PySide6 family. `requirements-build-lock.txt` pins installer/build tools. These
are tested constraints without hashes: they make selected versions repeatable,
but they are not a cryptographically locked wheel set and do not guarantee that
an upstream index will retain every platform wheel forever.

For a manual constrained installation:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --no-cache-dir -r requirements-build-lock.txt
python -m pip install --no-cache-dir --no-build-isolation \
  -c requirements-lock.txt -c requirements-gui-lock.txt -e ".[gui]"
python -m pip check
```

The double-click installers perform the same constrained installation. Python
bytecode created by their source check is redirected under `.venv`, and normal
runs disable bytecode writes, so package-source folders are not populated with
new cache files.

On Windows, keep the source checkout on a local drive. The `.bat` helpers
require their own folder as the working directory and stop with an error if
they cannot enter it, before installing, running tests, launching the GUI or
resetting an environment. UNC paths such as `\\server\share` are not supported
as source locations; copy the source to a local folder first. Paths containing
spaces are supported.

## Tests

The current two-laptop procedure starts at [`START_HERE.md`](../START_HERE.md).
In addition to the regression suite, a self-contained synthetic example has
known numerical answers and exercises the real command-line workflows. It
generates its resources outside the checkout; it does not contain biological
data. `test_demo_macos.command` and `test_demo_windows.bat` run both checks.

The compact regression suite exercises:

- Fisher table construction and zero/infinite odds-ratio boundaries;
- Benjamini-Hochberg adjustment, missing values, and the FDR boundary;
- protected output paths and manifest ownership;
- atomic directory/file publication;
- successful publication with non-fatal, auditable backup-cleanup failure;
- rollback when publication itself fails.

Run it with:

```bash
python -m unittest discover -s tests -v
```

Tests use synthetic temporary files only. They do not modify scientific inputs,
legacy results, or `refactored_outputs/`.

## Continuous integration

`.github/workflows/ci.yml` tests Python 3.10 and 3.12 on Linux and Python 3.11
on Windows and macOS, checks dependencies, exercises installed commands, and
builds source/wheel packages. The wheel build is a packaging smoke check, not a
standalone desktop distribution: current workflow dispatch deliberately uses
the checkout's top-level `scripts/` and `config/rbp_catalog.tsv`. The supported
distribution therefore remains an editable installation from the full
source checkout. A future standalone wheel or executable must first move those
resources behind package APIs and test every dispatched workflow.

The supported GitHub distribution makes the contents of `refactored_code/`
the repository root, so this workflow is discovered with the paths and working
directory exactly as tested. If `refactored_code/` remains nested inside the
larger research-project repository, GitHub will not discover this nested
workflow. Do not merely copy it upward: a parent-repository workflow must also
set `working-directory: refactored_code` and update dependency-cache paths.

`rna-mod-audit --skip-legacy-results` skips established result baselines but
still assumes the historical project data/four-callset layout. It is retained
for backward compatibility. `rna-mod-doctor --project-root ..` is the portable,
read-only diagnostic for a clean laptop: it checks the saved or auto-discovered
resource profile, workflow entry points, writable workspace, FASTA index,
catalogue panel sizes, binding-source mappings, POSTAR3 schema, and KnockRBP
metadata/file consistency without requiring historical bedMethyl datasets.

## Safe environment reset

`reset_environment_windows.bat` and `reset_environment_macos.command` remove
only the repository-local `.venv`, and only when its `pyvenv.cfg` marker is present.
They never target `data/`, legacy results, `refactored_outputs/`, or arbitrary
user-selected outputs. The removed environment is regenerable with the
installer; scientific results are unchanged.

## Release decisions still required

Before public distribution:

1. review the MIT licence selected during this preparation and the final
   software contributor list with the lab;
2. add further approved names/ORCIDs and publication metadata when available;
   Alexandru Zob and the lab repository URL are already recorded;
3. validate the constrained install on clean machines for each released
   operating system and architecture;
4. keep large third-party databases and restricted data outside the package.

`LICENSE` now contains the MIT text. It does not change any third-party data
terms. See `RELEASE_READINESS.md` for the checks completed in this preparation
and `CROSS_PLATFORM_ACCEPTANCE.md` for the remaining native laptop record.
