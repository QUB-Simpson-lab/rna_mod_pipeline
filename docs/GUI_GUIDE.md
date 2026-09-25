# Dataset-centred desktop launcher

The graphical launcher is an optional interface over the audited command-line
workflows. It does not reimplement filtering, overlap calculations, statistics,
or plotting.

For a first laptop trial, follow [`START_HERE.md`](../START_HERE.md). The
`test_demo_*` helpers run synthetic inputs and `run_demo_*` opens their separate
workspace. The ordinary `run_*` helpers continue to open the research workspace.

## What a dataset means

A dataset is one named bedMethyl input and its isolated output workspace. It is
not one of four hard-coded project slots. The launcher accepts any bedMethyl
that uses one of the currently supported chemistries:

- `m6a`;
- `m5c`;
- `pseu`.

The dataset name identifies the biological sample or replicate. The
modification value identifies the chemistry. For example,
`MCF7_m6A_rep1` and `MDA231_m6A_rep1` are two distinct datasets with the same
`m6a` chemistry.

By default, a dataset named `MCF7_m6A_rep1` writes below:

```text
refactored_outputs/datasets/MCF7_m6A_rep1/
```

The registry is stored in `refactored_outputs/dataset_registry.json`.
Registering a dataset never modifies the shipped RBP catalogue or Python
configuration.

New saved resource and dataset profiles use forward slashes for workspace-
relative paths, so those paths can move between Windows and macOS. Files
outside the workspace retain absolute paths and must be reselected when their
location changes. Older Windows backslash-relative profiles are migrated on
macOS only when the intended local target exists; a literal existing filename
with a backslash takes precedence. Ambiguous paths produce an explanatory
error rather than guessing. If a transferred old profile prevents startup,
use a new workspace or move its saved JSON aside and configure the local paths.

## Installation

Python 3.11 is recommended for a clean laptop installation. The tested and
supported range is Python 3.10-3.12.

### Double-click installation

Keep the helper files inside the `rna_mod_pipeline/` repository checkout.

On Windows:

1. Double-click `install_windows.bat` once.
2. After installation succeeds, double-click `run_windows.bat` whenever the
   GUI is needed.

On macOS:

1. Double-click `install_macos.command` once.
2. After installation succeeds, double-click `run_macos.command` whenever the
   GUI is needed.

Each installer creates `rna_mod_pipeline/.venv`, installs the core pipeline and
PySide6, checks dependency consistency, and verifies the installed commands.
The run helper always starts the GUI with
that environment. Internet access is normally required for the first
installation. Do not copy `.venv` between computers or operating systems;
each computer should create its own environment with its installer.

The code is installed in editable mode with the tested constraint files.
Changes to existing Python files under
`src/` or `scripts/` therefore take effect the next time the GUI or a command
is started; reinstalling is not necessary for ordinary code edits. The
installer compile-checks the Python source into a cache below `.venv`; normal
launches disable bytecode writes. Rerun the installer after changing
`pyproject.toml`, dependency files, or adding a new external package.

If macOS reports that a `.command` file is not executable after a non-Git file
transfer, open Terminal in `rna_mod_pipeline/` once and run:

```bash
chmod +x install_macos.command run_macos.command reset_environment_macos.command
```

### macOS or Linux

```bash
cd /path/to/rna_mod_workspace/rna_mod_pipeline
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --no-cache-dir -r requirements-build-lock.txt
python -m pip install --no-cache-dir --no-build-isolation \
  -c requirements-lock.txt -c requirements-gui-lock.txt -e ".[gui]"
rna-mod-gui --project-root ..
```

### Windows PowerShell

```powershell
cd C:\path\to\rna_mod_workspace\rna_mod_pipeline
py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --no-cache-dir -r requirements-build-lock.txt
python -m pip install --no-cache-dir --no-build-isolation -c requirements-lock.txt -c requirements-gui-lock.txt -e ".[gui]"
rna-mod-gui --project-root ..
```

The command-line pipeline remains usable without PySide6.

## Add and analyse a dataset

1. Open **Projects and datasets**.
2. Select **Resources…** and configure the reference, annotation, RBP-binding,
   expression, and KnockRBP files available on this computer. Leave genuinely
   unused optional resources blank. Saving checks that configured paths exist;
   paths inside the project workspace are stored relatively.
3. Select **Check resources…** to run the same portable integrity checks as
   `rna-mod-doctor`; this does not require a registered dataset.
4. Select **Add dataset…**.
5. Enter a unique dataset name.
6. Select its modification and bedMethyl file.
7. Select the matching FASTA and GTF when transcript or DRACH analysis is
   required, and confirm their genome-build and annotation-release labels.
8. Optionally record the sample, cell line, replicate, and a custom output
   root.
9. Open **Run workflow**, select the dataset and one workflow, inspect the
   exact command, and run validation before execution.

Shared binding, expression, KnockRBP, and STRING resources are requested only
by workflows that use them. No run requires all four historical bedMethyl
files.

For a new dataset, the usual order is:

1. **Phase 1 · Filter bedMethyl**;
2. for m6A, **Phase 1 · DRACH annotation** (this labels DRACH and non-DRACH
   calls; it does not discard the non-DRACH calls);
3. **Phase 1 · Transcript annotation**;
4. one or more loose or transcript-region overlap runs;
5. the relevant cross-database, expression, KnockRBP, STRING, or comparison
   workflow.

Each step remains independently runnable. The launcher validates that the
files needed by the selected step are present rather than silently rerunning
earlier steps.

The first interface runs one job at a time. It shows the exact argument list,
live output, run state, and a persistent job log. Existing results are not
overwritten unless **Advanced settings** is expanded and the overwrite option
is selected and the
underlying command accepts the existing manifest identity.

Advanced settings also expose analysis-facing FDR/plot thresholds,
database-specific source filters, expression model/column choices, KnockRBP
cell-line, dataset, and regulatory-heatmap-size controls, and STRING retry
controls. Deliberate
stress-test truncation controls are command-line-only and are not presented as
scientific GUI options. Loose-design downstream workflows are restricted to
all-transcript context; choose transcript-region design for 5′UTR, CDS, or
3′UTR analyses.

KnockRBP includes the cross-database regulatory network by default and expects
the corresponding cross-database workflow to be complete. Uncheck that option
to run target/DEG overlap alone, or select an explicit cross-database result in
Advanced settings.

## Compare two datasets

The **Compare two datasets** workflow requires:

- two distinct dataset names;
- the same biological modification;
- compatible reference and annotation identities;
- completed filtered and metagene tables for both datasets.

It compares exact site-set overlap and shared-site modification fraction and
coverage. If both profiles contain the same enrichment
`design × database` result, it also compares RBP effect sizes, directions, and
significance calls. An enrichment analysis present for only one dataset is
reported as unavailable for comparison, not as a null result.

When the same RBP resources are reused, agreement measures robustness to the
bedMethyl callset. It is not independent replication of CLIP experiments.

The comparison does not pool raw bedMethyl files. Technical-read pooling and
biological-replicate comparison require different assumptions, so pooling must
be implemented and audited as a separate workflow before it is offered in the
launcher.

## Results

The **Results** screen reads the run manifest, previews figures and bounded
parts of text tables, and checks that recorded outputs exist with their expected
sizes. Large site-level tables are never loaded in full by the interface.

Images use the whole preview area and **Fit** shows the complete figure without
stretching it. Resize the window or drag the divider beside the output list to
give the figure more space. For dense RBP labels, choose **100%** or **+**, then
drag the image or use its scrollbars. **−** zooms out and **Fit** restores the
whole-image view. **Open image** opens the unchanged original file in the
computer's normal image viewer. Zoom only changes the preview, not the PNG or
the analysis. Very tall plots necessarily have small labels in whole-image Fit
mode; use zoom to read them. The text preview is shown separately for tables.

## Cross-platform acceptance check

On a clean laptop:

1. Create a fresh virtual environment and install `.[gui]`.
2. Configure **Resources…**, then run `rna-mod-doctor --project-root ..` from
   `rna_mod_pipeline/`. Do not use the historical project audit unless the full
   four-callset Queen's project is present.
3. Start the launcher from a path containing a space to check path handling.
4. Add one small bedMethyl dataset and run Phase 1 into a new output root.
5. Close and restart the launcher and confirm the dataset registry reloads.
6. Test a deliberately missing FASTA or GTF and confirm preflight blocks the
   dependent workflow.
7. Cancel a test run and confirm an earlier completed output is unchanged.
8. Test STRING in offline mode with a populated cache, or separately test the
   network-enabled mode.
9. Run one representative full real-data workflow before relying on the new
   laptop for production analysis.

Scientific TSV outputs should agree with a command-line run on the same input.
PNG byte checksums can differ across operating systems because fonts and
rendering libraries differ.

The matched macOS/Windows checklist and sign-off table are in
[`CROSS_PLATFORM_ACCEPTANCE.md`](CROSS_PLATFORM_ACCEPTANCE.md).

Standalone `.app` or `.exe` packaging is deliberately deferred until this
source-and-virtual-environment pilot has passed. Desktop bundles must be built
and tested separately for each operating system and architecture; large
reference and database resources should remain external.
