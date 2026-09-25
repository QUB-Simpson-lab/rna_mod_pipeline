# RNA-modification RBP pipeline

This repository contains the maintained Queen's University Belfast RNA-
modification pipeline. Analyse any supported m6A, m5C or pseudouridine
bedMethyl file as a named dataset, or compare two datasets of the same
modification. Each workflow can be run from the GUI or command line.

For the two-laptop trial, start with [`START_HERE.md`](START_HERE.md).
The self-contained synthetic acceptance example needs no research databases
and checks known numerical answers. Its generated data are software fixtures,
not experimental findings. The current preparation status is recorded in
[`docs/RELEASE_READINESS.md`](docs/RELEASE_READINESS.md).

The legacy statistical engines and enrichment tables remain unchanged. The
legacy transcript-region plotting script and its overview plot artefacts were
updated in parallel with this refactored implementation so both use the same
boundary/non-estimable display policy. New analyses should write to a new
output directory.

The refactor itself does not require the established m6A, m5C, pseudouridine,
or replicate results to be rerun. It is a clean implementation for future
runs; parity and real-data checks are documented in
[`docs/REFACTOR_AUDIT.md`](docs/REFACTOR_AUDIT.md). The stabilisation review is
in [`docs/FINAL_HANDOFF_AUDIT.md`](docs/FINAL_HANDOFF_AUDIT.md).

## Included workflows

- Phase 1: bedMethyl filtering, m6A DRACH annotation, and transcript/metagene
  annotation.
- Loose overlap: oRNAment, ENCORI, or POSTAR3 enrichment against the complete
  coverage-qualified low-call background.
- Transcript-region overlap: enrichment stratified by exact transcript and
  5′UTR/CDS/3′UTR region.
- Cross-database validation with pairwise, triple-validation, membership, and
  summary plots.
- Expression integration.
- Cross-modification comparison.
- KnockRBP validation for loose or transcript-region results.
- STRING analysis for loose or transcript-region results.
- Reusable overview, modification-related, enriched/depleted, and regional
  plots.
- Same-modification dataset/callset robustness comparison.
- An optional dataset-centred desktop launcher for arbitrary supported
  bedMethyl inputs.

Loose and transcript-region all-RBP overview plots use the same visual
direction: finite estimates descend from highest to lowest on a `log2(OR)`
axis with `log2(OR) = 0` as the reference line. The loose plot shows Fisher
odds ratios; the transcript-region plot shows Mantel-Haenszel odds ratios.
Exact `OR = 0` results are shown as individual terminal purple
negative-infinity bars labelled by RBP. Positive-infinity and non-estimable tests are counted
in text on the figure and retained in the result/source tables.

Strict matched-site analysis is deliberately not included.

## Data flow

```text
bedMethyl + hg38 + GENCODE
        │
        └─ Phase 1 filtered/all-context site tables
                 ├─ loose overlap ───────────────┐
                 └─ transcript×region overlap ──┤
                                                ├─ cross-database validation
                                                ├─ expression integration
                                                ├─ KnockRBP target/DEG validation
                                                ├─ STRING functional context
                                                └─ cross-modification comparison
```

## Design principles

- Loose overlap runs exactly one modification and one database. Transcript
  overlap runs one modification with any selected database subset.
- Inputs, outputs, thresholds, windows, and source filters are command-line
  arguments.
- Existing outputs are protected unless `--overwrite` is supplied. Runs stage
  all expected files before publication, so a failed rerun preserves the
  previous complete result.
- An overwrite cannot silently strand plots or tables owned by the prior
  manifest. Use fresh output paths when changing the selected output set.
- Raw-data, legacy-result, and source-code trees are never valid output
  destinations, even with `--overwrite`.
- A custom loose-overlap output folder cannot be reused for another
  modification or database; use one directory per run.
- Within each design, all three binding-source adapters feed a shared overlap
  and enrichment engine.
- RBP aliases and tested RBP universes are versioned in
  `config/rbp_catalog.tsv`.
- Raw statistical calls and quality-qualified primary calls are retained in
  parallel. Boundary, partial-input, low-cluster, and leave-one-gene-out
  warning rows remain auditable but cannot silently enter primary consensus.
- Manifests store paths relative to the selected project root when possible;
  resources outside that root are recorded with absolute paths. They also
  record file sizes, SHA-256 checksums, package/library versions, and the Git
  commit/worktree state when available.
- Scientific corrections relative to the legacy code are listed in
  `docs/REFACTOR_AUDIT.md`.

## Installation

Python 3.10-3.12 is supported; Python 3.11 is recommended.

The supported distribution is the complete source checkout, including `scripts/`
and `config/rbp_catalog.tsv`, installed in editable mode. A wheel can be built
as a packaging check, but it is not a standalone workflow distribution and
must not be distributed without those checkout resources.

```bash
cd /path/to/rna_mod_workspace/rna_mod_pipeline
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --no-cache-dir -r requirements-build-lock.txt
python -m pip install --no-cache-dir --no-build-isolation \
  -c requirements-lock.txt -e .
```

Replace `python3.11` with Python 3.10 or 3.12 if necessary. After activation,
the examples consistently use the virtual environment's `python`.

The tested scientific package versions are recorded in
`requirements-lock.txt`; the tested optional desktop packages are in
`requirements-gui-lock.txt`, and installer tooling is in
`requirements-build-lock.txt`. These constrain tested versions but do not
provide package hashes. The looser compatibility ranges are in
`requirements.txt`, `requirements-gui.txt`, and `pyproject.toml`.

Install and start the optional graphical launcher with:

```bash
python -m pip install --no-cache-dir --no-build-isolation \
  -c requirements-lock.txt -c requirements-gui-lock.txt -e ".[gui]"
rna-mod-gui --project-root ..
```

For a double-click installation, use the platform helpers in this folder:

- Windows: `install_windows.bat`, followed by `run_windows.bat`.
- macOS: `install_macos.command`, followed by `run_macos.command`.

The guarded `reset_environment_windows.bat` and
`reset_environment_macos.command` helpers remove only the regenerable local
`.venv`; they never remove scientific inputs or outputs.

Each installer creates an isolated `.venv`, installs the scientific and GUI
dependencies, and checks both package imports. The run helper always uses that
environment and resolves the project directory relative to its own location,
so spaces in the path are supported.

The macOS installer also checks standard Homebrew and python.org installation
locations when the launch environment does not expose Python on `PATH`.

The launcher does not require the four historical bedMethyl files. Each
arbitrary m6A, m5C, or pseudouridine input is registered as a named dataset with
an isolated output workspace. See
[`docs/GUI_GUIDE.md`](docs/GUI_GUIDE.md).
Resource acquisition and exact layout are documented in
[`docs/RESOURCE_SETUP.md`](docs/RESOURCE_SETUP.md); environment constraints,
tests, and CI caveats are in
[`docs/INSTALLATION_AND_TESTING.md`](docs/INSTALLATION_AND_TESTING.md).
Use [`docs/CROSS_PLATFORM_ACCEPTANCE.md`](docs/CROSS_PLATFORM_ACCEPTANCE.md)
for the matched macOS and Windows release test.

After configuring resources in the GUI, or placing them in the established
project layout, run the portable check:

```bash
rna-mod-doctor --project-root ..
```

This validates the configured reference and binding panels without requiring
the four historical bedMethyl callsets or their result folders. The separate
`rna-mod-audit` command is the historical-project regression check.

## Independent execution

Each command processes only the requested unit. For example, a POSTAR3 m5C run
does not rerun m6A, oRNAment, or ENCORI. The main selectors are:

```text
--modification m6a|m5c|pseu|m6a_rep2
--database ornament|encori|postar3
--design loose|transcript-region
```

Run any script with `--help` for its complete, current argument list:

```bash
python scripts/loose_overlap.py --help
```

Complete commands for each script and a recommended end-to-end sequence are in
`docs/COMMAND_REFERENCE.md`.

Database-specific filters are validated against the selected database. For
example, ENCORI support thresholds cannot be silently applied to an
oRNAment-only run.

The downloaded oRNAment BED archive has an unresolved endpoint-convention
ambiguity. For reproducibility, loose and transcript-region workflows preserve
their established design-specific interpretations and record that policy in
their manifests. See [`docs/METHODS.md`](docs/METHODS.md) before comparing
oRNAment results between the two designs.

The overlap calculations themselves are independent. Some integrated
downstream questions necessarily need more than one source: cross-database
validation needs at least two databases (three by default), STRING selection
uses concordant ENCORI and POSTAR3 results, and transcript-region KnockRBP
target unions require both ENCORI and POSTAR3 case-binding matrices.

Cross-database figures are written below the selected comparison directory in
`plots/`. Their exact plotted values and selection flags are retained as TSV
source-data files beside the figures. Use `--skip-plots` for a table-only run.
KnockRBP likewise retains the complete regulatory network as TSV while its
heatmap uses a configurable top-effect subset with a separate plot-source TSV.

## Statistical designs

### Loose overlap

Cases are all filtered high-call sites. Controls are all positions in the same
bedMethyl file with coverage at least the chosen threshold and modification
fraction below the chosen high-call threshold. Sites are not matched by
transcript, region, coverage, motif, or distance. Binding means that an RBP
interval intersects the site-centred window. Each RBP is tested with a
two-sided Fisher exact test, followed by Benjamini-Hochberg correction within
one modification and database.

By default, the case table must be the complete set of raw calls meeting its
declared case thresholds. An intentional subset such as DRACH-only requires
`--allow-site-subset`; the manifest records how many qualifying calls were
omitted.

### Transcript-region overlap

Cases and controls are assigned to complete GENCODE Basic protein-coding
transcripts and 5′UTR, CDS, or 3′UTR. Only exact
`transcript_id × transcript_region` strata containing both cases and controls
contribute. All eligible controls are retained; there is no one-to-one matching
and no matching on coverage, sequence, or transcript distance. The primary
effect is a Mantel-Haenszel common odds ratio with gene-cluster inference.
Regional estimates are descriptive unless a formal region-interaction test is
performed separately.

### Binding-source interpretation

oRNAment contains sequence-based motif predictions. ENCORI and POSTAR3 aggregate
CLIP-derived binding records across available experiments and cell contexts.
Cross-database agreement therefore strengthens robustness but is not three
independent biological replications.

## Command reference

The compact command list covering all 14 entry points is in
`docs/ALL_SCRIPT_COMMANDS.md`. Full options, dependencies, and worked examples
are maintained in `docs/COMMAND_REFERENCE.md`.

## Installation validation

```bash
python scripts/validate_handoff.py --project-root ..
```

This validates the installation, catalogues, representative legacy outputs,
and portable manifests without rerunning multi-gigabyte source scans. It is
intended for the complete Queen's project. Use `rna-mod-doctor` for a clean
laptop or arbitrary-dataset installation. The compact synthetic regression
suite is retained in `tests/` and its final result is recorded in
`docs/FINAL_HANDOFF_AUDIT.md`.

## Operational notes

- Keep the reference FASTA, GENCODE GTF, bedMethyl files, and binding-database
  resources outside this code folder. Standard project-layout defaults locate
  them automatically; pass explicit paths when using another layout.
- A STRING run requires network access unless all required responses are already
  present in its cache.
- POSTAR3 is the largest input and should be stored uncompressed only when disk
  space permits; `.gz` is also supported by the source loader.
- Never interpret depletion as proof that an RBP is an “anti-reader,” or
  KnockRBP expression changes as proof of direct regulation.
- See `docs/SCOPE.md` for deliberately excluded legacy analyses, including the
  strict matched-site workflow and the RM2Target/WRE annotation layer.
- The optional desktop launcher remains a thin process layer over the command
  scripts. It does not contain scientific algorithms and is not required for
  command-line operation.

## Licence and credits

The software is supplied under the [MIT licence](LICENSE). External reference
genomes, binding databases and other research inputs retain their own terms;
they are not bundled or relicensed by this repository.

Alexandru Zob developed the studentship project with supervision from David
Simpson at Queen's University Belfast and day-to-day support from Stephen.
`CITATION.cff` currently identifies Alexandru as the software author; the lab
should review the final contributor list before publishing a formal release.
No unverified surnames, ORCIDs or publication identifiers have been added.
