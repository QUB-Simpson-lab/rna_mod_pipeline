# Command reference

The examples assume the shell is in the `rna_mod_pipeline/` repository inside
a dedicated workspace. `--project-root ..` therefore points to that workspace.
Explicit relative paths are resolved from the workspace root, not from the
current shell directory.

Every result-producing command refuses to replace existing results unless
`--overwrite` is supplied. Use a distinct output directory for analyses that
should coexist.
Outputs are staged and published only after successful completion, so a failed
overwrite retains the previous valid run. File-oriented workflows also verify
that existing files are owned by the matching run manifest.
An overwrite cannot omit or relocate files owned by the previous run. It may
add new destinations only when those paths do not already exist; this safely
allows an older result to acquire newly introduced plots. Rerunning a plotted
result with `--skip-plots`, moving its plot directory, or changing existing
output filenames requires fresh output paths.

## 1. `phase1_filter.py`

Filters one raw bedMethyl file. Defaults are coverage ≥20 and modification
percentage ≥20.

```bash
python scripts/phase1_filter.py \
  --project-root .. \
  --modification m6a
```

Default output:

```text
refactored_outputs/m6a/phase1/filtered_m6A.tsv
refactored_outputs/m6a/phase1/plots/01_coverage_distribution.png
refactored_outputs/m6a/phase1/plots/02_fraction_modified_distribution.png
refactored_outputs/m6a/phase1/plots/03_chromosome_distribution.png
refactored_outputs/m6a/phase1/plots/04_threshold_sensitivity.png
refactored_outputs/m6a/phase1/plots/04a_joint_threshold_heatmap.png
refactored_outputs/m6a/phase1/plots/05_filtered_chromosome_distribution.png
refactored_outputs/m6a/phase1/plots/06_filtered_fraction_distribution.png
refactored_outputs/m6a/phase1/plots/07_coverage_vs_fraction.png
```

Change the thresholds with `--min-coverage` and `--min-fraction`. Add
`--skip-plots` for a compute-only run.
The command rejects a bedMethyl whose modification codes conflict with
`--modification`.
Plot 04 changes coverage and fraction thresholds separately; plot 04a applies
their combinations jointly. These are per-run QC diagnostics, not complete
reruns of the downstream pipeline at alternative thresholds.

## 2. `phase1_drach.py`

Annotates m6A or its replicate with strand-corrected central five-mers and a
DRACH flag. It does not filter out non-DRACH calls.

```bash
python scripts/phase1_drach.py \
  --project-root .. \
  --modification m6a
```

This command applies only to `m6a` and `m6a_rep2`.
`is_DRACH` is missing, rather than false, where reference context cannot be
extracted.

## 3. `phase1_metagene.py`

Assigns sites to complete GENCODE Basic protein-coding transcripts and
5′UTR/CDS/3′UTR regions.

```bash
python scripts/phase1_metagene.py \
  --project-root .. \
  --modification m6a
```

For m6A, the default input is the DRACH-annotated all-context table. For m5C and
pseudouridine it is the filtered Phase 1 table.

## 4. `loose_overlap.py`

Runs exactly one modification and one binding database.

```bash
python scripts/loose_overlap.py \
  --project-root .. \
  --modification m6a \
  --database encori
```

Valid databases are `ornament`, `encori`, and `postar3`. Default output:

```text
refactored_outputs/m6a/loose/encori/
```

Each run writes its database-specific all-RBP overview as a descending
log2-Fisher-odds-ratio waterfall with `log2(OR) = 0` as the reference.
Observed `OR = 0` results are represented by individual purple
`log2(OR) = -infinity` boundary bars. Positive-infinity and non-estimable tests
are counted on the plot but omitted from the finite ranking; their exact
states remain in the enrichment TSV.

Equivalent independent commands are:

```bash
python scripts/loose_overlap.py --project-root .. --modification m6a --database ornament
python scripts/loose_overlap.py --project-root .. --modification m6a --database encori
python scripts/loose_overlap.py --project-root .. --modification m6a --database postar3
```

Useful options:

- `--window 10`
- `--background-min-coverage 20`
- `--background-max-fraction 20`
- `--case-min-coverage 20`
- `--case-min-fraction 20`
- `--minimum-clip-experiments 1` for ENCORI
- `--cell-types` and `--methods` for exact POSTAR3 source filtering
- `--site-annotation all|significant|none`
- `--skip-plots`

An unfiltered production run fails if a catalogue resource is missing or empty.
`--allow-empty-resources` is intended only for exploratory filtered subsets.
The complete case set is checked against raw bedMethyl by default. For an
intentional DRACH-only or other predeclared subset, add
`--allow-site-subset`; omitted qualifying calls are counted in the manifest.
Use a different `--output-dir` for every modification×database run. The command
rejects a folder whose manifest or result files belong to another loose run,
even with `--overwrite`.

## 5. `transcript_region_overlap.py`

Runs transcript×region-stratified overlap. It can process one database, a
selected subset, or all three.

Recommended all-database run:

```bash
python scripts/transcript_region_overlap.py \
  --project-root .. \
  --modification m6a
```

Default output:

```text
refactored_outputs/m6a/transcript_region/
```

Independent ENCORI-only run:

```bash
python scripts/transcript_region_overlap.py \
  --project-root .. \
  --modification m6a \
  --databases encori \
  --output-dir refactored_outputs/m6a/transcript_region_encori
```

Running all three together is more efficient because transcript opportunities
are constructed once. A one-database run is nevertheless a complete result for
that selected database and still enforces its full RBP catalogue.

The downstream transcript-region KnockRBP command is an integrated
ENCORI∪POSTAR3 target analysis and therefore needs a transcript run containing
both of those matrices. Cross-database and STRING steps likewise require the
database inputs implied by their scientific question.

Validate an existing run:

```bash
python scripts/transcript_region_overlap.py \
  --project-root .. \
  --validate-run refactored_outputs/m6a/transcript_region
```

The `--max-*`, `--rbps`, and `--chromosomes` arguments deliberately create
partial/test inputs and are recorded as such. Alternative complete thresholds,
windows, or strand policies are recorded as custom parameters rather than
mislabelled as truncated input.
ENCORI- and POSTAR3-specific filters are rejected when their corresponding
database is not selected.
An explicitly supplied oRNAment, ENCORI, or POSTAR3 resource path is likewise
rejected if that database is omitted from `--databases`. `--skip-input-hashes`
is valid only with `--validate-run`.

## 6. `plot_enrichment.py`

Regenerates transcript-region overview, top enriched/depleted,
modification-related, and optional regional plots. Confidence intervals remain
in the result tables but are not drawn. The all-RBP overview shows exact
`OR = 0` boundaries descriptively and omits other boundary and non-estimable
tests from its quantitative axis, reports their counts on the
PNG, and retains their rows and exclusion reasons in
`stratified_plot_source_data.tsv`. Finite estimates are shown in descending
order on the log2 odds-ratio scale, with `log2(OR) = 0` as the reference
line. Exact `OR = 0` boundaries are represented by individual terminal purple
negative-infinity bars labelled by RBP.
This command accepts only the transcript-region design and all-context overall
table. Supplying `--regional-results` adds separate 5′UTR, CDS, and 3′UTR
figures.

```bash
python scripts/plot_enrichment.py \
  --project-root .. \
  --modification m6a \
  --results refactored_outputs/m6a/transcript_region/stratified_enrichment_results.tsv \
  --regional-results refactored_outputs/m6a/transcript_region/region_specific_results.tsv \
  --output-dir refactored_outputs/m6a/transcript_region/replotted
```

Loose-overlap plots are generated directly by `loose_overlap.py`.

## 7. `cross_database.py`

Compares oRNAment, ENCORI, and POSTAR3 for one modification and design.

```bash
python scripts/cross_database.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --output-dir refactored_outputs/m6a/cross_database/loose
```

Transcript-region example:

```bash
python scripts/cross_database.py \
  --project-root .. \
  --modification m6a \
  --design transcript-region \
  --context all \
  --output-dir refactored_outputs/m6a/cross_database/transcript_region_all
```

For independently generated database folders, pass `--ornament`, `--encori`,
and `--postar3` explicitly. `--allow-missing` is required for a deliberately
incomplete comparison. A path supplied for a database that is not selected is
rejected instead of being silently ignored.

Figures and their TSV source data are generated by default under:

```text
<output-dir>/plots/
├── pairwise_scatter.png
├── pairwise_agreement.png
├── pairwise_shared_top.png
├── triple_validation.png
├── validation_summary.png
├── database_membership_counts.png
├── *_source.tsv
└── crossdb_plot_manifest.tsv
```

The membership plot reports exact database combinations and replaces an
area-proportional Venn diagram. Pairwise figures distinguish database absence,
tested non-significance, non-estimable results, boundary estimates, concordant
significance, and direction conflicts. No confidence intervals are drawn.
Use `--skip-plots` for tables only, `--plot-top-n` to change the top/shared and
triple panels, and `--plot-dpi` to change raster resolution.

## 8. `integrate_expression.py`

Integrates one enrichment result with either the in-house normalised expression
table or DepMap RNA-seq.

In-house profile:

```bash
python scripts/integrate_expression.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --database encori \
  --source nanopore \
  --expression data/Kate_231_0h_vs_6h_gene_counts_normalised.tsv \
  --enrichment refactored_outputs/m6a/loose/encori/encori_enrichment_results.tsv \
  --metagene refactored_outputs/m6a/phase1/filtered_m6A_metagene.tsv \
  --output-dir refactored_outputs/m6a/expression/nanopore
```

DepMap:

```bash
python scripts/integrate_expression.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --database encori \
  --source depmap \
  --expression data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --enrichment refactored_outputs/m6a/loose/encori/encori_enrichment_results.tsv \
  --metagene refactored_outputs/m6a/phase1/filtered_m6A_metagene.tsv \
  --model-id ACH-000768 \
  --tpm-threshold 1 \
  --output-dir refactored_outputs/m6a/expression/depmap
```

`--gene-lengths` optionally adds modification sites per kilobase alongside the
legacy raw site-count correlation.
For transcript-region results, use `--design transcript-region` and optionally
`--context 5UTR|CDS|3UTR`. The primary correlation subset is
detected/expressed modified genes; two clearly labelled exploratory subsets are
also written.

## 9. `knockrbp_validation.py`

Runs one modification, design, and transcript context. Multiple datasets for one
RBP remain separate; they cannot overwrite one another.

Loose design:

```bash
python scripts/knockrbp_validation.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --knockrbp-dir data/knockrbp \
  --metadata data/knockrbp/dataset_metadata.tsv \
  --encori-annotated refactored_outputs/m6a/loose/encori/filtered_m6A_encori_annotated.tsv \
  --postar3-annotated refactored_outputs/m6a/loose/postar3/filtered_m6A_postar3_annotated.tsv \
  --encori-enrichment refactored_outputs/m6a/loose/encori/encori_enrichment_results.tsv \
  --postar3-enrichment refactored_outputs/m6a/loose/postar3/postar3_enrichment_results.tsv \
  --cross-database refactored_outputs/m6a/cross_database/loose/cross_database_results.tsv \
  --regulatory-plot-top-n 40 \
  --output-dir refactored_outputs/m6a/knockrbp/loose
```

Transcript-region design:

```bash
python scripts/knockrbp_validation.py \
  --project-root .. \
  --modification m6a \
  --design transcript-region \
  --context 3UTR \
  --knockrbp-dir data/knockrbp \
  --metadata data/knockrbp/dataset_metadata.tsv \
  --transcript-run-dir refactored_outputs/m6a/transcript_region \
  --cross-database refactored_outputs/m6a/cross_database/transcript_region_3UTR/cross_database_results.tsv \
  --regulatory-plot-top-n 40 \
  --output-dir refactored_outputs/m6a/knockrbp/transcript_region_3UTR
```

Default eligible cell lines are MDA-MB-231 and MDA-MB-231-LM2. Override with
`--cell-lines`; select exact experiments with `--dataset-ids`.
The enrichment companions are inferred from the annotated-table directories
when they use standard names. Supplying them explicitly makes usable versus
empty/filtered resources auditable.
For transcript-region analysis, the `--context` value must match the supplied
cross-database result; the command rejects, for example, an all-context table
used with `--context 3UTR`.
The complete regulatory network remains in `rbp_regulatory_network.tsv`.
`--regulatory-plot-top-n` controls only the readable heatmap subset, ranked by
maximum absolute retained-DEG log2 fold change. The exact selected source rows
are written to `rbp_regulatory_network_plot_source.tsv`; grey heatmap cells mean
that no DEG passed the declared filters for that RBP/dataset pair, not zero
fold change.

## 10. `string_analysis.py`

Selects direction-concordant ENCORI+POSTAR3 RBP sets and queries STRING.

```bash
python scripts/string_analysis.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --cross-database refactored_outputs/m6a/cross_database/loose/cross_database_results.tsv \
  --cache-dir refactored_outputs/string_cache \
  --output-dir refactored_outputs/m6a/string/loose
```

Repeat a cached analysis without network access by adding `--offline`. STRING is
not version-pinned by its API, so retaining the cache and manifest is essential
for reproducibility.
The cache and output directories must be separate and non-nested.

## 11. `compare_modifications.py`

Compares selected modifications across selected databases. It writes the
RBP×modification matrix, shared/specific/direction-switching patterns, pairwise
effect correlations with FDR, overlap summaries, and Phase 1 site/metagene
context.

```bash
python scripts/compare_modifications.py \
  --project-root .. \
  --modifications m6a m5c pseu \
  --databases ornament encori postar3 \
  --design loose \
  --output-dir refactored_outputs/cross_modification/loose
```

Use repeated `--input MODIFICATION DATABASE PATH` arguments to compare custom or
independently generated result folders. Add `--skip-phase1-context` when
metagene/site-count plots are not required.
Without that switch, every requested modification must have a metagene input.
`--allow-missing` permits an intentionally incomplete comparison and records
every missing database and Phase-1 input.

## 12. `compare_datasets.py`

Compares two independently processed datasets with the same biological
modification. The profile JSON files identify the filtered and metagene tables,
reference/annotation identity, and any compatible enrichment tables.

```bash
python scripts/compare_datasets.py \
  --project-root .. \
  --left-profile refactored_outputs/datasets/MDA231_m6A_rep1/dataset_comparison_profile.json \
  --right-profile refactored_outputs/datasets/MDA231_m6A_rep2/dataset_comparison_profile.json \
  --output-dir refactored_outputs/datasets/comparison_MDA231_m6A_rep1_vs_MDA231_m6A_rep2
```

The desktop launcher creates these profile sidecars immediately before it
starts a dataset-comparison run, after preflight validation has passed.
Command-line users must provide already-existing profile JSON files;
`compare_datasets.py` does not create them from bedMethyl inputs.
The workflow compares exact site overlap and shared-site fraction/coverage
concordance. It also compares the intersection of declared
`design × database` enrichment results. A result available for only one dataset
is recorded as unavailable for comparison, not as a null finding.

This is a callset-robustness analysis. It neither pools raw replicates nor
constitutes independent replication of reused CLIP resources.

## 13. `pipeline_gui.py`

Starts the optional dataset-centred desktop launcher:

```bash
python -m pip install --no-build-isolation -c requirements-lock.txt -c requirements-gui-lock.txt -e ".[gui]"
rna-mod-gui --project-root ..
```

The launcher accepts arbitrary named m6A, m5C, or pseudouridine bedMethyl
inputs. It builds and displays the exact commands documented above; it does not
perform scientific calculations itself. See `docs/GUI_GUIDE.md`.

## 14. `validate_handoff.py`

Checks the required project inputs, RBP catalogues and source resources,
KnockRBP dataset mapping, command-wrapper structure, and current legacy
baselines:

```bash
python scripts/validate_handoff.py --project-root ..
```

`--skip-legacy-results` omits established result baselines but still requires
the historical project input layout. Use `rna-mod-doctor --project-root ..`
for a portable resource-profile check without the four callsets.

## Recommended complete order

For each modification:

1. Run `rna-mod-doctor --project-root ..` before the first analysis on a clean
   installation. Run `validate_handoff.py` only when the complete historical
   Queen's project is present.
2. Run `phase1_filter.py`.
3. Run `phase1_drach.py` for `m6a` and `m6a_rep2`; omit it for `m5c` and
   `pseu`.
4. Run `phase1_metagene.py`.
5. Run `loose_overlap.py` once per database.
6. Run `transcript_region_overlap.py`, then its `--validate-run` mode.
7. Run `plot_enrichment.py` only when regenerating transcript-region figures
   in a separate output directory.
8. Run `cross_database.py` once per design/context.
9. Run `integrate_expression.py`.
10. Run `knockrbp_validation.py`.
11. Run `string_analysis.py`.

After at least two modifications are complete, run
`compare_modifications.py`.

After two independently processed datasets of the same modification are
complete, use `compare_datasets.py` to assess callset and downstream effect
robustness.

The optional `pipeline_gui.py` provides the same workflows through the desktop
launcher and creates the profile sidecars needed by `compare_datasets.py`.

## Acceptance and resource verification utilities

`prepare_acceptance_demo.py` creates a new self-contained synthetic workspace
and can run the core workflows. It rejects a nonempty destination. The
example never changes the shipped biological RBP catalogue.

```bash
python scripts/prepare_acceptance_demo.py --workspace ../rna_mod_demo --run
python scripts/check_acceptance_demo.py --workspace ../rna_mod_demo
```

Add `--skip-plots` to both commands for a numerical-only test. For manual GUI
execution, omit `--run`, then start
`python scripts/pipeline_gui.py --project-root ../rna_mod_demo`.
See `examples/acceptance/README.md` for expected counts/effects and limitations.

`verify_resources.py` streams selected real resources and compares byte sizes
and SHA-256 hashes with `config/resource_snapshot_20260908.json`:

```bash
python scripts/verify_resources.py --workspace .. --groups reference ornament encori
```

Omit `--groups` to check all 435 files. Repeated `--map` arguments can relocate
recorded files/directories; they do not edit the GUI resource profile. New
inventories require `--create --output NEW_FILE` and never overwrite an
existing snapshot. See `docs/RESOURCE_SNAPSHOT.md` for exact commands.
