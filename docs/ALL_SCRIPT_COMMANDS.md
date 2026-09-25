# All pipeline script commands

This is the compact, copy-ready list for every runnable script in
`scripts/`. Run the commands from the repository root after
activating its virtual environment:

```bash
cd /path/to/rna_mod_workspace/rna_mod_pipeline
source .venv/bin/activate
```

The commands below use the current project layout and m6A as the worked
example. Run upstream commands before the downstream commands that consume
their outputs. Existing owned results are protected; add `--overwrite` only
when intentionally rerunning the same output set.

## Complete m6A command chain

### 1. Validate the installation and project resources

```bash
python scripts/validate_handoff.py --project-root ..
```

That command is the historical-project regression audit. On a clean laptop
without all four established callsets and result folders, use:

```bash
rna-mod-doctor --project-root ..
```

### 2. Filter and create all Phase 1 QC plots

```bash
python scripts/phase1_filter.py \
  --project-root .. \
  --modification m6a
```

### 3. Annotate DRACH context

```bash
python scripts/phase1_drach.py \
  --project-root .. \
  --modification m6a
```

Run this for `m6a` and `m6a_rep2`; omit it for `m5c` and `pseu`.

### 4. Annotate transcripts and metagene position

```bash
python scripts/phase1_metagene.py \
  --project-root .. \
  --modification m6a
```

### 5. Run loose enrichment independently for all three databases

```bash
python scripts/loose_overlap.py \
  --project-root .. \
  --modification m6a \
  --database ornament
```

```bash
python scripts/loose_overlap.py \
  --project-root .. \
  --modification m6a \
  --database encori
```

```bash
python scripts/loose_overlap.py \
  --project-root .. \
  --modification m6a \
  --database postar3
```

Each command generates its own descending log2-odds-ratio all-RBP overview and
top enriched/depleted plots.

### 6. Run transcript-region enrichment

```bash
python scripts/transcript_region_overlap.py \
  --project-root .. \
  --modification m6a
```

Validate the completed transcript-region run:

```bash
python scripts/transcript_region_overlap.py \
  --project-root .. \
  --validate-run refactored_outputs/m6a/transcript_region
```

### 7. Regenerate transcript-region plots separately

```bash
python scripts/plot_enrichment.py \
  --project-root .. \
  --modification m6a \
  --results refactored_outputs/m6a/transcript_region/stratified_enrichment_results.tsv \
  --regional-results refactored_outputs/m6a/transcript_region/region_specific_results.tsv \
  --output-dir refactored_outputs/m6a/transcript_region/replotted
```

This is optional because the transcript-region analysis already creates its
standard plots. Use it when a separate plot-only output is wanted.

### 8. Compare the three databases

Loose design:

```bash
python scripts/cross_database.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --context all \
  --output-dir refactored_outputs/m6a/cross_database/loose
```

Transcript-region all-context design:

```bash
python scripts/cross_database.py \
  --project-root .. \
  --modification m6a \
  --design transcript-region \
  --context all \
  --output-dir refactored_outputs/m6a/cross_database/transcript_region_all
```

### 9. Integrate expression

In-house normalised counts:

```bash
python scripts/integrate_expression.py \
  --project-root .. \
  --modification m6a \
  --source nanopore \
  --design loose \
  --context all \
  --database encori \
  --expression data/Kate_231_0h_vs_6h_gene_counts_normalised.tsv \
  --enrichment refactored_outputs/m6a/loose/encori/encori_enrichment_results.tsv \
  --metagene refactored_outputs/m6a/phase1/filtered_m6A_metagene.tsv \
  --output-dir refactored_outputs/m6a/expression/nanopore
```

DepMap/CCLE:

```bash
python scripts/integrate_expression.py \
  --project-root .. \
  --modification m6a \
  --source depmap \
  --design loose \
  --context all \
  --database encori \
  --expression data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
  --enrichment refactored_outputs/m6a/loose/encori/encori_enrichment_results.tsv \
  --metagene refactored_outputs/m6a/phase1/filtered_m6A_metagene.tsv \
  --model-id ACH-000768 \
  --tpm-threshold 1 \
  --output-dir refactored_outputs/m6a/expression/depmap
```

### 10. Run KnockRBP validation

Loose design:

```bash
python scripts/knockrbp_validation.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --context all \
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

Transcript-region example:

```bash
python scripts/knockrbp_validation.py \
  --project-root .. \
  --modification m6a \
  --design transcript-region \
  --context all \
  --knockrbp-dir data/knockrbp \
  --metadata data/knockrbp/dataset_metadata.tsv \
  --transcript-run-dir refactored_outputs/m6a/transcript_region \
  --cross-database refactored_outputs/m6a/cross_database/transcript_region_all/cross_database_results.tsv \
  --regulatory-plot-top-n 40 \
  --output-dir refactored_outputs/m6a/knockrbp/transcript_region_all
```

The full regulatory network is always retained in TSV form. The PNG uses the
top affected RBPs ranked by maximum absolute retained-DEG log2 fold change,
and `rbp_regulatory_network_plot_source.tsv` records the exact plotted rows.

### 11. Run STRING analysis

Loose design:

```bash
python scripts/string_analysis.py \
  --project-root .. \
  --modification m6a \
  --design loose \
  --context all \
  --cross-database refactored_outputs/m6a/cross_database/loose/cross_database_results.tsv \
  --cache-dir refactored_outputs/string_cache \
  --output-dir refactored_outputs/m6a/string/loose
```

Transcript-region design:

```bash
python scripts/string_analysis.py \
  --project-root .. \
  --modification m6a \
  --design transcript-region \
  --context all \
  --cross-database refactored_outputs/m6a/cross_database/transcript_region_all/cross_database_results.tsv \
  --cache-dir refactored_outputs/string_cache \
  --output-dir refactored_outputs/m6a/string/transcript_region_all
```

STRING needs network access on its first run. Add `--offline` only when the
required responses are already present in the cache.

### 12. Compare modifications

Run the relevant Phase 1 and overlap commands for m6A, m5C, and
pseudouridine first, then:

```bash
python scripts/compare_modifications.py \
  --project-root .. \
  --modifications m6a m5c pseu \
  --databases ornament encori postar3 \
  --design loose \
  --context all \
  --output-dir refactored_outputs/cross_modification/loose
```

### 13. Compare two datasets of the same modification

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
`compare_datasets.py` does not create them from bedMethyl inputs. This example
is therefore conditional on both named datasets having already been registered
and processed.

### 14. Start the graphical launcher

Double-click helpers are available for a clean installation:

- Windows: `install_windows.bat`, then `run_windows.bat`.
- macOS: `install_macos.command`, then `run_macos.command`.

The equivalent terminal commands are:

Install its optional dependency once:

```bash
python -m pip install --no-build-isolation -c requirements-lock.txt -c requirements-gui-lock.txt -e ".[gui]"
```

Then start it:

```bash
rna-mod-gui --project-root ..
```

## Modification substitutions

| Analysis | Raw bedMethyl used by default | Output root | DRACH step |
|---|---|---|---|
| `m6a` | `m6a/data/data.bedmethyl` | `refactored_outputs/m6a` | yes |
| `m5c` | `m5c/data/m5c_pseU__m5C.bedmethyl` | `refactored_outputs/m5c` | no |
| `pseu` | `pseudouridine/data/m5c_pseU__pseU_17802.bedmethyl` | `refactored_outputs/pseu` | no |
| `m6a_rep2` | `m6a_rep2/data/data_rep2_m6a.bedmethyl` | `refactored_outputs/m6a_rep2` | yes |

Replace `--modification m6a` and the matching `refactored_outputs/m6a` path
together. Use `--help` on any command for every optional selector. Detailed
explanations are in `docs/COMMAND_REFERENCE.md`.

## Laptop acceptance and resource utilities

These three preparation/checking scripts are additional to the 14 workflow
entry points above. Run them from the code folder with its installed Python.

```bash
python scripts/prepare_acceptance_demo.py --workspace ../rna_mod_demo --run
python scripts/check_acceptance_demo.py --workspace ../rna_mod_demo
python scripts/verify_resources.py --workspace .. --groups reference ornament encori
```

The first command creates explicitly synthetic inputs in a new external folder;
the second verifies their independently known answers. No real resources are
required for those two commands. The last command checks real external files
against the frozen project snapshot; select only resource groups you have.
Use `--help` for path remapping and other options. Double-click alternatives
for the demonstration are listed in `START_HERE.md`.
