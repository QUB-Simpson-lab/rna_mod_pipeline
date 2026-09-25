# Synthetic laptop acceptance example

This example creates artificial software-test inputs. It contains **no human
reference sequence, real RBP binding data, experimental cell line or research
results**. Every `DEMO_*` name and interaction is invented. Do not incorporate
these outputs into the scientific findings or manuscript.

The example needs no downloads or internet connection after installation. It
tests the actual public pipeline scripts and the GUI command builder. It does
not relax the production catalogue sizes, complete-case requirement, or evidence
quality policy.

## Run automatically

Install the pipeline first. From the source-code folder on macOS:

```bash
.venv/bin/python scripts/prepare_acceptance_demo.py --workspace ../rna_mod_demo --run
.venv/bin/python scripts/check_acceptance_demo.py --workspace ../rna_mod_demo
```

From Command Prompt on Windows:

```bat
.venv\Scripts\python.exe scripts\prepare_acceptance_demo.py --workspace ..\rna_mod_demo --run
.venv\Scripts\python.exe scripts\check_acceptance_demo.py --workspace ..\rna_mod_demo
```

`--workspace` must name a new or empty directory outside the source-code folder.
Preparation refuses to overwrite an existing workspace. Use a different folder
name for another automatic run. Keep paths quoted if they contain spaces.

The automatic run executes filtering, DRACH annotation, metagene annotation,
all three loose overlap databases, all three transcript-region databases,
transcript output validation, and both cross-database comparisons. Plotting is
enabled. Commands and individual logs are saved in `acceptance_logs/`.
The final checker writes `acceptance_check.json` after every check succeeds.

For a faster numerical test, supply `--skip-plots` to **both** commands. This does
not verify figure generation and is not a substitute for the plotted laptop test.

## Test through the GUI

Create a separate empty workspace without the `--run` option. Then launch the GUI
against that workspace:

```bash
# macOS, from the source folder
.venv/bin/python scripts/prepare_acceptance_demo.py --workspace ../rna_mod_gui_demo
.venv/bin/python scripts/pipeline_gui.py --project-root ../rna_mod_gui_demo
```

```bat
REM Windows Command Prompt, from the source folder
.venv\Scripts\python.exe scripts\prepare_acceptance_demo.py --workspace ..\rna_mod_gui_demo
.venv\Scripts\python.exe scripts\pipeline_gui.py --project-root ..\rna_mod_gui_demo
```

The resource paths and two datasets are already registered. Select
`SYNTHETIC_m6a`; keep all workflow thresholds at their defaults. Run:

1. Phase 1 filtering.
2. DRACH annotation.
3. Metagene annotation.
4. Loose overlap separately for oRNAment, ENCORI and POSTAR3.
5. Transcript-region overlap with all three databases selected.
6. Validate transcript results.
7. Cross-database comparison with loose design and all context.
8. Cross-database comparison with transcript-region design and all context.

Run the checker with `--workspace ../rna_mod_gui_demo` afterwards. The resource
profile uses the synthetic catalogue rather than the real biological catalogue.
Do not replace its resources with hg38 or real CLIP files. The second dataset,
`SYNTHETIC_identical_rep`, is an exact artificial copy for optional two-dataset GUI
tests; it is not a biological replicate. Its outputs are not required by the
core acceptance checker.

Also open several PNGs and TSVs through the Results page, close/reopen the app,
and verify that datasets and resources reload. Reload a completed output
manifest in Results; GUI job logs remain on disk, but this version has no
automatic job-history browser. Test cancellation
on a separate disposable run so you can inspect that completed outputs survive.
The numerical checker cannot establish these interactive behaviours.

## Independently known results

The artificial reference has 32 non-overlapping complete protein-coding
transcripts, 16 on each strand. Each has a 400-base 5′UTR, CDS and 3′UTR. Every
transcript-region stratum contains four high-call cases and four low-call
controls. Sites are 40 bases apart, so a ±10-base overlap window cannot reach a
neighbouring site's artificial interval.

| Quantity | Expected |
|---|---:|
| Raw bedMethyl rows | 800 |
| Qualifying cases, coverage ≥20 and fraction ≥20% | 384 |
| Eligible low-call controls | 384 |
| High-fraction rows excluded for coverage 10 | 32 |
| DRACH cases / non-DRACH cases | 288 / 96 |
| Cases in each of 5′UTR, CDS and 3′UTR | 128 |
| Cases on each strand | 192 |
| Contributing genes / transcript-region strata | 32 / 96 |
| Declared oRNAment / ENCORI / POSTAR3 RBPs | 133 / 281 / 216 |

Every database contains the following sentinels, followed by artificial neutral
RBPs that preserve the full declared panel size:

| Artificial RBP | Case overlaps | Control overlaps | Loose OR and MH OR |
|---|---:|---:|---:|
| DEMO_ENRICHED | 288 | 96 | 9 |
| DEMO_DEPLETED | 96 | 288 | 1/9 |
| DEMO_NULL | 192 | 192 | 1 |
| DEMO_ZERO | 0 | 192 | 0 |
| DEMO_INFINITY | 192 | 0 | ∞ |
| DEMO_ALL_BOUND | 384 | 384 | Not estimable |

The enriched pattern has two bound cases per region in eight genes, three in
sixteen genes, and four in eight genes, with one bound control throughout.
Depletion swaps case and control counts. Thus the global OR is exactly 9 or 1/9,
and the Mantel–Haenszel numerator/denominator are 108/12 or 12/108. Variation
between genes gives a nonzero gene-cluster standard error:
`sqrt(8/279)`, with 31 degrees of freedom. The checker evaluates the two-sided
Student-t tail independently through the incomplete beta function, and Fisher
probabilities directly from integer combinations; it does not compare outputs
against a saved output snapshot. Both finite directional sentinels must pass
the quality checks and FDR in all three databases and all transcript regions.

Neutral RBPs deliberately reverse their within-gene pattern between gene halves,
giving global OR=1 and p=1. Boundary and non-estimable sentinels must not become
primary triple consensus. The only primary triple-direction-concordant RBPs are
`DEMO_ENRICHED` and `DEMO_DEPLETED` for both designs.

Input hashes are deterministic across laptops. Run dates, absolute command paths,
plot metadata and compressed-output bytes may differ, so compare verified table
values and `acceptance_check.json`, not whole output-directory checksums.

## Scope and limitations

This compact example checks phase 1, database loaders, both overlap designs,
cross-database selection, expected boundary handling and PNG integrity. It does
not validate biological methods, performance on full human datasets, internet
connectivity, STRING API behaviour, expression integration, KnockRBP or
cross-modification analysis. Those workflows need their separate regression and
representative-data checks. The synthetic reference is intentionally labelled
`SYNTHETIC_NOT_HG38` and must never be used for an experimental bedMethyl file.

The machinery heatmaps are expected to be labelled empty placeholders because
invented `DEMO_*` proteins are not part of the biological machinery registry.
The all-RBP overview has a long flat centre because most artificial RBPs have
OR=1. These are intentional test properties, not missing outputs. The full check
requires all 66 expected PNGs and verifies that each can be decoded.
