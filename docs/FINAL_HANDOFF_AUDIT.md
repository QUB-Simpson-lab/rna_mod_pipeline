# Release stabilisation audit

> Historical 1.0.1 baseline. The 8 September 2026 source preparation and current
> verification are recorded in `RELEASE_READINESS.md`. Licence, citation,
> test counts and portability status below describe the August audit, not the
> later candidate. Native laptop sign-off remains separate.

**Audit date:** 6 August 2026  
**Audited distribution:** complete source checkout, version 1.0.1  
**Intended release status:** cross-platform pilot pending public-release approval

## Verdict

The complete source checkout is stable enough for native macOS and Windows
acceptance testing and for new analyses written to new output directories. No regression was found in the
scientific paths exercised below, and the historical m6A, m5C,
pseudouridine, and m6A-replicate results do **not** need to be rerun because of
this stabilisation work.

This is not a claim of absolute correctness. Native Windows execution has not
been run on this macOS workstation, the current wheel is a packaging smoke
artefact rather than a standalone workflow distribution, and public release is
blocked by the unresolved licensing and citation decisions described below.

## Stabilisation completed

- Resource paths are held in a portable, versioned profile and can be checked
  from the GUI or with `rna-mod-doctor` without requiring the four historical
  callsets.
- Dataset profiles can be added, edited, cloned, or removed without deleting
  their inputs or outputs. Paths inside the project are stored relatively.
- The GUI remains a thin launcher over the same 14 command-line workflows. It
  exposes relevant scientific options under a collapsed **Advanced settings**
  section, validates numeric values, rejects invalid design/context
  combinations, and shows the exact command.
- Output staging and manifest ownership protect complete previous results from
  failed or incompatible overwrites. The actual source checkout is protected
  regardless of the folder name used on another computer.
- Run manifests record input and output checksums, direct library versions,
  workspace/code Git state when available, and checksums of the three
  dependency-constraint files.
- A shared `quality_v1` policy now separates raw statistical calls from
  quality-qualified primary inference. Raw rows remain auditable. Partial
  inputs, boundary/non-estimable effects, invalid FDR values, low gene-cluster
  support, or observed leave-one-gene-out instability cannot silently enter a
  primary consensus.
- Cross-database, cross-modification, expression, KnockRBP, STRING, and dataset
  comparison outputs propagate those quality fields. KnockRBP recomputes
  consensus from the database-level evidence rather than trusting stale
  aggregate columns.
- The KnockRBP regulatory PNG no longer attempts to label every network edge.
  It shows a configurable, deterministic top-effect subset; the complete
  network and the exact plotted source rows are retained as TSV files.
- Scientific algorithms were separated from manifest and orchestration code so
  the maintained structural ceilings are restored.

## Current code inventory

| Measure | Verified value | Enforced ceiling |
|---|---:|---:|
| Runnable command wrappers | 14 | — |
| Package modules | 97 | — |
| Longest command wrapper | 298 lines | 300 |
| Longest package module | 500 lines | 500 |
| Longest production function | 236 lines | 250 outside transcript code |
| Distributed regression-test files | 8 | — |

The audit also caps top-level transcript functions and methods at 180 lines.
Strict matched-site analysis remains deliberately excluded from the release.

## Verification results

| Check | Final result |
|---|---|
| Unit and synthetic integration suite | **62/62 passed**, including GUI tests with headless Qt |
| Historical-project regression validator | **21/21 passed** |
| Portable full-resource doctor | **20 passed; 0 warnings; 0 skipped; 0 failed** |
| Command help smoke tests | **14/14 returned successfully** |
| Source GUI startup | 3 tabs, 14 workflows, checkout catalogue resolved |
| macOS helper syntax/permissions | all three scripts passed Bash parsing and are executable |
| Windows helper format | all three files are pure CRLF; control flow was statically reviewed |
| Wheel packaging smoke | version 1.0.1 wheel built and installed; all three console help commands loaded |
| Legacy transcript plotting self-test | 36/36 synthetic plots/source checks passed |
| Active overview PNG validation | 30/30 decoded successfully |
| Completed historical transcript-run validators | 4/4 passed |

This table records the 6 August stabilisation baseline. A 63rd regression test
was subsequently added for repository-name-independent transcript-region
catalogue discovery; it is part of the native laptop acceptance checklist.

The 62 tests cover Fisher tables and boundaries, Benjamini-Hochberg adjustment,
quality classification/propagation, output publication and rollback, protected
paths, manifest provenance, resource discovery, GUI forms and command
construction, dataset portability, cancellation cleanup, and the readable
KnockRBP plot selection.

The wheel check confirms package construction. It does not make the wheel a
standalone release: workflow dispatch still deliberately uses the checkout's
top-level `scripts/` and `config/rbp_catalog.tsv`. The supported distribution is
the complete checkout installed in editable mode.

## Real-project downstream smoke tests

All tests below read current project inputs and wrote to isolated temporary
directories. Their manifests validate with no checksum errors.

| Workflow | Verified result |
|---|---|
| Loose m6A cross-database comparison | 395 canonical RBPs; 27 primary/raw triple-significant; 17 triple direction-concordant |
| Transcript-region m6A cross-database comparison | 395 canonical RBPs; 23 primary/raw triple-significant; 15 triple direction-concordant |
| Transcript-region three-modification comparison | 1,890 rows: 1,636 primary robust, 177 sensitivity-only robustness warnings, 77 not inference-eligible |
| Transcript-region m6A KnockRBP | 10 datasets; 10 orthogonal rows; 508 regulatory edges; 81 primary-consensus edges |
| KnockRBP regulatory plot | 40 affected RBPs selected from the complete network; 64 plotted source rows retained |
| Transcript-region m6A STRING offline execution | 20 primary-selected RBPs; 12/12 requests served from cache; 20 outputs |
| Transcript-region m6A expression integration | 630 RBP rows, 77,114 gene rows, and 6 association rows |

The cross-modification refactor was checked against the immediately preceding
implementation: its four result tables and both plot-source tables were exactly
identical, including row order and values.

The loose cross-database audit illustrates why raw and primary calls are kept
separate. ENCORI IFIT2 and VIM have raw FDR below 0.05 but exact `OR = 0`;
they remain recorded as raw boundary calls and are excluded from primary
inference. Across that m6A ENCORI panel, 193 rows were raw significant and 191
were primary significant.

The offline STRING cache used for this software execution test contains
synthetic response fixtures. It validates selection, caching, parsing,
network/plot generation, and manifest handling; it must not be interpreted as
new biological STRING evidence.

## Legacy-result integrity

The refactor did not rewrite historical enrichment/statistical result tables.
A separate regression check found all 34 active loose/transcript enrichment
TSVs byte-identical to their tracked baseline. The 12 transcript plot-source
tables retain every pre-existing scientific column; only display/provenance
fields were added during the earlier overview-plot update.

The tracked legacy changes are confined to the intended all-RBP overview
plotting functions, regenerated overview PNGs, plot-source tables, and plot
manifests. No rerun of the historical enrichment calculations is required.

Two older downstream validator manifests under
`transcript_region_stratified_overlap/` still record the temporary recovery
   checkout path used when the project was restored. The corresponding current
files exist and match every recorded checksum: 44/44 KnockRBP inputs and 18/18
STRING inputs. This is stale path provenance, not detected data corruption, and
does not justify rerunning the analyses. It should be corrected only through a
deliberate provenance migration, not by editing scientific output tables.

## Scientific interpretation limits

- Loose enrichment is not matched for transcript, region, abundance, coverage,
  sequence, or RNA opportunity.
- Transcript-region analysis controls exact transcript and annotated region,
  but does not match coverage or abundance. Regional estimates remain
  descriptive without a formal interaction test.
- Missing robustness diagnostics are labelled and remain eligible under the
  core-only policy; their absence is not equivalent to demonstrated
  robustness.
- ENCORI and POSTAR3 aggregate experiments and cell contexts. Agreement among
  databases is not three independent MDA-MB-231 replications.
- oRNAment's downloaded genomic archive has an unresolved endpoint convention.
  The loose and transcript workflows preserve and record their established
  design-specific policies.
- Separate GENCODE v44 stop-codon bases follow the current region hierarchy and
  can be classified as 3′UTR. Stop-proximal interpretation must disclose this.
- KnockRBP changes may be indirect and depend on the DEG filter and gene
  universe. A network edge is not proof of direct RBP-to-RBP regulation.
- STRING connectivity does not prove a physical interaction in MDA-MB-231, and
  live API results are not release-pinned unless the response cache is retained.
- Cross-modification overlaps and correlations are descriptive; no permutation
  null or formal between-modification effect test is claimed.

## Distribution limitations and required decisions

1. `LICENSE` is intentionally a **no licence granted** placeholder. Obtain an
   authorised Queen's University Belfast/project-owner decision before any
   public GitHub release.
2. Replace the generic contributors and add approved ORCIDs, repository URL,
   and associated publication details in `CITATION.cff`.
3. Run the installer, doctor, GUI acceptance sequence, and one representative
   analysis on a clean native Windows machine before describing Windows as
   validated. Repeat for each macOS architecture intended for release.
4. Do not advertise a standalone `.exe`, `.app`, or wheel yet. Frozen builds
   must be created and tested separately for each operating system and
   architecture, while large reference/database resources remain external.
5. Preserve the exact prepared ENCORI/oRNAment panels and POSTAR3 snapshot used
   for this project. Their upstream portals do not provide one frozen archive
   that reconstructs every prepared resource file automatically.

## Release recommendation

Test a clean copy of the complete `refactored_code/` source on both the
maintainer's macOS and Windows laptops, following `README.md`,
`docs/RESOURCE_SETUP.md`, and `docs/GUI_GUIDE.md`. On each computer, create a
fresh environment with the provided installer, configure the external
resources, run `rna-mod-doctor --full`, and write the first analysis to a new
output directory. After both acceptance checks pass, publish the contents of
`refactored_code/` as the root of the lab repository at
`https://github.com/QUB-Simpson-lab/rna_mod_pipeline`. Use `rna-mod-audit` only
when the complete historical Queen's project layout is also present.
