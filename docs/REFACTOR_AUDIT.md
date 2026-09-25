# Refactor audit

> Historical refactor baseline. See `RELEASE_READINESS.md` for the current
> release preparation, updated code/test counts and licensing status. Numerical
> historical-project checks below are preserved as their original audit record.

## Audit conclusion

The handoff implementation covers the requested core workflow and is ready for
authorised internal use. It keeps each modification and binding database selectable at
the command line, separates scientific algorithms from command wrappers, and
does not import or alter the legacy scripts.

The legacy result folders do not need to be rerun because of this refactor.
New analyses should be written below `refactored_outputs/` or another new,
explicit output path. The raw-data trees, legacy-result trees, and
`refactored_code/` itself are protected from accidental output replacement.

## Audited scope

| Workflow | Command | Independently selectable unit |
|---|---|---|
| Phase 1 filtering | `phase1_filter.py` | one modification/replicate |
| m6A DRACH annotation | `phase1_drach.py` | m6A or m6A replicate 2 |
| transcript/metagene annotation | `phase1_metagene.py` | one modification/replicate |
| loose overlap | `loose_overlap.py` | one modification × one database |
| transcript-region overlap | `transcript_region_overlap.py` | one modification × any database subset |
| plotting | `plot_enrichment.py` | one overall result table, with an optional regional table |
| cross-database comparison and plots | `cross_database.py` | one modification/design/context × selected databases |
| expression integration | `integrate_expression.py` | one modification/design/context/database/source |
| KnockRBP | `knockrbp_validation.py` | one modification/design/context × selected datasets/cell lines |
| STRING | `string_analysis.py` | one modification/design/context × selected directions |
| cross-modification comparison | `compare_modifications.py` | any selected modifications × databases |
| same-modification dataset comparison | `compare_datasets.py` | two independently processed datasets with the same chemistry |
| desktop launcher | `pipeline_gui.py` | any registered bedMethyl dataset × selected workflow |

The strict matched-site workflow is not included. Other deliberate exclusions
are listed in `SCOPE.md`.

## Code-structure audit

- Fourteen small command wrappers expose the workflows; the longest is 298
  lines.
- Repeated interval, statistical, plotting, provenance, identity-validation,
  and alias logic is implemented once in the package.
- The 97 package modules are capped at 500 lines and command wrappers at 300.
  Transcript-region functions are additionally capped at 180 lines, and all
  other production functions at 250 lines. The pre-handoff verification suite
  checked these limits.
- Source files contain only short comments and focused docstrings; the detailed
  scientific explanation is kept in `METHODS.md`.
- The fixed RBP catalogue contains 133 oRNAment, 281 ENCORI, and 216 POSTAR3
  canonical tests.
- SF2 and SRSF1 source tracks are unioned and tested once as canonical SRSF1.
- Multiple KnockRBP datasets for the same RBP remain distinct by dataset ID and
  cell line.
- Cross-database plots live under each comparison output’s `plots/` directory,
  retain source-data TSVs, and use exact membership counts rather than a Venn
  dependency.
- Loose and transcript-region all-RBP overview plots use descending log2
  odds-ratio axes. Exact `OR = 0` results are represented by individual
  purple negative-infinity bars; other boundary and non-estimable tests are
  counted on the PNG and retained in the enrichment or plot-source table.
- Phase 1 filtering publishes the full eight-figure QC contract used by the
  established workflow, including raw/filtered chromosome distributions,
  marginal and joint threshold diagnostics, and filtered-site plots.

## Verification performed

### Automated checks

| Check | Result |
|---|---|
| Unit and synthetic integration suite | 62/62 passed; the compact suite remains in `tests/` |
| Python compilation | all package and command files passed |
| Command help smoke tests | all 14 commands returned valid help |
| Handoff validator | 21/21 checks passed |
| Legacy statistical engines/enrichment tables changed | none |
| Legacy transcript-region plotting changed | overview implementation and plot artefacts updated to match the refactor |

The pre-handoff verification covered threshold boundaries, modification-code mismatches,
interval-window boundaries, alias unions, missing/empty resources, complete
case-set enforcement, transcript assignment, Mantel-Haenszel calculations,
gene-cluster inference, FDR scopes, database absence versus null results,
replicate identity, repeated KnockRBP datasets, expression estimands, STRING
offline/cache behaviour, plotting edge cases, output-path protection, and
failed-overwrite preservation. It also covers arbitrary dataset registration,
bedMethyl inspection, shell-free GUI command construction, preflight failures,
process cancellation, manifest-gated completion, same-modification comparison,
and transcript-region enrichment discovery.

### Real-project checks

| Component | Observed result |
|---|---|
| Phase 1 default filtered counts | m6A 78,811; m5C 67,986; Ψ 22,930; m6A replicate 2 70,135 |
| m6A DRACH table | 78,811 retained; 58,199 DRACH-positive |
| m6A metagene table | 78,811 retained; 70,407 mapped |
| Loose m6A oRNAment parity | all 133 RBP overlap counts, odds ratios, p-values, and FDR values matched the established result to numerical tolerance |
| Transcript result panels | 630 rows per modification: 133 oRNAment + 281 ENCORI + 216 POSTAR3 |
| Loose m6A cross-database smoke run | 395 canonical RBPs; 27 significant in all three; 17 also direction-concordant |
| Transcript-region m6A cross-database smoke run | 395 canonical RBPs; 23 significant in all three; 15 also direction-concordant |
| Real loose m6A KnockRBP smoke run | 10 eligible datasets retained; 10 orthogonal rows; 508 perturbation-associated RBP-network edges |
| Transcript-region KnockRBP smoke run | m6A 3′UTR run completed for an explicitly selected YTHDF2 dataset with matching regional cross-database input |
| Real three-modification comparison | 1,890 database-specific matrix rows; 630 RBP-pattern rows; 18 effect-correlation rows; 9 significant-overlap rows |
| Primary in-house expression parity | 2,699 detected modified genes; Spearman ρ = -0.061833, p = 0.001309 |
| Primary DepMap expression parity | 8,399 expressed modified genes; Spearman ρ = -0.271182, p = 1.69 × 10^-141 |

The real smoke runs used existing project inputs but wrote only to temporary
directories outside the project.

## Intentional corrections relative to legacy code

1. ENCORI background coverage uses bedMethyl `Nvalid_cov`, not the BED score.
2. Interval lookup uses a complete prefix-maximum search, avoiding missed long
   intervals.
3. Modification codes are checked against the selected analysis.
4. m6A DRACH annotation retains all calls and distinguishes unavailable
   reference context from a genuine non-DRACH motif.
5. Transcript ties are deterministic: longest complete transcript, then
   lexicographically smallest transcript ID and region.
6. A production loose or transcript-region run verifies that its case table
   contains the complete qualifying raw-call set. Intentional subsets require
   an explicit flag.
7. Cross-database tables distinguish not covered, tested non-significant,
   non-estimable, significant concordant, and significant conflicting states.
8. Expression aliases and primary versus exploratory estimands are explicit.
9. KnockRBP resource availability is based on usable resource evidence, not
   merely the presence of a binary annotation column.
10. Repeated KnockRBP datasets cannot silently replace one another.
11. POSTAR3 comma- or semicolon-separated accessions are counted as individual
    experiments.
12. STRING plots retain resolved proteins with no qualifying edge as isolates.
13. All workflows stage a complete run before publication. Directory workflows
    require a matching directory manifest; file workflows replace only files
    owned by their matching run manifest and reject reruns that would strand
    previously owned plots or tables.
14. Run manifests contain checksums, package versions, and Git state, and store
    final rather than temporary output paths.
15. A custom loose-output directory cannot silently accumulate results from
    different modifications or databases.
16. Database-specific source filters fail when their database is not selected,
    rather than being silently ignored.

## Expected numerical differences

The refactor is not intended to preserve a known error merely for byte-for-byte
parity.

- Deterministic transcript tie-breaking changes 387 transcript assignments in
  the m6A table. Only 11 gene labels and 6 region labels change; the mapped-site
  count remains 70,407. There are 124 materially different metagene coordinates
  at an absolute tolerance of `1e-12`. These differences should be described as
  a reproducibility correction.
- A 2×2 table with no binding in either group is represented as non-estimable
  (`odds_ratio = NaN`, `p = 1`) rather than the legacy placeholder
  `odds_ratio = 0`. Significance and FDR interpretation are unchanged.
- Shared modification-related plot flags now use one registry across all
  database adapters. These flags affect display annotation only, not overlap
  statistics.
- The KnockRBP default retains every eligible dataset, so row/edge counts can
  exceed an older result that kept only one dataset per RBP.

## Remaining scientific limitations

- Loose enrichment is vulnerable to transcript, region, expression, coverage,
  and RNA-opportunity confounding.
- Transcript-region stratification controls exact transcript and region but
  does not match coverage or abundance. Its coverage-balance output must be
  reviewed.
- Region-specific estimates are descriptive comparisons; a formal
  RBP-by-region interaction model is not implemented.
- ENCORI and POSTAR3 aggregate experiments and cell contexts. Cross-database
  agreement is not equivalent to independent MDA-MB-231 replication.
- A legacy input without embedded analysis labels can be checked
  schema-wise but remains explicitly `unverified` for labels in the manifest.
- m6A and its second replicate share the same modification code, so replicate
  identity must come from the selected path/run manifest, not `mod_code`.
- KnockRBP target-overlap tests can reflect indirect perturbation effects and
  depend on the chosen gene universe.
- STRING results depend on the live API unless the recorded cache is retained;
  connectivity is not evidence of interaction in MDA-MB-231.
- Cross-modification set overlaps are descriptive; no permutation null is
  claimed.
- Expression-density analyses remain susceptible to transcript length and
  detection-depth effects unless a gene-length table is supplied.
- oRNAment’s downloaded genomic archive has an unresolved endpoint convention.
  Loose overlap retains `start == end` records under the historical
  interpretation, whereas transcript-region overlap excludes `end <= start`.
  A zero-span sensitivity analysis found only two borderline FDR transitions
  and no material change to the principal FXR2/FMR1 effects, but this does not
  resolve a possible general inclusive-end conversion. Cross-design oRNAment
  comparisons therefore require this caveat.
- Separate GENCODE v44 stop-codon bases are classified as 3′UTR by the current
  feature hierarchy (203 m6A, 11 pseudouridine, and 0 m5C sites). Cases and
  controls share the convention, but stop-proximal interpretations should
  disclose it.
- KnockRBP target and DEG symbols are compared exactly after limited RBP alias
  canonicalisation. Historical aliases outside that map may be missed.
- Dataset comparison uses exact genomic site keys and common available
  enrichment panels. It is a callset-robustness analysis, not raw-read pooling
  or independent replication of the reused RBP-binding resources.
- The optional launcher runs one job at a time and provides rerun, not
  checkpoint/resume. Its fast completion check verifies output presence and
  recorded size; the command manifests retain the stronger SHA-256 provenance.

## Handoff decision

No blocking code, configuration, or result-identity defect was found in the
exercised internal-handoff scope. The pipeline is suitable for postdoctoral
review and for new runs from the complete source checkout. Native Windows
acceptance, a licensing decision, approved citation metadata, and separate
standalone-app packaging remain required before a public release. Scientific
claims should continue to use the limitations above and the analysis-specific
manifests.

The final post-change verification, size inventory, and GUI pilot decision
are recorded in `FINAL_HANDOFF_AUDIT.md`.
