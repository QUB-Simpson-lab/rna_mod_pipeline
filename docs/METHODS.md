# Methods implemented by the refactored pipeline

## Site definition

ModKit bedMethyl files are read using the standard 18-column schema. A high-call
site satisfies:

```text
Nvalid_cov >= minimum coverage
fraction_modified >= minimum modification percentage
```

The established analysis uses coverage 20 and modification percentage 20.
`fraction_modified` is on a 0–100 percentage scale.
The selected analysis is checked against the bedMethyl `mod_code`; unexpected
or mixed modification codes fail rather than being silently discarded.

Each filtering run produces the complete Phase 1 quality-control suite:
coverage and fraction distributions before filtering, chromosome distributions
before and after filtering, separate marginal threshold-retention curves,
a joint coverage×fraction threshold heatmap, the filtered fraction
distribution, and the filtered coverage-versus-fraction scatterplot. The
marginal threshold plot changes one threshold at a time; the heatmap applies
both thresholds jointly.

For m6A, DRACH validation extracts a strand-corrected central five-mer from
hg38. It adds the motif and a Boolean DRACH annotation but does not remove
non-DRACH sites. Sites lacking an extractable five-base reference context are
labelled `unavailable_reference_context`, with nullable `is_DRACH`, and are not
counted as non-DRACH.

## Transcript annotation

GENCODE v44 Basic protein-coding transcripts with complete 5′UTR, CDS, and
3′UTR components are eligible. When several transcripts contain one site, the
longest eligible transcript is selected deterministically. Unmapped sites remain
available to the loose analysis but cannot enter transcript-region strata.
Equal-length ties use lexicographically smallest transcript ID and then region,
the same rule in Phase 1 and transcript-region analysis.

GENCODE v44 represents some stop-codon bases separately from CDS features while
also placing them inside UTR features. Under the implemented feature hierarchy,
those bases are therefore labelled 3′UTR. This affects 203 m6A sites, 11
pseudouridine sites, and no m5C sites in the current project. Cases and controls
use the same convention, but stop-proximal and 3′UTR interpretations should
acknowledge it.

## Loose RBP enrichment

For one modification and one binding database:

1. Cases are all filtered high-call sites.
2. Controls are all bedMethyl positions with adequate coverage and a
   modification percentage below the high-call threshold.
3. An RBP binds a site when any database interval overlaps the half-open window
   `[position-window, position+window+1)`.
4. Strand is ignored, matching the established analysis.
5. Binding is binary per site.
6. A two-sided Fisher exact test compares bound/unbound cases with
   bound/unbound controls.
7. Benjamini-Hochberg correction is applied across the declared RBP panel for
   that database and modification.

The default window is ±10 nucleotides.
The supplied case table must, by default, contain every raw bedMethyl call that
meets the case thresholds. This bidirectional completeness check detects a
wrong raw file or accidental truncation. `--allow-site-subset` permits an
intentional subset while retaining that provenance in the manifest. Cases must
still be genuine calls from the selected raw bedMethyl and meet the declared
thresholds. Case and background fraction definitions cannot overlap.

The loose all-RBP overview displays finite odds ratios in descending order on
a `log2(OR)` axis, with `log2(OR) = 0` as the reference. Exact `OR = 0`
results are shown as individual terminal purple negative-infinity bars rather
than assigned a finite effect. Positive-infinity and non-estimable tests are
counted in the figure note. Exact estimates and all result states remain in
the enrichment TSV.

For oRNAment, the downloaded
[genome-browser BED archive](https://rnabiology.ircm.qc.ca/oRNAment/downloads/)
is dominated by records with `end - start = 6`, consistent with its documented
[seven-nucleotide motifs](https://rnabiology.ircm.qc.ca/oRNAment/tutorial/) if
endpoints are inclusive, but the archive does not state the endpoint conversion
unambiguously. It also contains a small number of `start == end` records,
compatible with motif fragments split at exon junctions but not definitively
explained by the download documentation. The loose workflow preserves the
established interpretation: supplied coordinates are used directly and
zero-span records are retained as point-like records. This policy is recorded
in the run manifest.

## Transcript-region RBP enrichment

The same high- and low-call definitions are used, but every opportunity is
assigned to a selected transcript and region. Cases and controls contribute only
inside exact `transcript_id × region` strata that contain both groups.

The workflow reports:

- pooled Fisher odds ratio for description;
- Mantel-Haenszel common odds ratio;
- conventional CMH test;
- gene-cluster robust variance and Student-t inference;
- leave-one-gene-out influence diagnostics;
- overall and region-specific FDR.

The all-RBP overview displays finite, positive, non-boundary
Mantel-Haenszel estimates in descending order on a `log2(OR)` axis, with
`log2(OR) = 0` as the reference line. Exact `OR = 0` boundaries are shown as
individual terminal purple negative-infinity bars; the bars are descriptive
and do not give them a finite effect size. Positive-infinity
and non-estimable tests are omitted from the quantitative axis and counted on
the figure. All rows remain in the source TSV with an explicit display flag
and exclusion reason.

The design uses every eligible control. It does not match or subsample controls
by coverage, central sequence, 3′ distance, or case/control ratio. This is
therefore stratification, not the excluded strict matched-site design.

Transcript opportunities also require bidirectional agreement between the site
table and raw bedMethyl. Mixed modification codes, conflicting duplicates, or
an incomplete production case table fail. Coverage is not matched, so the
reported coverage-balance table and residual abundance/coverage confounding
must be considered when interpreting results.

The transcript-region resource loader preserves its established half-open BED
interpretation and excludes `end <= start` records, including oRNAment
zero-span rows. A project-wide sensitivity check that added those rows as
point-like intervals produced only two borderline FDR transitions (m6A RBM22
and m5C RBFOX2); the principal FXR2 and FMR1 effects were essentially
unchanged. This check does not resolve whether all oRNAment endpoints should be
converted from an inclusive convention. Consequently, loose and
transcript-region oRNAment opportunities are not strictly identical and should
not be treated as a clean design-only comparison. Harmonising them would
require an authoritative coordinate definition (or transcript-level source
data) followed by rerunning every oRNAment analysis.

## Cross-database validation

RBP aliases are canonicalised before database tables are joined. Pairwise
outputs distinguish absent coverage from a tested non-significant result.
Three-way outputs distinguish:

- significant in all three databases;
- significant and directionally concordant in all three;
- significant but directionally conflicting.

Arithmetic mean odds ratio is retained for legacy comparability. Mean log odds
ratio and its geometric transform are also reported because odds ratios are
multiplicative.

Cross-database figures use the same result states as the tables. Pairwise
scatterplots include only shared, finite, positive odds-ratio estimates, while
their source table retains excluded boundary and non-estimable rows. Agreement
plots report database coverage separately from inferential agreement among
shared RBPs. Exact database-membership counts are used instead of
area-proportional Venn diagrams. Top/shared and triple-validation panels show
odds ratios without confidence intervals. Every figure has a TSV source table
and a plot manifest.

Agreement is evidence of robustness across data resources, not three
independent experiments: oRNAment is predictive, ENCORI and POSTAR3 aggregate
CLIP records, and underlying experiments or biological contexts may overlap.

## Expression integration

The Nanopore profile classifies a gene as detected when its raw count is
greater than zero and reports normalised CPM. A zero count is treated as
low-depth/uninformative, not proof of absence. DepMap values are transformed
from `log2(TPM+1)` to TPM; the default MDA-MB-231 model is `ACH-000768`, and
TPM ≥1 is classified as expressed.

RBP names use the shared canonical alias map; historical construct labels with
no defensible human-gene mapping remain explicitly unmappable. The primary
legacy-compatible gene-level estimand includes genes with at least one
modification site that are detected/expressed. It correlates
`log2(expression+1)` with modification-site count using Spearman and Pearson
tests. Two explicitly exploratory estimands are also retained:

- all genes in the expression source, assigning zero sites to unmodified genes;
- all modified genes, irrespective of detection/expression status.

When gene lengths are supplied, sites per kilobase is an additional metric.
Benjamini-Hochberg correction is applied across all correlations in one run.
These are gene-level associations and do not demonstrate that expression
causes modification density.

## Evidence-quality policy

Every enrichment estimate retains its raw statistical direction and
significance. Primary significance additionally requires a complete input when
completeness is reported, a finite positive non-boundary odds ratio, a valid
FDR value, at least 20 gene clusters when that diagnostic is available, and no
observed leave-one-gene-out direction reversal or boundary. Missing robustness
diagnostics are labelled and treated as `primary_eligible_core_only`, not as
proof of robustness. Boundary, partial-input, or robustness-warning calls
remain available for sensitivity auditing rather than being deleted.

Cross-database aggregate fields are recomputed from these database-level raw
and primary calls. Missing database coverage remains `not_covered`; it is not
converted into a null or non-significant result.

## KnockRBP

Differential-expression rows pass `|log2FC| > 0.5`; datasets with adjusted
p-values additionally require `padj < 0.05`. Same-direction duplicate symbols
are resolved deterministically. Symbols reported in opposite directions within
one dataset are excluded as ambiguous.

Target overlap uses a one-sided Fisher exact test. Transcript-region analysis
uses the context-specific modified-gene universe as its primary background and
retains the historical 20,000-gene background as sensitivity analysis.
Multiple comparisons are corrected and small target/overlap counts are flagged.
Regulatory-network edges describe perturbation-associated expression changes,
not direct RBP-to-RBP regulation.

The complete regulatory edge table is retained. For readability, the heatmap
shows a configurable number of affected RBPs ranked by their maximum absolute
retained-DEG log2 fold change. Missing heatmap cells represent no DEG passing
the declared filters, not an observed zero. The exact plotted rows and ranks
are stored in a separate source-data TSV.

The binding-target union uses ENCORI and POSTAR3. Resource availability is
derived from enrichment/resource-loading tables when available; a declared row
whose usable interval count is zero is not treated as experimental coverage.
If a legacy annotated table has no resource-status companion, availability is
labelled annotation-only and unverified. Multiple KnockRBP datasets for one RBP
remain separate by dataset ID and cell line.

## STRING

RBP sets are selected by explicit significance, direction, database, and design
rules. STRING requests are cached by request hash. Medium-confidence
interactions use score 400 and high-confidence interactions use score 700.
Response payloads and software/input provenance are retained. STRING evidence is
functional context and does not demonstrate an interaction in MDA-MB-231.
Network plots include every resolved selected protein; resolved proteins with
no edge at the selected confidence threshold are shown as grey isolates.
STRING’s API does not expose a pinned database release here, so the cached
responses and their hashes define the reproducible result.

## Cross-modification comparison

For each database, the workflow preserves tested-significant,
tested-non-significant, tested-non-estimable, and not-covered states separately.
It classifies modification-specific, shared, and direction-switching RBP
patterns. Pearson and Spearman correlations use shared finite log2 odds ratios;
p-values receive BH correction both within database and across the complete
run. Significant-set intersections and Jaccard indices are descriptive and do
not by themselves provide a permutation p-value.

Phase-1 site-count and metagene plots are included when metagene inputs are
available. Use `--skip-phase1-context` to omit this layer entirely. Otherwise,
every requested modification requires a metagene table unless
`--allow-missing` is explicit. The manifest records included and missing
Phase-1 inputs; missing database coverage is never converted into a null
result.

## Multiple-testing scopes and interpretation

Loose FDR is within one modification×database panel. Transcript overall FDR is
within database; regional FDR is within database×region. Cross-modification
correlation FDR is reported both within database and across the run. Regional
effect comparisons are descriptive unless a formal region-interaction model
and a global region-aware correction are added.
