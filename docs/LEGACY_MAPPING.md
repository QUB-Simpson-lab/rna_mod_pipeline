# Legacy-to-refactored mapping

The refactor preserves the selected core scientific workflows while removing
duplicated configuration, plotting, interval, and statistical code.

| Legacy script or folder | Refactored command | Status |
|---|---|---|
| `scripts/filter.py` | `scripts/phase1_filter.py` | Included |
| `scripts/drach_validation.py` | `scripts/phase1_drach.py` | Included |
| `scripts/metagene_analysis.py` | `scripts/phase1_metagene.py` | Included |
| `scripts/rbp_overlap.py` | `scripts/loose_overlap.py --database ornament` | Included |
| `scripts/encori_overlap.py` | `scripts/loose_overlap.py --database encori` | Included |
| `scripts/postar3_overlap.py` | `scripts/loose_overlap.py --database postar3` | Included |
| `scripts/cross_database_validation.py` | `scripts/cross_database.py` | Included with clarified consensus outputs |
| `scripts/expression_integration.py` | `scripts/integrate_expression.py --source nanopore` | Included |
| `scripts/expression_integration_illumina.py` | `scripts/integrate_expression.py --source depmap` | Included |
| `comparison/scripts/compare_modifications.py` | `scripts/compare_modifications.py` | Included |
| `scripts/m6a_replicate_downstream_comparison.py` | `scripts/compare_datasets.py` | Generic site/enrichment robustness subset; monolithic legacy extras omitted |
| `scripts/knockrbp_validation.py` | `scripts/knockrbp_validation.py --design loose` | Included |
| `transcript_region_stratified_overlap/stratified_knockrbp_validation.py` | `scripts/knockrbp_validation.py --design transcript-region` | Included |
| `scripts/string_analysis_v2.py` | `scripts/string_analysis.py --design loose` | Included |
| `transcript_region_stratified_overlap/stratified_string_analysis.py` | `scripts/string_analysis.py --design transcript-region` | Included |
| `transcript_region_stratified_overlap/transcript_region_stratified_overlap.py` | `scripts/transcript_region_overlap.py` | Included |
| `transcript_region_stratified_overlap/stratified_plots.py` | `scripts/plot_enrichment.py` | Included |
| `transcript_region_stratified_overlap/regional_top_rbp_plots.py` | `scripts/plot_enrichment.py --regional-results PATH` | Included |
| `scripts/rbp_annotation.py` and `scripts/master_table.py` | none | RM2Target/WRE layer explicitly excluded |
| `scripts/rbp_enrichment_summary_plots.py` | core plot commands | Presentation-only duplication omitted |
| `scripts/generate_plot_4a.py` | Phase 1 plotting module | Removed as duplicate |
| `scripts/rebuild_ornament_alias_outputs.py` | Versioned RBP catalogue | Removed as one-off migration |
| source-specific database-comparison plot code | `scripts/cross_database.py` | Centralised |
| two separate expression implementations | one source-adapter workflow | Consolidated |
| two separate STRING implementations | one configurable workflow | Consolidated |
| loose and stratified KnockRBP data parsers | one hardened data layer | Consolidated |
| all strict/matched-site scripts and folders | none | Explicitly excluded |

The refactored loose and transcript-region all-RBP overviews both place finite
estimates from highest to lowest on a log2 odds-ratio scale. The loose plot
uses Fisher odds ratios; the transcript-region plot uses Mantel-Haenszel odds
ratios. Exact `OR = 0` results are shown as individual terminal purple
negative-infinity bars. Other boundary and non-estimable rows are counted on
the figure and retained in result/source tables.

The refactored Phase 1 filter also restores the complete eight-figure legacy
quality-control suite, including raw and filtered chromosome distributions,
separate marginal threshold-retention curves, the joint threshold heatmap, and
filtered-site diagnostic plots.

## Intentional scientific or engineering corrections

1. ENCORI background coverage is read from bedMethyl column 10
   (`Nvalid_cov`, zero-based index 9), rather than the BED score column.
2. Correct prefix-maximum interval search replaces bounded backward-search
   heuristics.
3. RBP aliases are represented once in a versioned catalogue.
4. Three-way significance and three-way direction agreement are separate
   outputs.
5. Non-significant, absent, and non-estimable results remain distinguishable.
6. STRING requests can be cached and replayed offline.
7. KnockRBP duplicate and direction-conflicting gene rows are handled
   deterministically.
8. Output manifests are portable across project locations.
9. Existing output directories are protected by default.
10. SF2 and SRSF1 source tracks are unioned and tested once as canonical SRSF1.
11. DRACH-unassessable reference contexts are missing values, not non-DRACH.
12. Input modification, design, database, and context labels are checked when
    present; unlabelled legacy inputs are marked unverified.
13. All-zero Fisher tables are reported as non-estimable (`OR=NaN`, `p=1`)
    rather than the legacy convention `OR=0`; FDR and significance are
    unchanged.
