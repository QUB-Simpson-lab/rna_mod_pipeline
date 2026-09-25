# Scope and exclusions

## Included

- Phase 1 filtering, DRACH annotation, and metagene annotation
- oRNAment, ENCORI, and POSTAR3 loose enrichment
- transcript-region-stratified enrichment
- loose per-site annotation TSVs, transcript-region case-binding matrices, and
  standard plots
- cross-database validation
- expression integration
- cross-modification comparison
- same-modification dataset/callset robustness comparison
- loose and transcript-region KnockRBP validation
- loose and transcript-region STRING analysis
- optional dataset-centred graphical launcher

Here, “Phase 2” means the three site-level binding-overlap workflows plus
cross-database validation. “Phase 3” means expression integration, KnockRBP,
and STRING. These are the reproducible core analyses selected for this release.

## Deliberately excluded

- strict/matched-site overlap and its context-design branches
- recovery and regeneration utilities
- historical one-off alias rebuilds
- duplicate plot-only scripts
- hypothesis and Word-document generators
- threshold-sensitivity orchestration
- the historical m6A-specific monolithic downstream-comparison script; the
  refactor instead includes a smaller generic site/enrichment robustness
  comparator
- target-gene and cancer-gene dossiers
- RM2Target writer/reader/eraser annotation, `rbp_annotation.py`, and
  `master_table.py`
- the legacy `rbp_enrichment_summary_plots.py` presentation-only summary
- TCGA survival analysis
- dedicated legacy POSTAR3 cross-cell-type stratification and comparison
  analyses (exact POSTAR3 cell-type and CLIP-method filtering within a loose or
  transcript-region run remains supported)
- extended legacy KnockRBP pathway, convergence, and formal cross-cell-line
  comparison analyses (the core workflow can still include datasets from
  multiple user-selected cell lines)
- data-download utilities

The excluded scientific analyses remain in the legacy project. They can be
refactored later without changing this core package.
