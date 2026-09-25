# Changelog

All notable release-code changes are recorded here. Scientific result changes
must also identify the affected workflow, inputs, parameters, and regenerated
manifests.

## 1.0.3 — Image-preview correction, 2026-09-08

- Gave figures the full Results preview area instead of sharing it with an
  unused text pane. Images fit proportionally on selection and window resize.
- Added Fit, 100%, zoom, drag/scroll and Open image controls without resaving
  or changing the original image. Missing/corrupt images show explicit errors;
  selecting another manifest clears stale previews.
- Added nine image-preview regression tests. No scientific engine, plot
  generation, dependency, installer or launcher changes in this patch.

## 1.0.2 — Laptop acceptance candidate, 2026-09-08

- Fixed cross-platform resource/dataset profile paths and added migration
  checks for older Windows-relative profiles.
- Added standard Homebrew/python.org fallbacks to macOS Python detection.
- Made Windows helpers stop before any environment operation if entering
  the source directory fails, with local-folder guidance and a regression check.
- Corrected loose-overlap reason labels and overview footnotes for all-bound
  and saturated contingency tables; counts, ORs, p-values and FDR are unchanged.
- Added a self-contained synthetic acceptance example with analytical expected
  counts/effects, verification scripts and double-click test/demo launchers.
- Recorded the exact 435-file external resource snapshot and a configurable
  standard-library checksum verifier; no research resources are bundled.
- Prepared MIT licensing, the documented software author's citation, a short
  start guide and explicit laptop/GitHub release instructions.
  Package metadata uses the supported SPDX licence format with setuptools>=77.
- Added reproducible demo checks to CI. Native Windows and manual laptop
  acceptance remain to be recorded; this is not a frozen `.exe` or `.app`.

### Included September preparation changes

- Made transcript-region catalogue discovery independent of the checkout name
  and added a regression test for a repository named `rna_mod_pipeline`.
- Added the lab repository metadata, cross-platform acceptance checklist, and
  standalone safeguards against committing external data or generated output.

## 1.0.1 — 2026-08-06

- Added bounded Python support and tested dependency constraints.
- Added installable `rna-mod-gui` and `rna-mod-audit` commands.
- Added the portable, resource-profile-aware `rna-mod-doctor` command.
- Added macOS, Windows, and Linux CI coverage for the compact regression suite.
- Made cleanup failures non-fatal after an output transaction has already been
  published; preserved artifacts are identified in a runtime warning.
- Added guarded environment-reset helpers that never target scientific data or
  output folders.
- Added contributor, citation, licensing-decision, testing, and resource-setup
  documentation.
- Added portable resource and dataset registries, GUI edit/clone/remove
  operations, side-effect-free validation, and safer cancellation handling.
- Added explicit evidence-quality fields so raw significance remains auditable
  while unstable or incomplete estimates are excluded from primary consensus
  selection.
- Recomputed KnockRBP consensus from database-level evidence and replaced the
  unreadable full-network heatmap with a configurable top-effect view whose
  exact source rows are retained separately.
- Extended manifests with code state, plotting/network/reference-library
  versions, and checksums of available dependency constraint files.

## 1.0.0 — Initial modular release

- Initial modular implementation of Phase 1, loose and transcript-region overlap,
  cross-database, expression, KnockRBP, STRING, cross-modification, plotting,
  dataset comparison, and optional GUI workflows.
