# Laptop-test release readiness — 1.0.3

## Image-preview follow-up, 8 September 2026

Alexandru reports that the Mac trial worked. The subsequent image review found
that the empty text-preview pane could reduce a figure to a thin scrollable
strip. Version 1.0.3 corrects the **Results viewer only**: proportional Fit by
default, full preview space for images, 100% and zoom controls, drag/scroll,
opening the original file, and clear missing/corrupt-image messages. Old image
content is cleared when the selected manifest changes. Original figures,
scientific calculations, dependencies and install/run helpers are unchanged.

Verification of this patch:

- Full suite: **109 tests, 107 passed, two native-Windows tests skipped**;
  includes nine new image-preview tests.
- **510 existing project PNGs** decoded using the same Qt image engine as the
  app; zero failures. This is loading verification, not a scientific figure audit.
- **40 real-image layout checks**, covering wide and tall images, small windows
  and normal/Retina-scale rendering; **60 image/text/clear transitions** passed.
- Full application screenshots of four real figures in Fit and 100% modes
  generated. Representative screenshots were visually inspected. Fit preserves
  the whole image; zoom and scrolling are necessary to read dense/tall figures.
- Original source PNG checksums unchanged. The maintained production code now
  contains 120 Python files (103 package modules and 17 scripts).

These checks used local macOS Python 3.10.8 / Qt 6.11.1 with offscreen rendering;
they do not constitute native Windows acceptance. The new preview controls
should be briefly checked when reopening the updated app. Use the **1.0.3 ZIP**
for the next transfer/upload, not the superseded 1.0.2 archive. The verification
below records the broader 1.0.2 preparation and remains its historical baseline.

## Version 1.0.2 preparation baseline

Prepared on 8 September 2026. This is a **source release candidate for laptop
testing**, not a claim of completed Windows acceptance or a standalone app.
Start with [START_HERE.md](../START_HERE.md).

## Changes in this preparation

| Area | Change | Scientific impact |
|---|---|---|
| GUI profiles | New relative paths use portable `/` separators. Existing Windows-relative paths can be migrated when unambiguous; foreign absolute paths require re-selection. | No change to analysis calculations. External resources still need valid paths on each laptop. |
| macOS installation | Checks standard Homebrew and python.org Python locations when Finder supplies a restricted command path. | Installation only; the source folder remains separate from outputs. |
| Windows helper safety | Every helper stops if it cannot enter its own folder, before installing, running or resetting an environment. | Prevents operations in an unintended working directory; use a local extracted source folder. |
| Acceptance example | Added a small synthetic reference, annotation, complete-sized binding panels and two registered example datasets, generated on demand. | Tests software with invented data; does not validate biological hypotheses. |
| Example helpers | Added `test_demo_*` and `run_demo_*` helpers for both platforms. | No changes to the normal research-workspace launchers. |
| Resource inventory | Added a frozen SHA-256 inventory and a read-only verifier with selectable groups and path mappings. | Identifies existing resource bytes; does not change database contents or prove provenance. |
| Degenerate loose results | Corrected the reason/plot-footer label for all-bound tables and other boundary cases. Previously an all-bound table could be described as having zero overlaps. | Metadata and explanatory text only: overlap counts, ORs, p-values, FDRs and annotations are unchanged. Previously saved tables are not rewritten. |
| Continuous integration | Added synthetic numerical checks to the existing Linux, macOS and Windows matrix. | These remote jobs still need to run after upload. |
| Release documentation | Added a short start guide, laptop checklist, GitHub instructions and resource-snapshot guide; refreshed older installation and command documentation. | Historical audit reports are explicitly labelled as historical. |
| Metadata | Version 1.0.2; MIT licence; machine-readable citation with Alexandru Zob's documented name. | The lab should review contributor attribution. Database licences are separate. |

No legacy analysis script was edited in this preparation. The parent workspace's
`docs/project_context.md` was updated to distinguish current state from its
historical notes. Original experimental inputs and scientific result directories
were not regenerated, moved or deleted.

## Verification actually performed

| Check | Observed result | Scope and limitation |
|---|---|---|
| Fresh macOS installation | Passed with Python 3.11.6 on macOS 15.6.1, Apple Silicon, in a disposable folder containing spaces | Installer was also exercised with a restricted Finder-like command path; not a manual Finder double-click test. |
| Dependency consistency | Passed | Installed constrained runtime and GUI dependencies in the disposable environment. |
| Full regression suite | 100 tests: 98 passed, 2 skipped | The two skipped cases require native Windows. Qt checks ran offscreen. |
| Synthetic end-to-end workflow | All 10 workflow commands passed | Phase 1, three loose databases, three transcript-region databases, validation and both cross-database designs. |
| Independent expected answers | Passed | Exact panel identities, counts, ORs, p-values, FDRs, strata and consensus checked against independently defined expectations. |
| Synthetic figures | 66 explicitly required PNGs passed file/integrity checks | The final ENCORI overview was also visually inspected; this is not a full manual review of every figure. |
| Python-version comparison | Recorded numeric sentinel results agreed between local Python 3.10 and fresh 3.11 runs | Not a Windows comparison and not bitwise identity of every output. |
| Script entry points | All 17 scripts accepted `--help` | 14 analysis/GUI/audit entry scripts plus 3 new acceptance/resource utilities. |
| Structural audit | Passed for 119 production Python files | 102 package modules + 17 entry scripts; largest module 500 lines, largest entry script 298 lines. |
| Historical-result audit | 21/21 checks passed | Read-only inspection of existing project outputs; no full scientific rerun. |
| Existing-resource doctor | 20 passed; 0 failed, warned or skipped | Checks this workspace's available reference/binding resources, not a fresh download. |
| Frozen resource verification | 435 files, 10,039,405,035 bytes; zero mismatches | Transfer identity only; details in [RESOURCE_SNAPSHOT.md](RESOURCE_SNAPSHOT.md). |
| Packaging | Source distribution and wheel built successfully | The supported GUI distribution remains the complete source checkout, not the wheel alone. |

The synthetic example contains 800 raw rows, 384 retained high-call sites,
384 eligible controls, 288 DRACH cases, 32 genes and 96 transcript-region strata.
Its enriched and depleted sentinel ORs are exactly 9 and 1/9. It also includes
null, zero, infinite and non-estimable cases. Its 133/281/216 canonical panel
sizes exercise the production quality rules without relaxing those rules.

## What remains for the laptop trial

1. Run the install, test and example-GUI helpers on the actual Mac and Windows
   laptops. Record versions, results and any error messages in
   [CROSS_PLATFORM_ACCEPTANCE.md](CROSS_PLATFORM_ACCEPTANCE.md).
2. Manually check resource selection, saved-dataset reload, command previews,
   result-table/plot opening and cancellation of a disposable job. Automated
   offscreen Qt tests do not replace these checks.
3. Run a representative small real-data analysis with the correct real
   reference and annotations. Never analyse experimental reads using the
   synthetic resource profile. Confirm memory/disk requirements for full inputs.
4. Review MIT/contributor metadata with the lab and upload only the source
   contents, following [GITHUB_RELEASE.md](GITHUB_RELEASE.md). Check the remote
   CI results before tagging a release.

Live STRING responses, availability of public downloads and new expression or
KnockRBP input formats were not end-to-end revalidated in the synthetic example.
The regression suite and historical audit cover selected downstream behaviour;
they cannot guarantee every future dataset or service response. No result here
establishes biological validity, database independence or direct RBP recognition.

## Safety and recovery

A pre-change source backup is retained in the parent research workspace at
`release_preparation_20260908/source_before_changes.tar.gz`. Its SHA-256 is
`3cf3eb2f21b71a2717429de31dd81f9dca3efbb7aad71c0941ae3c0edce0a1b7`.
This is a source backup, not a replacement for backing up raw data and results.

The laptop-test ZIP excludes virtual environments, experimental data, generated
outputs and caches. The original `refactored_code/` remains uninstalled; test
installation work was isolated in a temporary copy. Keep test reports outside
the source upload and preserve the exact source version tested on both laptops.
