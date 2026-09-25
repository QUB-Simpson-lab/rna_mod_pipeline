# Frozen resource inventory

`config/resource_snapshot_20260908.json` records the exact external files present
in the research workspace on **8 September 2026**. It contains workspace-relative
paths, byte sizes, SHA-256 checksums and source-track aliases. It contains no
reference sequence, binding records, raw bedMethyl data or patient-level data.
The large resources remain outside this code repository.

The measurement date is not a download date. Original download dates and
unrecorded database release labels are **not available in the inspected project
files**. Source URLs were copied from the existing resource instructions, not
used to infer a version. Empty URL lists mean that a confirmed URL was not
recorded for that group; they do not imply a resource is unavailable online.

| Group | Files | Exact bytes | Contents |
|---|---:|---:|---|
| `reference` | 3 | 4,842,338,003 | hg38 FASTA, its index, GENCODE v44 GTF |
| `ornament` | 134 | 346,747,671 | Catalogue-selected source BED tracks; 133 canonical RBPs |
| `encori` | 281 | 67,885,101 | Prepared project CLIP-track panel |
| `postar3` | 1 | 4,470,359,284 | Complete local `human.txt` snapshot |
| `knockrbp` | 14 | 4,801,614 | Metadata and all 13 metadata-resolved DEG JSON files |
| `expression` | 2 | 307,273,362 | In-house normalised counts and DepMap/CCLE expression table |
| **Total** | **435** | **10,039,405,035** | Approximately 9.35 GiB of external resources |

KnockRBP records identify cell line, dataset and whether a dataset is selected by
the default MDA-MB-231/MDA-MB-231-LM2 filter: 10 are selected by default; the
three other TNBC datasets remain available for explicit selection. This inventory
does not alter those filters. Splicing, editing, APA and Enrichr JSONs are not
inputs to the refactored expression-DEG validation, so they are not included.

The catalogue checksum and version are recorded. oRNAment SF2 and SRSF1 files
are both included and identify the same canonical SRSF1 result. Non-human
oRNAment tracks and `.tbi` indexes are not used by the pipeline and are omitted.
For BED sources, inventory creation follows the pipeline's preference for
`.bed.gz` over `.bed` when both exist.

## Verify after copying resources to another laptop

Run these commands from the code repository. The verifier needs only Python's
standard library, so pipeline dependencies need not be installed first.
`--workspace` means the directory containing `data/`.

macOS:

```bash
python3 scripts/verify_resources.py --workspace /path/to/rna_mod_workspace
```

Windows:

```powershell
py -3.11 scripts\verify_resources.py --workspace "D:\rna_mod_workspace"
```

A full verification streams about 10 GB from disk without copying resources.
It prints progress per group. The successful final line is:

```text
Verified 435 files; 0 failure(s). Only selected inventory files were checked.
```

You can verify just the resources required for your intended analysis:

```bash
python3 scripts/verify_resources.py --workspace /path/to/workspace --groups reference ornament
python3 scripts/verify_resources.py --workspace /path/to/workspace --groups encori
```

Other available groups are `postar3`, `knockrbp` and `expression`. Selecting one
does not require the others to exist. The report covers selected groups only;
it is not a complete installation acceptance test.

## Resources stored in different locations

Use repeated `--map` arguments to map a recorded file or directory to its actual
location. Quote the complete mapping when it contains spaces. Directory mappings
retain filenames under that directory; exact-file mappings can rename a file.

```powershell
py -3.11 scripts\verify_resources.py --workspace "D:\workspace" --groups ornament encori --map "data/HS=D:\RBP Resources\HS" --map "data/ENCORI=E:\ENCORI"
```

```bash
python3 scripts/verify_resources.py --workspace /path/to/workspace --groups reference --map "data/hg38.fa=/Volumes/References/hg38.fa" --map "data/hg38.fa.fai=/Volumes/References/hg38.fa.fai" --map "data/gencode.v44.annotation.gtf=/Volumes/References/gencode.v44.annotation.gtf"
```

An unrecognised mapping fails instead of silently checking a different directory.
These mappings affect this verification command only; select the same actual
locations in the GUI resource settings separately.

## Interpret the result

- Exit status `0`: all selected inventory files and the catalogue match.
- Exit status `1`: one or more files are missing or differ in size/checksum.
- Exit status `2`: malformed inventory, invalid options, or another input error.

A mismatch means **different bytes**, not automatically invalid biological data.
Recompressing a BED, rewrapping a FASTA, rebuilding an index, or downloading a new
database can change the checksum while retaining some or all biological content.
Do not silently overwrite the frozen inventory to make verification pass.
Compare the source and file format, run the resource doctor, and label a new
resource build if the change was deliberate.

Matching checksums establish transfer integrity; they do not establish source
authenticity, biological validity, database independence or MDA-MB-231 specificity.
ENCORI and POSTAR3 tracks aggregate multiple experimental contexts. A matching
panel is not evidence of binding in the project cell line itself.

The verifier checks listed files, not every unrelated file in the resource tree.
Continue to use `rna-mod-doctor --project-root /path/to/workspace` for structural
and catalogue checks, and the acceptance example for pipeline behaviour.

## Preserve the exact ENCORI panel

The snapshot identifies all 281 filenames and their checksums. It does not make
this prepared panel downloadable as a single provider archive. Transfer the
validated `data/ENCORI/` directory through the lab's approved storage, then run
the `--groups encori` verification. Keep download/source metadata alongside the
resources. Redistributing the code under MIT does not grant redistribution rights
for third-party databases or in-house expression data.

## Create a separate inventory for a future resource build

Creation expects the conventional workspace layout and the complete set of
resource groups above. It reads the catalogue and KnockRBP metadata to select
files. It never modifies resources or overwrites an existing snapshot.

```bash
python3 scripts/verify_resources.py --workspace /path/to/new_workspace --create --output /path/to/new_resource_snapshot.json
python3 scripts/verify_resources.py --workspace /path/to/new_workspace --snapshot /path/to/new_resource_snapshot.json
```

Use `--catalog /path/to/rbp_catalog.tsv` for a deliberately revised catalogue.
Creating a new inventory records what is present; it does not validate a new
database build. Preserve provenance and complete the appropriate audits first.
