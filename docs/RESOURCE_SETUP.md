# External resource setup

The code repository intentionally excludes large reference and binding files.
Keep them in the workspace `data/` tree; do not move them into the code
repository. Record the download date, source URL, release/build, file
size, and SHA-256 checksum before analysis.

For the first installation test, use the self-contained synthetic example in
`START_HERE.md`; none of the large downloads below is needed for that test.
For the existing research resources, the measured inventory and verification
commands are in [`RESOURCE_SNAPSHOT.md`](RESOURCE_SNAPSHOT.md). It covers 435
files and records measured hashes, not assumed download dates.

## Required layout

```text
rna_mod_workspace/
├── data/
│   ├── hg38.fa
│   ├── hg38.fa.fai
│   ├── gencode.v44.annotation.gtf
│   ├── human.txt
│   ├── HS/                       # oRNAment BED files
│   ├── ENCORI/                   # ENCORI BED files
│   ├── knockrbp/
│   │   ├── dataset_metadata.tsv
│   │   └── *_degs.json
├── refactored_outputs/          # generated results and saved GUI profiles
└── rna_mod_pipeline/            # this GitHub repository
```

The GUI can point to resources elsewhere, but the command examples use this
layout. The launch helpers treat the repository's parent as the workspace.
RM2Target, cancer-gene and survival resources are used by additional legacy
analyses and are not required by the maintained core described here.

## Reference genome and annotation

- Download UCSC `hg38.fa.gz` from
  <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz>.
- Download the GENCODE release 44 **CHR comprehensive gene annotation**,
  `gencode.v44.annotation.gtf.gz`, from
  <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz>.
  The CHR/ALL/PRI and Basic/Comprehensive alternatives are not interchangeable;
  this project used the file named above.
- Decompress both files. The expected names are `hg38.fa` and
  `gencode.v44.annotation.gtf`.
- Create the uncompressed FASTA index with `samtools faidx data/hg38.fa`, or
  open the FASTA once with pyfaidx. Confirm that `data/hg38.fa.fai` exists.

Do not mix `hg38`/GRCh38 coordinates with hg19/GRCh37 binding files. Record
whether chromosome names use `chr1` or `1`; adapters normalise only the formats
they explicitly support.

## RBP binding resources

### oRNAment

Use the official download page at
<https://rnabiology.ircm.qc.ca/oRNAment/downloads/> and select
`Homo_sapiens_oRNAment.bed.tar.gz`. The provider reports MD5
`58299ec277d8a0316324811892002f08` for the compressed archive. Extract the
human per-RBP BED resources into `data/HS/`. The frozen panel and aliases are
defined in `rna_mod_pipeline/config/rbp_catalog.tsv`; do not infer the tested
panel from whichever files happen to be present. The doctor currently expects
134 source BED files resolving to 133 canonical oRNAment RBPs because one
catalogue entry combines source names.

### ENCORI

Use the official ENCORI CLIP peak download interface at
<https://rnasysu.com/encori/download.php>, not the pan-cancer analysis page.
The pipeline input is the project's prepared 281-track panel: each
`source_names` value in `config/rbp_catalog.tsv` must resolve to
`<source_name>.bed` or `<source_name>.bed.gz` in `data/ENCORI/`. Each record
must contain at least chromosome, zero-based start, and end; when
`--minimum-clip-experiments` is greater than one, BED column 5 must contain the
support count used for filtering. Preserve the ENCORI dataset/cell-line/source
metadata separately because these tracks aggregate experiments from multiple
contexts. The standard overlap is not an MDA-MB-231-only analysis unless a
separate, explicitly filtered panel is constructed and labelled.

ENCORI does not distribute this exact frozen 281-file handoff panel as one
versioned archive. For internal transfer, preserve the already validated
`data/ENCORI/` directory with `config/resource_snapshot_20260908.json`. Check a
transferred copy from the code folder with
`python scripts/verify_resources.py --workspace .. --groups encori`.
A future re-download is
a new resource build and must be converted, audited against the catalogue, and
recorded rather than assumed identical to the current analysis.

### POSTAR3

Start from the official POSTAR portal at <https://postar.ncrnalab.org/> and
obtain the complete human binding table used for POSTAR3. Extract it as
`data/human.txt` (a `.gz` path is also accepted when selected explicitly).
The adapter expects at least ten tab-separated columns in this order:
chromosome, start, end, site identifier, strand, RBP, assay/method, cell type,
accession(s), and score. Keep the original compressed download, download date,
release label if supplied, byte size, and checksum outside version control.
The current 4.47-GB project table is a local frozen snapshot; do not substitute
a differently structured export without passing the doctor and a representative
overlap run.

Database content can change. The run manifest records the exact local input
checksum, while the separate resource record should retain the release or
download date. Cross-database agreement is not assumed to be independent
replication because source experiments can overlap.

## Other analysis resources

- KnockRBP: place each DEG JSON beside `dataset_metadata.tsv`; dataset IDs must
  be unique and metadata must identify RBP, cell line, and p-value availability.
- Expression: supply the relevant count/TPM table explicitly. Detection rules
  are workflow parameters and must be recorded with the output.
- STRING: network retrieval requires internet access unless a populated cache
  is supplied. Preserve the API response cache and manifest for reproducibility.
- RM2Target and cancer-gene lists are interpretation resources, not site-level
  binding evidence; keep their release/source provenance.

## Verification

With the complete project present:

```bash
cd /path/to/rna_mod_workspace/rna_mod_pipeline
rna-mod-audit --project-root ..
```

For a clean laptop or arbitrary dataset without historical results:

```bash
cd /path/to/rna_mod_workspace/rna_mod_pipeline
rna-mod-doctor --project-root ..
```

The audit validates expected files, catalogue membership, modification tables,
KnockRBP dataset mapping, and established baselines. A passing audit verifies
local integrity checks; it does not independently validate database biology.
