# RBP catalogue

`rbp_catalog.tsv` is the versioned source of truth for the RBP panels tested by
the refactored workflows.

| Column | Meaning |
|---|---|
| `catalog_version` | Version date for the frozen panel |
| `database` | `ornament`, `encori`, or `postar3` |
| `rbp_name` | Name retained for loose-output compatibility |
| `canonical_rbp` | Canonical identity used for cross-database comparisons |
| `source_names` | Semicolon-separated resource names whose intervals are unioned |

The frozen panel contains:

- 133 oRNAment RBPs;
- 281 ENCORI RBPs;
- 216 POSTAR3 RBPs.

For example, the oRNAment SF2 and SRSF1 resources are unioned and tested once as
SRSF1. The loose oRNAment output retains other historical names where necessary,
while cross-database and transcript-region workflows use `canonical_rbp`.

Do not add a new alias by editing an output table. Update this catalogue, bump
`catalog_version`, and rerun the catalogue tests.

The catalogue defines which RBPs are tested, not which are biologically related
to a modification. The single shared modification-related display registry is
`src/rna_mod_pipeline/config.py::KNOWN_RELATED`; loose and transcript plots use
that same registry. These flags are visual annotations, not additional
statistical evidence.
