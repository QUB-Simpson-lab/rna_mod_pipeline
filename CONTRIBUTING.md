# Contributing

This is research software. A code change and a scientific-method change are not
the same thing, and both must be reviewable.

## Before changing code

1. Create a short-lived branch from the reviewed release version.
2. State whether the change affects interfaces only, implementation only, or a
   scientific definition/statistic.
3. Do not edit raw inputs, legacy results, or existing run manifests in place.
4. Use a new output directory while validating a changed workflow.

## Required checks

From the repository root in a supported environment, run:

```bash
python -m unittest discover -s tests -v
python -m pip check
rna-mod-gui --help
rna-mod-audit --help
rna-mod-doctor --help
```

Run `rna-mod-audit --project-root ..` when the complete project data and legacy
baselines are available. Run `rna-mod-doctor --project-root ..` for a clean,
resource-profile-based installation. For a scientific change, rerun at least one compact
fixture plus a representative real dataset, compare the output tables, and
explain every intended numerical difference.

## Pull-request record

Include:

- the biological or software reason for the change;
- exact commands and environment versions;
- tests and real-data comparisons performed;
- workflows and columns that can change;
- backwards-compatibility and rerun requirements;
- new or updated resource provenance.

Do not commit `.venv`, caches, raw databases, experimental bedMethyl inputs,
generated scientific outputs, access tokens, or patient-identifiable data.
The acceptance example generates explicitly synthetic fixtures outside the
checkout. Update `CHANGELOG.md` for a release-facing change. The source uses
the MIT licence in `LICENSE`; external databases retain their own terms.
