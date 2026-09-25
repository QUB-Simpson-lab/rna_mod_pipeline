# Dataset-centred graphical launcher

## Implementation status

The source-and-virtual-environment pilot is implemented as an optional PySide6
layer. It accepts arbitrary named m6A, m5C, or pseudouridine bedMethyl datasets
and does not require the four historical project callsets. The command line
remains the scientific source of truth.

Start it with:

```text
python scripts/pipeline_gui.py --project-root ..
```

See `GUI_GUIDE.md` for installation and second-laptop testing.

## Architecture

The graphical interface is a separate, optional layer over the frozen command-
line pipeline. The command-line scripts remain the scientific source of truth.
The GUI constructs and runs those commands; it does not reimplement filtering,
overlap, statistics, or plotting. Its preflight checks provide early feedback,
while each command still performs its authoritative input and output
validation.

This keeps graphical and command-line runs scientifically identical and allows
the pipeline to remain usable on servers where a desktop interface is
unavailable.

## Implemented first version

The first version exposes:

- Phase 1 filtering, DRACH annotation, and metagene annotation;
- loose overlap;
- transcript-region overlap;
- cross-database validation;
- expression integration;
- KnockRBP validation;
- STRING analysis;
- cross-modification comparison;
- same-modification dataset robustness comparison;
- plot regeneration and historical-project validation.

For each workflow, the interface:

1. exposes the commonly used options while leaving final argument validation
   to the real command parser;
2. derives inputs and isolated outputs from the selected dataset profile;
3. displays the exact command before execution;
4. checks missing paths and incompatible options before starting;
5. runs one job at a time with live standard-output and error logs;
6. reports the current and final run state;
7. allows cancellation without deleting an earlier successful output;
8. opens and previews manifest-indexed results after validated completion.

Each bedMethyl is registered as a separate named dataset. A dataset name
identifies the sample or replicate; `m6a`, `m5c`, or `pseu` identifies the
chemistry. Outputs are isolated below the selected dataset workspace. Two
compatible same-modification workspaces can be compared without pooling their
raw calls.

Overwrite remains explicit. The interface does not delete files, silently chain
analyses, or infer missing scientific inputs.

## Implementation

PySide6 is used because it supports a responsive native interface and reliable
process control. It remains an optional
dependency, for example:

```text
pip install -e ".[gui]"
```

Tkinter is not recommended as the distribution default because its availability and
linked Tk version vary between the Python installations currently present on
this workstation.

A maintainable layout is used:

```text
scripts/pipeline_gui.py
src/rna_mod_pipeline/gui/app.py
src/rna_mod_pipeline/gui/workflow_registry.py
src/rna_mod_pipeline/gui/command_builder.py
src/rna_mod_pipeline/gui/runner.py
src/rna_mod_pipeline/gui/widgets.py
```

The GUI is split into focused state, form, resource, command, execution, and
results modules; scientific calculations remain outside the interface.
Commands are launched as an argument list with the active Python executable,
never through a shell string.

## Progress reporting

The current commands mainly report completed stages, so the GUI shows live logs
and an indeterminate activity indicator. It does not invent a percentage. If
detailed progress is later required, stable stage events should be added to the
command-line workflows first so both terminal and GUI consume the same events.

## Acceptance status

Completed locally:

- generated arguments for every registered workflow parse with the real
  command parser;
- tests cover command construction, paths with spaces, cancellation,
  manifest-gated completion, and window startup;
- an arbitrary registered bedMethyl completed a real Phase 1 run into an
  isolated output workspace;
- the GUI invokes the same scripts, so scientific output manifests retain the
  command-line schema;
- staged command publication protects an earlier completed output;
- the full command-line pipeline remains usable without PySide6;
- the main and workflow windows rendered successfully in a headless Qt
  environment and were visually inspected.

Still required before the lab GitHub release:

- install from a clean copy on both the maintainer's macOS and Windows laptops;
- complete the acceptance sequence in `GUI_GUIDE.md`;
- run at least one representative real-data downstream workflow before
  production use.

## Current distribution decision

Use the source checkout and a fresh virtual environment for the first
cross-laptop pilot. PySide6 is optional, and the complete command-line pipeline
continues to run without it.

Do not create a standalone `.app` or `.exe` until the pilot passes on the target
laptop. Frozen applications must be built separately for each operating system
and architecture, while large reference and binding resources should remain
external.
