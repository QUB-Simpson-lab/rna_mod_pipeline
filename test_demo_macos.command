#!/bin/bash
set -u
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1
VENV_PYTHON="$SCRIPT_DIR/.venv/bin/python"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUTF8=1
export QT_QPA_PLATFORM=offscreen

finish() {
    if [ -t 0 ]; then
        printf '\nPress Return to close this window...'
        read -r _ || true
    fi
    exit "$1"
}

if [ ! -x "$VENV_PYTHON" ]; then
    printf 'Run install_macos.command first.\n'
    finish 1
fi

printf 'Running the regression suite...\n'
"$VENV_PYTHON" -m unittest discover -s tests -v || finish "$?"
printf '\nPreparing and running synthetic data only...\n'
"$VENV_PYTHON" scripts/prepare_acceptance_demo.py \
    --workspace "$SCRIPT_DIR/../rna_mod_demo" --run || finish "$?"
"$VENV_PYTHON" scripts/check_acceptance_demo.py \
    --workspace "$SCRIPT_DIR/../rna_mod_demo" || finish "$?"
printf '\nAutomated checks passed. Open run_demo_macos.command to inspect the example.\n'
finish 0
