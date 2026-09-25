#!/bin/bash
set -u
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUTF8=1
status=1
if [ ! -x "$SCRIPT_DIR/.venv/bin/rna-mod-gui" ]; then
    printf 'Run install_macos.command first.\n'
elif [ ! -d "$SCRIPT_DIR/../rna_mod_demo/refactored_outputs" ]; then
    printf 'Run test_demo_macos.command first to create the synthetic example.\n'
else
    "$SCRIPT_DIR/.venv/bin/rna-mod-gui" --project-root "$SCRIPT_DIR/../rna_mod_demo"
    status="$?"
fi
if [ "$status" -ne 0 ] && [ -t 0 ]; then
    printf '\nPress Return to close this window...'
    read -r _ || true
fi
exit "$status"
