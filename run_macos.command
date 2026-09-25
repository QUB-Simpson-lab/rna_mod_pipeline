#!/bin/bash

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

pause_window() {
    if [ -t 0 ]; then
        printf "\nPress Return to close this window..."
        read -r _ || true
    fi
}

VENV_PYTHON="$SCRIPT_DIR/.venv/bin/python"
export PYTHONDONTWRITEBYTECODE=1

if [ ! -x "$VENV_PYTHON" ]; then
    printf "\nThe pipeline environment has not been installed.\n"
    printf "Run install_macos.command first.\n"
    pause_window
    exit 1
fi

if ! "$VENV_PYTHON" -c 'import rna_mod_pipeline, PySide6' >/dev/null 2>&1; then
    printf "\nThe graphical interface is not fully installed.\n"
    printf "Run install_macos.command again.\n"
    pause_window
    exit 1
fi

export PYTHONUTF8=1
"$SCRIPT_DIR/.venv/bin/rna-mod-gui" \
    --project-root "$SCRIPT_DIR/.."
status="$?"

if [ "$status" -ne 0 ]; then
    printf "\nThe graphical interface exited with an error.\n"
    printf "Review the message above or rerun install_macos.command.\n"
    pause_window
fi

exit "$status"
