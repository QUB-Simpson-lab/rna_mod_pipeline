#!/bin/bash

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="$SCRIPT_DIR/.venv"

pause_window() {
    if [ -t 0 ]; then
        printf "\nPress Return to close this window..."
        read -r _ || true
    fi
}

if [ ! -e "$VENV_DIR" ]; then
    printf "\nNo local .venv environment exists; nothing was removed.\n"
    pause_window
    exit 0
fi

if [ ! -f "$VENV_DIR/pyvenv.cfg" ]; then
    printf "\nERROR: Refusing to remove .venv because pyvenv.cfg is absent.\n"
    printf "The folder was left unchanged for manual inspection.\n"
    pause_window
    exit 1
fi

if ! rm -rf -- "$VENV_DIR" || [ -e "$VENV_DIR" ]; then
    printf "\nERROR: macOS could not fully remove the environment.\n"
    printf "Close programs using it, then run this reset helper again.\n"
    pause_window
    exit 1
fi
printf "\nThe regenerable local .venv environment was removed.\n"
printf "Scientific inputs and outputs were not touched.\n"
pause_window
exit 0
