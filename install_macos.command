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

fail() {
    status="$1"
    printf "\nERROR: Installation failed. Review the message above.\n"
    printf "Check the internet connection and confirm that Python 3.10-3.12 is installed.\n"
    pause_window
    exit "$status"
}

printf "\nRNA-modification pipeline installer\n"
printf "===================================\n\n"

VENV_PYTHON="$SCRIPT_DIR/.venv/bin/python"
AVAILABLE_KB="$(df -Pk "$SCRIPT_DIR" 2>/dev/null | awk 'NR == 2 {print $4}')"
export PYTHONDONTWRITEBYTECODE=1
export PIP_DISABLE_PIP_VERSION_CHECK=1

if [ -n "$AVAILABLE_KB" ] && [ "$AVAILABLE_KB" -lt 2097152 ]; then
    printf "ERROR: Less than 2 GB of free disk space is available.\n"
    printf "Free additional space before installing the GUI environment.\n"
    pause_window
    exit 1
elif [ -n "$AVAILABLE_KB" ] && [ "$AVAILABLE_KB" -lt 4194304 ]; then
    printf "WARNING: Less than 4 GB of free disk space is available.\n"
    printf "Installation can continue, but freeing more space is recommended.\n\n"
fi

if [ -x "$VENV_PYTHON" ]; then
    printf "Reusing the existing .venv environment.\n"
    "$VENV_PYTHON" -c \
        'import sys; raise SystemExit(0 if (3, 10) <= sys.version_info < (3, 13) else 1)' \
        || {
            printf "\nERROR: The existing .venv is not using supported Python 3.10-3.12.\n"
            printf "Run reset_environment_macos.command, then install again.\n"
            pause_window
            exit 1
        }
else
    PYTHON_BIN=""
    for candidate in python3.11 python3.12 python3.10 python3 \
        /opt/homebrew/bin/python3.11 /usr/local/bin/python3.11 \
        /Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11 \
        /opt/homebrew/bin/python3.12 /usr/local/bin/python3.12 \
        /Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12 \
        /opt/homebrew/bin/python3.10 /usr/local/bin/python3.10 \
        /Library/Frameworks/Python.framework/Versions/3.10/bin/python3.10; do
        if command -v "$candidate" >/dev/null 2>&1 &&
            "$candidate" -c \
                'import sys; raise SystemExit(0 if (3, 10) <= sys.version_info < (3, 13) else 1)' \
                >/dev/null 2>&1; then
            PYTHON_BIN="$candidate"
            break
        fi
    done

    if [ -z "$PYTHON_BIN" ]; then
        printf "ERROR: Supported Python 3.10, 3.11, or 3.12 was not found.\n"
        printf "Install Python from https://www.python.org/downloads/macos/\n"
        printf "Then run this installer again.\n"
        pause_window
        exit 1
    fi

    printf "Creating .venv with %s...\n" "$PYTHON_BIN"
    "$PYTHON_BIN" -m venv "$SCRIPT_DIR/.venv" || fail "$?"
fi

printf "Installing the tested build tools...\n"
"$VENV_PYTHON" -m pip install --no-cache-dir \
    -r "$SCRIPT_DIR/requirements-build-lock.txt" || fail "$?"

printf "Installing the constrained pipeline and graphical interface...\n"
"$VENV_PYTHON" -m pip install --no-cache-dir --no-build-isolation \
    -c "$SCRIPT_DIR/requirements-lock.txt" \
    -c "$SCRIPT_DIR/requirements-gui-lock.txt" \
    -e ".[gui]" || fail "$?"

printf "Checking dependency consistency...\n"
"$VENV_PYTHON" -m pip check || fail "$?"

printf "Checking all Python source files...\n"
PYTHONPYCACHEPREFIX="$SCRIPT_DIR/.venv/pycache" \
    "$VENV_PYTHON" -m compileall -q "$SCRIPT_DIR/src" "$SCRIPT_DIR/scripts" \
    || fail "$?"

printf "Checking the graphical interface imports...\n"
"$VENV_PYTHON" -c \
    'import sys, PySide6, rna_mod_pipeline, rna_mod_pipeline.gui; print("Python", sys.version.split()[0]); print("PySide6", PySide6.__version__)' \
    || fail "$?"
"$SCRIPT_DIR/.venv/bin/rna-mod-gui" --help >/dev/null \
    || fail "$?"
"$SCRIPT_DIR/.venv/bin/rna-mod-doctor" --help >/dev/null \
    || fail "$?"

printf "\nInstallation completed successfully.\n"
printf "Double-click run_macos.command to open the graphical interface.\n"
pause_window
exit 0
