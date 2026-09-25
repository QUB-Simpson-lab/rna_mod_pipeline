@echo off
setlocal
cd /d "%~dp0" || (
    echo ERROR: Could not open the source folder.
    echo Copy the source folder to a local drive, then run this helper again.
    pause
    exit /b 1
)

echo.
echo RNA-modification pipeline installer
echo ===================================
echo.

set "VENV_PYTHON=%CD%\.venv\Scripts\python.exe"
set "PYTHONDONTWRITEBYTECODE=1"
set "PIP_DISABLE_PIP_VERSION_CHECK=1"

if exist "%VENV_PYTHON%" (
    echo Reusing the existing .venv environment.
    "%VENV_PYTHON%" -c "import sys; raise SystemExit(0 if sys.version_info[:2] in ((3,10),(3,11),(3,12)) else 1)"
    if errorlevel 1 goto bad_venv
    goto install
)

set "PYTHON_LAUNCHER="

where py >nul 2>nul
if not errorlevel 1 (
    py -3.11 -c "import sys" >nul 2>nul
    if not errorlevel 1 set "PYTHON_LAUNCHER=py -3.11"
    if not defined PYTHON_LAUNCHER (
        py -3.12 -c "import sys" >nul 2>nul
        if not errorlevel 1 set "PYTHON_LAUNCHER=py -3.12"
    )
    if not defined PYTHON_LAUNCHER (
        py -3.10 -c "import sys" >nul 2>nul
        if not errorlevel 1 set "PYTHON_LAUNCHER=py -3.10"
    )
)

if not defined PYTHON_LAUNCHER (
    where python >nul 2>nul
    if not errorlevel 1 (
        python -c "import sys; raise SystemExit(0 if sys.version_info[:2] in ((3,10),(3,11),(3,12)) else 1)" >nul 2>nul
        if not errorlevel 1 set "PYTHON_LAUNCHER=python"
    )
)

if not defined PYTHON_LAUNCHER (
    where python3 >nul 2>nul
    if not errorlevel 1 (
        python3 -c "import sys; raise SystemExit(0 if sys.version_info[:2] in ((3,10),(3,11),(3,12)) else 1)" >nul 2>nul
        if not errorlevel 1 set "PYTHON_LAUNCHER=python3"
    )
)

if not defined PYTHON_LAUNCHER goto no_python

echo Creating .venv with %PYTHON_LAUNCHER%...
%PYTHON_LAUNCHER% -m venv ".venv"
if errorlevel 1 goto failed

:install
echo Installing the tested build tools...
"%VENV_PYTHON%" -m pip install --no-cache-dir -r "requirements-build-lock.txt"
if errorlevel 1 goto failed

echo Installing the constrained pipeline and graphical interface...
"%VENV_PYTHON%" -m pip install --no-cache-dir --no-build-isolation -c "requirements-lock.txt" -c "requirements-gui-lock.txt" -e ".[gui]"
if errorlevel 1 goto failed

echo Checking dependency consistency...
"%VENV_PYTHON%" -m pip check
if errorlevel 1 goto failed

echo Checking all Python source files...
set "PYTHONPYCACHEPREFIX=%CD%\.venv\pycache"
"%VENV_PYTHON%" -m compileall -q "src" "scripts"
if errorlevel 1 goto failed

echo Checking the graphical interface imports...
"%VENV_PYTHON%" -c "import rna_mod_pipeline, rna_mod_pipeline.gui, PySide6; print('Python', __import__('sys').version.split()[0]); print('PySide6', PySide6.__version__)"
if errorlevel 1 goto failed
"%CD%\.venv\Scripts\rna-mod-gui.exe" --help >nul
if errorlevel 1 goto failed
"%CD%\.venv\Scripts\rna-mod-doctor.exe" --help >nul
if errorlevel 1 goto failed

echo.
echo Installation completed successfully.
echo Double-click run_windows.bat to open the graphical interface.
echo.
pause
exit /b 0

:no_python
echo.
echo ERROR: Supported Python 3.10, 3.11, or 3.12 was not found.
echo Install Python from https://www.python.org/downloads/windows/
echo Select "Add python.exe to PATH" during installation, then run this file again.
echo.
pause
exit /b 1

:bad_venv
echo.
echo ERROR: The existing .venv is not using supported Python 3.10-3.12.
echo Run reset_environment_windows.bat, then install again.
echo.
pause
exit /b 1

:failed
echo.
echo ERROR: Installation failed. Review the message above.
echo Check the internet connection and confirm that Python 3.10-3.12 is installed.
echo.
pause
exit /b 1
