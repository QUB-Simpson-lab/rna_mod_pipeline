@echo off
setlocal
cd /d "%~dp0" || (
    echo ERROR: Could not open the source folder.
    echo Copy the source folder to a local drive, then run this helper again.
    pause
    exit /b 1
)
set "PYTHONDONTWRITEBYTECODE=1"
set "PYTHONUTF8=1"
if not exist "%CD%\.venv\Scripts\rna-mod-gui.exe" goto not_installed
if not exist "%CD%\..\rna_mod_demo\refactored_outputs" goto no_demo
"%CD%\.venv\Scripts\rna-mod-gui.exe" --project-root "%CD%\..\rna_mod_demo"
if errorlevel 1 goto failed
exit /b 0
:not_installed
echo Run install_windows.bat first.
pause
exit /b 1
:no_demo
echo Run test_demo_windows.bat first to create the synthetic example.
pause
exit /b 1
:failed
echo The demonstration GUI exited with an error. Keep the output for review.
pause
exit /b 1
