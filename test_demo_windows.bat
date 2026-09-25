@echo off
setlocal
cd /d "%~dp0" || (
    echo ERROR: Could not open the source folder.
    echo Copy the source folder to a local drive, then run this helper again.
    pause
    exit /b 1
)
set "VENV_PYTHON=%CD%\.venv\Scripts\python.exe"
set "PYTHONDONTWRITEBYTECODE=1"
set "PYTHONUTF8=1"
set "QT_QPA_PLATFORM=offscreen"
if not exist "%VENV_PYTHON%" goto not_installed
echo Running the regression suite...
"%VENV_PYTHON%" -m unittest discover -s tests -v
if errorlevel 1 goto failed
echo.
echo Preparing and running synthetic data only...
"%VENV_PYTHON%" scripts\prepare_acceptance_demo.py --workspace "%CD%\..\rna_mod_demo" --run
if errorlevel 1 goto failed
"%VENV_PYTHON%" scripts\check_acceptance_demo.py --workspace "%CD%\..\rna_mod_demo"
if errorlevel 1 goto failed
echo.
echo Automated checks passed. Open run_demo_windows.bat to inspect the example.
pause
exit /b 0
:not_installed
echo Run install_windows.bat first.
pause
exit /b 1
:failed
echo.
echo A test did not complete successfully. Keep the output for review.
pause
exit /b 1
