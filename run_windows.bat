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

if not exist "%VENV_PYTHON%" (
    echo.
    echo The pipeline environment has not been installed.
    echo Run install_windows.bat first.
    echo.
    pause
    exit /b 1
)

"%VENV_PYTHON%" -c "import rna_mod_pipeline, PySide6" >nul 2>nul
if errorlevel 1 (
    echo.
    echo The graphical interface is not fully installed.
    echo Run install_windows.bat again.
    echo.
    pause
    exit /b 1
)

set "PYTHONUTF8=1"
"%CD%\.venv\Scripts\rna-mod-gui.exe" --project-root "%CD%\.."
set "GUI_EXIT=%ERRORLEVEL%"

if not "%GUI_EXIT%"=="0" (
    echo.
    echo The graphical interface exited with an error.
    echo Review the message above or rerun install_windows.bat.
    echo.
    pause
)

exit /b %GUI_EXIT%
