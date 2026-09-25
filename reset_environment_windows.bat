@echo off
setlocal
cd /d "%~dp0" || (
    echo ERROR: Could not open the source folder.
    echo Copy the source folder to a local drive, then run this helper again.
    pause
    exit /b 1
)

set "VENV_DIR=%CD%\.venv"

if not exist "%VENV_DIR%" (
    echo.
    echo No local .venv environment exists; nothing was removed.
    echo.
    pause
    exit /b 0
)

if not exist "%VENV_DIR%\pyvenv.cfg" (
    echo.
    echo ERROR: Refusing to remove .venv because pyvenv.cfg is absent.
    echo The folder was left unchanged for manual inspection.
    echo.
    pause
    exit /b 1
)

rmdir /s /q "%VENV_DIR%"
if exist "%VENV_DIR%" (
    echo.
    echo ERROR: Windows could not fully remove the environment.
    echo Close programs using it, then run this reset helper again.
    echo.
    pause
    exit /b 1
)

echo.
echo The regenerable local .venv environment was removed.
echo Scientific inputs and outputs were not touched.
echo.
pause
exit /b 0
