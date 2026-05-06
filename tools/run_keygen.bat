@echo off
cd /d "%~dp0"

REM ── Find Python ──────────────────────────────────────────────────────────────
set PYTHON=
for %%P in (python3 python py) do (
    if not defined PYTHON (
        where %%P >nul 2>&1 && set PYTHON=%%P
    )
)

if not defined PYTHON (
    echo.
    echo  ERROR: Python not found in PATH.
    echo  Please install Python 3.8+ from https://python.org
    echo  and make sure to tick "Add Python to PATH" during install.
    echo.
    pause
    exit /b 1
)

REM ── Run the GUI ───────────────────────────────────────────────────────────────
echo Starting NextGen-Amplicon Key Generator...
%PYTHON% keygen.py

if errorlevel 1 (
    echo.
    echo  === keygen.py exited with an error (see above) ===
    pause
)
