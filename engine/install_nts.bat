@echo off
setlocal
cd /d "%~dp0"
echo === HRMS NTS toolkit installer ===
set "BASEPY="
where py >nul 2>nul && set "BASEPY=py -3"
if not defined BASEPY where python >nul 2>nul && set "BASEPY=python"
if not defined BASEPY (
  echo Python 3.10+ not found. Install from https://www.python.org/downloads/ and re-run.
  pause & exit /b 1
)
echo Creating virtual environment in .venv ...
%BASEPY% -m venv .venv || (echo venv creation failed & pause & exit /b 1)
call ".venv\Scripts\activate.bat"
python -m pip install --upgrade pip
python -m pip install -e .[dev] || (echo pip install failed & pause & exit /b 1)
echo.
python -m hrms_nts doctor
echo.
echo Installation finished.  Next:  copy config.local.example.yaml to config.local.yaml, edit the data paths, then run_nts.bat doctor
pause
endlocal
