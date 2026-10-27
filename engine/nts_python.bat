@echo off
rem Helper: sets PY to a Python that can run the toolkit (the .venv made by install_nts.bat first, then any Python with rdkit + pyteomics).
set "PY="
if exist "%~dp0.venv\Scripts\python.exe" set "PY=%~dp0.venv\Scripts\python.exe"
rem inside pyMSscreen-NTS the environment made by START_pyMSscreen-NTS.bat sits one folder up
if not defined PY if exist "%~dp0..\.venv\Scripts\python.exe" set "PY=%~dp0..\.venv\Scripts\python.exe"
if not defined PY (
  where py >nul 2>nul && (py -3 -c "import rdkit,pyteomics" >nul 2>nul && set "PY=py -3")
)
if not defined PY (
  where python >nul 2>nul && (python -c "import rdkit,pyteomics" >nul 2>nul && set "PY=python")
)
if not defined PY (
  echo No suitable Python found. Run install_nts.bat first.
  exit /b 1
)
exit /b 0
