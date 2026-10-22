@echo off
rem ==========================================================================================
rem  pyMSscreen-NTS - double-click to start.  First start installs everything (5-10 min);
rem  later starts take a few seconds.  Opens http://127.0.0.1:8765 in your browser.
rem  Close the black window to stop the application.
rem ==========================================================================================
setlocal
cd /d "%~dp0"
title pyMSscreen-NTS
if not defined PORT set "PORT=8765"
if not defined WORKING_DIR set "WORKING_DIR=%USERPROFILE%\Documents\pyMSscreen-NTS_data"

rem ---- already running? just open the browser --------------------------------------------
netstat -ano | findstr /r /c:"127.0.0.1:%PORT% .*LISTENING" >nul 2>nul && (
  echo pyMSscreen-NTS is already running - opening the browser.
  start "" "http://127.0.0.1:%PORT%/"
  exit /b 0
)

rem ---- 1. Python environment ---------------------------------------------------------------
if not exist ".venv\Scripts\python.exe" (
  echo === First start: creating the Python environment in .venv ===
  set "BASEPY="
  where py >nul 2>nul && set "BASEPY=py -3"
  if not defined BASEPY where python >nul 2>nul && set "BASEPY=python"
  if not defined BASEPY (
    echo Python 3.10 or newer was not found. Install it from https://www.python.org/downloads/ ^(tick "Add to PATH"^) and start again.
    pause & exit /b 1
  )
  call %%BASEPY%% -m venv .venv || (echo Could not create .venv & pause & exit /b 1)
  ".venv\Scripts\python.exe" -m pip install --upgrade pip
  ".venv\Scripts\python.exe" -m pip install -r requirements-windows.txt || (echo Package installation failed - see messages above. & pause & exit /b 1)
)

rem ---- 2. Web interface (built once with Node.js) -----------------------------------------
if not exist "frontend\dist\index.html" (
  echo === First start: building the web interface ===
  where npm >nul 2>nul || (
    echo Node.js is needed once to build the interface. Install the LTS version from https://nodejs.org and start again.
    pause & exit /b 1
  )
  pushd frontend
  call npm install --no-audit --no-fund || (popd & echo npm install failed & pause & exit /b 1)
  call npm run build:windows || (popd & echo Interface build failed & pause & exit /b 1)
  popd
)

rem ---- 3. Run --------------------------------------------------------------------------------
if not exist "%WORKING_DIR%" mkdir "%WORKING_DIR%"
set "PYTHONIOENCODING=utf-8"
".venv\Scripts\python.exe" backend\run_windows.py
if errorlevel 1 pause
endlocal
