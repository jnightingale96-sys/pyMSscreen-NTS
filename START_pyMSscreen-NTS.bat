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
rem .venv\installed.ok is written only after every package installed, so an interrupted first start is resumed next time
if not exist ".venv\installed.ok" (
  echo === First start: setting up the Python environment in .venv ===
  set "HERE=%~dp0"
  call set "TAIL=%%HERE:~60,1%%"
  if defined TAIL (
    echo This folder's path is too long for Windows ^(some packages exceed the 260-character limit^):
    echo   %~dp0
    echo Move the folder to a short path such as C:\NTS\ and start again.
    pause & exit /b 1
  )
  if not exist ".venv\Scripts\python.exe" (
    set "BASEPY="
    where py >nul 2>nul && set "BASEPY=py -3"
    if not defined BASEPY where python >nul 2>nul && set "BASEPY=python"
    if not defined BASEPY (
      echo Python 3.10 or newer was not found. Install it from https://www.python.org/downloads/ ^(tick "Add to PATH"^) and start again.
      pause & exit /b 1
    )
    call %%BASEPY%% -m venv --without-pip .venv || (echo Could not create .venv & pause & exit /b 1)
  )
  ".venv\Scripts\python.exe" -m pip --version >nul 2>nul || ".venv\Scripts\python.exe" -m ensurepip --upgrade --default-pip || (echo Could not install pip into .venv - start again, or delete the .venv folder and start again. & pause & exit /b 1)
  ".venv\Scripts\python.exe" -m pip install --upgrade pip
  ".venv\Scripts\python.exe" -m pip install -r requirements-windows.txt || (echo Package installation failed - see messages above. Start again to retry. & pause & exit /b 1)
  echo ok> ".venv\installed.ok"
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

rem "START_pyMSscreen-NTS.bat install" = set everything up without starting the app (used by Claude Code / scripts)
if /i "%~1"=="install" (
  echo Installed. Command line: .venv\Scripts\python -m hrms_nts doctor   ^(run from the engine folder^)
  exit /b 0
)

rem ---- 3. Run --------------------------------------------------------------------------------
if not exist "%WORKING_DIR%" mkdir "%WORKING_DIR%"
set "PYTHONIOENCODING=utf-8"
".venv\Scripts\python.exe" backend\run_windows.py
if errorlevel 1 pause
endlocal
