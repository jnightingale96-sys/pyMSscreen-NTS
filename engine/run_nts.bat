@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"
call "%~dp0nts_python.bat" || (pause & exit /b 1)

rem ---- direct mode: run_nts.bat <any hrms_nts arguments> ---------------------------------------------------
if not "%~1"=="" goto direct
goto menu

:direct
%PY% -m hrms_nts %*
exit /b %errorlevel%

:menu
echo.
echo ============ HRMS non-target screening (Level 5/4/3/2b/2a/1 + semi-quant) ============
echo   1  Check installation / paths (doctor)
echo   2  Build MS/MS spectral library (MassBank/GNPS/in-house MSP -^> SQLite^)   [one-off]
echo   3  RUN full screening on mzML folder/files
echo   4  Targeted prescreen from compounds.csv (pyMSscreen-style QA)
echo   5  Predict product ions / MRM transitions for a SMILES
echo   6  Measure calibrant areas from a standard mzML (semi-quant)
echo   7  Convert Thermo .raw folder to centroided mzML
echo   8  Run the automated tests (pytest)
echo   9  Validate a results folder against expected.csv
echo   0  Exit
set /p CH=Choose:
if "%CH%"=="1" %PY% -m hrms_nts doctor & goto menu
if "%CH%"=="2" %PY% -m hrms_nts build-library & goto menu
if "%CH%"=="3" goto run
if "%CH%"=="4" goto pre
if "%CH%"=="5" goto pred
if "%CH%"=="6" goto cal
if "%CH%"=="7" goto conv
if "%CH%"=="8" %PY% -m pytest tests & goto menu
if "%CH%"=="9" goto val
if "%CH%"=="0" exit /b 0
goto menu

:run
set /p IN=mzML folder or file (wildcards ok):
set /p MX=Matrix [water/grain/soil/extract] (default water):
if "%MX%"=="" set MX=water
set /p SQ=Semi-quant mode [single/logie/curve/none] (default single):
if "%SQ%"=="" set SQ=single
set /p IS=In-silico engine [builtin/metfrag/both/none] (default builtin):
if "%IS%"=="" set IS=builtin
%PY% -m hrms_nts run "%IN%" --matrix %MX% --semiquant %SQ% --insilico %IS%
pause & goto menu

:pre
set /p CP=compounds.csv:
set /p IN=mzML folder or file:
%PY% -m hrms_nts prescreen "%CP%" "%IN%"
pause & goto menu

:pred
set /p SM=SMILES:
set /p AD=Adduct (default [M+H]+):
if "%AD%"=="" set AD=[M+H]+
%PY% -m hrms_nts predict-ions --smiles "%SM%" --adduct "%AD%"
pause & goto menu

:cal
set /p MZ=Standard mzML:
%PY% -m hrms_nts calibrate "%MZ%" --out calibrant_areas.csv
pause & goto menu

:conv
set /p RD=Folder with .raw files:
set /p OD=Output folder for mzML:
%PY% -m hrms_nts convert "%RD%" "%OD%"
pause & goto menu

:val
set /p RS=Results folder: 
set /p EX=expected.csv: 
set /p AB=absent.csv (optional, Enter to skip): 
if "%AB%"=="" (%PY% -m hrms_nts validate "%RS%" "%EX%") else (%PY% -m hrms_nts validate "%RS%" "%EX%" --absent "%AB%")
pause & goto menu
