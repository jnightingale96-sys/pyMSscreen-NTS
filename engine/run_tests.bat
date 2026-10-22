@echo off
cd /d "%~dp0"
call "%~dp0nts_python.bat" || (pause & exit /b 1)
%PY% -m pytest tests
pause
