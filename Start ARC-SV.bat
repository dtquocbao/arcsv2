@echo off
cd /d "%~dp0"
REM Developer build from source: Start ARC-SV.bat --build
if /I "%~1"=="--build" (
  powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\start-arcsv.ps1" -Build
) else (
  powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\start-arcsv.ps1"
)
if errorlevel 1 pause
