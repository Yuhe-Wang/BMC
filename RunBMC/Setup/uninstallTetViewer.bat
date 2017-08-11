@echo off 
Title uninstalling
REM test if we have administrator privilege
bcdedit>nul
if errorlevel 1 (echo Please right-click to run with administrator privilege
pause
goto end)

echo cleaning the registered association. . . 
reg delete "HKCR\MDOSE" /f
reg delete "HKCR\.mdose" /f
cd /d "%~dp0"
refresh_icon.exe
pause Finished!
:end
