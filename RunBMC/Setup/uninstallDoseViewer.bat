@echo off 
Title uninstalling
REM test if we have administrator privilege
bcdedit>nul
if errorlevel 1 (echo Please right-click to run with administrator privilege
pause
goto end)

echo cleaning the registered association. . . 
reg delete "HKCR\DOSE" /f
reg delete "HKCR\.dose" /f
reg delete "HKCR\PHTM" /f
reg delete "HKCR\.phtm" /f
reg delete "HKCR\VIEWRAYRED" /f
reg delete "HKCR\.red" /f
cd /d "%~dp0"
refresh_icon.exe
pause Finished!
:end
