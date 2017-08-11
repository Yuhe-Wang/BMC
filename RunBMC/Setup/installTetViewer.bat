@echo off 
Title Associating TetViewer.exe¡­ 
REM test if we have administrator privilege
bcdedit>nul
if errorlevel 1 (echo Please right-click to run with administrator privilege
pause
goto end)

cd /d "%~dp0"
pushd ..

echo associating files. . . 
echo Some anti-virus software may prevent this associating process. Please allow it manually.
reg add "HKCR\.mdose" /f /ve /t REG_SZ /d "MDOSE">nul 2>nul 
reg add "HKCR\MDOSE\DefaultIcon" /f /ve /t REG_SZ /d "%CD%\TetViewer.exe,4">nul 2>nul 
reg add "HKCR\MDOSE\shell\open(&O)\command" /f /ve /t REG_SZ /d "%CD%\TetViewer.exe \"%%1\"">nul 2>nul 
popd
refresh_icon.exe
echo Association finished!
pause 

:end
