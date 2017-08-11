@echo off
cd /d "%~dp0"
copy /y curl.exe ..\
copy /y 7z.exe ..\

Title Associating DoseViewer.exe¡­ 
REM test if we have administrator privilege
bcdedit>nul
if errorlevel 1 (echo Please right-click to run with administrator privilege
    pause
    goto end
)

reg import "TDR_20Seconds_import&reboot.reg"
pushd ..

echo associating files...
echo Some anti-virus software may prevent this associating process. Please allow it manually.
reg add "HKCR\.phtm" /f /ve /t REG_SZ /d "PHTM">nul 2>nul 
reg add "HKCR\.dose" /f /ve /t REG_SZ /d "DOSE">nul 2>nul 
reg add "HKCR\.red" /f /ve /t REG_SZ /d "VIEWRAYRED">nul 2>nul
reg add "HKCR\PHTM\DefaultIcon" /f /ve /t REG_SZ /d "%CD%\DoseViewer.exe,2">nul 2>nul 
reg add "HKCR\PHTM\shell\open(&O)\command" /f /ve /t REG_SZ /d "%CD%\DoseViewer.exe \"%%1\"">nul 2>nul 
reg add "HKCR\DOSE\DefaultIcon" /f /ve /t REG_SZ /d "%CD%\DoseViewer.exe,1">nul 2>nul 
reg add "HKCR\DOSE\shell\open(&O)\command" /f /ve /t REG_SZ /d "%CD%\DoseViewer.exe \"%%1\"">nul 2>nul 
reg add "HKCR\VIEWRAYRED\DefaultIcon" /f /ve /t REG_SZ /d "%CD%\DoseViewer.exe,3">nul 2>nul 
reg add "HKCR\VIEWRAYRED\shell\open(&O)\command" /f /ve /t REG_SZ /d "%CD%\DoseViewer.exe \"%%1\"">nul 2>nul 
popd
refresh_icon.exe
echo Association and TDR setting finished, please reboot to take effect.
pause 

:end
