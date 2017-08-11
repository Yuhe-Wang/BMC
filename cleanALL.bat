@echo off
del /f /q RunBMC\*.dll
del /f /q RunBMC\*.exe
del /f /q RunBMC\*.encrypt
rmdir /s /q RunBMC\Materials

cd ..
rmdir /s /q Build
cd Dependencies
rmdir /s /q wxWidgets-3.1\build\msw\vc_mswudll
rmdir /s /q wxWidgets-3.1\lib\vc_dll




