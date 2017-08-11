@echo off
echo Going to clean *.dose *.log *.skip in each job directory

for /f "delims=" %%a in ('dir /ad /b') do (
del "%%a\*.dose"
del "%%a\*.log"
del "%%a\*.skip"
)

:end
echo.
echo Cleaning Finished!
echo.