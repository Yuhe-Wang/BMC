@echo off
echo drag any file to this batch script
echo to copy it to each job directory

if "%1" == "" (
echo no file name specified, please don't double click
pause
goto end)

for /f "delims=" %%a in ('dir /ad /b') do (
copy "%1" "%%a"
)

:end