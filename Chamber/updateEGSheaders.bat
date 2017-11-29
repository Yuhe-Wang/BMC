@echo off
xcopy /y %HEN_HOUSE%\egs++\*.h ..\include\egs++\
xcopy /y %HEN_HOUSE%\lib\win2k-cl\egs_config1.h ..\include\egs++\
xcopy /y %HEN_HOUSE%\egs++\dso\win2k-cl\egspp.lib ..\lib\egspp.lib