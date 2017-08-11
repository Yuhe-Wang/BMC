Exact the package to anywhere you like and check the "Install" directory:

Install.bat will set the environment path for 3DMatViewer.exe and copy showMat.m to %USERPROFILE%\Documents\MATLAB
you need to restart MATLAB after the installation. Try showMat(rand(100,100,100)) to see if it works.

Uninstall.bat will delete the environment path for 3DMatViewer.exe and delete showMat.m as well.