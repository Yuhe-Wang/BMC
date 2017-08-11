The DoseViewer.exe is designed to view dose in my format or ViewRay's format. It also provides a GUI to control Monte Carlo simulation engines (currently gZeus and gMeshMM are available).

Installation:
(1) Extract these files to any folder you like.
(2) Navigate to the directory SetDoseViewer, and run install.bat with administrator privilege. This batch script associates several file types with this program so you can open them by just double-clicking. It also sets the TDR time to avoid GPU driver crashing due to possible long kernel time. Reboot is required. 
(3) Attention: You may have many copies of DoseViewer at different directories on one PC, but only one installation is required. You can execute the scripts in setDoseViewer\ to exact independent dose calculation packages with console UI.

Uninstallation:
(1) Run uninstall.bat with administrator privilege to cancel those associations.
(2) Delete any file belong to DoseViewer.

Brief description:
(1) DoseViewer.exe  The main program to control MC simulation and view the dose distribution.
(2) 7z.exe, curl.exe, trash.exe   Auxiliary programs that shouldn't be runned by user.
(3) ConfigEmail.exe  Help you create an encripted file for your gmail account.
(4) ConvertBeams.exe  Convert the customized beams in config file to ViewRay's format (*.beams)
(5) ConvertPhantom.exe  Convert the customized phantom in config file to ViewRay's format (*.red)
(6) gZeus.exe gMeshMM.exe  Console programs to launch MC simulations with corresponding engines. Don't rename them because they launch the engines based on their own file names. 


If you have any question, please address to yuhe.wang@wustl.edu