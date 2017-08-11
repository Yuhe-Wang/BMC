The program gMesh.exe is a GPU adaptation of the original PENELOPE code on ViewRay MRIDian platform. 
It produces identical results and runs much faster. If you use faster GPU, better performance is expected.

Installation:
(1) Exact all files to any directory you like
(2) Double click "TDR_20Seconds_import&reboot.reg" to modify the "Timeout Detection & Recovery"(TDR) from 2 seconds to 20 seconds. It will avoid GPU driver crashing problem because some times the GPU kernel may take more than TDR time to finish.
Reboot is required!
(3) Additionally, you need two data files: ViewRayCo.phsp and PENELOPE_DataBase.7z. They're not included in this package because they're shared by gZeus.exe (another faster GPU MC code). Put them in any directory. Remember you need to specify this path in your configuration file when launching a simualtion.
(4) Read the comments in the example *.py file in BatchMode\example to decide the config parameters. Double click gMeshMM.exe to launcing the simulation.

Please install DoseViewer to view the dose result

File Discription:
(1) gMesh.exe, gMeshMM.exe -> The main GPU programs. They read jobList.py and search the job configration files to execute simulation jobs in sequence. gMesh.exe is optimized only for single material, while gMeshMM.exe can handle multiple materials.
(2) ConfigEmail.exe -> A program to help you create encrypted email account file. If you're worried about the password in jobList.py being seen by others, this would be a better solution.
(3) ConvertPhantom.exe -> A program to convert the customized phantom in config file to ViewRay's format.(Drag the config file to this program)
(4) ConvertBeams.exe -> A program to convert the customized beams in config file to ViewRay's format.(Drag the config file to this program)
(5) 7z.exe -> Third party program being called to exact necessary data from PENELOPE_DataBase.7z (Don't run by yourself)
(6) curl.exe -> Third party program being called to send emails. Don't execute it by your own.(Don't run by yourself)

If you have any question, please address to yuhe.wang@wustl.edu


