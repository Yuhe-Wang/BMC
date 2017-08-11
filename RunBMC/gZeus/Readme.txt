The program gZeus.exe is a GPU adaptation of the original Zeus code on windows platform. 
It produces almost identical results and runs 4.7 times faster. 
(GT650M vs i7 3630M 8 threads). If you use faster GPU, better performance is expected.

Installation:
Double click "TDR_20Seconds_import&reboot.reg" to modify the "Timeout Detection & Recovery"(TDR) from 2 seconds to 20 seconds.
This will avoid GPU driver crashing problem because some times the GPU kernel may take over TDR time to finish.
Reboot and installation is done!

Double click gZeus.exe to give a test run, and you should get a dose file with the path "BatchMode\example\gZeus.dose". The lowest GPU compute capability version required is 2.0 (Fermi), and you may need to update your GPU driver according to the promt.

Please install DoseViewer to view the generated dose.

Brief usage:
(1) gZeus.exe -> The main GPU program. It reads jobList.py and search the job configration files to execute jobs in sequence.
(2) ConfigEmail.exe -> A program to help you create encrypted email account file. If you're worried about the password in jobList.py being seen by others, this would be a better solution.
(3) ConvertBeams.exe -> A program to convert text beam configuration(.py or PlanOverview.txt) to ViewRay's binary format(.beams). If you wana do benchmark against the original Zeus, this can help provides identical inputs. 
(4) ConvertPhantom.exe -> A program to convert text phantom configuration(.py) to ViewRay's binary format(.red). If you wana do benchmark against the original Zeus, this can help provides identical inputs. 
(5) curl.exe -> Third party program being called to send emails. Don't execute it by your own.


If you have any question, please address to yuhe.wang@wustl.edu



