# BMC
GPU accelerated Monte Carlo simulation platform for Photon Therapy

# Dependency list:
PINIT
mathgl-2.3.2 
dcmtk-3.6.1 
zlib-1.2.11 
szip-2.1.1 
hdf5-1.8.19 
matio-1.5.10

# Build guide:
(1) Download or git clone this BMC project
(2) Download or git clone my another porject named Dependencies, and place it in the same folder of BMC. It should look like this:
MyProject
|__BMC
|__Dependencies
(3) Install the required tools for Windows or Linux platform, 
	Linux (Ubuntu 16.04.3 x64)
	0) sudo apt install clang libomp-dev
	   sudo update-alternatives --config c++
	   sudo update-alternatives --config cc
	1) gfortran
	2) cmake
	3) CUDA Toolkit 8 or 9 (8 tested on 14.04, 9 tested on 16.04)

		`sudo dpkg -i cuda-repo-ubuntu1604-9-0-local-rc_9.0.103-1_amd64.deb`
		`sudo apt-key add /var/cuda-repo-<version>/7fa2af80.pub`
		`sudo apt-get update`
		`sudo apt-get install cuda`

	4) libgtk3-dev
	5) libmatio-dev
	6) p7zip-full (for gMeshMM to call 7z) & curl (to send email when job is finished)

	Windows (Win10 x64) (PINIT precompiled by ICC)
	1) Visual Studio 2015/2017 professional
	2) CUDA toolkit 8 (VS 2015) /9 (VS2017)
	3) cmake (>= 3.8.2), add to %path%
	4) MinGW gfortran x64 (PINIT.dll was precompiled by ICC. You can also use MinGW gfortran to regenerate it)
	5) 7z.exe (for gMeshMM to call 7z) & trash.exe have been preincluded(for DoseViewer to call) & curl (to send email when job is finished)
	6) If you are using VS 2017, you will run into problems compiling cuda programs. You can create a bin directory in C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\
	, and create a symbol link named cl.exe in it pointing to C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\VC\Tools\MSVC\14.10.25017\bin\HostX64\x64\cl.exe


(4) For Windows, open the Visual studio promt(your may need to create it by yourself due to version and installation variance. please change working directory to BMC) in the folder BMC, and type in build.bat
	For Linux, open the terminal in the folder BMCm and type in ./build.sh (you may need to do chmod +x build.sh first)

# Run BMC
The executable files and running configurations are all in folder BMC/RunBMC

Since github restricts the largest uploading file to 100M, I have to compress ViewRayCo.phsp by 7zip.

You should extract it by 7zip and put it at the same place of ViewRayCo.7z before running BMC simulations.

In Windows, you should also run installDoseViewer.bat and installTetViewer.bat in RunBMC\Setup to setup necessary association and copy 7z.exe and trash.exe to parent directory

In both platform, the folder RunBMC is portable. You can move it to wherever you want.

# Disclaim

BMC is maily develeped on Windows platform, and Linux platform hasn't been thoroughly tested. Basically, most of the functions run OK. If you run into any bug, please report to yuhewang.ustc@gmail.com.
I'll handle it as soon as possible. Thank you!



