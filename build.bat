@echo off
cl 1>nul 2>nul
if errorlevel 1 (
    echo The build script must run under Visual Studio promt!
	echo.
	pause
	goto end
)

echo Weclome to building process of BMC!
echo.
echo You may want to modify the script to customize build of PINIT
echo Or change the cmake generator
echo.

set GEN="Visual Studio 14 2015 Win64"
set GEN="Visual Studio 15 2017 Win64"
set GEN="NMake Makefiles"
Rem it can be Debug or Release
set BuildType=Release
set usePreCompile=yes

echo Building dependencies with nmake:
echo. 

cd ..
if not exist Build mkdir Build
cd Build

set proj="PINIT"
echo Building %proj%...
if not exist %proj% (
    mkdir %proj%
	if %usePreCompile%==yes (
		mkdir bin
		copy ..\Dependencies\%proj%\%proj%.dll bin\
	) else gfortran ..\Dependencies\%proj%\PInterface.f ..\Dependencies\%proj%\penelope.f -shared -fPIC -m64 -o bin\%proj%.dll -O3 -fno-underscoring -static-libgcc -static-libgfortran 
)
echo %proj% building done!

rem compile cmake dependencies
set dependsList=mathgl-2.3.2 dcmtk-3.6.1 zlib-1.2.11 szip-2.1.1 hdf5-1.8.19 matio-1.5.10
for %%i in (%dependsList%) do (
	if not exist %%i/xxx (
		mkdir %%i
		cd %%i
		cmake -G %GEN% -D CMAKE_BUILD_TYPE=%BuildType% ../../Dependencies/%%i
		nmake install
		rem for %%j in (*.sln) do Devenv %%j /Build "%BuildType%" /Project INSTALL
		cd ..
	)
)


set proj="wxWidgets-3.1"
if not exist %proj% (
    mkdir %proj%
    cd ../Dependencies/wxWidgets-3.1/build/msw
	if %BuildType%==Release (
		nmake BUILD=release /f makefile.vc
	)else (nmake BUILD=debug /f makefile.vc)
	cd ../..
	rem cd = wxWidgets-3.1
    xcopy /e /y include ..\..\Build\include\wx-3.1\
	xcopy /y lib\vc_dll\*.dll ..\..\Build\bin\
	xcopy /y lib\vc_dll\*.lib ..\..\Build\lib\
	xcopy /e /y lib\vc_dll\mswu ..\..\Build\lib\mswu\
    cd ../../Build
)
echo %proj% building done!

set GEN="Visual Studio 15 2017 Win64"
set proj="BMC"
rem if not exist %proj% (
    mkdir %proj%
    cd %proj% 
	if %BuildType%==Debug set BuildType=MyDebug
    cmake -G %GEN% -D CMAKE_BUILD_TYPE=%BuildType% ../../%proj%
	rem for compile purpose only
    rem nmake install
	for %%j in (*.sln) do Devenv %%j /Build "%BuildType%" /Project INSTALL
    cd ..
rem )
echo %proj% building done!

cd ..
copy Build\bin\*.dll BMC\RunBMC

:end


