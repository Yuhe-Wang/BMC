#!/bin/sh
echo welcome to the building process of BMC
echo You may want to modify the build type
echo.
cd ..
BuildType="Release"
if ! test -d ./Build
then
    mkdir Build
    mkdir Build/lib
fi
cd Build

proj="PINIT"
if ! test -d ./${proj}
then
    mkdir $proj
    gfortran ../Dependencies/${proj}/PInterface.f ../Dependencies/${proj}/penelope.f -shared -fPIC -m64 -o lib/PINIT.so -O3 -fno-underscoring
fi


proj="mathgl-2.3.2"
if ! test -d ./${proj}
then
    mkdir $proj 
    cd $proj 
    cmake -D CMAKE_BUILD_TYPE=$BuildType ../../Dependencies/${proj}
    make  -j 8 install
    cd ..
fi

proj="dcmtk-3.6.1"
if ! test -d ./${proj}
then
    mkdir $proj 
    cd $proj 
    cmake -D CMAKE_BUILD_TYPE=$BuildType ../../Dependencies/${proj}
    make -j 8 install
    cd ..
fi

proj="wxWidgets-3.1"
if ! test -d ./${proj}
then
    mkdir $proj 
    cd $proj
    chmod +x ../../Dependencies/${proj}/configure
    installpath=$(readlink -m ../../Build)
	if ${BuildType}==Debug
	then
		wxflag='--enable-debug'
	else
		wxflag=''
	fi
    ../../Dependencies/${proj}/configure --with-opengl --prefix=$installpath --enable-debug ${wxflag}
    make -j 8 install
    cd ..
fi

proj="BMC"
#if ! test -d ./${proj}
#then
    mkdir $proj 
    cd $proj
	if ${BuildType}==Debug 
	then BuildType=MyDebug
	fi
    cmake -D CMAKE_BUILD_TYPE=$BuildType ../../${proj}
    make -j 8
    cd ..
#fi

cd lib
cp $(readlink *.so) ../../BMC/RunBMC
cd ../../BMC/RunBMC
chmod +x *.app
cd .. # cd = BMC
