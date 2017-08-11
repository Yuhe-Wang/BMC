echo Please make sure boost, 7z and curl were installed!
#sudo apt-get install p7zip-full #install 7z:  
#sudo apt-get install curl #install curl:  
#sudo apt-get install libboost-all-dev #install boost: 

echo compiling ConfigEmail
g++ -Wno-unused-result ../ConfigEmail/ConfigEmail.cpp -Wno-unused-result -std=c++11 -O3 -o ConfigEmail

echo compiling PINIT
gfortran ../PINIT/PInterface.f ../PINIT/penelope.f -fPIC -shared -o PINIT.so -O3 -fno-underscoring

echo compiling SourceHead
g++ -Wno-unused-result -fopenmp ../SourceHead/SourceHead.cpp -std=c++11 -O3 -fPIC -shared -o libSourceHead.so

echo compiling Tools
g++ -Wno-unused-result -fopenmp ../Tools/lz4.cpp ../Tools/quicklz.cpp ../Tools/Tools.cpp -std=c++11 -O3 -fPIC -shared -o libTools.so

nvccpath=/usr/local/cuda-7.5/bin/
outputName=${PWD##*/}

for icppf in *.cpp
do
  name=$(ls $icppf | cut -d. -f1)
  echo "Compling ${icppf}"
  #${nvccpath}nvcc -O3 -Xcompiler -fPIC -Xcompiler -Wno-unused-result -Xcompiler -fopenmp -std=c++11 -gencode arch=compute_20,code=sm_20  -odir "." -M -o "${name}.d" "${icppf}"
  ${nvccpath}nvcc -O3 -Xcompiler -fPIC -Xcompiler -Wno-unused-result -Xcompiler -fopenmp -std=c++11 --compile  -x c++ -odir "." -o "${name}.o" "${icppf}"
  linkList="${linkList} ${name}.o"
  #echo $linkList
done

for icuf in *.cu
do
  name=$(ls $icuf | cut -d. -f1)
  echo "Compling ${icuf}"
  #${nvccpath}nvcc -O3 -Xcompiler -fPIC -Xcompiler -Wno-unused-result -Xcompiler -fopenmp -std=c++11 -gencode arch=compute_20,code=sm_20  -odir "." -M -o "${name}.d" "${icuf}"
  ${nvccpath}nvcc -O3 -Xcompiler -fPIC -Xcompiler -Wno-unused-result -Xcompiler -fopenmp -std=c++11 --compile --relocatable-device-code=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20  -x cu -odir "." -o  "${name}.o" "${icuf}"
  linkList="${linkList} ${name}.o"
done

echo "Linking to generate the target file: ${outputName}"
${nvccpath}nvcc --cudart static -L"./" --relocatable-device-code=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20 -link -o "${outputName}"  ${linkList} -lTools -lboost_system -lSourceHead -lboost_filesystem -lgomp

#clean the intermediate files
#rm *.d
rm *.o

chmod 764 run.sh
echo Rember to change the path of PENELOPE_DataBase.7z and ViewRayCo.phsp