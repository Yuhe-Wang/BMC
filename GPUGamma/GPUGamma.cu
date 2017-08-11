
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "GPUGamma.h"

#define at1(x,y,z) ((x)+NX1*((y)+NY1*(z))) //this should be the same as ArrayManager


__constant__ int NX, NY, NZ;
__constant__ int NSX, NSY, NSZ;
__constant__ double dx, dy, dz;
__constant__ double DTA2;
__constant__ double percent;
__constant__ double maxDose;
__constant__ int SK, SB;
__constant__ SFloat* gd1;
__constant__ SFloat* gd2;
__constant__ SFloat* gg;
__constant__ bool localDelta;
__device__ int nPass;

__global__ void gamma(size_t off_set, size_t up_bound)
{
	size_t it = off_set + blockIdx.x * blockDim.x + threadIdx.x;
	if (it >= up_bound) return;
	//let's decompose the it to ix,iy,iz; // it should have the same memory layout order as ArrayMgr
	int iz = it / (NX*NY);
	int irest = it - iz*NX*NY;
	int iy = irest / NX;
	int ix = irest- iy*NX;
	
	double deltaDose2 = percent*maxDose;
	deltaDose2 *= deltaDose2; //default relative to maxDose
	SFloat min_g = 1e9;

	//search around this voxel
	for (int jx = -NSX; jx <= NSX; ++jx)
		for (int jy = -NSY; jy <= NSY; ++jy)
			for (int jz = -NSZ; jz <= NSZ; ++jz)
			{
				//get the index of the searched voxel
				int six = SK*ix + jx;
				int siy = SK*iy + jy;
				int siz = SK*iz + jz;
				int NX1 = SK*NX + SB;
				int NY1 = SK*NY + SB;
				int NZ1 = SK*NZ + SB;
				if (0 <= six&&six < NX1 && 0 <= siy&&siy < NY1 && 0 <= siz&&siz < NZ1)   //index inside the matrix
				{
					double dist2 = (jx*dx*jx*dx + jy*dy*jy*dy + jz*dz*jz*dz) / (SK*SK);
					if (dist2 < DTA2) //inside the sphere, calculate gamma
					{
						double diff2 = gd1[at1(six, siy, siz)] - gd2[it];
						diff2 = diff2*diff2;
						if (localDelta)
						{
							deltaDose2 = percent*gd2[it];
							deltaDose2 *= deltaDose2;
						}
						SFloat gi = (SFloat)sqrt(dist2 / DTA2 + diff2 / deltaDose2);
						if (gi < min_g) min_g = gi;
					}
				}
			}
	//found the min_g  => count pass number
	if (min_g <= 1) atomicAdd(&nPass, 1);
	gg[it] = min_g;
}

int GPUGammaAnalysis(ArrayMgr<SFloat>& d1, ArrayMgr<SFloat>& d2, ArrayMgr<SFloat>& g, double hdx, double hdy, double hdz, double hpercent, double DTA, bool local, int NInterp)
{
	if (!d1.equalDims(d2) || !d1.equalDims(g)) return -1; //dimension check
	ArrayMgr<SFloat> di1;
	if (!interp3N(d1, di1, NInterp)) return -1;

	g = 0;//make sure the initial content is zero

	//d2 is treated as basement, and scan di1 to find matches
	int hnx = d2.getWidth(1);
	int hny = d2.getWidth(2);
	int hnz = d2.getWidth(3);
	if (0 == DTA) DTA = 0.01; //as denominator, it shouldn't be zero

	int hNSX = int(DTA / hdx * NInterp);
	int hNSY = int(DTA / hdy * NInterp);
	int hNSZ = int(DTA / hdz * NInterp);
	double hDTA2 = DTA*DTA;
	int arrayLength = hnx*hny*hnz;

	SFloat maxv, minv;
	d2.getMaxMin(maxv, minv);
	double hmaxv = maxv;

	int hSK = NInterp;
	int hSB = 1 - NInterp;

	int NT = GPUAvailable();
	int hNPass = 0;
	for (int it = 0; it < NT; ++it)
	{
		cudaSetDevice(it);
		//copy these constants to GPU
		cudaMemcpyToSymbol(NX, &hnx, sizeof(int));
		cudaMemcpyToSymbol(NY, &hny, sizeof(int));
		cudaMemcpyToSymbol(NZ, &hnz, sizeof(int));
		cudaMemcpyToSymbol(NSX, &hNSX, sizeof(int));
		cudaMemcpyToSymbol(NSY, &hNSY, sizeof(int));
		cudaMemcpyToSymbol(NSZ, &hNSZ, sizeof(int));
		cudaMemcpyToSymbol(dx, &hdx, sizeof(double));
		cudaMemcpyToSymbol(dy, &hdy, sizeof(double));
		cudaMemcpyToSymbol(dz, &hdz, sizeof(double));
		cudaMemcpyToSymbol(SK, &hSK, sizeof(int));
		cudaMemcpyToSymbol(SB, &hSB, sizeof(int));
		cudaMemcpyToSymbol(percent, &hpercent, sizeof(double));
		cudaMemcpyToSymbol(DTA2, &hDTA2, sizeof(double));
		cudaMemcpyToSymbol(maxDose, &hmaxv, sizeof(double));
		cudaMemcpyToSymbol(nPass, &hNPass, sizeof(int));//initiate it to zero
		cudaMemcpyToSymbol(localDelta, &local, sizeof(bool));
		//resize memory space for di1, d2, g
		
		SFloat* hd1 = NULL;
		SFloat* hd2 = NULL;
		SFloat* hg = NULL;
		if (cudaSuccess != cudaMalloc(&hd1, di1.getInnerLength()* sizeof(SFloat))) return -1;
		cudaMalloc(&hd2, arrayLength* sizeof(SFloat));
		cudaMalloc(&hg, arrayLength* sizeof(SFloat));
		cudaMemcpy(hd1, di1.getP(), di1.getInnerLength()* sizeof(SFloat), cudaMemcpyHostToDevice);//copy the content to GPU
		cudaMemcpy(hd2, d2.getP(), arrayLength* sizeof(SFloat), cudaMemcpyHostToDevice);
		cudaMemcpy(hg, g.getP(), arrayLength* sizeof(SFloat), cudaMemcpyHostToDevice); //make sure the initial value are zero
		cudaMemcpyToSymbol(gd1, &hd1, sizeof(SFloat*)); //copy the pointers to GPU
		cudaMemcpyToSymbol(gd2, &hd2, sizeof(SFloat*));
		cudaMemcpyToSymbol(gg, &hg, sizeof(SFloat*));
	}

	//call GPU kernels with openMP
	const int blockSize = 256;
	const int gridSize = 128;
	const int NAverageLoad = arrayLength / NT;
#pragma omp parallel num_threads(NT)
	{
		int it = omp_get_thread_num();
		cudaSetDevice(it); //set which GPU this thread will operate on
		int NWorkLoad = NAverageLoad;
		if (it == NT - 1) NWorkLoad = arrayLength - (NT - 1)*NAverageLoad;//the last thread may have more workload
		int NBatch = NWorkLoad / (blockSize*gridSize) + 1;
		int up_bound = (it+1)*NAverageLoad;
		if (it == NT - 1) up_bound = arrayLength;
		for (int ib = 0; ib < NBatch; ++ib)
		{
			gamma << <gridSize, blockSize >> >(it*NAverageLoad + ib*blockSize*gridSize, up_bound);
			cudaDeviceSynchronize();
		}
	}
	
	
// 	string err = cudaGetErrorString(cudaGetLastError());
// 	printf("cuda error: %s\n", err.c_str());
	
	int totNPass = 0;
	SFloat* hd1 = NULL;
	SFloat* hd2 = NULL;
	SFloat* hg = NULL;
	ArrayMgr<SFloat> gtemp = g;
	for (int it = 0; it < NT; ++it)
	{
		cudaSetDevice(it);
		//get the pointer address in order to free them
		cudaMemcpyFromSymbol(&hd1, gd1, sizeof(SFloat*));
		cudaMemcpyFromSymbol(&hd2, gd2, sizeof(SFloat*));
		cudaMemcpyFromSymbol(&hg, gg, sizeof(SFloat*));
		//copy g back to CPU memory
		cudaMemcpy(gtemp.getP(), hg, arrayLength* sizeof(SFloat), cudaMemcpyDeviceToHost);
		cudaMemcpyFromSymbol(&hNPass, nPass, sizeof(int));
		g += gtemp;
		totNPass += hNPass;
		
		//free the allocated memory
		cudaFree(hd1);
		cudaFree(hd2);
		cudaFree(hg);
	}
	return totNPass;
}

int GPUAvailable()
{
	int devCount = 0;
	cudaGetDeviceCount(&devCount);
	return devCount;
}
