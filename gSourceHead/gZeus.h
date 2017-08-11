#pragma once
#ifndef _PENELOPE_H_
#define _PENELOPE_H_

#include "../Tools/Tools.h"
#include "config.h"
#include "Phantom.h"
#include "curand_kernel.h"

//resolve the compiler's warning
#ifndef __CUDACC__
#undef __device__ 
#undef __host__
#undef __forceinline__
#endif
#include "cuda_runtime.h"
//#define USE_MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

//#define DEPOSIT_ENERGY_IN_PHOTON_INTERACTION 1
//#define USE_INTERPOLATED_LOG 1


/************************Start: Constant Definitions **************************/
#define PI GlueF(3.141592653589793238463)
#define Es GlueF(5.10998902e5)
#define TEs GlueF(1.021997804e6)
#define alpha (7.297352536e-3)
#define INV_ELECTRON_MASS GlueF(1.956951337e-6)
#define FLTEPSILON GlueF(1.192092896e-07)


/************************End: Constant Definitions **************************/

// class GRNG //simple version, need test to verify it's uniform if force calling fillbuffer
// {
// public:
// 	__device__ void init(int seed) //init the generator by its index in an array or threads
// 	{
// 		ISEED1 = 1234 + 2*seed;
// 		ISEED2 = 1234 + 2 * seed + 1;
// 	}
// 	__device__ ZFloat operator () () //default return random number in (0,1)
// 	{
// 		const double USCALE = 1.0 / 2.147483563e9;
// 		long I1 = ISEED1 / 53668;
// 		ISEED1 = 40014 * (ISEED1 - I1 * 53668) - I1 * 12211;
// 		if (ISEED1 < 0) ISEED1 += 2147483563;
// 
// 		long I2 = ISEED2 / 52774;
// 		ISEED2 = 40692 * (ISEED2 - I2 * 52774) - I2 * 3791;
// 		if (ISEED2 < 0) ISEED2 += 2147483399;
// 
// 		long IZ = ISEED1 - ISEED2;
// 		if (IZ < 1) IZ += 2147483562;
// 		return ZFloat(IZ*USCALE);
// 	}
// 	__device__ void printStatus()
// 	{
// 		printf("IS1 = %d, IS2 = %d\n", ISEED1, ISEED2);
// 	}
// private:
// 	long ISEED1;
// 	long ISEED2;
// };

class GRNG
{
public:
	__device__ void init(int seed) //init the generator by its index in an array or threads
	{
		curand_init(1234, seed, 0, &state);
	}
	__device__ __forceinline__ ZFloat operator () () //default return random number in (0,1)
	{
#ifdef USE_SINGLE_PRECISION
		return curand_uniform(&state);
#else
		return curand_uniform_double(&state);
#endif	
	}
private:
	curandState state;
};

/*************************Start: type/class definition *********************************/
struct GPUConfig
{
	void addGPUPointer(void* p)
	{
		gpu2free.push_back(p);
	}
	void freeGPU()
	{
		cudaSetDevice(id);
		int n = (int)gpu2free.size();
		for (int i = 0; i < n; ++i) cudaFree(gpu2free[i]);
		gpu2free.resize(0);
	}

	int id; //GPU index

	vector<void*> gpu2free;

	//some pointers in the GPU. Will be used to release the resource after simulation is finished
	ParticleR *d_stackBuff;
	GRNG* d_RNGState;
	//pointers for phantom
	SFloat* d_ph;
	SFloat* d_doseScore;

	//these are most likely the same
	int NBlock, BlockSize, NBatch;
	int SourceReuseTimes;
	int refillPeriod;

	double hist; //accumulate how many histories has be simulated in this GPU

	int getNOneBatch(){ return NBlock*BlockSize*NBatch; }
};



class GZEUS //class wrapper of the simulation kernel
{
public:
	vector<GPUConfig> gc;
	int getGPUConfig(ConfigFile* gcf);
	void init(ConfigFile *cf);

	void GRNGStat();

	GZEUS()
	{
		_phant = NULL;
	}
	~GZEUS()
	{
		if(_phant) delete _phant;
		freeGPU();
	}

	Phantom* _phant;

private:

	void phantom2GPU();

	void initGPU();

	void freeGPU(); 

	ConfigFile* _cf;
};

//create a object of name(type *p, int ns) in each GPU thread, and access the data by calling get##name(int i). type is float or double
#define Texture1Ddata(type, name) \
texture<type, cudaTextureType1D, cudaReadModeElementType> name##TextureData;\
class name \
{ \
	type* cuArray; \
public: \
	name(const type* p, int ns) \
	{ \
		cuArray = NULL; \
		cudaErrorCheck(cudaMalloc(&cuArray, ns* sizeof(type))); \
		cudaErrorCheck(cudaMemcpy(cuArray, p, ns * sizeof(type), cudaMemcpyHostToDevice)); \
		cudaErrorCheck(cudaBindTexture(NULL, name##TextureData, cuArray, name##TextureData.channelDesc, ns * sizeof(type))); \
		name##TextureData.normalized = false; \
		name##TextureData.addressMode[0] = cudaAddressModeClamp; \
		name##TextureData.addressMode[1] = cudaAddressModeClamp; \
	} \
	~name() \
	{ \
		cudaErrorCheck(cudaUnbindTexture(&name##TextureData)); \
		if (cuArray) cudaErrorCheck(cudaFree(cuArray)); \
	} \
}; \
__device__ __forceinline__ type get##name(int i) \
{ \
	return tex1Dfetch(name##TextureData, i); \
} 

//create a object of name(type *p, int ns) in each GPU thread, and access the data by calling get##name(int i). type is float or double
#define Texture2DFloat(name) \
texture<float, cudaTextureType2D, cudaReadModeElementType> name##TextureData;\
class name \
{ \
	cudaArray* cuArray; \
public: \
	name(const float* p, int width, int height) \
	{ \
		cuArray = NULL; \
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0,cudaChannelFormatKindFloat); \
		cudaErrorCheck(cudaMallocArray(&cuArray, &channelDesc, width, height)); \
		cudaErrorCheck(cudaMemcpyToArray(cuArray, 0, 0, p, width * height * sizeof(float), cudaMemcpyHostToDevice)); \
		cudaErrorCheck(cudaBindTextureToArray(name##TextureData, cuArray, channelDesc)); \
	} \
	~name() \
	{ \
		cudaErrorCheck(cudaUnbindTexture(&name##TextureData)); \
		if (cuArray) cudaErrorCheck(cudaFreeArray(cuArray)); \
	} \
}; \
__device__ __forceinline__ float get##name(float i, float j) \
{ \
	return tex2D(name##TextureData, i , j); \
} 
/*************************End: type/class definition *********************************/

void cudaErrorCheck(cudaError_t cudaStatus, const char* func = NULL);

void inline cudaKernelCheck(int i, int it = -1);

int ConvertSMVer2Cores(int major, int minor);

void printGPUProperties(int i);

#endif