#pragma once
#ifndef _PENELOPE_H_
#define _PENELOPE_H_

#include "../Tools/Tools.h"
#include "Phantom.h"
#include "../SourceHead/SourceHead.h"
#include "../ZeusData/ZeusData.h"
//#include "device_functions.h"
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
#if __CUDA_ARCH__ < 350
#define USE_TEXTURE_CACHE 1
#endif

#define USE_SINGLE_PRECISION

//notice GlueF will operate on float constant or float function like sin, sqrt
#ifdef USE_SINGLE_PRECISION
typedef float ZFloat;
#define GlueF(a) a##f
#else
typedef double ZFloat;
#define GlueF(a) a
#endif
/************************Start: Constant Definitions **************************/
#define PI GlueF(3.141592653589793238463)
#define Es GlueF(5.10998902e5)
#define TEs GlueF(1.021997804e6)
#define alpha (7.297352536e-3)
#define INV_ELECTRON_MASS GlueF(1.956951337e-6)
#define FLTEPSILON GlueF(1.192092896e-07)

#define NMatMax 5 //the maximum material number
#define NWOODCOCK 4096
#define NLAMDAPHOTON 2048
#define NLAMDACOMPTON 256
#define NLAMDAPAIR 1024
#define NQSURFACE_E 256
#define NQSURFACE_Q 128
#define NSCREENINGPARAMETER 512
#define NRANGE 1024
#define NINVERSERANGE 1024
/************************End: Constant Definitions **************************/

void cudaErrorCheck(cudaError_t cudaStatus, const char* func = NULL);

void inline cudaKernelCheck(int i, int it = -1);

int ConvertSMVer2Cores(int major, int minor);

void printGPUProperties(int i);

/*************************Start: type/class definition *********************************/

class ParticleR //store status information of a particle
{
public:
	__device__ void changeByCos(ZFloat cosTheta, ZFloat phi)
	{
		ZFloat cosPhi = GlueF(cos)(phi);
		ZFloat sinPhi = GlueF(sqrt)(1 - cosPhi*cosPhi);
		if (phi > PI) sinPhi = -sinPhi;
		ZFloat dxy = u*u + v*v;
		ZFloat dxyz = dxy + w*w;
		if (GlueF(fabs)(dxyz - GlueF(1.0)) > GlueF(1e-14))
		{
			ZFloat norm = 1 / GlueF(sqrt)(dxyz);
			u *= norm;
			v *= norm;
			w *= norm;
			dxy = u*u + v*v;
		}

		if (dxy > GlueF(1.0e-28))
		{
			ZFloat sinTheta = GlueF(sqrt)((1 - cosTheta*cosTheta) / dxy);
			ZFloat up = u;
			u = u*cosTheta + sinTheta*(up*w*cosPhi - v*sinPhi);
			v = v*cosTheta + sinTheta*(v*w*cosPhi + up*sinPhi);
			w = w*cosTheta - dxy*sinTheta*cosPhi;
		}
		else
		{
			ZFloat sinTheta = GlueF(sqrt)(1 - cosTheta*cosTheta);
			if (w > 0)
			{
				u = sinTheta*cosPhi;
				v = sinTheta*sinPhi;
				w = cosTheta;
			}
			else
			{
				u = -sinTheta*cosPhi;
				v = -sinTheta*sinPhi;
				w = -cosTheta;
			}
		}
	}
	//data
	ZFloat x, y, z; //position, unit cm
	ZFloat u, v, w; //direction vector
	int ivx, ivy, ivz, iabsv; //voxel index for current particle

	ParticleType type;

	ZFloat E; //energy, unit eV
	ZFloat weight;
};

class ParticleStack
{
	//method
public:
	__device__ void push(ParticleR &par)
	{
		if (cur < 10)
		{
			ss[cur] = par; //copy and store the particle
			++cur;
		}
		else
		{
			printf("exceed the max depth of the GPU stack!!!!!!!!!!!!!!!\n");
		}
	}
	__device__ __forceinline__ bool empty()
	{
		if (cur > 0) return false;
		else return true;
	}
	__device__ __forceinline__ void pop()
	{
		if (cur <= 0) printf("no element to pop in the stack!");
		--cur;
	}
	__device__ __forceinline__ ParticleR& top(){ return ss[cur - 1]; }
	__device__ void init(ParticleR* pp){ cur = 0; ss = pp; }
	//data
private:
	ParticleR* ss;
	int cur; //point to next free space, 0 at initial time
};

// class GRNG //simple version, need test to verify it's uniform if force calling fillbuffer
// {
// public:
// 	__device__ void init(int seed) //init the generator by its index in an array or threads
// 	{
// 		long In = 1234 + seed;
// 
// 		In = 16807 * (In % 127773) - 2836 * (In / 127773);
// 		if (In < 0) In += 2147483647;
// 		ISEED1 = In;
// 		In = 16807 * (In % 127773) - 2836 * (In / 127773);
// 		if (In < 0) In += 2147483647;
// 		ISEED2 = In;
// 	}
// 	__device__ __forceinline__ ZFloat operator () () //default return random number in (0,1)
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
// 	int ISEED1;
// 	int ISEED2;
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

struct GPUConfig
{
	void createStream()
	{
		cudaSetDevice(id);
		cudaStreamCreate(&kernelstream);
		cudaStreamCreate(&copystream);
	}
	void destroyStream()
	{  
		cudaSetDevice(id);
		cudaStreamDestroy(kernelstream);
		cudaStreamDestroy(copystream);
	}
	cudaStream_t kernelstream;
	cudaStream_t copystream;
	int id; //GPU index

	//some pointers in the GPU. Will be used to release the resource after simulation is finished
	ParticleR* d_InitParsA;
	ParticleR* d_InitParsB;
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

class SourcePool //a wrapper of SourceHead
{
public:
	SourcePool(vector<GPUConfig>* gcp,int NT, int NOneBatch, int NSampleStack = 400)
	{
		_gcp = gcp;
		_b_initBuff = true;
		if (0 == NT) NT = get_thread_num();
		_NGenerateThread = NT;
		_PoolCapacity = NOneBatch / NT;

		_np = new int[NT]; //how many particle was actually held in each thread
		_hist = new int[NT]; //how many histories have been sampled for each thread
		_rng = new PRNG[NT]; // create rng for each thread

		_PPool = new ParticleR*[NT]; //create pool pointer
		_PPoolTemp = new Particle*[NT];
		for (int i = 0; i < NT; ++i)
		{
			_PPool[i] = new ParticleR[_PoolCapacity + NSampleStack]; //resize memory for each thread
			_PPoolTemp[i] = new Particle[NSampleStack];
			_np[i] = 0;
			_hist[i] = 0;
			_rng[i].init(1234 + i);//init the rng for each thread
		}
	}
	~SourcePool()
	{
		//delete particle pool
		for (int i = 0; i < _NGenerateThread; ++i)
		{
			delete[] _PPool[i];
			delete[] _PPoolTemp[i];
		}
		delete[] _PPool;
		delete[] _PPoolTemp;
		delete[] _np;
		delete[] _rng;
		delete[] _hist;

#ifdef STAT_HEAD
		statOutput();
#endif	
	}
	int prepareCopy()
	{
		for (int i = 0; i < _NGenerateThread; ++i) _hist[i] = 0;
		//Let's fill each pool with multi-thread
#ifdef USE_OPENMP
#pragma omp parallel num_threads(_NGenerateThread)
#endif
		{
			int it = omp_get_thread_num(); //thread index starting from 0
			int rest = _np[it] - _PoolCapacity;//rest particles left in the buffer part
			if (rest > 0)
			{
				for (int j = 0; j < rest; ++j) //need move the buffered particles to the left side
					_PPool[it][j] = _PPool[it][j + _PoolCapacity];
				_np[it] = rest;
			}
			else _np[it] = 0;

			while (_np[it] < _PoolCapacity)
			{
				int nadd = SourceHead_Sample(&_rng[it], _PPoolTemp[it]);
				ParticleR *pbase = _PPool[it] + _np[it];
				for (int i = 0; i < nadd; ++i)
				{
					ParticleR* p = pbase + i;
					p->E = (ZFloat)_PPoolTemp[it][i].E;
					p->u = (ZFloat)_PPoolTemp[it][i].u;
					p->v = (ZFloat)_PPoolTemp[it][i].v;
					p->w = (ZFloat)_PPoolTemp[it][i].w;
					p->x = (ZFloat)_PPoolTemp[it][i].x;
					p->y = (ZFloat)_PPoolTemp[it][i].y;
					p->z = (ZFloat)_PPoolTemp[it][i].z;
					p->type = photon;
					p->weight = (ZFloat)_PPoolTemp[it][i].weight;
				}
				_np[it] += nadd;
				++_hist[it];
			}
		}
		int sum = 0;
		for (int i = 0; i < _NGenerateThread; ++i) sum += _hist[i];
#ifdef STAT_HEAD
		statHead();
#endif
		copy2GPU2();
		return sum; //return how many histories have generated
	}
	void copy2GPU(int ig) //should be called by each GPU thread
	{
		vector<GPUConfig>& gc = *_gcp;
		ParticleR* pp = gc[ig].d_InitParsB;
		if (_b_initBuff) pp = gc[ig].d_InitParsA;

		for (int i = 0; i < _NGenerateThread; ++i)
		{
			cudaMemcpyAsync(pp + i*_PoolCapacity, _PPool[i], _PoolCapacity * sizeof(ParticleR), cudaMemcpyHostToDevice, gc[ig].copystream);
			//cudaMemcpy(pp + i*_PoolCapacity, _PPool[i], _PoolCapacity * sizeof(ParticleR), cudaMemcpyHostToDevice);
			//cudaStreamSynchronize(gc[ig].copystream);
		}
		cudaStreamSynchronize(gc[ig].copystream);
		_b_initBuff = !_b_initBuff; //use the other array next time
	}
	void copy2GPU2() //called by one thread only
	{
		vector<GPUConfig>& gc = *_gcp;
		int ng = (int)gc.size();
#ifdef USE_OPENMP
#pragma omp parallel num_threads(ng)
#endif
		{
			int ig = omp_get_thread_num();
			cudaSetDevice(gc[ig].id);
			ParticleR* pp = _b_initBuff ? gc[ig].d_InitParsA : gc[ig].d_InitParsB;

			for (int i = 0; i < _NGenerateThread; ++i)
			{
				cudaMemcpyAsync(pp + i*_PoolCapacity, _PPool[i], _PoolCapacity * sizeof(ParticleR), cudaMemcpyHostToDevice, gc[ig].copystream);
			}
			cudaStreamSynchronize(gc[ig].copystream);
		}
		_b_initBuff = !_b_initBuff; //use the other array next time
	}
	ParticleR* getAP(int ig) //must call after copy2GPU finished and before another perpareCopy is called
	{
		vector<GPUConfig>& gc = *_gcp;
		if (_b_initBuff == true) return gc[ig].d_InitParsB;
		else return gc[ig].d_InitParsA;
	}
private:
#ifdef STAT_HEAD
	void statHead();
	void statOutput();
#endif
	int _NGenerateThread;
	int _PoolCapacity; //how many particles it can hold in each line
	ParticleR** _PPool;
	Particle** _PPoolTemp;
	int *_np;   //how many particle was actually held in each line
	int *_hist; //how many histories have been sampled for each line
	PRNG *_rng; //rng array
	bool _b_initBuff; //true means pointing to InitPars_A, false means pointing to InitPars_B
	vector<GPUConfig>* _gcp; //pointer to GPU config
};

class GFREE //class wrapper of the simulation kernel
{
public:
	vector<GPUConfig> gc;
	int getGPUConfig(ConfigFile* gcf);
	void init(ConfigFile *cf);

	void GRNGStat();

	GFREE()
	{
		_phant = NULL;
	}
	~GFREE()
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
		name##TextureData.normalized = false; \
		name##TextureData.addressMode[0] = cudaAddressModeClamp; \
		name##TextureData.addressMode[1] = cudaAddressModeClamp; \
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

#endif