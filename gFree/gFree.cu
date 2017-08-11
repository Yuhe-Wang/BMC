#include "gFree.h"

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: variables in device memory <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
__constant__ int NBatch;
__constant__ int NStackDepth;
__constant__ int FixedSplit = 1;
__constant__ int NMaxSplit = 50;
__constant__ int SIMU_ELECTRON = 0;

__constant__ ZFloat EAbsPhoton = GlueF(50e3);
__constant__ ZFloat EAbsElectron = GlueF(50e3);
__constant__ ZFloat ERangeCut = GlueF(10e3);
__constant__ ZFloat EMaxCSDA = GlueF(200e3);

#include "WaterCS.h"
//}}

//{{ pointers to accelerating access in the device
__constant__ ParticleR* InitPars;
__constant__ ParticleR* InitParsA;
__constant__ ParticleR* InitParsB;
__constant__ ParticleR* StackBuff;//memory pointer for stacks
__constant__ GRNG* RNGState;
//}}

//{{ data for phantom
__constant__ int NX, NY, NZ; //voxel number
__constant__ ZFloat DX, DY, DZ; // voxel size, unit cm
__constant__ ZFloat InvDX, InvDY, InvDZ;
__constant__ ZFloat LX, LY, LZ; // side length Lx=DX*NX
__constant__ ZFloat xo, yo, zo;
__constant__ ZFloat MaxDensity;
__constant__ ZFloat Bx, By, Bz; //unit magnetic field direction
__constant__ ZFloat rf;
__constant__ int uniform;
__constant__ SFloat* doseScore; //pointer to dose counter
__constant__ SFloat* ph; //pointer to phantom

Texture1Ddata(SFloat, Phant)
//}}
/*>>>>>>>>>>>>>>>>>>>>>>>>> end: variables in device memory >>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#if __CUDACC_VER_MAJOR__ < 8
//just provide double float version of atomicAdd in case we need it.
__device__ __forceinline__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
			__double_as_longlong(val +
			__longlong_as_double(assumed)));
		// Note: uses integer comparison to avoid hang in case of NaN (since NaN !=NaN)
	} while (assumed != old);
	return __longlong_as_double(old);
}
#endif

__device__ void exitKernel(const char inf[])
{
	printf("error: %s\n\n", inf);
	asm("trap;");
}
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: phantom method definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

__device__ __forceinline__ ZFloat getDensity(int iabsv)
{
	if (uniform) return MaxDensity;
	else return getPhant(iabsv);
}

__device__ __forceinline__ ZFloat getDensity(ParticleR& p)
{
	if (uniform) return MaxDensity;
	else return getPhant(p.iabsv);
}

__device__ __forceinline__ void deposit(ParticleR& p, ZFloat DE)
{
	atomicAdd(doseScore + p.iabsv, DE);
}

__device__ __forceinline__ bool lineInPhantom(ParticleR& p)
{
	//first translate the coordinate system
	p.x -= xo;
	p.y -= yo;
	p.z -= zo;
	//assuming the incident particle go straight line, judge if it will enter the phantom.
	//if true, give the interaction position
	const ZFloat Delta = GlueF(1e-5);
	if (p.x < 0 || p.x >= LX || p.y < 0 || p.y >= LY || p.z < 0 || p.z >= LZ) //initial position lays outside the phantom
	{
		if (p.x < 0)
		{
			//now it's outside the phantom	
			if (p.u > 0)
			{
				ZFloat t = p.x / p.u;
				p.x = 0;
				p.y -= t*p.v;
				p.z -= t*p.w;
				if (0 <= p.y &&p.y < LY && 0 <= p.z&&p.z < LZ) return true;
			}
			else return false;
		}
		else if (p.x >= LX)
		{
			if (p.u < 0)
			{
				ZFloat t = (LX - Delta - p.x) / p.u;
				p.x = LX - Delta;
				p.y += t*p.v;
				p.z += t*p.w;
				if (0 <= p.y &&p.y < LY && 0 <= p.z&&p.z < LZ) return true;
			}
			else return false;
		}

		if (p.y < 0)
		{
			//now it's outside the phantom
			if (p.v > 0)
			{
				ZFloat t = p.y / p.v;
				p.y = 0;
				p.x -= p.u*t;
				p.z -= p.w*t;
				if (0 <= p.x &&p.x < LX && 0 <= p.z&&p.z < LZ) return true;
			}
			else return false;
		}
		else if (p.y >= LY)
		{
			if (p.v < 0)
			{
				ZFloat t = (LY - Delta - p.y) / p.v;
				p.y = LY - Delta;
				p.x += t*p.u;
				p.z += t*p.w;
				if (0 <= p.x &&p.x < LX && 0 <= p.z&&p.z < LZ) return true;
			}
			else return false;
		}

		if (p.z < 0)
		{
			//now it's outside the phantom
			if (p.w > 0)
			{
				ZFloat t = p.z / p.w;
				p.z = 0;
				p.y -= t*p.v;
				p.x -= t*p.u;
				if (0 <= p.y &&p.y < LY && 0 <= p.x&&p.x < LX) return true;
			}
			else return false;
		}
		else if (p.z >= LZ)
		{
			if (p.w < 0)
			{
				ZFloat t = (LZ - Delta - p.z) / p.w;
				p.z = LZ - Delta;
				p.y += t*p.v;
				p.x += t*p.u;
				if (0 <= p.y &&p.y < LY && 0 <= p.x&&p.x < LX) return true;
			}
			else return false;
		}
	}
	else return true;


	return false;
}
__device__ bool lineIn(ParticleR& p)
{
	if (lineInPhantom(p)) //forward detect the voxel density to skip the air. Is it necessary?
	{
		//prepare the voxel index
		p.ivx = int(p.x*InvDX);
		p.ivy = int(p.y*InvDY);
		p.ivz = int(p.z*InvDZ);
		p.iabsv = at(p.ivx, p.ivy, p.ivz); //this will be recalculated anyway
		p.x -= p.ivx*DX;
		p.y -= p.ivy*DY;
		p.z -= p.ivz*DZ;

		return true; //now it's ready to transport
	}
	else return false;

	
}
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end: phantom method definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/



/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: Kernel definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
__global__ void initThreads(int cardOffset)//init the random number generator, call it before any run
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	RNGState[it].init(cardOffset + it);
}

__device__  __forceinline__ bool refill(size_t it, ParticleR* pInit, ParticleR& p, int& cur)
{
	while (cur < NBatch)
	{
		p = pInit[it*NBatch + cur];
		++cur;
		if (!lineIn(p)) continue;
		else return true;
	}
	return false;
}


__global__ void
//__launch_bounds__(128, 16)
gFreeRun(ParticleR* pInit) //let's implement in a simple way first
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	ParticleR p; //current particle to process
	//pStack.init(StackBuff + it*NStackDepth); //assign the stack pointer
	int initCur = 0;//reset the current particle index
	const ZFloat INF = GlueF(1e9);
	while (refill(it, pInit, p, initCur)) // keep on feeding in the particles
	{
		while (true) //iterate until current particle exits the phantom
		{
			deposit(p, GlueF(1.0));//deposit to this voxel first

			ZFloat sx, sy, sz;
			if (p.u > 0) sx = (DX - p.x) / p.u;
			else if (p.u < 0) sx = -p.x / p.u;
			else sx = INF;

			if (p.v > 0) sy = (DY - p.y) / p.v;
			else if (p.v < 0) sy = -p.y / p.v;
			else sy = INF;

			if (p.w > 0) sz = (DZ - p.z) / p.w;
			else if (p.w < 0) sz = -p.z / p.w;
			else sz = INF;

			//find the shortest distance, and update the voxel position
			if (sx <= sy && sx <= sz)
			{
				if (p.u > 0)
				{
					p.ivx += 1;
					if (p.ivx > NX - 1) break;
					p.x = 0;
					
				}
				else
				{
					p.ivx += -1;
					if (p.ivx < 0) break;
					p.x = DX;
				}
				p.y += sx*p.v;
				p.z += sx*p.z;
				p.iabsv = at(p.ivx, p.ivy, p.ivz);
			}
			else if (sy <= sx&&sy <= sz)
			{
				if (p.v > 0)
				{
					p.ivy += 1;
					if (p.ivy > NY - 1) break;
					p.y = 0;
				}
				else
				{
					p.ivy += -1;
					if (p.ivy < 0) break;
					p.y = DY;
				}
				p.x += sy*p.u;
				p.z += sy*p.w;
				p.iabsv = at(p.ivx, p.ivy, p.ivz);
			}
			else //sz the shortest
			{
				if (p.w > 0)
				{
					p.ivz += 1;
					if (p.ivz > NZ - 1) break;
					p.z = 0;
				}
				else
				{
					p.ivz += -1;
					if (p.ivz < 0) break;
					p.z = DZ;
				}
				p.x += sz*p.u;
				p.y += sz*p.v;
				p.iabsv = at(p.ivx, p.ivy, p.ivz);
			}
		}
	}
}

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end: Kernel definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


/*<<<<<<<<<<<<<<<<<<<<<<<<< start: tool functions of cuda <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
void cudaErrorCheck(cudaError_t cudaStatus,const char* func)
{
	if (cudaSuccess != cudaStatus)
	{
		if (func) Log("cuda error: %s in function module-> %s\n", cudaGetErrorString(cudaStatus), func);
		else Log("cuda error: %s\n", cudaGetErrorString(cudaStatus));
		exitApp("cuda function call failed!");
	}
}

void inline cudaKernelCheck(int i, int it)
{
#ifdef CUDA_KERNEL_CHECK
#ifdef DEBUG
	cudaDeviceSynchronize();//make sure all kernels are finished
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaSuccess != cudaStatus)
	{
		if (it != -1) Log("thread id %d with cuda kernel number %d error: %s", it, i, cudaGetErrorString(cudaStatus));
		else Log("cuda kernel number %d error: %s", i, cudaGetErrorString(cudaStatus));
		exitApp("cuda execuation error in executeJob()");
	}
#endif
#endif
}

int ConvertSMVer2Cores(int major, int minor)
{
	// Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
	typedef struct
	{
		int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
		int Cores;
	} sSMtoCores;

	sSMtoCores nGpuArchCoresPerSM[] =
	{
		{ 0x10, 8 }, // Tesla Generation (SM 1.0) G80 class
		{ 0x11, 8 }, // Tesla Generation (SM 1.1) G8x class
		{ 0x12, 8 }, // Tesla Generation (SM 1.2) G9x class
		{ 0x13, 8 }, // Tesla Generation (SM 1.3) GT200 class
		{ 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
		{ 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
		{ 0x30, 192 }, // Kepler Generation (SM 3.0) GK10x class
		{ 0x32, 192 }, // Kepler Generation (SM 3.2) GK10x class
		{ 0x35, 192 }, // Kepler Generation (SM 3.5) GK11x class
		{ 0x37, 192 }, // Kepler Generation (SM 3.7) GK21x class
		{ 0x50, 128 }, // Maxwell Generation (SM 5.0) GM10x class
		{ -1, -1 }
	};

	int index = 0;

	while (nGpuArchCoresPerSM[index].SM != -1)
	{
		if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
		{
			return nGpuArchCoresPerSM[index].Cores;
		}

		index++;
	}

	// If we don't find the values, we default use the previous one to run properly
	printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[index - 1].Cores);
	return nGpuArchCoresPerSM[index - 1].Cores;
}

int cudaGetMaxGflopsDeviceID(vector<GPUConfig>& gc) // This function returns the best GPU (with maximum GFLOPS)
{
	int current_device = 0, sm_per_multiproc = 0;
	int best_SM_arch = 0;
	int devices_prohibited = 0;

	unsigned long long max_compute_perf = 0;
	cudaDeviceProp deviceProp;

	// Find the best major SM Architecture GPU device
	for (unsigned int i = 0; i < gc.size(); ++i)
	{
		current_device = gc[i].id;
		cudaGetDeviceProperties(&deviceProp, current_device);

		// If this GPU is not running on Compute Mode prohibited, then we can add it to the list
		if (deviceProp.computeMode != cudaComputeModeProhibited)
		{
			if (deviceProp.major > 0 && deviceProp.major < 9999)
			{
				best_SM_arch = max(best_SM_arch, deviceProp.major);
			}
		}
		else
		{
			devices_prohibited++;
		}
	}

	// Find the best CUDA capable GPU device
	int gc_i = 0;
	for (unsigned int i = 0; i < gc.size(); ++i)
	{
		current_device = gc[i].id;

		cudaGetDeviceProperties(&deviceProp, current_device);

		// If this GPU is not running on Compute Mode prohibited, then we can add it to the list
		if (deviceProp.computeMode != cudaComputeModeProhibited)
		{
			if (deviceProp.major == 9999 && deviceProp.minor == 9999)
			{
				sm_per_multiproc = 1;
			}
			else
			{
				sm_per_multiproc = ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
			}

			unsigned long long compute_perf = (unsigned long long) deviceProp.multiProcessorCount * sm_per_multiproc * deviceProp.clockRate;

			if (compute_perf > max_compute_perf)
			{
				// If we find GPU with SM major > 2, search only these
				if (best_SM_arch > 2)
				{
					// If our device==dest_SM_arch, choose this, or else pass
					if (deviceProp.major == best_SM_arch)
					{
						max_compute_perf = compute_perf;
						gc_i = i;
					}
				}
				else
				{
					max_compute_perf = compute_perf;
					gc_i = i;
				}
			}
			else if (compute_perf = max_compute_perf)
			{
				if (gc[gc_i].id == 0) gc_i = i; //if GPUs' flops are identical, don't choose 0 device because the OS would use it in priority
			}
		}
	}
	return gc_i; // max_perf_device;
}  

bool speedCompareFunc(pair<unsigned long long, int> lhs, pair<unsigned long long, int> rhs)
{
	return lhs.first > rhs.first;
}
void cudaGetGflopsList(vector<int>& speedList) // This function returns the best GPU (with maximum GFLOPS)
{
	cudaDeviceProp deviceProp;

	int devCount = 0;
	cudaErrorCheck(cudaGetDeviceCount(&devCount));
	if (devCount < 1) exitApp("Cannot find any CUDA-capable GPU on your computer!");
	vector<pair<unsigned long long, int>> speed;
	// Find the best major SM Architecture GPU device
	for (int i = 0; i < devCount; ++i)
	{
		cudaGetDeviceProperties(&deviceProp, i);

		// If this GPU is not running on Compute Mode prohibited, then we can add it to the list
		if (deviceProp.computeMode != cudaComputeModeProhibited)
		{
			int sm_per_multiproc = 1;
			if (deviceProp.major == 9999 && deviceProp.minor == 9999)
			{
				sm_per_multiproc = 1;
			}
			else
			{
				sm_per_multiproc = ConvertSMVer2Cores(deviceProp.major, deviceProp.minor);
			}
			unsigned long long compute_perf = (unsigned long long) deviceProp.multiProcessorCount * sm_per_multiproc * deviceProp.clockRate;
			speed.push_back(make_pair(compute_perf, i));
		}
	}
	sort(speed.begin(), speed.end(), speedCompareFunc);
	int ng = (int)speed.size();
	speedList.resize(ng);
	for (int i = 0; i < ng; ++i) speedList[i] = speed[i].second;
}

void printGPUProperties(int i)
{
	// Get device properties
	Log("\nCUDA Device #%d\n", i);
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, i);

	Log("Major revision number:         %d\n", devProp.major);
	Log("Minor revision number:         %d\n", devProp.minor);
	Log("Name:                          %s\n", devProp.name);
	Log("Total global memory:           %lu\n", devProp.totalGlobalMem);
	Log("Total shared memory per block: %lu\n", devProp.sharedMemPerBlock);
	Log("Total registers per block:     %d\n", devProp.regsPerBlock);
	Log("Warp size:                     %d\n", devProp.warpSize);
	Log("Maximum memory pitch:          %lu\n", devProp.memPitch);
	Log("Maximum threads per block:     %d\n", devProp.maxThreadsPerBlock);
	for (int i = 0; i < 3; ++i)
		Log("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
	for (int i = 0; i < 3; ++i)
		Log("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
	Log("Clock rate:                    %d\n", devProp.clockRate);
	Log("Total constant memory:         %lu\n", devProp.totalConstMem);
	Log("Texture alignment:             %lu\n", devProp.textureAlignment);
	Log("Concurrent copy and execution: %s\n", (devProp.deviceOverlap ? "Yes" : "No"));
	Log("Number of multiprocessors:     %d\n", devProp.multiProcessorCount);
	Log("Number of total cores:     %d\n", ConvertSMVer2Cores(devProp.major, devProp.minor)*devProp.multiProcessorCount);
	Log("Kernel execution timeout:      %s\n", (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
}
/*>>>>>>>>>>>>>>>>>>>>>>>>>> end: tool functions of cuda >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/



/*<<<<<<<<<<<<<<<<<<<<<<<<< start: GFREE method definitions <<<<<<<<<<<<<<<<<<<<<<<<<<*/
void GFREE::phantom2GPU() //copy phantom info to GPU. Make sure the phantom has loaded data
{
	int NVoxel = _phant->getNVoxel();
	for (unsigned int i = 0; i < gc.size(); ++i)
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));
		//resize data in the GPU
		cudaErrorCheck(cudaMalloc(&gc[i].d_ph, NVoxel*sizeof(SFloat)),"resize ph");
		cudaErrorCheck(cudaMemcpyToSymbol(ph, &gc[i].d_ph, sizeof(SFloat *)),"copy ph pointer to GPU constant"); //set the array pointer
		cudaErrorCheck(cudaMalloc(&gc[i].d_doseScore, NVoxel*sizeof(SFloat)),"resize doseScore");
		cudaErrorCheck(cudaMemcpyToSymbol(doseScore, &gc[i].d_doseScore, sizeof(SFloat *)), "copy doseScore pointer to GPU constant"); //set the array pointer
		//set the initial value of phantom and dose counter
		cudaErrorCheck(cudaMemcpy(gc[i].d_ph, _phant->ph.getP(), sizeof(SFloat)*NVoxel, cudaMemcpyHostToDevice), "init the value of ph in GPU");
		cudaErrorCheck(cudaMemcpy(gc[i].d_doseScore, _phant->dose.getP(), sizeof(SFloat)*NVoxel, cudaMemcpyHostToDevice), "init the value of doseScore in GPU");

		//copy the rest constant
		cudaErrorCheck(cudaMemcpyToSymbol(NX, &_phant->NX, sizeof(int)), "copy NX to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(NY, &_phant->NY, sizeof(int)), "copy NY to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(NZ, &_phant->NZ, sizeof(int)), "copy NZ to GPU");
		ZFloat temp = (ZFloat)_phant->DX;
		cudaErrorCheck(cudaMemcpyToSymbol(DX, &temp, sizeof(ZFloat)), "copy DX to GPU");
		temp = (ZFloat)(1/_phant->DX);
		cudaErrorCheck(cudaMemcpyToSymbol(InvDX, &temp, sizeof(ZFloat)), "copy InvDX to GPU");
		temp = (ZFloat)_phant->DY;
		cudaErrorCheck(cudaMemcpyToSymbol(DY, &temp, sizeof(ZFloat)), "copy DY to GPU");
		temp = (ZFloat)(1 / _phant->DY);
		cudaErrorCheck(cudaMemcpyToSymbol(InvDY, &temp, sizeof(ZFloat)), "copy InvDY to GPU");
		temp = (ZFloat)_phant->DZ;
		cudaErrorCheck(cudaMemcpyToSymbol(DZ, &temp, sizeof(ZFloat)), "copy DZ to GPU");
		temp = (ZFloat)(1 / _phant->DZ);
		cudaErrorCheck(cudaMemcpyToSymbol(InvDZ, &temp, sizeof(ZFloat)), "copy InvDZ to GPU");
		temp = (ZFloat)_phant->LX;
		cudaErrorCheck(cudaMemcpyToSymbol(LX, &temp, sizeof(ZFloat)), "copy LX to GPU");
		temp = (ZFloat)_phant->LY;
		cudaErrorCheck(cudaMemcpyToSymbol(LY, &temp, sizeof(ZFloat)));
		temp = (ZFloat)_phant->LZ;
		cudaErrorCheck(cudaMemcpyToSymbol(LZ, &temp, sizeof(ZFloat)));
		temp = (ZFloat)_phant->xo;
		cudaErrorCheck(cudaMemcpyToSymbol(xo, &temp, sizeof(ZFloat)), "copy xo to GPU");
		temp = (ZFloat)_phant->yo;
		cudaErrorCheck(cudaMemcpyToSymbol(yo, &temp, sizeof(ZFloat)));
		temp = (ZFloat)_phant->zo;
		cudaErrorCheck(cudaMemcpyToSymbol(zo, &temp, sizeof(ZFloat)));
		temp = (ZFloat)_phant->Bx;
		cudaErrorCheck(cudaMemcpyToSymbol(Bx, &temp, sizeof(ZFloat)), "copy Bx to GPU");
		temp = (ZFloat)_phant->By;
		cudaErrorCheck(cudaMemcpyToSymbol(By, &temp, sizeof(ZFloat)));
		temp = (ZFloat)_phant->Bz;
		cudaErrorCheck(cudaMemcpyToSymbol(Bz, &temp, sizeof(ZFloat)));
		temp = (ZFloat)_phant->MaxDensity;
		cudaErrorCheck(cudaMemcpyToSymbol(MaxDensity, &temp, sizeof(ZFloat)));
		temp = (ZFloat)_phant->rf;
		cudaErrorCheck(cudaMemcpyToSymbol(rf, &temp, sizeof(ZFloat)));
		cudaErrorCheck(cudaMemcpyToSymbol(uniform, &_phant->uniform, sizeof(int)), "copy uniform to GPU");
	}
}

void GFREE::initGPU()
{
	RunTimeCounter rc;
	if (!initWaterCS("ZeusCrossSections.binary", gc)) exitApp("Cannot initiate the water cross sections");
	int NGPU = (int)gc.size();
	
	int split = 50, h_NStackDepth = 40;
	_cf->getValue("NMaxSplit", split);
	_cf->getValue("Particle stack depth", h_NStackDepth);
	string str = "true";
	int h_fixedSplit = 1;
	_cf->getValue("Fixed Split", str);
	if (str.compare("yes") != 0) h_fixedSplit = 0;
	int h_simuElectron = 1;
	_cf->getValue("Simulate electron",str);
	if (str.compare("yes") != 0) h_simuElectron = 0;

	double fv = 0;
	ZFloat h_EAbsPhoton = 50e3; //unit eV
	if (_cf->getValue("EAbsPhoton", fv)) h_EAbsPhoton = (ZFloat)fv;
	ZFloat h_EAbsElectron = 50e3; //unit eV
	if (_cf->getValue("EAbsElectron", fv)) h_EAbsElectron = (ZFloat)fv;
	ZFloat h_EMaxCSDA = 200e3; //unit eV
	if (_cf->getValue("EMaxCSDA", fv)) h_EMaxCSDA = (ZFloat)fv;

	for (int i = 0; i < NGPU; ++i) //for each GPU
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));

		int NGPUThread = gc[i].NBlock*gc[i].BlockSize;

		//resize initial particle memory in GPU
		cudaErrorCheck(cudaMalloc(&gc[i].d_InitParsA, NGPUThread * gc[i].NBatch * sizeof(ParticleR)));
		cudaErrorCheck(cudaMemcpyToSymbol(InitParsA, &gc[i].d_InitParsA, sizeof(ParticleR *)));
		cudaErrorCheck(cudaMalloc(&gc[i].d_InitParsB, NGPUThread * gc[i].NBatch * sizeof(ParticleR)));
		cudaErrorCheck(cudaMemcpyToSymbol(InitParsB, &gc[i].d_InitParsB, sizeof(ParticleR *)));

		//resize memory for the particle stack in GPU
		cudaErrorCheck(cudaMalloc(&gc[i].d_stackBuff, NGPUThread * h_NStackDepth* sizeof(ParticleR)));
		cudaErrorCheck(cudaMemcpyToSymbol(StackBuff, &gc[i].d_stackBuff, sizeof(ParticleR*)));
		cudaErrorCheck(cudaMemcpyToSymbol(NStackDepth, &h_NStackDepth, sizeof(int)));

		//resize memory for the GRNG status in GPU
		cudaErrorCheck(cudaMalloc(&gc[i].d_RNGState, NGPUThread* sizeof(GRNG)));
		cudaErrorCheck(cudaMemcpyToSymbol(RNGState, &gc[i].d_RNGState, sizeof(GRNG*)));

		cudaErrorCheck(cudaMemcpyToSymbol(NBatch, &gc[i].NBatch, sizeof(int)));
		cudaErrorCheck(cudaMemcpyToSymbol(NMaxSplit, &split, sizeof(int)));
		cudaErrorCheck(cudaMemcpyToSymbol(FixedSplit, &h_fixedSplit, sizeof(int)));
		cudaErrorCheck(cudaMemcpyToSymbol(SIMU_ELECTRON, &h_simuElectron, sizeof(int)));

		cudaErrorCheck(cudaMemcpyToSymbol(EAbsPhoton, &h_EAbsPhoton, sizeof(ZFloat)));
		cudaErrorCheck(cudaMemcpyToSymbol(EAbsElectron, &h_EAbsElectron, sizeof(ZFloat)));
		cudaErrorCheck(cudaMemcpyToSymbol(EMaxCSDA, &h_EMaxCSDA, sizeof(ZFloat)));
	}
	Log("\nIt costs %f seconds to init GPU ", rc.stop());
}

void GFREE::freeGPU()//do some clean up
{
	for (unsigned int i = 0; i < gc.size(); ++i)
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));
		cudaErrorCheck(cudaFree(gc[i].d_InitParsA));
		cudaErrorCheck(cudaFree(gc[i].d_InitParsB));
		cudaErrorCheck(cudaFree(gc[i].d_stackBuff));
		cudaErrorCheck(cudaFree(gc[i].d_RNGState));
		cudaErrorCheck(cudaFree(gc[i].d_ph));
		cudaErrorCheck(cudaFree(gc[i].d_doseScore));
	}
}

void GFREE::init(ConfigFile* cf)
{
	_cf = cf;
	initGPU();

	_phant = new Phantom;
	string lastDoseFile;
	bool loadLastDose = false;
	if (cf->getValue("last dose file", lastDoseFile))
	{
		loadLastDose = _phant->readDose(lastDoseFile.c_str());
	}
	ConfigFile* ph_cf = cf->getBlock("PHANTOM");
	if (ph_cf == NULL) exitApp("Cannot find the phantom configuration");
	if (!loadLastDose) _phant->loadPhantom(ph_cf);
	else Log("load last dose file successfully with %d existing histories", _phant->getHist());
	SourceHead_GetPrescrition(&(_phant->prescriptionDose), &(_phant->treatmentFraction));

	phantom2GPU();
}

int GFREE::getGPUConfig(ConfigFile* gcf)
{
	if (NULL == gcf) exitApp("cannot find GPU configuration!");
	int devCount = 0;
	cudaErrorCheck(cudaGetDeviceCount(&devCount));
	if (devCount < 1) exitApp("Cannot find any CUDA-capable GPU on your computer!");
	string GPU_Query;
	gcf->getValue("GPU Query", GPU_Query);
	if (GPU_Query.compare("yes") == 0)
	{
		int devCount;
		cudaErrorCheck(cudaGetDeviceCount(&devCount));
		Log("There are %d CUDA devices listed as follow:\n", devCount);
		for (int i = 0; i < devCount; ++i) printGPUProperties(i);
		printf("\nDo you want to continue executing GPU computation? y/n\n");
		if (getchar() != 'y') exit(0);
	}
	
	int NBlock = 128, BlockSize = 256, NBatch = 100, GRNG_Refill_Period = 70, Source_Reuse_Times = 10;
	string rngStat;
	gcf->getValue("GPU Block Num", NBlock);
	gcf->getValue("GPU Block Dim", BlockSize);
	gcf->getValue("GPU Batch Num", NBatch);
	gcf->getValue("GPU RNG Statistic", rngStat);
	gcf->getValue("GRNG Refill Period", GRNG_Refill_Period);
	gcf->getValue("Source Reuse Times", Source_Reuse_Times);
	//double GPU_Weight = 0;
	GPUConfig gpuc;
	//gpuc.id = GPU_Index;
	gpuc.NBlock = NBlock;
	gpuc.BlockSize = BlockSize;
	gpuc.NBatch = NBatch;
	gpuc.refillPeriod = GRNG_Refill_Period;
	gpuc.SourceReuseTimes = Source_Reuse_Times;

	
	vector<int> GPU_in_speed;
	cudaGetGflopsList(GPU_in_speed);

	vector<int> GPU_Index;
	if (!gcf->getValue("GPU Index", GPU_Index)) //no specific GPU index
	{
		int NGPU = 0;
		gcf->getValue("GPU Num", NGPU);
		if (NGPU <= 0) exitApp("Invalid GPU index configuration!");
		int NGPUAvailable = (int)GPU_in_speed.size();
		for (int i = 0; i < NGPU; ++i)
		{
			if (i < NGPUAvailable) GPU_Index.push_back(GPU_in_speed[i]);
			else break;
		}
	}
	
	for (unsigned int i = 0; i < GPU_Index.size(); ++i)
	{
		if (GPU_Index[i] >= 0 && GPU_Index[i] < devCount)
		{
			gpuc.id = GPU_Index[i];
			gc.push_back(gpuc);
		}
		else exitApp("Invalid GPU index");
	}

	//find the best GPU as the main thread, and optimize the work load
	int main_id = 0;
	if (gc.size() > 1)	main_id = cudaGetMaxGflopsDeviceID(gc);


	Log("/******************* The following GPU will be used ***************************/");
	for (unsigned int i = 0; i < gc.size(); ++i)	printGPUProperties(gc[i].id);
	Log("/************************ End GPU description *********************************/\n\n");

	//create streams of GPU control
	for (unsigned int i = 0; i < gc.size(); ++i) gc[i].createStream();
	return main_id;
}
/*>>>>>>>>>>>>>>>>>>>>>>>>>> end: GFREE method definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


RunTimeCounter sourceCounter;

void getSource(SourcePool* sp, volatile int* hist)
{
	sourceCounter.start();
	*hist = sp->prepareCopy();//prepare one batch before any run

	Log("time cost to generate particle = %f s", sourceCounter.stop());
}

void executeJob(const char* configFileName, MPS& configMacro) //execute one job according to the config file
{
	RunTimeCounter totalTime;
	ConfigFile cf(configFileName); //parse the total config file
	//find out where is the config file located
	string cdir(configFileName);
	size_t pos = cdir.find_last_of("\\/");
	cdir.erase(++pos, string::npos);
	cf.macroReplace(string("$cdir$"), cdir);
	cf.macroReplace(configMacro);

	string logDir, logAppend, logDescription;
	cf.getValue("log file directory", logDir);
	cf.getValue("log file append", logAppend);
	cf.getValue("log description", logDescription);
	if (0 == logAppend.compare("no")) logAppend = "w";
	else logAppend = "a";
	Log.setLogName(0, logAppend, logDir);//start log recording for this job
	Log("The job config file name = %s", configFileName);
	if (logDescription.compare("NA") != 0) Log("Short description: %s", logDescription.c_str());
	Log("Start log time = %s\n\n", Log.timeNow());

	double fNSIMU = 1e7;
	cf.getValue("NSIMU", fNSIMU);

	//search and config the GPU part ->> get gc
	GFREE zeus;
	vector<GPUConfig>& gc = zeus.gc;//makes the name shorter
	ConfigFile *gcf = cf.getBlock("GPU");
	int main_id = zeus.getGPUConfig(gcf);

	//initialize GFREE and SourceHead by configurations
	ConfigFile *scf = cf.getBlock("SOURCEHEAD");
	SourceHead_Init(scf);
	//config the source particle pool
	int NGT = 1, NGStack = 400;
	scf->getValue("NThread", NGT);
	scf->getValue("Sample Stack Depth", NGStack);
	int NOneFetch = gc[0].getNOneBatch();
	SourcePool sp(&gc, NGT, NOneFetch, NGStack);
// 	string dataDir;
// 	scf->getValue("DataDir",dataDir);
// 	if (!ZeusData_load(dataDir.c_str())) exitApp("Cannot load Zeus cross-sections correctly!");
	zeus.init(&cf); //prepare GPU data and phantom

	Log("\nCalculating dose, please wait patiently...\n\n");
	RunTimeCounter rc; //count the calculating time

	//note history number != particle number in the source pool
	const int NGPU = (int)gc.size();
	//history number generated once by the source pool, only modified by one thread,
	//but accessed by multiple threads. It can be read only after the generating thread finished.
	volatile int histNew = 0; //histNew isn't always the same as hisAdd because histNew is modified in an isolated thread
	volatile int histAdd = 0; //this variable is shared by all threads, so add the key word "volatile" for safe
	std::thread sthread; //source generating thread, unattached
	sthread = std::thread(&getSource, &sp, &histNew); //start a fetch prepare before any run
	sthread.join(); //wait for the source generating thread to end
	histAdd = histNew; //means the main GPU has histAdd new histories

// #ifdef USE_OPENMP
// #pragma omp parallel num_threads(NGPU)
// #endif
// 	{
// 		int it = omp_get_thread_num(); //thread index starting from 0
// 		//cudaErrorCheck(cudaSetDevice(gc[it].id));
// 		sp.copy2GPU(it); //copy particles to GPU
// 	}//here all GPU has got initial particles

	RunTimeCounter kernelCounter;
	RunTimeCounter copyCounter;

	

// 	omp_lock_t mainGPUMutex;
// 	omp_init_lock(&mainGPUMutex);

	vector<SFloat*> dose(NGPU);//to store dose from all GPU cards
	vector<SFloat*> uncertainty(NGPU);
	int NVoxel = zeus._phant->getNVoxel();
	SFloat* hist = new SFloat[NVoxel];
	int nBatch = 1; //count the batch number that has been done

	//each thread takes care of one GPU
#ifdef USE_OPENMP
#pragma omp parallel num_threads(NGPU)
#endif
	{
		int it = omp_get_thread_num();
		cudaErrorCheck(cudaSetDevice(gc[it].id)); //set which GPU this thread will operate on
#ifdef USE_TEXTURE_CACHE
		Phant aPhant(zeus._phant->ph.getP(), NVoxel);
		WaterQS aWaterQS((float*)WaterQSData, NQSURFACE_E, NQSURFACE_Q);
#endif
		const int NBlock = gc[it].NBlock;
		const int BlockSize = gc[it].BlockSize;
		const int Source_Reuse_Times = gc[it].SourceReuseTimes;
		double fHMax = fNSIMU / NGPU;
		double hist = 0;

		//std::mt19937 mtrng((unsigned)std::chrono::system_clock::now().time_since_epoch().count());
		//generate a random thread for each working thread
		//int seed = int((it + 1) * 12345789 * (double(mtrng()) / mtrng.max() + 0.5)); //???

		dose[it] = new SFloat[NVoxel]; //to store the final dose of this thread
		SFloat* gdose = new SFloat[NVoxel]; //to fetch temporary dose from GPU
		uncertainty[it] = new SFloat[NVoxel];
		//initialize the dose score in the CPU end
		memset(dose[it], 0, sizeof(SFloat)*NVoxel);
		memset(gdose, 0, sizeof(SFloat)*NVoxel);
		memset(uncertainty[it], 0, sizeof(SFloat)*NVoxel);

		initThreads << <NBlock, BlockSize >> >(it*NBlock*BlockSize);
		cudaKernelCheck(0);

		int source_reuse = 0; //counter about reuse source
		ParticleR* pInit = sp.getAP(it); //get the source particle pointer
#pragma omp barrier //make sure all cards get their source pointers before another generating thread is launched
		if (it == main_id) sthread = std::thread(&getSource, &sp, &histNew); //generating particles in another thread
		while(true) //calculating main loop, end when hist >= fHMax
		{
			hist += histAdd; //no mater whether we update the source, this loop will accumulate histAdd more histories

			/****************** Begin a batch run on GPU ********************/
			if (it == main_id) kernelCounter.start(); //only count the kernel time for the main GPU

			//cudaErrorCheck(cudaSetDevice(gc[it].id));//for safe
			gFreeRun << <NBlock, BlockSize,0, gc[it].kernelstream>> >(pInit);
			cudaStreamSynchronize(gc[it].kernelstream);//wait for the kernel to finish

			//print the speed information
			if (it == main_id)
			{
				Log("time cost to execute kernels = %f s", kernelCounter.stop());
				double time = rc.stop(true);
				double speed = hist / time;
				double rest = 0;
				if (fHMax > hist) rest = (fHMax - hist) / speed;
				else rest = 0;
				Log("GPU processed ------------------------ %3.1f%%,   speed = %d h/s\n", hist*100.0 / fHMax, int(speed));
				Log("Time escaped = %.1f min, left time expected = %.1f min", time / 60.0, rest / 60.0);
				++nBatch;
			}
			/****************** End a batch run on GPU *********************/

			//after one batch, we need to fetch dose from GPU to calculate the uncertainty
			cudaErrorCheck(cudaMemcpy(gdose, gc[it].d_doseScore, sizeof(SFloat)*NVoxel, cudaMemcpyDeviceToHost)); //fetch the batch dose from GPU
			SFloat minv = gdose[0];
			SFloat maxv = minv;
			for (int i = 0; i < NVoxel; ++i)
			{
				minv = min(minv, gdose[i]);
				maxv = max(maxv, gdose[i]);
				dose[it][i] += gdose[i];
				uncertainty[it][i] += gdose[i] * gdose[i];
				gdose[i] = 0;
			}
			if (it == main_id) Log("max dose = %g, min dose = %g", maxv, minv);

			if (hist >= fHMax)
			{
				if (it == main_id) sthread.join(); //make sure the source generating thread ends safely
#pragma omp barrier 
				break;
			}
			cudaErrorCheck(cudaMemcpy(gc[it].d_doseScore, gdose, sizeof(SFloat)*NVoxel, cudaMemcpyHostToDevice)); //reset the dose counter in GPU
			++source_reuse; //count how many time the source has been used
			//try to reuse the source particle to save CPU end time
			if (source_reuse >= Source_Reuse_Times) //need to regenerate initial particles
			{
				if (it == main_id)
				{
					printf("REUSE detected!\n\n\n\n\n\n");
					sthread.join(); //wait for the particles have been copied to GPU
					histAdd = histNew;
				}
#pragma omp barrier //wait until all cards got new particles
				pInit = sp.getAP(it); //update the particle array pointer
#pragma omp barrier //wait until all card threads updated the incident photon array pointer
				source_reuse = 0; //reset the reuse counter
				if (it == main_id) //prepare next batch of incident particles simultaneously
				{
					sthread = std::thread(&getSource, &sp, &histNew);
				}
			}
		}

		//finish in this GPU thread
		gc[it].hist = hist;
		if (it == main_id) Log("GPU processed ------------------------ 100%%,   speed = %d h/s\n", int(hist / rc.stop(true)));
		if (it == main_id) Log("\nWait all GPUs to finish their job...\n");
		delete[] gdose;
		gc[it].destroyStream(); //remember to release the GPU stream resource 
	} //end openMP

	Log("All GPUs have finished their simulation job! Collecting dose...\n\n");

	//merge dose and dose^2 from all GPU cards
	double totHist = gc[0].hist;
	for (int i = 1; i < NGPU; ++i)
	{
		totHist += gc[i].hist;
		//add up dose and uncertainty in different GPU card
		for (int j = 0; j < NVoxel; ++j)
		{
			dose[0][j] += dose[i][j];
			uncertainty[0][j] += uncertainty[i][j];
		}
		delete[] dose[i];
		delete[] uncertainty[i];
	}

	// Get the dose. The dose comes back as Gray per photon times normalization (the constant passed to the 
	// getDose() function below). The purpose of normalization is to converted to real dose 
	// by multiplying with twice the current source activity (because there are two photons per decay)
	// and the total beam on time. For a fresh 15,000 Ci source and 1 minute beam on time, the factor is 
	// 2 * 1.5e4 * 3.7e10 * 60 = 6.66e16. 
	// Here we simply use a normalization that will produce numerical values typically in the range 1-10.

	double norm = 1.0889e15 * SourceHead_BeamOnTime();
	double err50 = zeus._phant->addDose(dose[0], uncertainty[0], NGPU * nBatch, totHist, norm);
	//zeus._phant->addDose(dose[0],totHist, norm);

	delete[] dose[0];
	delete[] uncertainty[0];
	string outname,viewrayFormat;
	cf.getValue("output file name", outname);
	outname += ".dose";

	cf.getValue("ViewRay's format", viewrayFormat);
	if (viewrayFormat.compare("yes")==0) zeus._phant->output(outname.c_str(), 1);
	else zeus._phant->output(outname.c_str(), 0);

	SourceHead_Delete(); //release the source head resource safely
	//ZeusData_unload();

	Log("Time statistics for main GPU:");
	Log("Total copy time =%.2f", copyCounter.getStoredTime());
	Log("Total kernel time = %.1f s", kernelCounter.getStoredTime());
	Log("Total SourceHead time = %.1f s\n\n", sourceCounter.getStoredTime());
	Log("Mixed running time = %.1f minutes, total history number = %.0f", rc.stop(true) / 60.0, totHist);
	Log("The overall simulating speed = %d hist/sec\n\n", int(totHist / rc.stop(true)));
	Log("End log time = %s\n\n", Log.timeNow());
	Log("/##############################################################################/\n\n");
	Log.closeFile();
} //end executeJob