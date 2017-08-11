#pragma once
#ifndef _PENELOPE_H_
#define _PENELOPE_H_

#include "../Tools/Tools.h"
#include "Phantom.h"
#include "../SourceHead/SourceHead.h"
#include "octree_gpu.h"

#define CallStyle //_stdcall //if using CVF, please define _stdcall


#ifdef WIN32
#include <io.h>
#endif

//resolve the compiler's warning
#ifndef __CUDACC__
#undef __device__ 
#undef __host__
#undef __forceinline__
#undef __global__
#undef __constant__
#endif
#include "cuda_runtime.h"
//#define USE_MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: macro definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
//#define STEP_DEBUG_STATUS //for single thread debug purpose
#define STEP_NUM 170000
//#define ACCU_DEBUG
#define INCIDENT_NUM 1000000

// 0 means only compiling kernel-by-kernel method, 1 means only compiling one-thread-one-history method, -1 means both
#define EXECUTE_TYPE -1
#define TRANSPORT_OCTREE 0 // whether to include the octree transportation ( 1 yes, 0 no )
#define MULTI_MATERIAL 0 // whether to transport through multiple materials ( 1 yes, 0 no )
/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end: macro definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/************************Start: Constant Definitions **************************/
#define PI 3.141592653589793238463
#define Es 5.10998902e5
#define TEs 1.021997804e6
#define alpha 7.297352536e-3
#define TRUNC 1.01538698  //used in Gaussian sampling

#define NEGrid 200
#define NE 96
#define NA 606
#define NP 128
#define NO 64
#define NOCO 64
#define NRX 15000
#define NTP 8000
#define NBW 32

#define NMatMax 5 //the maximum material number
/************************End: Constant Definitions **************************/

__device__ void exitKernel(const char inf[]);

/*************************Start: type/class definition *********************************/
enum KnockMode { SOFT, HARD };
enum KonckEvent
{
	GRA_EVT=0,
	GCO_EVT=1,
	GPH_EVT=2,
	GPP_EVT=3,

	ESOFT_EVT=4,
	EEL_EVT=5,
	EIN_EVT=6,
	EBR_EVT=7,
	ESI_EVT=8,

	PSOFT_EVT=9,
	PEL_EVT=10,
	PIN_EVT=11,
	PSI_EVT=12,
	PAN_EVT=13,

	DELTA_EVT=14,
	EXIT_EVT=15, // the particle leave the interaction region
	INIT_EVT=16,
	NEWMAT_EVT=17 //moved into a new material
};

struct CommonConst
{
	double Etab_h[NEGrid];
};

class Material
{
	//interface:
public:
	void load(int mati);
	__device__ void start(); //called before initializing a new primary particle simulation
	__device__ void jump(size_t it); //It includes moving as well.
	__device__ double knock(size_t it);

	__device__ double SimpleJump(size_t it);
	__device__ bool SimpleKnock(size_t it);
	__device__ bool SimpleCompton(size_t it);

	__device__ void GJump(size_t);
	__device__ void EJump(size_t);
	__device__ void PJump(size_t);
	__device__ double ESoftKnock(size_t it);
	__device__ double PSoftKnock(size_t it);

	//photon knock
	__device__ double GRA(size_t);
	__device__ double GCO(size_t);
	__device__ double GPH(size_t);
	__device__ double GPP(size_t);

	__device__ double EEL(size_t);
	__device__ double EIN(size_t);
	__device__ double EBR(size_t); //suitable for both electron and positron
	__device__ double ESI(size_t);

	//positron knock
	__device__ double PEL(size_t);
	__device__ double PIN(size_t);
	__device__ double PSI(size_t);
	__device__ double PAN(size_t);
	__device__ double PAN(size_t, Particle& p);
private:

	//other process
	__device__ double EELt(double TA, double TB, double pc, size_t it);//used for both electron and positron
	__device__ double EELd(double pc, size_t it);
	__device__ double PELd(double pc, size_t it);
	
	__device__ void PANR();
	__device__ double relax(int IZZ, int ISH);
	//data:
public:
	int M; //material index number, starting from 1
	double massDensity;
	double totalSigmaG(int IE);

	//simulation parameters
	double Eabs[3]; //absorption cutoff energy
	double C1; //average angular deflection, 
	double C2; //maximum average fractional energy loss
	double Wcc; //cutoff energy loss for hard inelastic collision, eV
	double Wcr; //cutoff energy loss for hard bremsstrahlung emission, eV
	//
	double Ecutr;
	//composition data
	double weight[30], Zt, At, numDensity;
	int elemNum;
	int IZ[30];

	//photo simulation data 
	double SGRA[NEGrid], SGCO[NEGrid], SGPH[NEGrid], SGPP[NEGrid];
	//CGRA
	double X2COH[241], PDCOH[241];
	// Compton Scattering
	double FCO[NOCO], UICO[NOCO], FJ0[NOCO];
	int KZCO[NOCO], KSCO[NOCO], NOSCCO;

	//CGPP00, pair production, 
	double ZEQPP, Fpp, Bpp;// they're F0(MAXMAT,2), BCB(MAXMAT)

	// electron simulation table
	double SEHEL[NEGrid], SEHIN[NEGrid], SEISI[NEGrid], SEHBR[NEGrid],
		SETOT[NEGrid], CSTPE[NEGrid], RSTPE[NEGrid],
		DEL[NEGrid], W1E[NEGrid], W2E[NEGrid], T1E[NEGrid], T2E[NEGrid],
		RNDCE[NEGrid], AE[NEGrid], BE[NEGrid];
	double SPHEL[NEGrid], SPHIN[NEGrid], SPISI[NEGrid], SPHBR[NEGrid], SPAN[NEGrid],
		SPTOT[NEGrid], CSTPP[NEGrid], RSTPP[NEGrid],
		W1P[NEGrid], W2P[NEGrid], T1P[NEGrid], T2P[NEGrid],
		RNDCP[NEGrid], AP[NEGrid], BP[NEGrid];

	//CELSEP, electron/position elastic collision
	double EELMAX, PELMAX;
	//CEELDB
	double XSE[NP][NEGrid], PSE[NP][NEGrid], ASE[NP][NEGrid], BSE[NP][NEGrid];
	int ITLE[NP][NEGrid], ITUE[NP][NEGrid];
	//CPELDB
	double XSP[NP][NEGrid], PSP[NP][NEGrid], ASP[NP][NEGrid], BSP[NP][NEGrid];
	int ITLP[NP][NEGrid], ITUP[NP][NEGrid];

	// E/P inelastic collision
	double meanExcitE, OP2, F[NO], Ui[NO], WRI[NO];
	int KZ[NO], KS[NO], NOSC;
	//CEINAC
	double EINAC[NEGrid][NO], PINAC[NEGrid][NO];

	//CEBR, Bremsstrahlung emission
	double PBcut[NEGrid], WBcut[NEGrid], PDFB[NEGrid][NBW], PACB[NEGrid][NBW];
	//double Q1B[NEGrid][21], Q2B[NEGrid][21]; //need for special init, derived from BP1 and BP2

	//auxiliary arrays
	double BP1[6][21][4], BP2[6][21][4];
};

class ParticleStack
{
	//method
public:
	__device__ void push(Particle &par);
	__device__ __forceinline__ bool empty()
	{
		if (cur > 0) return false;
		else return true;
	}
	__device__ __forceinline__ void pop()
	{
		if (cur <= 0) exitKernel("no element to pop in the stack!");
		--cur;
	}
	__device__ __forceinline__ Particle& top(){ return ss[cur - 1]; }
	__device__ void init(Particle* pp){ cur = 0; ss = pp; }
	//data
private:
	Particle* ss;
	unsigned short cur; //point to next free space, 0 at initial time
};

class GRNG
{
public:
	__device__ void init(int aSeed, int threadID)
	{
		ISEED1 = aSeed + 2 * threadID;
		ISEED2 = ISEED1 + 1;
	}
	__device__ __forceinline__ double operator() ()
	{
		const double USCALE = 1.0 / 2.147483563e9;
		int I1, I2, IZ;

		I1 = ISEED1 / 53668;
		ISEED1 = 40014 * (ISEED1 - I1 * 53668) - I1 * 12211;
		if (ISEED1 < 0) ISEED1 += 2147483563;

		I2 = ISEED2 / 52774;
		ISEED2 = 40692 * (ISEED2 - I2 * 52774) - I2 * 3791;
		if (ISEED2 < 0) ISEED2 += 2147483399;

		IZ = ISEED1 - ISEED2;
		if (IZ < 1) IZ += 2147483562;
		return IZ*USCALE;
	}

private:
	int ISEED1;
	int ISEED2;
};

class ThreadVars //parameters that vary with different thread
{
public:
	GRNG rng;
	Particle p;
	ParticleStack pStack;
	int IE; //left index of interval
	double XIE; //portion in that interval
	double st, W1, W2, T1, T2, DS1, DST;
	bool KDelta;
	KnockMode  Mode; //0 artificial soft event, 1 hard event
	KonckEvent evt;
	int InitCur;
#if	MULTI_MATERIAL == 1
	char cMatID; // current material id
#endif
};

struct GPUConfig // should not copy this variable
{
	GPUConfig()
	{
		d_InitPars = NULL;
		d_doseScore = NULL;
		d_gOctreeRoot = NULL;
	}

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

		if (d_gOctreeRoot) deleteGPUOctree(d_gOctreeRoot);
	}

	int id; //GPU index

	//some pointers in the GPU that will be used in simulation
	Particle* d_InitPars;
	SFloat* d_doseScore;

	vector<void*> gpu2free;

	//pointer for the Octree
	goctree_node* d_gOctreeRoot;

	//these are most likely the same
	int NBlock, BlockSize, NBatch;
	int SourceReuseTimes;

	double hist; //accumulate how many histories has be simulated in this GPU

	int getNOneBatch(){ return NBlock*BlockSize*NBatch; }
};

class SourcePool //a wrapper of SourceHead
{
public:
	SourcePool(int NT, int NOneBatch, int NSampleStack=400)
	{
		_NGenerateThread = NT;
		if (NOneBatch%NT != 0) exitApp("Unequal thread load in SourcePool. Please make sure NOneBatch/NThread == 0");
		_PoolCapacity = NOneBatch / NT;

		_np = new int[NT]; //how many particle was actually held in each thread
		_hist = new int[NT]; //how many histories have been sampled for each thread
		_rng = new PRNG[NT]; // create rng for each thread

		_PPool = new Particle*[NT]; //create pool pointer
		for (int i = 0; i < NT; ++i)
		{
			_PPool[i] = new Particle[_PoolCapacity + NSampleStack]; //resize memory for each thread
			_np[i] = 0;
			_hist[i] = 0;
			_rng[i].init(1234 + i); //init the rng for each thread
		}
	}
	~SourcePool()
	{
		//delete particle pool
		for (int i = 0; i < _NGenerateThread; ++i)
		{
			delete[] _PPool[i];
		}
		delete[] _PPool;
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
			if (rest >= 0)
			{
				for (int j = 0; j < rest; ++j) //need move the buffered particles to the left side
					_PPool[it][j] = _PPool[it][j + _PoolCapacity];
				_np[it] = rest;
			}
			else _np[it] = 0;

			while (_np[it] < _PoolCapacity)
			{
				_np[it] += SourceHead_Sample(&_rng[it], _PPool[it] + _np[it]);
				++_hist[it];
			}
		}
		int sum = 0;
		for (int i = 0; i < _NGenerateThread; ++i) sum += _hist[i];
#ifdef STAT_HEAD
		statHead();
#endif
		return sum; //return how many histories have generated
	}
	void copy2GPU(Particle* d_InitPars)
	{
		for (int i = 0; i < _NGenerateThread; ++i)
		{
			cudaMemcpy(d_InitPars + i*_PoolCapacity, _PPool[i], _PoolCapacity * sizeof(Particle), cudaMemcpyHostToDevice);
		}
	}
private:
#ifdef STAT_HEAD
	void statHead();
	void statOutput();
#endif
	int _NGenerateThread;
	int _PoolCapacity; //how many particles it can hold in each line
	Particle** _PPool;
	int *_np;   //how many particle was actually held in each line
	int *_hist; //how many histories have been sampled for each line
	PRNG *_rng; //rng array
};

class PENELOPE //class wrapper of the simulation kernel
{
public:
	vector<GPUConfig> gc;
	int getGPUConfig(ConfigFile* gcf);
	void init(ConfigFile *cf);

	PENELOPE() //inner default value
	{
		simu_SEC = 0;
		simu_POS = 0;
		simu_GRA = 1;
		simu_GCO = 1;
		simu_GPH = 1;
		simu_GPP = 1;
		_UseOctree = false;
	}
	~PENELOPE()
	{
		if(_phant) delete _phant;
		if(_pMat) delete [] _pMat;
		int ngpu = (int)gc.size();
		for (int i = 0; i < ngpu; ++i) gc[i].freeGPU();
	}

	//simulation switches
	int simu_GRA;
	int simu_GCO;
	int simu_GPH;
	int simu_GPP;

	int simu_EEL;
	int simu_EIN;
	int simu_EBR;
	int simu_ESI;

	int simu_SEC;
	int simu_POS;

	int simu_PEL;
	int simu_PIN;
	int simu_PSI;
	int simu_PAN;

	Phantom* _phant; //pointer to the phantom in CPU
	Material* _pMat; //pointer to the material array in CPU
private:

	void loadCommonData();

	void phantom2GPU();

	void initMaterial(ConfigFile *cf);

	void initOctree(ConfigFile *cf);

	void initGPU();

	ConfigFile* _cf;
	int _NMAT; // how many materials specified in PENELOPE config block, doesn't have to equal the material number in phantom.NMAT
	double _DSMax;
	int _NStackDepth;
	vector<int> _MatList; //PENELOPE material id list
	bool _UseOctree;
};

/*************************End: type/class definition *********************************/

void executeJob(const char* configFileName, double processWeight, MPS& configMacro); //execute one job according to the config file

#endif