#ifndef _GSOURCEHEAD_H_
#define _GSOURCEHEAD_H_
#include "config.h"

class GRNG;
class ParticleR;



class CoPhaseSpace4GPU
{
public:
	void init(void* vps);
	__device__ int sample(GRNG &rng, bool isPrimary, int bin, ZFloat& E, Vector& x, Vector& u);
private:

	ZFloat         _prob1[57600];         //!< Probability for a 1.17 MeV photon as a function of x,y position in a plane just above the MLC
	ZFloat         _pScoarse[225];        //!< Probability for a scattered photon to fall into the first angular table
	ZFloat         _pScoarse1[225];       //!< Probability for a scattered photon to fall into the second angular table

	unsigned short _scatSpectrum[225][2 * 254];    //!< Scatter particle energy spectra in _Nx1*_Ny1 tiles prepared for convenient sampling
	unsigned short _thetaPhiPflu[225][2 * 65536];  //!< Primary particle fluence as a function of sin(theta),phi in _Nx1*_Ny1 tiles	
	unsigned short _thetaPhiSflu[225][2 * 65536];  //!< Scatter particle fluence as a function of sin(theta),phi in _Nx1*_Ny1 tiles
	unsigned short _thetaPhiSflu1[225][2 * 16384]; //!< Scatter particle fluence as a function of sin(theta),phi in _Nx1*_Ny1 tiles
	unsigned short _thetaPhiSfluC[225][2 * 4096];  //!< Used for sampling sin(theta),phi of primary particles
};

class SEGMENT4GPU
{
public:
	void init(void* beam, void* sg);
	
	__device__ int sample(GRNG &rng, int iseg, ParticleR pars[]);

	// to determine sampling weight, must be initialized manually
	ZFloat wfirst;
	int wsecond;

	ZFloat nSample;
	ZFloat pPrim;
	unsigned short primAT[2 * 240 * 240];
	unsigned short scatAT[2 * 240 * 240];
	char primMask[240 * 240];
	char scatMask[240 * 240];

	ZFloat gantryAngle;
	ZFloat cosPhi, sinPhi;
};

__global__ void gSample(); // expose interface for C++ host code to call

__device__ int sample(GRNG &rng, ParticleR pars[]); // expose interface for other CUDA function to call

void SourceHead4GPU_Init(ConfigFile* cf, void* gc);// expose interface for initialization

double SourceHead4GPU_BeamOnTime();

void SourceHead4GPU_GetIsoCenter(double* px, double* py, double* pz);

void SourceHeadGPU_GetPrescription(double* pDose, int* fraction);

void getHostData(ConfigFile* cf, CoPhaseSpace4GPU*&hps, SEGMENT4GPU*&hseg, int& hnSeg, int& nMaxPhoton, vector<vector<pair<double, double> > >& aListofSegments); //expose interface to SourceHead.cpp
#endif