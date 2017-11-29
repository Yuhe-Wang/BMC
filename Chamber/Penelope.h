#pragma once
#ifndef _PENELOPE_H_
#define _PENELOPE_H_

#include "../Tools/Tools.h"
using MonteCarlo::Vector;

/************************Start: Constant Definitions **************************/
const int NEGrid = 200;
const int LNEGrid = 200;
const int NMatMax = 10; //the maximum material number
const int NE = 96;
const int NA = 606;
const int NP = 128;
const int NO = 64;
const int NOCO = 64;
const int NRX = 15000;
const int NTP = 8000;
const int NBW = 32;
const int NSECDEPTH = 400;//the depth of the secondary particle stack
/************************End: Constant Definitions **************************/

class Material
{
	//interface:
public:
	void load(int mati);
	void start(int it); //called before initializing a new primary particle simulation
	double jump(int it);
	double knock(int it);

private:
	//photon knock
	double GRA(int it);
	double GCO(int it);
	double GPH(int it);
	double GPP(int it);
	//electron knock
	double EEL(double TA, double TB, double pc, int it);//used for both electron and positron
	double EELd(double pc, int it);
	double EIN(int it);
	double EBR(int it); //suitable for both electron and positron
	double ESI(int it);
	//positron knock
	double PELd(double pc, int it);
	double PIN(int it);
	double PSI(int it);
	double PAN(int it);

	//other process
	void PANR(int it);
	double relax(int IZZ, int ISH, int it);
	//data:
public:
	int M; //material index number, starting from 0
	double invXCSE; // short for inverse photon cross-section enhancement
private:
	//simulation parameters
	double Eabs[3]; //absorption cutoff energy
	double C1; //average angular deflection, 
	double C2; //maximum average fractional energy loss
	double Wcc; //cutoff energy loss for hard inelastic collision, eV
	double Wcr; //cutoff energy loss for hard bremsstrahlung emission, eV
	//
	double Ecutr;
	//composition data
	double weight[30], Zt, At, massDensity, numDensity;
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

class PENParticle //basic particle information
{
public:
	double E, weight;
	Vector x, v;
	char RegionID;
	ParticleType type;
};

class SecondaryStack
{
	//method
public:
	SecondaryStack() :cur(0){};
	bool empty()
	{
		if (cur > 0) return false;
		else return true;
	}
	void push(double E, double weight, char RegionID, Vector& x, Vector& v, ParticleType type);
	void pop(double& E, double& weight, char& RegionID, Vector& x, Vector& v, ParticleType& type);
	//data
private:
	PENParticle ss[NSECDEPTH];
	int cur; //point to next free space, 0 at initial time
};

class PENRNG
{
public:
	void init(int aSeed, int threadID)
	{
		ISEED1 = aSeed + 2 * threadID;
		ISEED2 = ISEED1 + 1;
	}
	inline double operator () ()
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

class ThreadParameters //parameters that vary with different thread
{
	//only particle p and secStack should be accessed by the public
public:
	//current particle properties to track
	double E, weight;
	char RegionID; //current media id
	Vector x, v;
	ParticleType type;

	PENRNG rng;
	SecondaryStack secStack;

	//auxiliary variables during one particle's simulation
	int IE; //left index of interval
	double XIE; //portion in that interval
	//inverse mean free path related parameters
	double p1, p2, p3, p4, p5, p6, p7, st; //store the total cross section for different interaction events
	double W1, W2, T1, T2, DS1, DST;
	bool Mode; //false artificial soft event, true hard event
	bool KSoftI; //if soft stopping interaction is active
	bool KSoftE; //if soft stopping scattering is active
	bool KDelta; //false, a hard interaction follows; true, a delta interaction follows
};

class PENELOPE
{
public:
	PENELOPE() :mat(NULL){};
	void init(ConfigFile *cf, int Thread = 1, double Emax = 0, vector<int> ids = vector<int>());
	~PENELOPE();
	
	Material* getMat(){ return mat; }
	ThreadParameters* getTP();

private:
	Material* mat;
};

// C type function interface
// variable exposure
#endif