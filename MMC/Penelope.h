#pragma once
#ifndef _PENELOPE_H_
#define _PENELOPE_H_

#define NSECDEPTH 40
#include "../Tools/Tools.h"
#include "../SourceHead/SourceHead.h"
#include "Phantom.h"
/************************Start: Constant Definitions **************************/

const double RE = 2.817940285e-13; //classical electron radius
const double PIRE2 = PI*RE*RE;
const double TRUNC = 1.01538698;//used in Gaussian sampling
const double alpha = 1.0 / 137.03599976; //fine-structure constant alpha = 7.29735257e-3
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

const double WB[] = //the discrete value of the reduced energy of photon
{
	1.0e-12, 0.025e0, 0.05e0, 0.075e0, 0.1e0, 0.15e0, 0.2e0, 0.25e0,
	0.3e0, 0.35e0, 0.4e0, 0.45e0, 0.5e0, 0.55e0, 0.6e0, 0.65e0, 0.7e0,
	0.75e0, 0.8e0, 0.85e0, 0.9e0, 0.925e0, 0.95e0, 0.97e0, 0.99e0,
	0.995e0, 0.999e0, 0.9995e0, 0.9999e0, 0.99995e0, 0.99999e0, 1.0e0
};

const int NBW = 32;
/************************End: Constant Definitions **************************/
enum KnockMode { SOFT, HARD };

struct ProcessData //data shared by all processes. This wrap makes it easier transport through MPI
{
	double ELtab[NEGrid]; //log(E) grid table
	double Etab[NEGrid];
	double ELow, ELowL, EUp; //energy up and low bound
	double DLE; //delta log(E)

	double FLCOH;
	double EPH[8000], XPH[8000][10];
	int IPHF[99], IPHL[99], NPHS[99];
	double ET[15000], FR[15000];
	int IAL[15000], IS1[15000], IS2[15000], IFIRST[99][9], ILAST[99][9];
	double XESI[6000][9], XPSI[6000][9];
	int IESIF[99], NSESI[99], IPSIF[99], NSPSI[99];
	double EB[99][30];
	double BET[6];
};


class Material
{
	//interface:
public:
	void load(int mati);
	void start(); //called before initializing a new primary particle simulation
	double jump(double DSMax, int it);
	double knock(int it);

private:
	void EIMFP(KnockMode type);
	void PIMFP(KnockMode type);

	//photon knock
	double GRA();
	double GCO(int it);
	double GPH();
	double GPP();
	//electron knock
	double EEL(double TA, double TB, double pc);//used for both electron and positron
	double EELd(double pc, int it);
	double EIN();
	double EBR(); //suitable for both electron and positron
	double ESI();
	//positron knock
	double PELd(double pc);
	double PIN();
	double PSI();
	double PAN();

	//other process
	void PANR();
	void relax(int IZZ, int ISH);
	//data:
public:
	int M; //material index number, starting from 1
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
	void push(Particle &par)
	{
		if (cur == NSECDEPTH) exitApp("exceed the max depth of the stack!");
		ss[cur] = par; //copy and store the particle
		++cur;
	}
	void pop()
	{
		if (cur <= 0) exitApp("no element to pop in the stack!");
		--cur;
	}
	Particle& top(){ return ss[cur - 1]; }
	//data
private:
	Particle ss[NSECDEPTH];
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
	Particle p;
	PENRNG rng;
	SecondaryStack secStack;
	int IE; //left index of interval
	double XIE; //portion in that interval
	//inverse mean free path related parameters
	double p1, p2, p3, p4, p5, p6, p7, st; //store the total cross section for different interaction events
	double W1, W2, T1, T2, DS1, DST;
	KnockMode  Mode; //0 artificial soft event, 1 hard event
	bool KSoftI; //if soft stopping interaction is active
	bool KSoftE; //if soft stopping scattering is active
	bool KDelta; //false, a hard interaction follows; true, a delta interaction follows
};



void initPenelope(ConfigFile *cf, Material*& mat, int NThread);

void executeJob(const char* configFileName,  MPS& configMacro); //execute one job according to the config file

#endif