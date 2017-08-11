#pragma once
#ifndef _PHANTOM_H_
#define _PHANTOM_H_
#include "../Tools/Tools.h"

#define PhantomBufferSize 256*256*256

/**** this class depends highly on your specific problem ****/
class Phantom //define the geometry of the phantom
{
	//methods
public:
	Phantom()
	{
		Hist = 0;
		uniform = 0;
		Bs = 0;
		CT_kVp = 110;
		nBatch = 0;
		prescriptionDose = 0;
		treatmentFraction = 1;
		trimDensityThreshold = 0.02;
	}
	void loadPhantom(ConfigFile *cf);
	void output(const char* fname, int FormatType);
	bool readDose(const char* fname);
	void mergeDose(Phantom* rp);
	
	double getMaxDensity() 
	{
		MaxDensity = ph[0];
		uniform = 1;
		int len = getNVoxel();
		for (int i = 0; i < len; ++i)
		{
			if (ph[i] != ph[0]) uniform = 0;
			if (MaxDensity < ph[i]) MaxDensity = ph[i];
		}
		return MaxDensity;
	}
	double getHist(){ return Hist; } //this is for logging purpose
	int getNVoxel(){ return NX*NY*NZ; }
	double addDose(SFloat* d, SFloat* u, int NB, double fH, double aNorm) //return average uncertainty of 50% dose region
	{
		SFloat norm = SFloat(1.602e-16*aNorm / (DX*DY*DZ));
		int NVoxel = getNVoxel();
		int CB = nBatch + NB; //current total batch number
		double CH = Hist + fH; //current history number
		SFloat D = 0, U = 0;
		SFloat maxDose = 0;
		for (int i = 0; i < NVoxel; ++i)
		{
			D = SFloat(Hist*dose[i]);//previous total dose
			if(nBatch) U = (nBatch*uncertainty[i]*uncertainty[i] + 1) * D * D / nBatch; //previous sum dose^2
			else U = 0;

			//normalize the additional dose and uncertainty
			d[i] *= norm;
			u[i] *= norm*norm;
			//if (ph[i]>0)
			{
				SFloat dNorm = 1.0f;// / ph[i];
				d[i] *= dNorm;
				u[i] *= dNorm*dNorm;
				D += d[i]; //current total dose
				U += u[i]; //current sum dose^2
			}
// 			else
// 			{
// 				D = 0;
// 				U = 0;
// 			}
			
			
			dose.a(i) = SFloat(D / CH); //current average dose per history
			if (D) uncertainty.a(i) = SFloat(sqrt(U / (D*D) - 1.0 / CB)); //current standard deviation % based on batches
			else uncertainty.a(i) = -1;

			if (dose.a(i) > maxDose) maxDose = dose.a(i);
		}
		//update batch number and history number
		nBatch = CB;
		Hist = CH;

		//find out the uncertainty of 50% dose region
		int n50 = 0;
		double errN50 = 0;
		for (int i = 0; i < NVoxel; ++i)
		{
			if (dose[i]>0.5*maxDose)
			{
				++n50;
				errN50 += uncertainty[i]*uncertainty[i];
			}
		}
		double ave_err = sqrt(errN50 / n50);
		Log("Total batch number = %d", CB);
		Log("D>50%%*Dmax region occupies %.1f%% all voxels, with average uncertainty = %f",100.0*double(n50)/NVoxel,ave_err);
		return ave_err;
	}
	void addDose(SFloat* d, double fhist, double aNorm) //return average uncertainty of 50% dose region
	{
		SFloat norm = SFloat(1.602e-16*aNorm / (DX*DY*DZ) / fhist);
		int NVoxel = getNVoxel();
		for (int i = 0; i < NVoxel; ++i)
		{
			//normalize the additional dose and uncertainty
			dose.a(i) = d[i]*(norm / ph.a(i));
		}
		Hist = fhist;
	}
	//data
	ArrayMgr<SFloat> ph;
	ArrayMgr<SFloat> dose; //store average dose per history (not per batch)
	ArrayMgr<SFloat> uncertainty; // store standard deviation based on batch
	
	int NX, NY, NZ; //voxel number
	double DX, DY, DZ; // voxel size, unit cm
	double Bs; //magnetic field strength
	double Bx, By, Bz; //unit magnetic field direction
	double rf; //radius factor= 100/C/B

	double Hist; 
	int nBatch;
	int uniform;
	double MaxDensity;
	double xo, yo, zo; // position of the cuboid corner in isocenter coordinate system;
	double LX, LY, LZ; // side length Lx=DX*NX, not need to output
	double COPX, COPY, COPZ;//DICOM coordinates info, keep it just for DoseViewer
	double CT_kVp; //keep it just for DoseViewer
	double prescriptionDose;
	int treatmentFraction;
	double trimDensityThreshold;
};

#endif