#pragma once
#ifndef _PHANTOM_H_
#define _PHANTOM_H_
#include "../Tools/Tools.h"

#define PhantomBufferSize 256*256*256

#include "gZeus.h"

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
	void output(const char* fname);
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
	int getiabs(int ix, int iy, int iz) { return at(ix, iy, iz); }
	double getHist(){ return Hist; } //this is for logging purpose
	int getNVoxel(){ return NX*NY*NZ; }
	void addDose(double fhist, double aNorm) //return average uncertainty of 50% dose region
	{
		SFloat norm = SFloat(1.602e-16*aNorm / (DX*DY*DZ) / fhist);
		int NVoxel = getNVoxel();
		for (int i = 0; i < NVoxel; ++i)
		{
			//normalize the additional dose and uncertainty
			dose.a(i) *= (norm / ph.a(i));
		}
		Hist = fhist;
	}
	void deposit(int iabsv, ZFloat DE)
	{
#ifdef _DEBUG
		const int LEN = int(1e2);
		static int* dbg_i = new int[LEN];
		static double* dbg_E = new double[LEN];
		static int ind = 0;
#endif

#ifdef USE_OPENMP
#pragma omp critical
#endif
		{
			dose.a(iabsv) += SFloat(DE);
#ifdef _DEBUG
			dbg_i[ind] = iabsv;
			dbg_E[ind] = DE;
			++ind;
			if (ind == LEN)
			{
				FILE* fp = fopen("debug.dat", "w");
				for (int i = 0; i < LEN; ++i)
				{
// 				fwrite(&dbg_i[i], sizeof(int), 1, fp);
// 				fwrite(&dbg_E[i], sizeof(double), 1, fp);
					fprintf(fp, "%10d\t%.15g\n", dbg_i[i], dbg_E[i]);
				}
					
				fclose(fp);
				delete[] dbg_i;
				delete[] dbg_E;
				exitApp("Exported the debug data!");
			}
#endif
		}
	}

	bool lineInPhantom(ParticleR& p)
	{
		p.x -= (xo - 2e-15);
		p.y -= (yo - 2e-15);
		p.z -= (zo - 2e-15);
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

	bool lineIn(ParticleR& p)
	{
		if (lineInPhantom(p)) //forward detect the voxel density to skip the air. Is it necessary?
		{
			const ZFloat step = min(DX, min(DY, DZ));
			const int preNum = 1; //how many step it will forwardly detect
			const ZFloat dmx = step*p.u;
			const ZFloat dmy = step*p.v;
			const ZFloat dmz = step*p.w;
			while (true)
			{
				int ix = int((p.x + preNum*dmx) /DX);
				int iy = int((p.y + preNum*dmy) /DY);
				int iz = int((p.z + preNum*dmz) /DZ);
				if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) return false;//it will leave the phantom without scattering
				if (ph.a(ix, iy, iz) > GlueF(0.04)) break; //stop when it get close to the target
				//advance the particle
				p.x += dmx;
				p.y += dmy;
				p.z += dmz;
			}
		}
		else return false;

		//prepare the voxel index
		p.ivx = int(p.x/DX);
		p.ivy = int(p.y/DY);
		p.ivz = int(p.z/DZ);
		//p.iabsv = at(p.ivx, p.ivy, p.ivz); //this will be recalculated anyway, because the first step is to move
		p.x -= p.ivx*DX;
		p.y -= p.ivy*DY;
		p.z -= p.ivz*DZ;

		return true; //now it's ready to transport
	}

	bool photonFlight(ParticleR & p, ZFloat ds) //return whether particle leaves the phantom
	{
		p.x += ds*p.u + p.ivx*DX;
		p.y += ds*p.v + p.ivy*DY;
		p.z += ds*p.w + p.ivz*DZ;
		if (p.x < 0 || p.x >= LX || p.y < 0 || p.y >= LY || p.z < 0 || p.z >= LZ) return true;
		//calculate the voxel index
		p.ivx = int(p.x / DX);
		p.ivy = int(p.y / DY);
		p.ivz = int(p.z / DZ);

		p.x -= DX*p.ivx;
		p.y -= DY*p.ivy;
		p.z -= DZ*p.ivz;
		p.iabsv = at(p.ivx, p.ivy, p.ivz);
		return false;
	}
	int electronFlight(ParticleR& p, ZFloat Eend);
	SFloat getDensity(ParticleR& p)
	{
		return ph.a(p.iabsv);
	}
	//data
	ArrayMgr<SFloat> ph;
	ArrayMgr<char> mids; //material ids starting from 1
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