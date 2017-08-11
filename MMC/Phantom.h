#pragma once
#ifndef _PHANTOM_H_
#define _PHANTOM_H_
#include "../Tools/Tools.h"

#define at(x,y,z) ((x) + NX*((y)+ NY*(z))) //Column first memory layout, this should be the same as ArrayManager
#define TEs 1.021997804e6

extern SFloat *density;
extern SFloat maxDensity;
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
	double getHist(){ return Hist; } //this is for logging purpose
	int getNVoxel(){ return NX*NY*NZ; }
	void addDose(double fhist, double aNorm) //return average uncertainty of 50% dose region
	{
		SFloat norm = SFloat(1.602e-16*aNorm / (DX*DY*DZ)/fhist);
		int NVoxel = getNVoxel();
		for (int i = 0; i < NVoxel; ++i)
		{
			//normalize the additional dose and uncertainty
			dose.a(i) *= (norm / ph.a(i));
		}
		Hist = fhist;
	}
	
	bool lineIn(Particle& p)
	{
		double px = p.x;
		double py = p.y;
		double pz = p.z;
		//first translate the coordinate system
		px -= xo;
		py -= yo;
		pz -= zo;
		//assuming the incident particle go straight line, judge if it will enter the phantom.
		//if true, give the interaction position
		const double Delta = 1e-5;
		double sx, sy, sz;
		double x = -1.0, y = -1.0, z = -1.0;
		if (px < 0 || px >= LX || py < 0 || py >= LY || pz<0 || pz>LZ) //initial position lays outside the phantom
		{
			double pu = p.u;
			double pv = p.v;
			double pw = p.w;
			if (px < 0)
			{

				//now it's outside the phantom	
				if (pu > 0)
				{
					sx = 0;
					sy = py - px*pv / pu;
					sz = pz - px*pw / pu;
					if (0 <= sy &&sy < LY && 0 <= sz&&sz < LZ)
					{
						x = sx;
						y = sy;
						z = sz;
					}
				}
				else return false;
			}
			else if (px >= LX)
			{
				if (pu < 0)
				{
					sx = LX - Delta;
					sy = py + (sx - px)*pv / pu;
					sz = pz + (sx - px)*pw / pu;
					if (0 <= sy &&sy < LY && 0 <= sz&&sz < LZ)
					{
						x = sx;
						y = sy;
						z = sz;
					}
				}
				else return false;
			}

			if (py < 0)
			{
				//now it's outside the phantom
				if (pv > 0)
				{
					sy = 0;
					sx = px - py*pu / pv;
					sz = pz - py*pw / pv;
					if (0 <= sx &&sx < LX && 0 <= sz&&sz < LZ)
					{
						x = sx;
						y = sy;
						z = sz;
					}
				}
				else return false;
			}
			else if (py >= LY)
			{
				if (pv < 0)
				{
					sy = LY - Delta;
					sx = px + (sy - py)*pu / pv;
					sz = pz + (sy - py)*pw / pv;
					if (0 <= sx &&sx < LX && 0 <= sz&&sz < LZ)
					{
						x = sx;
						y = sy;
						z = sz;
					}
				}
				else return false;
			}

			if (pz < 0)
			{
				//now it's outside the phantom
				if (pw > 0)
				{
					sz = 0;
					sy = py - pz*pv / pw;
					sx = px - pz*pu / pw;
					if (0 <= sy &&sy < LY && 0 <= sx&&sx < LX)
					{
						x = sx;
						y = sy;
						z = sz;
					}
				}
				else return false;
			}
			else if (pz >= LZ)
			{
				if (pw < 0)
				{
					sz = LZ - Delta;
					sy = py + (sz - pz)*pv / pw;
					sx = px + (sz - pz)*pu / pw;
					if (0 <= sy &&sy < LY && 0 <= sx&&sx < LX)
					{
						x = sx;
						y = sy;
						z = sz;
					}
				}
				else return false;
			}

			if (x != -1.0) //got effective intersection point (x,y,z)
			{
				const double step = min(DX, min(DY, DZ));
				const int preNum = 1; //how many step it will forwardly detect
				const double dmx = step*pu;
				const double dmy = step*pv;
				const double dmz = step*pw;
				while (true)
				{

					int ix = int((x + preNum*dmx) / DX);
					int iy = int((y + preNum*dmy) / DY);
					int iz = int((z + preNum*dmz) / DZ);
					if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) return false;//it will leave the phantom without scattering
					if (ph.a(ix, iy, iz) > 0.04) break; //stop when it get close to the target
					x += dmx;
					y += dmy;
					z += dmz;
				}
				px = x; //store the position
				py = y;
				pz = z;
			}
			else return false;
		}

		//prepare the voxel index
		int ix = int(px / DX);
		int iy = int(py / DY);
		int iz = int(pz / DZ);
		px -= ix*DX;
		py -= iy*DY;
		pz -= iz*DZ;

		//write back
		p.x = px;
		p.y = py;
		p.z = pz;
		p.ivx = ix;
		p.ivy = iy;
		p.ivz = iz;
		p.iabsv = at(ix, iy, iz);

		return true; //now it's ready to transport
	}

	bool move(Particle & p, double ds) //return whether particle leaves the phantom
	{
		//move for electron/photon. Note that coordinates are relative to each voxel
		int pType = p.type;
		if (pType == photon)
		{
			//note ds is the water equivalent length
			ds /= MaxDensity; //this is the step in media with max density(woodcock tracking)
			double px = p.x + ds*p.u + p.ivx*DX;
			double py = p.y + ds*p.v + p.ivy*DY;
			double pz = p.z + ds*p.w + p.ivz*DZ;
			if (px < 0 || px >= LX || py < 0 || py >= LY || pz < 0 || pz >= LZ) return true;
			//calculate the voxel index
			int ix = int(px / DX);
			int iy = int(py / DY);
			int iz = int(pz / DZ);

			p.x = px - DX*ix;
			p.y = py - DY*iy;
			p.z = pz - DZ*iz;
			p.ivx = ix;
			p.ivy = iy;
			p.ivz = iz;
			p.iabsv = at(ix, iy, iz);
			return false;
		}

		//it's for photon and positron
		double pE = p.E;
		//fetch data from the slow global memory
		double px = p.x;
		double py = p.y;
		double pz = p.z;
		double pu = p.u;
		double pv = p.v;
		double pw = p.w;
		int ivx = p.ivx;
		int ivy = p.ivy;
		int ivz = p.ivz;
		int iabsv = p.iabsv;

		if (rf > 0) //we need to consider magnetic field effect
		{
			const double deltax = 0.01;
			double voxden = ph[iabsv];
			double invDen = 1 / voxden;
			double q = -1;
			if (pType == positron) q = 1;
			double uwx = -q*Bx;
			double uwy = -q*By;
			double uwz = -q*Bz;
			double Rb = rf*sqrt(pE*(pE + TEs));
			double maxStep = sqrt(2 * Rb*deltax); //max allowed distance to move to ensure accuracy

			while (true)
			{
				double step = ds*invDen;
				bool finalStep = true;
				if (step > maxStep)
				{
					step = maxStep;
					finalStep = false;
				}
				//check if it intersect with the boundary of current voxel
				short idex, dvox; //idex = 1,2,3 means x,y,z direction; dvox = +1, -1 means moving in positive or negative direction
				bool intersect = false;

				if (pv > 0)
				{
					double next = DY - py;
					if (pv*step > next)
					{
						step = next / pv; idex = 2; dvox = 1; intersect = true;
					}
				}
				else if (pv < 0)
				{
					double next = -py;
					if (pv*step < next)
					{
						step = next / pv; idex = 2; dvox = -1; intersect = true;
					}
				}

				if (pw > 0)
				{
					double next = DZ - pz;
					if (pw*step > next)
					{
						step = next / pw; idex = 3; dvox = 1; intersect = true;
					}
				}
				else if (pw < 0)
				{
					double next = -pz;
					if (pw*step < next)
					{
						step = next / pw; idex = 3; dvox = -1; intersect = true;
					}
				}

				if (pu > 0)
				{
					double next = DX - px;
					if (pu*step > next)
					{
						step = next / pu; idex = 1; dvox = 1; intersect = true;
					}
				}
				else if (pu < 0)
				{
					double next = -px;
					if (pu*step < next)
					{
						step = next / pu; idex = 1; dvox = -1; intersect = true;
					}
				}

				if (intersect) finalStep = false;

				if (!finalStep)
				{
					ds -= step*voxden; //record the left distance to go
				}

				//move the electron/positron
				px += pu*step;
				py += pv*step;
				pz += pw*step;

				double vuw = pu*uwx + pv*uwy + pw*uwz;
				double vperpx = pu - vuw * uwx,
					vperpy = pu - vuw * uwy,
					vperpz = pu - vuw * uwz;
				double vxwx = vperpy*uwz - vperpz*uwy,
					vxwy = vperpz*uwx - vperpx*uwz,
					vxwz = vperpx*uwy - vperpy*uwx;
				// The step-length dependent variables f1 & f2
				double f1, f2;
				double arg = step / Rb;
				if (arg < 0.2)
				{
					// arg is small, so use power series expansion of sine and cosine
					double arg2 = arg*arg;
					f1 = -0.5*arg2 + 0.0416666667*arg2*arg2;  // for 0.2, relative error is 2.2e-6
					f2 = arg - 0.16666667*arg*arg2;           // for 0.2, relative error is 1.3e-5, absolute error is 2.6e-6
				}
				else
				{
					f1 = cos(arg) - 1;
					f2 = sin(arg);
				}
				// Direction change 
				double dvx = f1*vperpx - f2*vxwx;  // would simplify to f1*_v.x - f2*_v.y;
				double dvy = f1*vperpy - f2*vxwy;  // would simplify to f1*_v.y + f2*_v.x;
				double dvz = f1*vperpz - f2*vxwz;  // would simplify to 0 (i.e., component along the magnetic field remains constant).

				if (finalStep)
				{
					//update direction and break the loop
					pu += dvx;
					pv += dvy;
					pw += dvz;

					break;
				}

				//not the final step, we may need to check the direction
				if (intersect)
				{
					//
					// We are entering a new voxel. But because we are also changing the direction, we need to verify that the direction has not changed 
					// in a way that we are actually coming back into the voxel from which we came. The condition for coming back into the same voxel 
					// is that v*(v + dv) < 0, where v and dv are the direction and the direction change of the component crossing the voxel boundary
					//
					switch (idex)
					{
					case 3:
						if (pw*pw + pw*dvz >= 0) //enter a new voxel
						{
							ivz += dvox;
							if (dvox > 0) {
								if (ivz >= NZ) return true;
								pz = 0;
							}
							else {
								if (ivz < 0) return true;
								pz = DZ;
							}
						}
						else intersect = false;
						break;
					case 2:
						if (pv*pv + pv*dvy >= 0)
						{
							ivy += dvox;
							if (dvox > 0) {
								if (ivy >= NY) return true;
								py = 0;
							}
							else {
								if (ivy < 0) return true;
								py = DY;
							}
						}
						else intersect = false;
						break;
					default:
						if (pu*pu + pu*dvx >= 0)
						{
							ivx += dvox;
							if (dvox > 0) {
								if (ivx >= NX) return true;
								px = 0;
							}
							else {
								if (ivx < 0) return true;
								px = DX;
							}
						}
						else intersect = false;
					}
				}// end if (intersect)

				//still intersect after the direction check, need to update the voxel density and index
				if (intersect)
				{
					iabsv = at(ivx, ivy, ivz);
					voxden = ph[iabsv];
					invDen = 1 / voxden;
				}

				//update direction and keep going
				pu += dvx;
				pv += dvy;
				pw += dvz;
			}

			//write back direction and position
			p.u = pu;
			p.v = pv;
			p.w = pw;
			p.x = px;
			p.y = py;
			p.z = pz;

			p.ivx = ivx;
			p.ivy = ivy;
			p.ivz = ivz;
			p.iabsv = iabsv;
			return false;
		}
		else //no magnetic field curving
		{
			while (true)
			{
				double voxden = ph[iabsv];
				double step = ds / voxden;
				//check if it intersect with the boundary of current voxel
				short idex, dvox; //idex = 1,2,3 means x,y,z direction; dvox = +1, -1 means moving in positive or negative direction
				bool intersect = false;
				if (pv > 0)
				{
					double next = DY - py;
					if (pv*step > next)
					{
						step = next / pv; idex = 2; dvox = 1; intersect = true;
					}
				}
				else if (pv < 0)
				{
					double next = -py;
					if (pv*step < next)
					{
						step = next / pv; idex = 2; dvox = -1; intersect = true;
					}
				}

				if (pw > 0)
				{
					double next = DZ - pz;
					if (pw*step > next)
					{
						step = next / pw; idex = 3; dvox = 1; intersect = true;
					}
				}
				else if (pw < 0)
				{
					double next = -pz;
					if (pw*step < next)
					{
						step = next / pw; idex = 3; dvox = -1; intersect = true;
					}
				}

				if (pu > 0)
				{
					double next = DX - px;
					if (pu*step > next)
					{
						step = next / pu; idex = 1; dvox = 1; intersect = true;
					}
				}
				else if (pu < 0)
				{
					double next = -px;
					if (pu*step < next)
					{
						step = next / pu; idex = 1; dvox = -1; intersect = true;
					}
				}

				if (!intersect) //end in current voxel
				{
					//update the position, but don't need to update ivx, ivy, ivz, iabsv
					px += pu*step;
					py += pv*step;
					pz += pw*step;

					break;
				}

				//it enters another voxel, so we need to subtract the previous displacement
				ds -= step*voxden;

				//update the voxel index
				if (idex == 3)
				{
					ivz += dvox;
					if (dvox > 0)
					{
						if (ivz >= NZ) return true;
						pz = 0;
					}
					else {
						if (ivz < 0) return true;
						pz = DZ;
					}
				}
				else if (idex == 2)
				{
					ivy += dvox;
					if (dvox > 0)
					{
						if (ivy >= NY) return true;
						py = 0;
					}
					else
					{
						if (ivy < 0) return true;
						py = DY;
					}
				}
				else
				{
					ivx += dvox;
					if (dvox > 0)
					{
						if (ivx >= NX) return true;
						px = 0;
					}
					else
					{
						if (ivx < 0) return true;
						px = DX;
					}
				}

				iabsv = at(ivx, ivy, ivz);
			}//end of moving forward loop

			//write back the position and index. Note that the direction doesn't change.
			p.x = px;
			p.y = py;
			p.z = pz;

			p.ivx = ivx;
			p.ivy = ivy;
			p.ivz = ivz;
			p.iabsv = iabsv;
			return false;
		}
	}

	void deposit(Particle & p, double DE)
	{	
#ifdef USE_OPENMP
#pragma omp critical
#endif
		{
			dose.a(p.iabsv) += SFloat(DE*p.weight);
		}
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