#include "gZeus.h"

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: variables in device memory <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
//extern __shared__ int sharedMem[];
__constant__ int NBatch;
__constant__ int NStackDepth;
__constant__ int FixedSplit = 1;
__constant__ int ForwardDetect = 1; //whether do forward advancing movement when guidling in the photon
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

__device__ ZFloat getDensity(int iabsv)
{
	if (uniform) return MaxDensity;
	else return ph[iabsv]; 
	//else return getPhant(iabsv);
}

__device__ __forceinline__ ZFloat getDensity(ParticleR& p)
{
	return getDensity(p.iabsv);
}

__device__ __forceinline__ void deposit(int iabsv, ZFloat DE)
{
	atomicAdd(doseScore + iabsv, DE);
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
	//if (!lineInPhantom(p)) return false;

	if(ForwardDetect) //forward detect the voxel density to skip the air. Is it necessary?
	{
		const ZFloat step = min(DX, min(DY, DZ));
		const int preNum = 1; //how many step it will forwardly detect
		const ZFloat dmx = step*p.u;
		const ZFloat dmy = step*p.v;
		const ZFloat dmz = step*p.w;
		while (true)
		{
			int ix = int((p.x + preNum*dmx) *InvDX);
			int iy = int((p.y + preNum*dmy) *InvDY);
			int iz = int((p.z + preNum*dmz) *InvDZ);
			if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) return false;//it will leave the phantom without scattering
			if (getDensity(at(ix, iy, iz)) > GlueF(0.04)) break; //stop when it get close to the target
			//advance the particle
			p.x += dmx;
			p.y += dmy;
			p.z += dmz;
		}
	}

	//prepare the voxel index
	p.ivx = int(p.x*InvDX);
	p.ivy = int(p.y*InvDY);
	p.ivz = int(p.z*InvDZ);
	//p.iabsv = at(p.ivx, p.ivy, p.ivz); //this will be recalculated anyway
	p.x -= p.ivx*DX;
	p.y -= p.ivy*DY;
	p.z -= p.ivz*DZ;
	
	return true; //now it's ready to transport
}

__device__ ZFloat tperp(ParticleR& p)
{
	ZFloat tx = DX - p.x;
	if (p.x < tx) tx = p.x;
	ZFloat ty = DY - p.y;
	if (p.y < ty) ty = p.y;
	ZFloat tz = DZ - p.z;
	if (p.z < tz) tz = p.z;
	return tx < ty && tx < tz ? tx : ty < tz ? ty : tz; //min of (tx, ty, tz)
}

__device__ bool photonFlight(ParticleR & p, ZFloat ds) //return whether particle leaves the phantom
{
	p.x += ds*p.u + p.ivx*DX;
	p.y += ds*p.v + p.ivy*DY;
	p.z += ds*p.w + p.ivz*DZ;
	if (p.x < 0 || p.x >= LX || p.y < 0 || p.y >= LY || p.z < 0 || p.z >= LZ) return true;
	//calculate the voxel index
	p.ivx = int(p.x*InvDX);
	p.ivy = int(p.y*InvDY);
	p.ivz = int(p.z*InvDZ);

	p.x -= DX*p.ivx;
	p.y -= DY*p.ivy;
	p.z -= DZ*p.ivz;

	p.iabsv = at(p.ivx, p.ivy, p.ivz);
	return false;
}

// __device__ __forceinline__ bool intersect(ParticleR &p, ZFloat &step, short& idex, short& dvox)
// {
// 	bool b_intersect = false;
// 	if (p.v > 0)
// 	{
// 		ZFloat next = DY - p.y;
// 		if (p.v*step > next)
// 		{
// 			step = next / p.v;
// 			idex = 2;
// 			dvox = 1;
// 			b_intersect = true;
// 		}
// 	}
// 	else if (p.v < 0)
// 	{
// 		ZFloat next = -p.y;
// 		if (p.v*step < next)
// 		{
// 			step = next / p.v;
// 			idex = 2;
// 			dvox = -1;
// 			b_intersect = true;
// 		}
// 	}
// 
// 	if (p.w > 0)
// 	{
// 		ZFloat next = DZ - p.z;
// 		if (p.w*step > next)
// 		{
// 			step = next / p.w;
// 			idex = 3;
// 			dvox = 1;
// 			b_intersect = true;
// 		}
// 	}
// 	else if (p.w < 0)
// 	{
// 		ZFloat next = -p.z;
// 		if (p.w*step < next)
// 		{
// 			step = next / p.w;
// 			idex = 3;
// 			dvox = -1;
// 			b_intersect = true;
// 		}
// 	}
// 
// 	if (p.u > 0)
// 	{
// 		ZFloat next = DX - p.x;
// 		if (p.u*step > next)
// 		{
// 			step = next / p.u;
// 			idex = 1;
// 			dvox = 1;
// 			b_intersect = true;
// 		}
// 	}
// 	else if (p.u < 0)
// 	{
// 		ZFloat next = -p.x;
// 		if (p.u*step < next)
// 		{
// 			step = next / p.u;
// 			idex = 1;
// 			dvox = -1;
// 			b_intersect = true;
// 		}
// 	}
// 
// 	return b_intersect;
// }
// 
// __device__ __forceinline__ bool intersect(ParticleR &p, ZFloat &step)
// {
// 	bool b_intersect = false;
// 	if (p.v > 0)
// 	{
// 		ZFloat next = DY - p.y;
// 		if (p.v*step > next)
// 		{
// 			step = next / p.v;
// 			b_intersect = true;
// 		}
// 	}
// 	else if (p.v < 0)
// 	{
// 		ZFloat next = -p.y;
// 		if (p.v*step < next)
// 		{
// 			step = next / p.v;
// 			b_intersect = true;
// 		}
// 	}
// 
// 	if (p.w > 0)
// 	{
// 		ZFloat next = DZ - p.z;
// 		if (p.w*step > next)
// 		{
// 			step = next / p.w;
// 			b_intersect = true;
// 		}
// 	}
// 	else if (p.w < 0)
// 	{
// 		ZFloat next = -p.z;
// 		if (p.w*step < next)
// 		{
// 			step = next / p.w;
// 			b_intersect = true;
// 		}
// 	}
// 
// 	if (p.u > 0)
// 	{
// 		ZFloat next = DX - p.x;
// 		if (p.u*step > next)
// 		{
// 			step = next / p.u;
// 			b_intersect = true;
// 		}
// 	}
// 	else if (p.u < 0)
// 	{
// 		ZFloat next = -p.x;
// 		if (p.u*step < next)
// 		{
// 			step = next / p.u;
// 			b_intersect = true;
// 		}
// 	}
// 
// 	return b_intersect;
// }

__device__ __forceinline__ bool chvox(ParticleR &p, ZFloat &step, short& idex, short& dvox)
{
	switch (idex)
	{
	case 3:
		if (dvox > 0)
		{
			++p.ivz;
			if (p.ivz >= NZ) return false;
			p.z = 0;
		}
		else
		{
			--p.ivz;
			if (p.ivz < 0) return false;
			p.z = DZ;
		}
		p.x += p.u*step;
		p.y += p.v*step;
		break;
	case 2:
		if (dvox > 0)
		{
			++p.ivy;
			if (p.ivy >= NY) return false;
			p.y = 0;
		}
		else
		{
			--p.ivy;
			if (p.ivy < 0) return false;
			p.y = DY;
		}
		p.x += p.u*step;
		p.z += p.w*step;
		break;
	default:
		if (dvox > 0)
		{
			++p.ivx;
			if (p.ivx >= NX) return false;
			p.x = 0;
		}
		else
		{
			--p.ivx;
			if (p.ivx < 0) return false;
			p.x = DX;
		}
		p.y += p.v*step;
		p.z += p.w*step;
	}
	p.iabsv = at(p.ivx, p.ivy, p.ivz);
	return true;
}

__device__ int electronFreeFlight(ParticleR& p, ZFloat Eend)
{
	ZFloat voxden = getDensity(p);
	ZFloat range = (p.E - ERangeCut)*WaterRangeCS(p.E);
	ZFloat finalRange = Eend > ERangeCut ? (Eend - ERangeCut)*WaterRangeCS(Eend) : 0;
	
	while (true)
	{
		ZFloat step = (range - finalRange) / voxden;
		//check if it intersect with the boundary of current voxel
		short idex, dvox; //idex = 1,2,3 means x,y,z direction; dvox = +1, -1 means moving in positive or negative direction
		bool intersect = false;

		if (p.v > 0)
		{
			ZFloat next = DY - p.y;
			if (p.v*step > next)
			{
				step = next / p.v; idex = 2; dvox = 1; intersect = true;
			}
		}
		else if (p.v < 0)
		{
			ZFloat next = -p.y;
			if (p.v*step < next)
			{
				step = next / p.v; idex = 2; dvox = -1; intersect = true;
			}
		}

		if (p.w > 0)
		{
			ZFloat next = DZ - p.z;
			if (p.w*step > next)
			{
				step = next / p.w; idex = 3; dvox = 1; intersect = true;
			}
		}
		else if (p.w < 0)
		{
			ZFloat next = -p.z;
			if (p.w*step < next)
			{
				step = next / p.w; idex = 3; dvox = -1; intersect = true;
			}
		}

		if (p.u > 0)
		{
			ZFloat next = DX - p.x;
			if (p.u*step > next)
			{
				step = next / p.u; idex = 1; dvox = 1; intersect = true;
			}
		}
		else if (p.u < 0)
		{
			ZFloat next = -p.x;
			if (p.u*step < next)
			{
				step = next / p.u; idex = 1; dvox = -1; intersect = true;
			}
		}

		ZFloat newEnergy = Eend;
		if (intersect)
		{
			range -= step*voxden;
			newEnergy = WaterInverseRangeCS(range);
			if (newEnergy < ERangeCut) newEnergy = 0;
		}
		deposit(p.iabsv, (p.E - newEnergy)*p.weight);
		//up the energy to the new one
		p.E = newEnergy;
		if (p.E < ERangeCut) return 1; // if this is the final step to local absorption, we don't need to update the position and the direction

		//move the electron
		p.x += p.u*step;
		p.y += p.v*step;
		p.z += p.w*step;

		if (!intersect) break;

		if (2 == idex) //y direction
		{
			p.ivy += dvox;
			if (dvox > 0)
			{
				if (p.ivy >= NY) return 1;
				p.y = 0;
			}
			else
			{
				if (p.ivy < 0) return 1;
				p.y = DY;
			}
		}
		else if (3 == idex) //z direction
		{
			p.ivz += dvox;
			if (dvox > 0)
			{
				if (p.ivz >= NZ) return 1;
				p.z = 0;
			}
			else
			{
				if (p.ivz < 0) return 1;
				p.z = DZ;
			}
		}
		else //x direction
		{
			p.ivx += dvox;
			if (dvox > 0)
			{
				if (p.ivx >= NX) return 1;
				p.x = 0;
			}
			else
			{
				if (p.ivx < 0) return 1;
				p.x = DX;
			}
		}

		p.iabsv = at(p.ivx, p.ivy, p.ivz);
		voxden = getDensity(p);
	}
	return 0;
}

__device__ int electronFlight(ParticleR& p, ZFloat Eend)
{
	if (rf == 0) return electronFreeFlight(p, Eend);

	//move for electron/photon. Note that coordinates are relative to each voxel
	const ZFloat deltax = GlueF(0.01);
	ZFloat voxden = getDensity(p);
	ZFloat range = (p.E - ERangeCut)*WaterRangeCS(p.E);
	ZFloat finalRange = Eend > ERangeCut ? (Eend - ERangeCut)*WaterRangeCS(Eend) : 0;
	ZFloat e = GlueF(0.5)*(p.E + Eend);
	if (e < EAbsElectron) e = EAbsElectron; // limit to Eabs

	//ZFloat uwx = Bx;
	//ZFloat uwy = By;
	//ZFloat uwz = Bz;
	ZFloat Rb = rf*GlueF(sqrt)(e*(e + TEs));
	ZFloat Rbi = 1 / Rb; //in case this value being used many times
	ZFloat maxStep = GlueF(sqrt)(2 * Rb*deltax); //max allowed distance to move to ensure accuracy

	while (true)
	{
		ZFloat step = (range - finalRange) / voxden;

		bool finalStep = true;
		if (step > maxStep)
		{
			step = maxStep;
			finalStep = false;
		}

		//check if it intersect with the boundary of current voxel
		int idex, dvox; //idex = 1,2,3 means x,y,z direction; dvox = +1, -1 means moving in positive or negative direction
		bool intersect = false;

		if (p.v > 0)
		{
			ZFloat next = DY - p.y;
			if (p.v*step > next)
			{
				step = next / p.v; idex = 2; dvox = 1; intersect = true;
			}
		}
		else if (p.v < 0)
		{
			ZFloat next = -p.y;
			if (p.v*step < next)
			{
				step = next / p.v; idex = 2; dvox = -1; intersect = true;
			}
		}

		if (p.w > 0)
		{
			ZFloat next = DZ - p.z;
			if (p.w*step > next)
			{
				step = next / p.w; idex = 3; dvox = 1; intersect = true;
			}
		}
		else if (p.w < 0)
		{
			ZFloat next = -p.z;
			if (p.w*step < next)
			{
				step = next / p.w; idex = 3; dvox = -1; intersect = true;
			}
		}

		if (p.u > 0)
		{
			ZFloat next = DX - p.x;
			if (p.u*step > next)
			{
				step = next / p.u; idex = 1; dvox = 1; intersect = true;
			}
		}
		else if (p.u < 0)
		{
			ZFloat next = -p.x;
			if (p.u*step < next)
			{
				step = next / p.u; idex = 1; dvox = -1; intersect = true;
			}
		}

		if (intersect) finalStep = false;

		ZFloat newEnergy = Eend;
		if (!finalStep)
		{
			range -= step*voxden;
			newEnergy = WaterInverseRangeCS(range);
			if (newEnergy < ERangeCut) newEnergy = 0;
		}
		deposit(p.iabsv, (p.E - newEnergy)*p.weight);
		//up the energy to the new one
		p.E = newEnergy;
		if (p.E < ERangeCut) return 1; // if this is the final step to local absorption, we don't need to update the position and the direction

		//move the electron/positron
		p.x += p.u*step;
		p.y += p.v*step;
		p.z += p.w*step;

		ZFloat vuw = p.u*Bx + p.v*By + p.w*Bz;
		ZFloat vperpx = p.u - vuw * Bx,
			vperpy = p.v - vuw * By,
			vperpz = p.w - vuw * Bz;
		ZFloat vxwx = vperpy*Bz - vperpz*By,
			vxwy = vperpz*Bx - vperpx*Bz,
			vxwz = vperpx*By - vperpy*Bx;
		// The step-length dependent variables f1 & f2
		ZFloat f1, f2;
		ZFloat arg = step * Rbi;
// 		if (arg < GlueF(0.2))
// 		{
// 			// arg is small, so use power series expansion of sine and cosine
// 			ZFloat arg2 = arg*arg;
// 			f1 = -GlueF(0.5)*arg2 + GlueF(0.0416666667)*arg2*arg2;  // for 0.2, relative error is 2.2e-6
// 			f2 = arg - GlueF(0.16666667)*arg*arg2;           // for 0.2, relative error is 1.3e-5, absolute error is 2.6e-6
// 		}
// 		else
/*		{*/
			f1 = GlueF(cos)(arg)-1;
			f2 = GlueF(sin)(arg);
/*		}*/

		// Direction change 
		ZFloat dvx = f1*vperpx - f2*vxwx;  // would simplify to f1*_v.x - f2*_v.y;
		ZFloat dvy = f1*vperpy - f2*vxwy;  // would simplify to f1*_v.y + f2*_v.x;
		ZFloat dvz = f1*vperpz - f2*vxwz;  // would simplify to 0 (i.e., component along the magnetic field remains constant).

		//update the direction
		p.u += dvx;
		p.v += dvy;
		p.w += dvz;

		if (finalStep)
		{
			ZFloat tp = tperp(p);
			ZFloat range = (p.E - ERangeCut)*WaterRangeCS(p.E);
			if (range < tp*voxden) // can deposit without further simulation
			{
				deposit(p.iabsv, p.E*p.weight);
				return 1; //end the simulation of this electron
			}
			return 0;
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
				if (p.w*(p.w - dvz) >= 0) //enter a new voxel
				{
					p.ivz += dvox;
					if (dvox > 0) {
						if (p.ivz >= NZ) return 1;
						p.z = 0;
					}
					else {
						if (p.ivz < 0) return 1;
						p.z = DZ;
					}
				}
				else intersect = false;
				break;
			case 2:
				if (p.v*(p.v - dvy) >= 0)
				{
					p.ivy += dvox;
					if (dvox > 0) {
						if (p.ivy >= NY) return 1;
						p.y = 0;
					}
					else {
						if (p.ivy < 0) return 1;
						p.y = DY;
					}
				}
				else intersect = false;
				break;
			default:
				if (p.u*(p.u - dvx) >= 0)
				{
					p.ivx += dvox;
					if (dvox > 0) {
						if (p.ivx >= NX) return 1;
						p.x = 0;
					}
					else {
						if (p.ivx < 0) return 1;
						p.x = DX;
					}
				}
				else intersect = false;
			}
		}// end if (intersect)

		//still intersect after the direction check, need to update the voxel density and index
// 		if (intersect) // in GPU, just update all even if some are not necessary
// 		{
			p.iabsv = at(p.ivx, p.ivy, p.ivz);
			voxden = getDensity(p);
/*		}*/

// 		//update direction and keep going
// 		p.u += dvx;
// 		p.v += dvy;
// 		p.w += dvz;
	}
}


/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end: phantom method definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/



/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: Kernel definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
__global__ void initThreads(int cardOffset)//init the random number generator, call it before any run
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	RNGState[it].init(cardOffset + it);
}

__device__  __forceinline__ bool refill(size_t it, ParticleR* pInit, ParticleStack& pStack, ParticleR& p, int& cur)
{
	if (pStack.empty())
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
	else
	{
		p = pStack.top(); //it must be already in the phantom
		pStack.pop();
		return true;
	}	
}

__device__ void Rotate(ZFloat& ux, ZFloat& uy, ZFloat& uz, ZFloat costh, ZFloat cosph, ZFloat sinph) {
	ZFloat costh2 = costh*costh;
	ZFloat rho2 = ux * ux + uy * uy;
	if (rho2 > 0 && costh2 < 1) {
		ZFloat a = GlueF(sqrt)((1 - costh2) / rho2);
		ZFloat xrho = ux * a;
		ZFloat yrho = uy * a;
		ZFloat ztmp = uz * cosph;
		ux = ux * costh - yrho * sinph + ztmp * xrho;
		uy = uy * costh + xrho * sinph + ztmp * yrho;
		uz = uz * costh - rho2 * a * cosph;
	}
	else {
		if (costh2 >= 1) {
			if (costh < 0) {
				ux = -ux; uy = -uy; uz = -uz;
			}
			return;
		}
		ZFloat b = GlueF(sqrt)(1 - costh2);
		uy = b * sinph;
		if (uz > 0) {
			ux = b * cosph;
			uz = costh;
		}
		else {
			ux = -b * cosph;
			uz = -costh;
		}
	}
};

// __device__ void Rotate(ZFloat& ux, ZFloat& uy, ZFloat& uz, ZFloat costh, GRNG& rng) //rotate with a random phi
// {
// 	ZFloat phi = 2 * PI*rng();
// 	ZFloat cphiCompt = GlueF(cos)(phi);
// 	ZFloat sphiCompt = (phi < PI ? 1 : -1)*GlueF(sqrt)(1 - cphiCompt*cphiCompt);
// }
__device__ void samcom(ZFloat energy, ZFloat &efrac, ZFloat &costh, GRNG& rng)
{
	ZFloat ko = energy*INV_ELECTRON_MASS;
	ZFloat broi = 1 + 2 * ko; ZFloat bro = 1 / broi;
	ZFloat br, temp;
	if (ko < 10) {
		// "low" energy case: uniformly between bro and bro1.
		ZFloat bro1 = 1 - bro;
		ZFloat ko2 = ko*ko;
		ZFloat rejmax = ko2*(broi + bro);
		ZFloat br2;
		do {
			br = bro + bro1*rng(); br2 = br*br;
		} while (rng()*br2*rejmax > ko2*br*(br2 + 1) - (1 - br)*(br*broi - 1));
		temp = (1 - br) / (ko*br);
	}
	else {
		// "high" energy case: the usual way 
		ZFloat broi2 = broi*broi;
		ZFloat alpha1 = GlueF(log)(broi);
		ZFloat alpha2 = ko*(broi + 1)*bro*bro;
		ZFloat alphaS = alpha1 + alpha2;
		ZFloat sint;
		do {
			br = rng()*alphaS < alpha1 ? GlueF(exp)(alpha1*rng())*bro : GlueF(sqrt)(rng()*(broi2 - 1) + 1)*bro;
			temp = (1 - br) / (ko*br); sint = temp*(2 - temp);
		} while (rng()*(1 + br*br) < br*sint);
	}
	efrac = br;
	costh = 1 - temp;
}



__device__ ZFloat samsca(ZFloat e, GRNG& rng)
{
	// The screening parameter at this energy.
	ZFloat ie = -1 / e;
	ZFloat b = WaterScreenCS(ie);
	ZFloat oneb = GlueF(1.0) + b;

	// Determine energy bin of the pre-computed q surface.
	// Instead of linearly interpolating between energies in each qsurf() evaluation below, we 
	// use the lower or higher energy of the bin in which the energy falls with the correponding probability
	int je; ZFloat pe;
	WaterQSGetEnergyIndex1(ie, je, pe);
	if (pe > 0)
	{
		if (rng() < pe) ++je;
	}

	ZFloat u;
	while (1)
	{
		u = rng();//u is in [0, 1)
		if (rng() < WaterQSurface(je, u)) break;
	}
	ZFloat mu = (oneb - u*(oneb + b)) / (oneb - u); // oneb>=1.0, so (oneb-u) is always positive
	return mu;
}

__device__ void electronRun(ParticleR& p, GRNG& rng)
{
	ZFloat fuelxt = 0, ebefor = 0, nextEnergy = 0;

	while (true)
	{
		ebefor = p.E - fuelxt;
		if (ebefor > EAbsElectron)
		{
			ZFloat Eloss = min(ebefor, EMaxCSDA);
			ZFloat fuel = Eloss*rng();
			nextEnergy = p.E - fuelxt - fuel;
			if (nextEnergy <= ERangeCut) nextEnergy = 0;
			fuelxt = Eloss - fuel;
		}
		else
		{
			if (p.E < ERangeCut) return;
			nextEnergy = 0;
		}

		int status = electronFlight(p, nextEnergy);



		if (status == 1) return; //below the cut-off energy or exit the phantom

		// Check if we can discontinue transport because the electron cannot escape the current voxel
// 		ZFloat tp = tperp(p);
// 		ZFloat range = (p.E - ERangeCut)*WaterRangeCS(p.E);
// 		if (range < tp*getDensity(p))
// 		{
// 			deposit(p.iabsv, p.E*p.weight);
// 			return; //end the simulation of this electron
// 		}

		//do elastic multi-scattering
		ZFloat costhe = samsca(ebefor, rng);
		if (costhe < 1)
		{
			ZFloat cphi, sphi;
			randomAzimuth(rng(), cphi, sphi);
			Rotate(p.u, p.v, p.w, costhe, cphi, sphi);
		}
	}
}

__device__ void photonRun(ParticleR& p, GRNG& rng, int& hasSec)
{
	//this version put electron in stack and simulate later

	while (true) //break until this photon is finished
	{
		ZFloat lamph = WaterPhotonCS(p.E);
		ZFloat lammin = 1 / (lamph*MaxDensity);
		
		ZFloat s = -lammin*GlueF(log)(rng());
		if (photonFlight(p, s)) return; //finish simulating this photon
		ZFloat lamden = lammin * getDensity(p);

		ZFloat rn = rng();
		if (rn < lamden*lamph)
		{
			ZFloat lamco = WaterComptonCS(p.E);
			if (rn < lamden * lamco)// It's a Compton interaction
			{
				ParticleR pe = p; // back the status of this photon

				ZFloat efracCompt = 1, costheCompt, cphiCompt, sphiCompt;
				samcom(p.E, efracCompt, costheCompt, rng);
				randomAzimuth(rng(), cphiCompt, sphiCompt);

				// Electron energy
				efracCompt = 1 - efracCompt;
				pe.E = p.E*efracCompt;
				// Compute Compton electron direction
				ZFloat e0 = p.E * INV_ELECTRON_MASS;
				ZFloat cost;
				if (efracCompt > FLTEPSILON)
				{
					cost = (1 + e0) * GlueF(sqrt)(efracCompt / (e0*(2.0 + e0*efracCompt)));
					if (cost > 1) cost = 1;
				}
				else cost = 0;
				Rotate(pe.u, pe.v, pe.w, cost, -cphiCompt, -sphiCompt);

				p.E -= pe.E; // the scattered photon's energy

				electronRun(pe, rng);

				if (p.E > EAbsPhoton)
				{
					Rotate(p.u, p.v, p.w, costheCompt, cphiCompt, sphiCompt);
				}
				else return;
			}
			else //should be pair production or photoelectric absorption
			{
				bool doPair = p.E > TEs; //only possible to do pair production
				if (doPair)
				{
					ZFloat lampair = WaterPairCS(p.E);
					if (rn > lamden*(lamco + lampair)) doPair = false;
				}

				if (doPair) //it's a pair production
				{
					ParticleR pe(p);
					pe.type = electron;
					ZFloat Epair1 = rng() * (p.E - TEs);
					pe.E = Epair1;
					electronRun(pe, rng);
					pe.E = p.E - TEs - Epair1;
					pe.copyDirection(p);
					pe.copyPosition(p);
					electronRun(pe, rng);

					p.w = 2 * rng() - 1;
					ZFloat sinthe = GlueF(sqrt)(1.0 - p.w*p.w);
					ZFloat phi = 2 * PI * rng();
					p.u = sinthe * GlueF(cos)(phi);
					p.v = sinthe * GlueF(sin)(phi);
					p.E = Es; //for debug, it should be Es

					size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
					StackBuff[NStackDepth*it] = p; //since the event is rare, we can use slow global memory
					hasSec = 1;

					//The other photon has opposite direction
					p.u = -p.u;
					p.v = -p.v;
					p.w = -p.w;
				}
				else //it's a photoelectric absorption
				{
					p.type = electron;
					electronRun(p, rng);
					return;//finish simulating this photon
				}
			}
		}
		//else //reject the scattering and continue next move
	}
}

__device__ void smartPhotonRun(ParticleR& p, GRNG& rng, int& hasSec)
{
	int nsplit = NMaxSplit;
	while (true)
	{
		if (!FixedSplit)
		{
			//splitting number
			ZFloat e = p.E*GlueF(1e-6);
			ZFloat e2 = e*e;
			ZFloat pInitial = e > GlueF(0.3) ? GlueF(-0.053139) + GlueF(1.0695) * e - GlueF(0.24783) * e2 + GlueF(0.028566) * e*e2 : GlueF(0.25);
			nsplit = (int)(pInitial*NMaxSplit);
		}

		ZFloat rnno1 = rng();
		ZFloat rnno2 = rng();
		ZFloat delta = 1 / ZFloat(nsplit);
		ZFloat eweight = p.weight*delta;

		ZFloat lamph = WaterPhotonCS(p.E);
		ZFloat lammin = 1 / (lamph*MaxDensity);
		ZFloat lamco = -1;
		ZFloat lampair = -1;

		int keepID = (int)(rng()*nsplit);
		// These will be used to remember the relevant Compton scattering variables
		ZFloat eCompElec = -1;
		ZFloat eu = p.u, ev = p.v, ew = p.w;
		ZFloat costheCompt, cphiCompt, sphiCompt; //used to calculate the direction of the scattered photon

		ZFloat Epair1 = -1; // remember the pair production energy

		ParticleR pOld = p; // remember the initial status
		p.E = 0; //default to exit
		for (int isplit = 0; isplit < nsplit; ++isplit)
		{
			ParticleR pi = pOld; //start with the initial status
			ZFloat rnnoi = 1 - delta*(rnno1 + isplit);
			if (rnnoi <= 0) break;
			//ZFloat s = lammin*lambdaCalculator(rnnoi);
			ZFloat s = -lammin*GlueF(log)(rnnoi);

			if (photonFlight(pi, s)) break;

			// update lammin with density in the voxel and get material id.
			ZFloat lamden = lammin * getDensity(pi);

			if (lamph < 0) lamph = WaterPhotonCS(pi.E); // it's so weird that this statement improved the performance

			// Check if a real interaction
			if (rnno2 < lamden*lamph) // yes. 
			{
				if (lamco < 0) lamco = WaterComptonCS(pi.E);

				if (rnno2 < lamden * lamco) // It's a Compton interaction
				{
					if (eCompElec < 0) // Haven't sampled a Compton interaction yet, so do it now.
					{
						//keepCompton = (int)(rng()*nsplit);
						// Sample the interaction
						ZFloat efracCompt = 1;
						samcom(pi.E, efracCompt, costheCompt, rng);
						randomAzimuth(rng(), cphiCompt, sphiCompt);

						// Electron energy
						efracCompt = 1 - efracCompt;
						eCompElec = pi.E*efracCompt;
						// Compute Compton electron direction
						ZFloat e0 = pi.E * INV_ELECTRON_MASS;
						//ZFloat efrac1 = 1 - efracCompt;
						ZFloat cost;
						if (efracCompt > FLTEPSILON)
						{
							cost = (1 + e0) * GlueF(sqrt)(efracCompt / (e0*(2.0 + e0*efracCompt)));
							if (cost > 1) cost = 1;
						}
						else cost = 0;
						Rotate(eu, ev, ew, cost, -cphiCompt, -sphiCompt);
					}

					//if (isplit == keepCompton && eCompGamma > EAbsPhoton)
					if (isplit == keepID)
					{
						ZFloat eCompGamma = pi.E - eCompElec;
						p.E = eCompGamma;
						if (eCompGamma > EAbsPhoton)
						{
							Rotate(p.u, p.v, p.w, costheCompt, cphiCompt, sphiCompt);
							p.x = pi.x;
							p.y = pi.y;
							p.z = pi.z;
							p.ivx = pi.ivx;
							p.ivy = pi.ivy;
							p.ivz = pi.ivz;
							p.iabsv = pi.iabsv;
						}
						
					}

					// Now, instead of first pushing the Compton electron onto the stack and later getting it back, we simply transport it here.
					pi.type = electron;
					pi.E = eCompElec;
					pi.weight = eweight;
					pi.u = eu;
					pi.v = ev;
					pi.w = ew;
					electronRun(pi, rng);
				}
				else // Not a Compton, so check if pair or photo.
				{
					bool doPair = pi.E > TEs;
					if (doPair)
					{
						if (lampair < 0) lampair = WaterPairCS(pi.E);
						if (rnno2 > lamden*(lamco + lampair)) doPair = false;
					}
					if (doPair) // It's a pair production -> the photon disappears (but we add annihilation photons below as needed).
					{
						if (Epair1 < 0) // Haven't sampled a pair event yet, so do it now. 
						{
							Epair1 = rng() * (pi.E - TEs);
						}

						if (isplit == keepID)
						{
							p.w = 2 * rng() - 1;
							ZFloat sinthe = GlueF(sqrt)(1.0 - p.w*p.w);
							ZFloat phi = 2 * PI * rng();
							p.u = sinthe * GlueF(cos)(phi);
							p.v = sinthe * GlueF(sin)(phi);
							p.E = Es; //for debug, it should be Es
							//copy the position
							p.x = pi.x; p.y = pi.y; p.z = pi.z;
							p.ivx = pi.ivx; p.ivy = pi.ivy; p.ivz = pi.ivz; p.iabsv = pi.iabsv;

							size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
							StackBuff[NStackDepth*it] = p; //since the event is rare, we can use slow global memory
							hasSec = 1;

							//The other photon has opposite direction
							p.u = -p.u;
							p.v = -p.v;
							p.w = -p.w;
						}

// 						ParticleR* pars = (ParticleR*)sharedMem;
// 						ParticleR& pb = pars[threadIdx.x]; //current particle to process
// 						pb = pi;

						// Put an e+/e- pair on the stack. We do not distinguish between electrons and positrons at this point.
						ZFloat bx = pi.x; ZFloat by = pi.y; ZFloat bz = pi.z; //backup position
						int bix = pi.ivx; int biy = pi.ivy; int biz = pi.ivz; int biabs = pi.iabsv;

						pi.type = electron;
						pi.E = Epair1;
						pi.weight = eweight;
						electronRun(pi, rng);
						//restore the position, direction
						pi.x = bx; pi.y = by; pi.z = bz;
						pi.ivx = bix; pi.ivy = biy; pi.ivz = biz; pi.iabsv = biabs;
						pi.u = pOld.u; pi.v = pOld.v; pi.w = pOld.w;
						pi.E = pOld.E - TEs - Epair1;
// 						pi = pb;
// 						pi.E = pOld.E - TEs - Epair1;
						electronRun(pi, rng);
					}
					else
					{
						// It's a photo absorption -> the photon disappears
						pi.type = electron;
						pi.weight = eweight;
						electronRun(pi, rng);
					}
				}
			}
			else // The interaction was rejected. 
			{
				if (isplit == keepID) //copy the position, and energy
				{
					p.x = pi.x; p.y = pi.y; p.z = pi.z;
					p.ivx = pi.ivx; p.ivy = pi.ivy; p.ivz = pi.ivz; p.iabsv = pi.iabsv;
					p.E = pi.E;
				}
			}
		}
		if (p.E < EAbsPhoton) return;
	}
}

__device__ void ComptonPhotonRun(ParticleR& p, GRNG &rng)
{
	while (true)
	{
		ZFloat lamph = WaterPhotonCS(p.E);
		ZFloat lammin = 1 / (lamph*MaxDensity);
		ZFloat lamco = WaterComptonCS(p.E);
		ZFloat s = -lammin*GlueF(log)(rng());

		if (photonFlight(p, s)) return;

		// update lammin with density in the voxel and get material id.
		ZFloat lamden = lammin * getDensity(p);
		ZFloat rnno = rng();
		// Check if a real interaction
		if (rnno < lamden*lamph) // yes. 
		{
			if (rnno < lamden * lamco) // It's a Compton interaction
			{
				ZFloat efracCompt, costheCompt;
				samcom(p.E, efracCompt, costheCompt, rng);
	
				//deposit the energy of electron
				ZFloat eCompGamma = efracCompt*p.E;
				deposit(p.iabsv, (p.E - eCompGamma)*p.weight);
				
				if (eCompGamma > EAbsPhoton)
				{
					ZFloat cphiCompt, sphiCompt;
					randomAzimuth(rng(), cphiCompt, sphiCompt);

					Rotate(p.u, p.v, p.w, costheCompt, cphiCompt, sphiCompt);
					p.E = eCompGamma;
				}
				else
				{
					//deposit(p.iabsv, eCompGamma*p.weight);
					return;
				}

			}
			else // Not a Compton, deposit all energy here
			{
				deposit(p.iabsv, p.E*p.weight);
				return;
			}
		}
	}
}
__device__ void SmartComptonPhotonRun(ParticleR& p, GRNG &rng)
{
	int nsplit = NMaxSplit;
	while (true)
	{
		if (!FixedSplit)
		{
			//splitting number
			ZFloat e = p.E*GlueF(1e-6);
			ZFloat e2 = e*e;
			ZFloat pInitial = e > GlueF(0.3) ? GlueF(-0.053139) + GlueF(1.0695) * e - GlueF(0.24783) * e2 + GlueF(0.028566) * e*e2 : GlueF(0.25);
			nsplit = (int)(pInitial*NMaxSplit);
		}

		ZFloat rnno1 = rng();
		ZFloat rnno2 = rng();
		ZFloat delta = 1 / ZFloat(nsplit);
		ZFloat eweight = p.weight*delta;

		ZFloat lamph = WaterPhotonCS(p.E);
		ZFloat lammin = 1 / (lamph*MaxDensity);
		ZFloat lamco = WaterComptonCS(p.E);

		int keepID = (int)(rng()*nsplit);
		// These will be used to remember the relevant Compton scattering variables
		ZFloat costheCompt = 0, eCompElec = -1;

		//remember initial position
		ZFloat px = p.x, py = p.y, pz = p.z;
		int pivx = p.ivx, pivy = p.ivy, pivz = p.ivz;
		ParticleR pi = p; //the direction of pi will not change
		p.E = 0; //default to lose all energy and exit
		for (int isplit = 0; isplit < nsplit; ++isplit)
		{
			//reset the position
			pi.x = px;
			pi.y = py;
			pi.z = pz;
			pi.ivx = pivx;
			pi.ivy = pivy;
			pi.ivz = pivz;

			ZFloat rnnoi = 1 - delta*(1 - rnno1 + isplit);
			
			//ZFloat s = lammin*lambdaCalculator(rnnoi); //this implementation is slower than calling log directly
			ZFloat s = -lammin*GlueF(log)(rnnoi);

			if (photonFlight(pi, s)) break;

			ZFloat lamden = lammin*getDensity(pi);

			// Check if a real interaction
			if (rnno2 < lamden*lamph) // yes. 
			{
				if (rnno2 < lamden * lamco) // It's a Compton interaction
				{
					if (eCompElec < 0) // Haven't sampled a Compton interaction yet, so do it now.
					{
						ZFloat efracCompt = 1;
						samcom(pi.E, efracCompt, costheCompt, rng);
						eCompElec = pi.E*(1 - efracCompt);
					}
					if (isplit == keepID)
					{
						p.E = pi.E - eCompElec;
						if (p.E > EAbsPhoton)
						{
							ZFloat cphiCompt, sphiCompt;
							randomAzimuth(rng(), cphiCompt, sphiCompt);

							Rotate(p.u, p.v, p.w, costheCompt, cphiCompt, sphiCompt);
							//copy position
							p.x = pi.x;
							p.y = pi.y;
							p.z = pi.z;
							p.ivx = pi.ivx;
							p.ivy = pi.ivy;
							p.ivz = pi.ivz;
							p.iabsv = pi.iabsv;
						}
					}

					deposit(pi.iabsv, eCompElec*eweight); //deposition the energy of the electron
				}
				else // Not a Compton, deposit all energy here
				{
					deposit(pi.iabsv, pi.E*eweight);
				}
			}
			else // The interaction was rejected. 
			{
				if (isplit == keepID)
				{
					//copy energy and position. Direction is actually unchanged
					p.E = pi.E;
					p.x = pi.x;
					p.y = pi.y;
					p.z = pi.z;
					p.ivx = pi.ivx;
					p.ivy = pi.ivy;
					p.ivz = pi.ivz;
					p.iabsv = pi.iabsv;
				}
			}
		}
		if (p.E < EAbsPhoton)
		{
			//if (p.E > 0) deposit(p.iabsv, p.E*p.weight);
			return;
		}
	}
}


__global__ void
//__launch_bounds__(128, 16)
gZeusSmartRun(ParticleR* pInit) //let's implement in a simple way first
{
	unsigned int it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	GRNG rng = RNGState[it]; //get the RNG status
	//ParticleR* pars = (ParticleR*)sharedMem;
	//ParticleR& p = pars[threadIdx.x]; //current particle to process
	ParticleR p; //current particle to process
	int hasSec = 0;
	int cur = 0;//reset the current particle index
	int NT = blockDim.x * gridDim.x;
	while (true)
	{
		//first try to fetch particle from particle stack, and then try to fetch from the particle buffer
		if (hasSec)
		{
			p = StackBuff[NStackDepth*it];
			hasSec = 0;
		}
		else //fetch from the buffer
		{
			while (true)
			{
				if (cur == NBatch) // exhausted the buffer
				{
					RNGState[it] = rng; //record the status of current RNG
					return; // exit this thread
				}
				p = pInit[cur*NT + it];
				//p = pInit[it*NBatch + cur];
				++cur;
				if (lineIn(p)) break;
			}
		}

		smartPhotonRun(p, rng, hasSec);
	}

}

__global__ void gZeusRun(ParticleR* pInit) //let's implement in a simple way first
{
	unsigned int it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	GRNG rng = RNGState[it]; //get the RNG status
	//ParticleR* pars = (ParticleR*)sharedMem;
	//ParticleR& p = pars[threadIdx.x]; //current particle to process
	ParticleR p; //current particle to process
	int hasSec = 0;
	int cur = 0;//reset the current particle index
	while (true)
	{
		//first try to fetch particle from particle stack, and then try to fetch from the particle buffer
		if (hasSec)
		{
			p = StackBuff[NStackDepth*it];
			hasSec = 0;
		}
		else //fetch from the buffer
		{
			while (true)
			{
				if (cur == NBatch) // exhausted the buffer
				{
					RNGState[it] = rng; //record the status of current RNG
					return; // exit this thread
				}
				p = pInit[it*NBatch + cur];
				++cur;
				if (lineIn(p)) break;
			}
		}

		photonRun(p, rng, hasSec);
	}

}

__global__ void gZeusComptonRun(ParticleR* pInit) //let's implement in a simple way first
{
	unsigned int it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	GRNG rng = RNGState[it]; //get the RNG status
	ParticleR p; //current particle to process
	int cur = 0;//reset the current particle index
	while (true)
	{
		while (true)
		{
			if (cur == NBatch) // exhausted the buffer
			{
				RNGState[it] = rng; //record the status of current RNG
				return; // exit this thread
			}
			p = pInit[it*NBatch + cur];
			++cur;
			if (lineIn(p)) break;
		}

		ComptonPhotonRun(p, rng);
	}

}

__global__ void gZeusSmartComptonRun(ParticleR* pInit) //let's implement in a simple way first
{
	unsigned int it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	GRNG rng = RNGState[it]; //get the RNG status
	ParticleR p; //current particle to process
	int cur = 0;//reset the current particle index
	int NT = blockDim.x * gridDim.x;
	while (true)
	{
		while (true)
		{
			if (cur == NBatch) // exhausted the buffer
			{
				RNGState[it] = rng; //record the status of current RNG
				return; // exit this thread
			}
			p = pInit[cur*NT + it];
			++cur;
			if (lineIn(p)) break;
		}

		SmartComptonPhotonRun(p, rng);
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



/*<<<<<<<<<<<<<<<<<<<<<<<<< start: GZEUS method definitions <<<<<<<<<<<<<<<<<<<<<<<<<<*/
void GZEUS::phantom2GPU() //copy phantom info to GPU. Make sure the phantom has loaded data
{
	int NVoxel = _phant->getVoxelNum();
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

void GZEUS::initGPU()
{
	RunTimeCounter rc;
	if (!initWaterCS("ZeusCrossSections.binary", gc)) exitApp("Cannot initiate the water cross sections");
	int NGPU = (int)gc.size();
	ConfigFile* zcf = _cf->getBlock("ZEUS");
	if (zcf == NULL) zcf = _cf->getBlock("PHANTOM"); //in case cannot find the Zeus block; use default values instead

	int split = 50, h_NStackDepth = 40;
	zcf->getValue("NMaxSplit", split);
	zcf->getValue("Particle stack depth", h_NStackDepth);
	string str = "yes";
	int h_fixedSplit = 1; //default yes
	if (zcf->getValue("Fixed Split", str) && str.compare("yes") != 0) h_fixedSplit = 0;
	int h_simuElectron = 1; //default yes
	if (zcf->getValue("Simulate electron", str) && str.compare("yes") != 0) h_simuElectron = 0;
	int h_forwardDetection = 1;//default yes
	if (zcf->getValue("Forward detection", str) && str.compare("yes") != 0) h_forwardDetection = 0;

	double fv = 0;
	ZFloat h_EAbsPhoton = 50e3; //unit eV
	if (zcf->getValue("EAbsPhoton", fv)) h_EAbsPhoton = (ZFloat)fv;
	ZFloat h_EAbsElectron = 50e3; //unit eV
	if (zcf->getValue("EAbsElectron", fv)) h_EAbsElectron = (ZFloat)fv;
	ZFloat h_EMaxCSDA = 200e3; //unit eV
	if (zcf->getValue("EMaxCSDA", fv)) h_EMaxCSDA = (ZFloat)fv;

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
		cudaErrorCheck(cudaMemcpyToSymbol(ForwardDetect, &h_forwardDetection, sizeof(int)));

		cudaErrorCheck(cudaMemcpyToSymbol(EAbsPhoton, &h_EAbsPhoton, sizeof(ZFloat)));
		cudaErrorCheck(cudaMemcpyToSymbol(EAbsElectron, &h_EAbsElectron, sizeof(ZFloat)));
		cudaErrorCheck(cudaMemcpyToSymbol(EMaxCSDA, &h_EMaxCSDA, sizeof(ZFloat)));
	}
	Log("\nIt costs %f seconds to init GPU ", rc.stop());
}

void GZEUS::freeGPU()//do some clean up
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
		gc[i].destroyStream();
	}
}

void GZEUS::init(ConfigFile* cf)
{
	_cf = cf;
	initGPU();

	ConfigFile* ph_cf = cf->getBlock("PHANTOM");
	_phant = new Phantom;
	_phant->loadPhantom(ph_cf);
	string lastDoseFile;
	if (cf->getValue("proceed last simulation", lastDoseFile) && lastDoseFile.compare("yes") == 0)
	{
		cf->getValue("output file name", lastDoseFile);
		lastDoseFile += ".dose";
		if (!_phant->previousDose(lastDoseFile.c_str())) exitApp("Cannot load last dose file to continue the simulation!");
		else Log("load last dose file successfully with %.0f existing histories", _phant->getHist());
	}
	SourceHead_GetPrescrition(&(_phant->prescriptionDose), &(_phant->treatmentFraction));

	phantom2GPU();
}

int GZEUS::getGPUConfig(ConfigFile* gcf)
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
/*>>>>>>>>>>>>>>>>>>>>>>>>>> end: GZEUS method definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


RunTimeCounter sourceCounter;

void getSource(SourcePool* sp, volatile int* hist)
{
	sourceCounter.start();
	*hist = sp->prepareCopy();//prepare one batch before any run

	Log("time cost to generate particle = %f s", sourceCounter.stop());
}

int NThread, NProcess, pID;
bool b_peek = false;

BinaryFile gBF; //to store the raw dose data

ProgressCallBack progressCB = NULL;
PeekDoseCallBack peekDoseCB = NULL;
JobFinishedCallBack jobFinishedCB = NULL;
LogCallBack exitCB = NULL;

// Begin: essential for command line call
bool b_thread_active = false;
bool b_abort = false;

#ifdef WIN32
#define DEXPORT __declspec(dllexport)
#else
#define DEXPORT __attribute__ ((visibility ("default")))
#endif


extern "C" DEXPORT void startSimulation(const char* configFileName, MPS& configMacro, bool bWait) //launch a thread to do the simulation
{
	b_thread_active = true;
	std::thread thd(executeJob, configFileName, 1, std::ref(configMacro));
	if (bWait) thd.join();
	else thd.detach();
}

extern "C" DEXPORT void stopSimulation()
{
	if (b_thread_active) b_abort = true; //only if the thread is on
}
// End: essential for command line call

extern "C" DEXPORT void peekDose() //mark to generate intermediate dose
{
	if (b_thread_active) b_peek = true; //only if the thread is on
}

//interfaces to register the call-back function
extern "C" DEXPORT void setProgressCallBack(ProgressCallBack pcb){ progressCB = pcb; }
extern "C" DEXPORT void setPeekDoseCallBack(PeekDoseCallBack pcb){ peekDoseCB = pcb; }
extern "C" DEXPORT void setJobFinishedCallBack(JobFinishedCallBack jcb){ jobFinishedCB = jcb; }
extern "C" DEXPORT void setLogCallBack(LogCallBack lcb){ Log.setCallBack(lcb); }
extern "C" DEXPORT void setExitCallBack(LogCallBack ecb){ exitCB = ecb; }

void exitApp(const char *inf)
{
	Log("fatal error: %s", inf);
	Log.flush(); //flush the log file's buffer 
	if (exitCB) exitCB(inf); // let the call back function display the error
	else
	{
		Log("\nPress enter key to exit...");
		getchar();
	}
	exit(-1);
}

void executeJob(const char* configFileName, double processWeight, MPS& configMacro) //execute one job according to the config file
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
	Log.setLogName(pID, logAppend, logDir);//start log recording for this job
	Log("The job config file name = %s", configFileName);
	if (logDescription.compare("NA") != 0) Log("Short description: %s", logDescription.c_str());
	Log("Start log time = %s\n\n", Log.timeNow());

	double fNSIMU = 1e7;
	cf.getValue("NSIMU", fNSIMU);
	fNSIMU *= processWeight; //real workload on this node
	double targetErr = -1;
	cf.getValue("target uncertainty", targetErr);
	double targetErr2 = targetErr*targetErr;
	double thresholdRegion = 0.5; //default 50% of the max dose as threshold for uncertainty calculation
	cf.getValue("threshold region", thresholdRegion);
	string outname;
	cf.getValue("output file name", outname);

	ConfigFile* zcf = cf.getBlock("ZEUS");
	if (zcf == NULL) zcf = cf.getBlock("PHANTOM"); //in case cannot find the Zeus block; use default values instead
	int nsplit = 50;
	zcf->getValue("NMaxSplit", nsplit);
	bool simuElectron = true;
	string str;
	if (zcf->getValue("Simulate electron", str) && str.compare("yes") != 0) simuElectron = false;

	//search and config the GPU part ->> get gc
	GZEUS zeus;
	vector<GPUConfig>& gc = zeus.gc;//makes the name shorter
	ConfigFile *gcf = cf.getBlock("GPU");
	int main_id = zeus.getGPUConfig(gcf);

	//initialize GZEUS and SourceHead by configurations
	ConfigFile *scf = cf.getBlock("SOURCEHEAD");
	SourceHead_Init(scf);
	
#ifdef SOURCE_STATIStICS
	//for the source energy statistics
	PRNG _rng;
	_rng.init(1234);
	Particle pars[100];
	int NES = int(1e6);
	double* Ens = new double[NES];
	int nsam = 0;
	while (nsam < NES)
	{
		int np = SourceHead_Sample(&_rng, pars);
		for (int i = 0; i < np; ++i)
		{
			Ens[nsam] = pars[i].E;
			++nsam;
			if (nsam >= NES) break;
		}
	}
	FILE* fps = fopen("E.txt", "wb");
	fwrite(Ens, NES, 1, fps);
	fclose(fps);
	delete[] Ens;
#endif

	//string dataDir;
	//scf->getValue("DataDir",dataDir);
	//if (!ZeusData_load(dataDir.c_str())) exitApp("Cannot load Zeus cross-sections correctly!");
	zeus.init(&cf); //prepare GPU data and phantom

	fNSIMU -= zeus._phant->Hist; //subtract previous histories
	if (fNSIMU <= 0)
	{
		Log("Don't need to run any more history! Skip this task...\n\n");
		return;
	}

	//config the source particle pool
	int NGT = 1, NGStack = 400;
	scf->getValue("NThread", NGT);
	scf->getValue("Sample Stack Depth", NGStack);
	int NOneFetch = gc[0].getNOneBatch();
	SourcePool sp(&gc, NGT, NOneFetch, zeus._phant, NGStack);

	Log("\nCalculating dose, please wait patiently...\n\n");
	cf.getValue("log short mode", logDescription);
	if (logDescription.compare("yes") == 0) Log.shortMode(true);
	RunTimeCounter rc; //count the calculating time

	//note history number != particle number in the source pool
	const int NGPU = (int)gc.size();
	//history number generated once by the source pool, only modified by one thread,
	//but accessed by multiple threads. It can be read only after the generating thread finished.
	volatile int histNew = 0; //histNew isn't always the same as hisAdd because histNew is modified in an isolated thread
	volatile int histAdd = 0; //this variable is shared by all threads, so add the key word "volatile" for safe
	const int firstSeed = zeus._phant->seedBegin( NGPU *gc[0].NBlock*gc[0].BlockSize);// note it's different from gPENELOPE
	std::thread sthread; //source generating thread, unattached

	RunTimeCounter kernelCounter;
	RunTimeCounter copyCounter;

	sthread = std::thread(&getSource, &sp, &histNew);

	vector<SFloat*> dose(NGPU);//to store dose from all GPU cards
	vector<SFloat*> uncertainty(NGPU);
	int NVoxel = zeus._phant->getVoxelNum();
	int nBatch = 0; //count the batch number that has been done
	const int NBlock = gc[0].NBlock;
	const int BlockSize = gc[0].BlockSize;

	// module of auto reuse source particles
	volatile int Source_Reuse_Times = gc[0].SourceReuseTimes; //maybe modified by one thread, so use volatile for safe
	bool autoSourceReuse = false; //used for auto reuse
	if (Source_Reuse_Times < 1)
	{
		autoSourceReuse = true;
		Source_Reuse_Times = 1;
	}

	bool b_targetErrReached = false;
	//each thread takes care of one GPU
#ifdef USE_OPENMP
#pragma omp parallel num_threads(NGPU)
#endif
	{
		int it = omp_get_thread_num();
		cudaErrorCheck(cudaSetDevice(gc[it].id)); //set which GPU this thread will operate on

#ifdef USE_SINGLE_PRECISION
		WaterQS aWaterQS((float*)h_WaterQSData, NQSURFACE_Q, NQSURFACE_E);
#endif
		cudaErrorCheck(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
		double fHMax = fNSIMU / NGPU;
		double hist = 0;

		//std::mt19937 mtrng((unsigned)std::chrono::system_clock::now().time_since_epoch().count());
		//generate a random thread for each working thread
		//int seed = int((it + 1) * 12345789 * (double(mtrng()) / mtrng.max() + 0.5)); //???

		SFloat* gdose = new SFloat[NVoxel]; //to fetch temporary dose from GPU
		memset(gdose, 0, sizeof(SFloat)*NVoxel);
		// resize memory for CPU end storage
		dose[it] = new SFloat[NVoxel]; //to store the final dose of this thread
		uncertainty[it] = new SFloat[NVoxel]; //to store the uncertainty
		// initialize the dose score
		for (int i = 0; i < NVoxel; ++i)
		{
			dose[it][i] = zeus._phant->dose[i] / NGPU;
			uncertainty[it][i] = zeus._phant->uncertainty[i] / NGPU;
		}

		int seed = firstSeed + it*NBlock*BlockSize; //make sure all seeds are unique; note it's different from gPENELOPE
		initThreads << <NBlock, BlockSize >> >(seed);
		cudaKernelCheck(0);

		int source_reuse = Source_Reuse_Times; //count how many times has been reused; force to initially generate incident particles
		ParticleR* pInit = NULL;
		while(true) //calculating main loop, end when hist >= fHMax
		{
			//need to regenerate initial particles
			if (source_reuse >= Source_Reuse_Times)
			{
				if (it == main_id) //it's the main thread
				{
					//wait until the source is ready
					sthread.join();
					histAdd = histNew; //means the main GPU has histAdd new histories
				}
#pragma omp barrier //wait until all GPU threads arrive here
				pInit = sp.getAP(it); //update the particle array pointer
				source_reuse = 0; //reset the reuse counter
#pragma omp barrier //wait until all GPUs received the data
				//if (it == main_id) sthread = std::thread(&getSource, &sp, &histNew);
				if (it == main_id) sthread = std::thread(&getSource, &sp, &histNew); //start a fetch prepare
			}
			++source_reuse; //count how many time the source has been used
			hist += histAdd; // histAdd more histories will be simulated
			
			/****************** Begin a batch run on GPU ********************/
			if (it == main_id) kernelCounter.start(); //only count the kernel time for the main GPU
			int sharedSize = 0; // sizeof(float)*BlockSize * 16;
			if (simuElectron)
			{
				if (nsplit > 1) gZeusSmartRun << <NBlock, BlockSize, sharedSize, gc[it].kernelstream >> >(pInit);
				else gZeusRun << <NBlock, BlockSize, sharedSize, gc[it].kernelstream >> >(pInit);
			}
			else
			{
				if (nsplit>1) gZeusSmartComptonRun << <NBlock, BlockSize, sharedSize, gc[it].kernelstream >> >(pInit);
				else gZeusComptonRun << <NBlock, BlockSize, sharedSize, gc[it].kernelstream >> >(pInit);
			}
			
			cudaStreamSynchronize(gc[it].kernelstream);//wait for the kernel to finish

			//print the speed information
			if (it == main_id)
			{
				Log("time cost to execute kernels = %f s", kernelCounter.stop());
				if (targetErr <= 0)
				{
					double time = rc.stop(true);
					double speed = hist / time;
					double rest = 0;
					if (fHMax > hist) rest = (fHMax - hist) / speed;
					else rest = 0;
					Log("GPU processed ------------------------ %3.1f%%,   speed = %d h/s\n", hist*100.0 / fHMax, int(speed));
					Log("Time escaped = %.1f min, left time expected = %.1f min", time / 60.0, rest / 60.0);
					if (progressCB) progressCB(hist*100.0 / fHMax, speed, time, rest, 0);
				}
				++nBatch;
				if (nBatch == 10 && autoSourceReuse) // may change Source_Reuse_Times based on the performance of the first 10 batches
				{
					int rt = (int)round(sourceCounter.getStoredTime() / kernelCounter.getStoredTime());
					if (rt > 1) Source_Reuse_Times = rt;
				}
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
			cudaErrorCheck(cudaMemcpy(gc[it].d_doseScore, gdose, sizeof(SFloat)*NVoxel, cudaMemcpyHostToDevice)); //reset the dose counter in GPU
			
			if (it == main_id && b_peek)
			{
				SFloat* d = new SFloat[NVoxel];
				SFloat* u = new SFloat[NVoxel];
				//add up dose and uncertainty in different GPU card
				for (int j = 0; j < NVoxel; ++j)
				{
					d[j] = u[j] = 0;
					for (int i = 0; i < NGPU; ++i)
					{
						d[j] += dose[i][j];
						u[j] += uncertainty[i][j];
					}
				}
				//Log("\nOutputing the intermediate dose:");
				Phantom tempPhant(*zeus._phant);
				double norm = 1.0889e15 * SourceHead_BeamOnTime();
				tempPhant.addDose(d, u, NGPU * nBatch, NGPU * hist, norm, thresholdRegion);
				delete[] d;
				delete[] u;
				tempPhant.getBinaryFile(gBF);

				b_peek = false;
				if (peekDoseCB) peekDoseCB(gBF);
			}

			if (targetErr > 0) //calculate the uncertainty
			{
				if (it == main_id)
				{
					double err2 = zeus._phant->peekUncertainty(dose[it], uncertainty[it], nBatch + zeus._phant->nBatch / NGPU, thresholdRegion);
					double histEstimate = hist*err2 / (targetErr2*NGPU);
					double time = rc.stop(true);
					double speed = hist / time;
					double rest = (histEstimate - hist) / speed;
					Log("GPU processed ------------------------ %3.1f%%,   speed = %d h/s\n", hist*100.0 / histEstimate, int(speed));
					Log("Time escaped = %.1f min, left time expected = %.1f min", time / 60.0, rest / 60.0);
					if (err2 < targetErr2*NGPU) b_targetErrReached = true;
					if (progressCB) progressCB(hist*100.0 / histEstimate, speed, time, rest, sqrt(err2 / NGPU)*100.0);
				}
#pragma omp barrier //make sure all work threads got the break signal
				if (b_targetErrReached) break; //all work threads will break the while loop
			}

			if (targetErr <= 0 && hist >= fHMax || b_abort) break;
		}

		//finish in this GPU thread
		gc[it].hist = hist;
		if (!b_abort && it == main_id) //if it's aborted, print nothing
		{
			Log("GPU processed ------------------------ 100%%,   speed = %d h/s\n", int(hist / rc.stop(true)));
			Log.shortMode(false);
			Log("\nWait all GPUs to finish their job...\n");
		}
		delete[] gdose;
	} //end openMP

	Log("All GPUs have finished their simulation job! Collecting dose...\n\n");

	double totHist = NGPU*gc[0].hist;
	SFloat* d = new SFloat[NVoxel];
	SFloat* u = new SFloat[NVoxel];
	//add up dose and uncertainty in different GPU card
	for (int j = 0; j < NVoxel; ++j)
	{
		d[j] = u[j] = 0;
		for (int i = 0; i < NGPU; ++i)
		{
			d[j] += dose[i][j];
			u[j] += uncertainty[i][j];
		}
	}
	for (int i = 0; i < NGPU; ++i)
	{
		delete[] dose[i];
		delete[] uncertainty[i];
	}
	Log.shortMode(false);
	double norm = 1.0889e15 * SourceHead_BeamOnTime();
	zeus._phant->addDose(d, u, NGPU * nBatch, totHist, norm, thresholdRegion);
	zeus._phant->getBinaryFile(gBF);
	delete[] d;
	delete[] u;

	if (b_abort) outname += "_abort";// make sure it wouldn't overwrite the last file
	outname += ".dose";
	string vrFormat;
	cf.getValue("ViewRay format", vrFormat);
	if (vrFormat.compare("yes") == 0) zeus._phant->output(outname.c_str(), 1);
	else zeus._phant->output(outname.c_str());
	Log("Wait for the source generating thread to finish...");
	sthread.join();
	SourceHead_Delete(); //release the source head resource safely
	Log("Source reuse times = %d", Source_Reuse_Times);
	Log("\nTime statistics for main GPU:");
	//Log("Total copy time =%.2f s", copyCounter.getStoredTime());
	Log("Total kernel time = %.2f s", kernelCounter.getStoredTime());
	Log("Total SourceHead time = %.2f s\n\n", sourceCounter.getStoredTime());
	Log("Mixed running time = %.2f minutes, total history number = %g", rc.stop(true) / 60.0, totHist);
	Log("The overall simulating speed = %d hist/sec\n\n", int(totHist / rc.stop(true)));
	Log("End log time = %s\n\n", Log.timeNow());
	Log("/##############################################################################/\n\n");
	Log.closeFile();

	if (jobFinishedCB) jobFinishedCB(b_abort, gBF);
	if (b_abort)
	{
		b_abort = false;
		if (jobFinishedCB == NULL) exit(0); //in command line execution mode, exit the program directly
	}
	b_thread_active = false;
} //end executeJob