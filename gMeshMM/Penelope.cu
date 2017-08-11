#include "Penelope.h"


/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: variables in device memory <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
//{{common constant data, need initialization from the CPU
__constant__ double Etab[NEGrid];
__constant__ double WoodcockLogSigma[NEGrid];
__constant__ double ELow, ELowL, EUp; //energy up and low bound
__constant__ double invDLE; // 1.0 / delta log(E)
__constant__ double FLCOH;
__constant__ double WB[32] =
{
	1.0e-12, 0.025e0, 0.05e0, 0.075e0, 0.1e0, 0.15e0, 0.2e0, 0.25e0,
	0.3e0, 0.35e0, 0.4e0, 0.45e0, 0.5e0, 0.55e0, 0.6e0, 0.65e0, 0.7e0,
	0.75e0, 0.8e0, 0.85e0, 0.9e0, 0.925e0, 0.95e0, 0.97e0, 0.99e0,
	0.995e0, 0.999e0, 0.9995e0, 0.9999e0, 0.99995e0, 0.99999e0, 1.0e0
};
__device__ double EPH[8000], XPH[8000][10];
__device__ int IPHF[99], IPHL[99], NPHS[99];
__device__ double ET[15000], FR[15000];
__device__ int IAL[15000], IS1[15000], IS2[15000], IFIRST[99][9], ILAST[99][9];
__device__ double XESI[6000][9], XPSI[6000][9];
__constant__ int IESIF[99], NSESI[99], IPSIF[99], NSPSI[99];
__device__ double EB[99][30];
__constant__ double BET[6];

__constant__ double DSMax;
__constant__ int NBatch;
__constant__ int NStackDepth;
__constant__ int simuPositron;
__constant__ int simuSecondary;
//}}

//{{ pointers to accelerating access in the device
__constant__ ThreadVars* tv;
__constant__ Particle* InitPars;
__constant__ Material* mat;
__constant__ Particle* stackBuff;
__constant__ goctree_node* gOctreeRoot;
//}}

__device__ int ContinueSignal;

//{{ data for phantom
__constant__ int NX, NY, NZ; //voxel number
__constant__ double DX, DY, DZ; // voxel size, unit cm
__constant__ double invDX, invDY, invDZ;
__constant__ double LX, LY, LZ; // side length Lx=DX*NX
__constant__ double xo, yo, zo;
__constant__ double MaxDensity; // only effective if there's only one material
__constant__ double Bx, By, Bz; //unit magnetic field direction
__constant__ double rf;
__constant__ int uniform;
__constant__ SFloat* doseScore;//pointer to dose counter
__constant__ SFloat* ph; //pointer to phantom
__constant__ char* matid;
__constant__ bool ThroughOctree = false; // whether to go through the octree
//}}
#ifdef STEP_DEBUG_STATUS
__device__ double DebugTack[2*STEP_NUM+20];
__device__ int DebugCur = 0;
#endif
#ifdef ACCU_DEBUG
__device__ int NIncident = 0;
#endif
/*>>>>>>>>>>>>>>>>>>>>>>>>> end: variables in device memory >>>>>>>>>>>>>>>>>>>>>>>>>>>*/

CommonConst comConst;

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: phantom method definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

__device__ __forceinline__ void deposit(Particle & p, double DE)
{

	atomicAdd(doseScore + p.iabsv, float(DE*p.weight));
#ifdef STEP_DEBUG_STATUS
	DebugTack[2*DebugCur] = p.iabsv;
	DebugTack[2* DebugCur + 1] = DE;
	++DebugCur;
#endif
}

#if(TRANSPORT_OCTREE == 1)
__device__ bool goThroughOctree(Particle& p, size_t it) //return if the particle's energy is lower than threshold after passing through the octree
{
	goctree_node* cNode = gOctreeRoot->lineIn(p.x, p.y, p.z, p.u, p.v, p.w);
	
	if (cNode != NULL)
	{
		int oldMat = -1, cMat = -1;
		double ds = INFINIT_LENGTH;
		MatInfo matIDs[MAX_OBJ_OVERLAP];
		for (int i = 0; i < MAX_OBJ_OVERLAP; ++i) matIDs->id = -1;
		
		cNode = cNode->step_octree(ds, p.x, p.y, p.z, p.u, p.v, p.w, matIDs, cMat);
		while (cNode != NULL) //still in the octree
		{
			if (cMat == oldMat) //moving in the same material, and then scatter
			{
				if (!mat[cMat].SimpleKnock(it)) return false; // energy too low
			}
			//else need to re-jump
			ds = mat[cMat].SimpleJump(it);
			cNode = cNode->step_octree(ds, p.x, p.y, p.z, p.u, p.v, p.w, matIDs, cMat);
		}
	}

	return true; //must still have energy higher than threshold
}
#endif

__device__ bool lineIn(Particle & p, size_t it)
{
#if(TRANSPORT_OCTREE == 1)
	if (ThroughOctree)
	{
		if (!goThroughOctree(p, it)) return false;
	}
#endif

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

				int ix = int((x + preNum*dmx) * invDX);
				int iy = int((y + preNum*dmy) * invDY);
				int iz = int((z + preNum*dmz) * invDZ);
				if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) return false;//it will leave the phantom without scattering
				if (ph[at(ix, iy, iz)] > 0.04) break; //stop when it get close to the target
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
	int ix = int(px * invDX);
	int iy = int(py * invDY);
	int iz = int(pz * invDZ);
	px -= ix*DX;
	py -= iy*DY;
	pz -= iz*DZ;

	//write back to global memory
	p.x = px;
	p.y = py;
	p.z = pz;
	p.ivx = ix;
	p.ivy = iy;
	p.ivz = iz;
	p.iabsv = at(ix, iy, iz);

#if	MULTI_MATERIAL == 1
	tv[it].cMatID = matid[p.iabsv]; // find the current material index
#endif
	return true; //now it's ready to transport
}

__device__ bool movePhoton(Particle & p, double ds) //return whether particle leaves the phantom
{
	double px = p.x + ds*p.u + p.ivx*DX;
	double py = p.y + ds*p.v + p.ivy*DY;
	double pz = p.z + ds*p.w + p.ivz*DZ;
	if (px < 0 || px >= LX || py < 0 || py >= LY || pz < 0 || pz >= LZ) return true;
	//calculate the voxel index
	int ix = int(px * invDX);
	int iy = int(py * invDY);
	int iz = int(pz * invDZ);

	p.x = px - DX*ix;
	p.y = py - DY*iy;
	p.z = pz - DZ*iz;
	p.ivx = ix;
	p.ivy = iy;
	p.ivz = iz;
	p.iabsv = at(ix, iy, iz);
	return false;
}

// Note ds is the distance for that material, so it should scale with the density in each voxel
__device__ int moveEP(Particle & p, double ds) //return 0 if particle remains in the phantom, 1 if particle leaves the phantom, and 2 if it enters a new material
{
#if	MULTI_MATERIAL == 1
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index, to check current mat id
#endif

	//move for electron/photon. Note that coordinates are relative to each voxel
	int pType = p.type;
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

		double q = -1;
		if (pType == positron) q = 1;
		double uwx = -q*Bx;
		double uwy = -q*By;
		double uwz = -q*Bz;
		double Rb = rf*sqrt(pE*(pE + TEs));
		double maxStep = sqrt(2 * Rb*deltax); //max allowed distance to move to ensure accuracy

		while (true)
		{
			double step = ds / voxden;
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
				vperpy = pv - vuw * uwy,
				vperpz = pw - vuw * uwz;
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
							if (ivz >= NZ) return 1;
							pz = 0;
						}
						else {
							if (ivz < 0) return 1;
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
							if (ivy >= NY) return 1;
							py = 0;
						}
						else {
							if (ivy < 0) return 1;
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
							if (ivx >= NX) return 1;
							px = 0;
						}
						else {
							if (ivx < 0) return 1;
							px = DX;
						}
					}
					else intersect = false;
				}
			}// end if (intersect)

			//still intersect after the direction check, need to update the voxel density and index
			if (intersect)
			{
				//check the matid first, if it's a new material, break
				iabsv = at(ivx, ivy, ivz);
#if	MULTI_MATERIAL == 1
				char mid = matid[iabsv];
				if (mid != tv[it].cMatID)
				{
					tv[it].cMatID = mid; // update the mat id
					//update direction and break
					pu += dvx;
					pv += dvy;
					pw += dvz;
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
					return 2;
				}
#endif
				voxden = ph[iabsv];

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
		return 0;
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
					if (ivz >= NZ) return 1;
					pz = 0;
				}
				else {
					if (ivz < 0) return 1;
					pz = DZ;
				}
			}
			else if (idex == 2)
			{
				ivy += dvox;
				if (dvox > 0)
				{
					if (ivy >= NY) return 1;
					py = 0;
				}
				else
				{
					if (ivy < 0) return 1;
					py = DY;
				}
			}
			else
			{
				ivx += dvox;
				if (dvox > 0)
				{
					if (ivx >= NX) return 1;
					px = 0;
				}
				else
				{
					if (ivx < 0) return 1;
					px = DX;
				}
			}

			iabsv = at(ivx, ivy, ivz);
#if	MULTI_MATERIAL == 1
			char mid = matid[iabsv];
			if (mid != tv[it].cMatID)
			{
				tv[it].cMatID = mid; // update the mat id
				//write back position and index
				p.x = px;
				p.y = py;
				p.z = pz;

				p.ivx = ivx;
				p.ivy = ivy;
				p.ivz = ivz;
				p.iabsv = iabsv;
				return 2;
			}
#endif
		}//end of moving forward loop

		//write back the position and index. Note that the direction doesn't change.
		p.x = px;
		p.y = py;
		p.z = pz;

		p.ivx = ivx;
		p.ivy = ivy;
		p.ivz = ivz;
		p.iabsv = iabsv;
		return 0;
	}
}

#ifdef STEP_DEBUG_STATUS
__device__ void trackStatus(Particle & p)
{
	DebugTack[7 * DebugCur] = p.u;
	DebugTack[7 * DebugCur + 1] = p.v;
	DebugTack[7 * DebugCur + 2] = p.w;
	DebugTack[7 * DebugCur + 3] = p.x + DX*p.ivx;
	DebugTack[7 * DebugCur + 4] = p.y + DY*p.ivy;
	DebugTack[7 * DebugCur + 5] = p.z + DZ*p.ivz;
	DebugTack[7 * DebugCur + 6] = p.E;
	++DebugCur;
}
#endif
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end: phantom method definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

__device__ void ParticleStack::push(Particle &par)
{
	if (cur == NStackDepth)
	{
		printf("!!!!!!!!!!!!!!!exceed the max depth of the GPU stack!!!!!!!!!!!!!!!\n\
			   		!!!!!!!!!!!!!!!will abandon this particle now!!!!!!!!!!!!!!!!!!!!!!\n\
							!!!!!!!!!!!!!!!current maximum stack depth is %d !!!!!!!!!!!!!!!!\n\
									!!!!!!!!!!!!!!!you should increase NStackDepth later!!!!!!!!!!!!!!!\n", NStackDepth);
		asm("trap;");
	}
	else
	{

		ss[cur] = par; //copy and store the particle
		++cur;
	}
}

__device__ void exitKernel(const char inf[])
{
	printf("error: %s\n\n", inf);
	asm("trap;");
}

__device__ __forceinline__ void calcG12(double b, double &g1, double &g2) //used in pair production
{
	double b2 = b*b;
	double barctan = 0;
	if (b > 1e-10) barctan = b*atan(1.0 / b);
	else barctan = b* PI / 2.0;
	double lnb = 2 * log(1 + b2);
	double A = b2*(4 - 4 * barctan - 3 * log(1.0 + 1.0 / b2));
	g1 = 7.0 / 3.0 - lnb - 6 * barctan - A;
	g2 = 11.0 / 6.0 - lnb - 3 * barctan + 0.5*A;
}

#if __CUDACC_VER_MAJOR__ < 8
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

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: Kernel definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
__global__ void initThreads(int seed)//init the random number generator and particle stack, call it before any run
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	tv[it].rng.init(seed, it); //use thread id + seed as new seed ??? reliable ???
	tv[it].pStack.init(stackBuff + it*NStackDepth);
}

#if (EXECUTE_TYPE == 1 || EXECUTE_TYPE < 0)
//{{{{{{ This block implements the traditional "one-thread-one-history" approach
__global__ void PENELOPE_Kernel()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	Particle &p = tv[it].p;
	ParticleStack &pStack = tv[it].pStack;
	int Cur = 0;
	do
	{
		//fetch the primary particles
		p = InitPars[it*NBatch + Cur];
		++Cur;
		tv[it].Mode = SOFT;
		//move the primary particle to the phantom
		if (lineIn(tv[it].p, it))
		{
			while (true) //keep jumping and knocking
			{
				char mid = 0;
#if	MULTI_MATERIAL == 1
				do 
				{
					mid = tv[it].cMatID; // update the current material index
					mat[mid].jump(it);
				} while (tv[it].evt == NEWMAT_EVT);
#else
				mat[mid].jump(it);
#endif
				if (tv[it].evt != EXIT_EVT) //still inside the object
				{
					double DE = mat[mid].knock(it); //knock and loss some energy. it may update tv[it].evt

					if (DE != 0)
					{
						deposit(p, DE); //record the energy loss in the phantom lattice
#ifdef STEP_DEBUG_STATUS
						if (DebugCur >= STEP_NUM)
						{
							ContinueSignal = -1;
							return;
						}
#endif
					}
				}

				if (tv[it].evt == EXIT_EVT) //we need to stop simulating current particle
				{
					if (pStack.empty()) break; //got no particle at all, so we stop the loop
					//still have secondary particles, so fetch the top one, and continue the while loop (starting a new jump)
					p = pStack.top();
					pStack.pop();
					deposit(p, -p.E);
#ifdef STEP_DEBUG_STATUS
					if (DebugCur >= STEP_NUM)
					{
						ContinueSignal = -1;
						return;
					}
#endif
					tv[it].Mode = SOFT;
					//tv[it].evt = INIT_EVT; // we don't need INIT_EVT in this mode
#if	MULTI_MATERIAL == 1
					tv[it].cMatID = matid[p.iabsv]; // update the current material index
#endif
					continue;
				}
				//else do nothing, continue current particle's simulation
			}
		}
		//else do nothing, continue next fetching and simulating process
	} while (Cur < NBatch);
}

__device__ void Material::jump(size_t it)
{
	Particle &p = tv[it].p;
	if (p.type == photon) GJump(it);
	else if (p.type == electron) EJump(it);
	else PJump(it);
}

__device__ double Material::knock(size_t it)
{
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	if (p.type == photon)
	{
		if (evt == GRA_EVT) return GRA(it);
		else if (evt == GCO_EVT) return GCO(it);
		else if (evt == GPH_EVT) return GPH(it);
		else if (evt == GPP_EVT) return GPP(it);
		else return 0; //DELTA_EVT
	}
	else if (p.type == electron)
	{
		if (evt == ESOFT_EVT) return ESoftKnock(it);
		else if (evt == EEL_EVT) return EEL(it);
		else if (evt == EIN_EVT) return EIN(it);
		else if (evt == EBR_EVT) return EBR(it);
		else if (evt == ESI_EVT) return ESI(it);
		else return 0; //DELTA_EVT
	}
	else // p.type == positron
	{
		if (evt == PSOFT_EVT) return PSoftKnock(it);
		else if (evt == PEL_EVT) return PEL(it);
		else if (evt == PIN_EVT) return PIN(it);
		else if (evt == EBR_EVT) return EBR(it);
		else if (evt == PSI_EVT) return PSI(it);
		else if (evt == PAN_EVT) return PAN(it);
		else return 0; //DELTA_EVT
	}
}
//}}}}}}}
#endif


#if (EXECUTE_TYPE == 0 || EXECUTE_TYPE < 0)
__global__ void start() //start the GPU simulation
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	tv[it].p = InitPars[it*NBatch];//fetch one particle from the cache
	tv[it].InitCur = 1;
	tv[it].Mode = SOFT;
	tv[it].evt = INIT_EVT;

	if (!lineIn(tv[it].p, it)) tv[it].evt = EXIT_EVT;

#ifdef ACCU_DEBUG
	++NIncident;
#endif
}

__global__ void photonJump()//kernel: init the random number generator
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].p.type == photon && tv[it].evt != EXIT_EVT) mat[mid].GJump(it);

}
__global__ void photonRA()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == GRA_EVT) mat[mid].GRA(it);
}
__global__ void photonCO()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == GCO_EVT)
	{
		double E = mat[mid].GCO(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void photonPH()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == GPH_EVT)
	{
		double E = mat[mid].GPH(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void photonPP()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == GPP_EVT)
	{
		double E = mat[mid].GPP(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}

__global__ void electronJump()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].p.type == electron && tv[it].evt != EXIT_EVT) mat[mid].EJump(it);
}
__global__ void electronSoft()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == ESOFT_EVT)
	{
		double E = mat[mid].ESoftKnock(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void electronEL()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == EEL_EVT) mat[mid].EEL(it);
}
__global__ void electronIN()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == EIN_EVT)
	{
		double E = mat[mid].EIN(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void electronBR()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == EBR_EVT)
	{
		double E = mat[mid].EBR(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void electronSI()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == ESI_EVT)
	{
		double E = mat[mid].ESI(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}

__global__ void positronJump()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].p.type == positron && tv[it].evt != EXIT_EVT) mat[mid].PJump(it);
}
__global__ void positronSoft()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == PSOFT_EVT)
	{
		double E = mat[mid].PSoftKnock(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void positronEL()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == PEL_EVT) mat[mid].PEL(it);
}
__global__ void positronIN()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == PIN_EVT)
	{
		double E = mat[mid].PIN(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void positronSI()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == PSI_EVT)
	{
		double E = mat[mid].PSI(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}
__global__ void positronAN()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	char mid = 0;
#if	MULTI_MATERIAL == 1
	mid = tv[it].cMatID; // find the current material index
#endif
	if (tv[it].evt == PAN_EVT)
	{
		double E = mat[mid].PAN(it);
		if (E != 0) deposit(tv[it].p, E);
	}
}

__global__ void fillParticle() //handle threads which are labeled with EXIT_EVT
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	Particle& p = tv[it].p;
	int& InitCur = tv[it].InitCur;
	ParticleStack &pStack = tv[it].pStack;
	bool hasSource = true;

#ifdef STEP_DEBUG_STATUS
	if (DebugCur >= STEP_NUM)
	{
		ContinueSignal = -1;
		return;
	}
#endif


	if (tv[it].evt == EXIT_EVT)
	{
		if (!simuSecondary || pStack.empty())
		{
#ifdef ACCU_DEBUG
			if (NIncident >= INCIDENT_NUM)
			{
				ContinueSignal = -2;
				return;
			}
#endif
			//try to load a photon from the source pool
			if (InitCur < NBatch)
			{
				//fetch the primary particles
				p = InitPars[it*NBatch + InitCur];
				++InitCur;
				tv[it].Mode = SOFT;
				tv[it].evt = INIT_EVT;
				//move the primary particle to the phantom
				if (!lineIn(p, it)) tv[it].evt = EXIT_EVT;

#ifdef ACCU_DEBUG
				++NIncident;
#endif
			}
			else hasSource = false;//cannot get any incident photon for this thread
		}
		else //need to simulate the secondary particles stored in the stack
		{
			p = pStack.top();
			pStack.pop();
			deposit(p, -p.E);
			tv[it].Mode = SOFT;
			tv[it].evt = INIT_EVT;
#if	MULTI_MATERIAL == 1
			tv[it].cMatID = matid[p.iabsv]; // find the current material index
#endif
		}
	}

	if (hasSource) ++ContinueSignal; //if ContinueSignal!=0, then the loop continues
}
#endif
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end: Kernel definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


/*<<<<<<<<<<<<<<<<<<<<<<<<< start: Material method definitions <<<<<<<<<<<<<<<<<<<<<<<<*/
void Material::load(int mati)
{
	M = mati;
	char name[20];
	sprintf(name, "Materials/mat%d.dat", M);

	FILE* fp = fopen(name, "rb");
	if (fp == NULL) exitApp("Cannot open Material file!");
	fread(&M, sizeof(int), 1, fp);
	fread(&Eabs[electron], sizeof(double), 1, fp);
	fread(&Eabs[photon], sizeof(double), 1, fp);
	fread(&Eabs[positron], sizeof(double), 1, fp);
	fread(&C1, sizeof(double), 1, fp);
	fread(&C2, sizeof(double), 1, fp);
	fread(&Wcc, sizeof(double), 1, fp);
	fread(&Wcr, sizeof(double), 1, fp);
	fread(&Ecutr, sizeof(double), 1, fp);
	fread(&Zt, sizeof(double), 1, fp);
	fread(&At, sizeof(double), 1, fp);
	fread(&massDensity, sizeof(double), 1, fp);
	fread(&numDensity, sizeof(double), 1, fp);
	fread(&elemNum, sizeof(int), 1, fp);
	fread(weight, sizeof(double), 30, fp);
	fread(IZ, sizeof(int), 30, fp);
	fread(SGRA, sizeof(double), 200, fp);
	fread(SGCO, sizeof(double), 200, fp);
	fread(SGPH, sizeof(double), 200, fp);
	fread(SGPP, sizeof(double), 200, fp);
	fread(X2COH, sizeof(double), 241, fp);
	fread(PDCOH, sizeof(double), 241, fp);
	fread(FCO, sizeof(double), 64, fp);
	fread(UICO, sizeof(double), 64, fp);
	fread(FJ0, sizeof(double), 64, fp);
	fread(KZCO, sizeof(int), 64, fp);
	fread(KSCO, sizeof(int), 64, fp);
	fread(&NOSCCO, sizeof(int), 1, fp);
	fread(&ZEQPP, sizeof(double), 1, fp);
	fread(&Fpp, sizeof(double), 1, fp);
	fread(&Bpp, sizeof(double), 1, fp);
	fread(SEHEL, sizeof(double), 200, fp);
	fread(SEHIN, sizeof(double), 200, fp);
	fread(SEISI, sizeof(double), 200, fp);
	fread(SEHBR, sizeof(double), 200, fp);
	fread(SETOT, sizeof(double), 200, fp);
	fread(CSTPE, sizeof(double), 200, fp);
	fread(RSTPE, sizeof(double), 200, fp);
	fread(DEL, sizeof(double), 200, fp);
	fread(W1E, sizeof(double), 200, fp);
	fread(W2E, sizeof(double), 200, fp);
	fread(T1E, sizeof(double), 200, fp);
	fread(T2E, sizeof(double), 200, fp);
	fread(RNDCE, sizeof(double), 200, fp);
	fread(AE, sizeof(double), 200, fp);
	fread(BE, sizeof(double), 200, fp);

	fread(SPHEL, sizeof(double), 200, fp);
	fread(SPHIN, sizeof(double), 200, fp);
	fread(SPISI, sizeof(double), 200, fp);
	fread(SPHBR, sizeof(double), 200, fp);
	fread(SPAN, sizeof(double), 200, fp);
	fread(SPTOT, sizeof(double), 200, fp);
	fread(CSTPP, sizeof(double), 200, fp);
	fread(RSTPP, sizeof(double), 200, fp);
	fread(W1P, sizeof(double), 200, fp);
	fread(W2P, sizeof(double), 200, fp);
	fread(T1P, sizeof(double), 200, fp);
	fread(T2P, sizeof(double), 200, fp);
	fread(RNDCP, sizeof(double), 200, fp);
	fread(AP, sizeof(double), 200, fp);
	fread(BP, sizeof(double), 200, fp);

	fread(&EELMAX, sizeof(double), 1, fp);
	fread(&PELMAX, sizeof(double), 1, fp);

	fread(XSE, sizeof(double), 128 * 200, fp);
	fread(PSE, sizeof(double), 128 * 200, fp);
	fread(ASE, sizeof(double), 128 * 200, fp);
	fread(BSE, sizeof(double), 128 * 200, fp);
	fread(ITLE, sizeof(int), 128 * 200, fp);
	fread(ITUE, sizeof(int), 128 * 200, fp);

	fread(XSP, sizeof(double), 128 * 200, fp);
	fread(PSP, sizeof(double), 128 * 200, fp);
	fread(ASP, sizeof(double), 128 * 200, fp);
	fread(BSP, sizeof(double), 128 * 200, fp);
	fread(ITLP, sizeof(int), 128 * 200, fp);
	fread(ITUP, sizeof(int), 128 * 200, fp);

	fread(&meanExcitE, sizeof(double), 1, fp);
	fread(&OP2, sizeof(double), 1, fp);

	fread(F, sizeof(double), 64, fp);
	fread(Ui, sizeof(double), 64, fp);
	fread(WRI, sizeof(double), 64, fp);
	fread(KZ, sizeof(int), 64, fp);
	fread(KS, sizeof(int), 64, fp);

	fread(&NOSC, sizeof(int), 1, fp);

	for (int i = 0; i < 64; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&EINAC[j][i], sizeof(double), 1, fp);
	}
	for (int i = 0; i < 64; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&PINAC[j][i], sizeof(double), 1, fp);
	}

	fread(PBcut, sizeof(double), 200, fp);
	fread(WBcut, sizeof(double), 200, fp);

	for (int i = 0; i < 32; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&PDFB[j][i], sizeof(double), 1, fp);
	}
	for (int i = 0; i < 32; ++i)
	{
		for (int j = 0; j < 200; ++j) fread(&PACB[j][i], sizeof(double), 1, fp);
	}

	for (int i = 0; i < 6; ++i)
	{
		for (int j = 0; j < 21; ++j)
		{
			for (int k = 0; k < 4; ++k) fread(&BP1[i][j][k], sizeof(double), 1, fp);
		}
	}
	for (int i = 0; i < 6; ++i)
	{
		for (int j = 0; j < 21; ++j)
		{
			for (int k = 0; k < 4; ++k) fread(&BP2[i][j][k], sizeof(double), 1, fp);
		}
	}

	fclose(fp);
}

double Material::totalSigmaG(int IE)
{
	double p1 = exp(SGRA[IE]); //N*sigma for Rayleigh scattering
	double p2 = exp(SGCO[IE]); //N*sigma for Compton scattering
	double p3 = SGPH[IE];
	double p4 = 0;
	if (comConst.Etab_h[IE] > TEs) p4 = exp(SGPP[IE]);
	return p1 + p2 + p3 + p4; //sum of N*sigma for 4 mechanism, including woodcock tracing
}

#if (TRANSPORT_OCTREE == 1)
__device__ double Material::SimpleJump(size_t it)
{
	/*this function move this particle and label what interaction will happen next */
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	double pE = p.E; //cache the variable 

	double XE = (log(pE) - ELowL) * invDLE;
	int IE = (int)XE;
	double XIE = XE - IE;
	//write back the energy index
	tv[it].IE = IE;
	tv[it].XIE = XIE;

	double p1 = exp((1 - XIE)*SGRA[IE] + XIE*SGRA[IE + 1]); //N*sigma for Rayleigh scattering
	double p2 = exp((1 - XIE)*SGCO[IE] + XIE*SGCO[IE + 1]); //N*sigma for Compton scattering
	double p3 = SGPH[IE];
	double p4 = 0;
	if (pE > TEs) p4 = exp((1 - XIE)*SGPP[IE] + XIE*SGPP[IE + 1]);
	double st = p1 + p2 + p3 + p4; //sum of N*sigma for 4 mechanism, including woodcock tracing
	double ds = -log(rng()) / st;//sample the jump distance

	//determine what interaction happens next
	double sst = st*rng();
	if (sst < p1) evt = GRA_EVT; //Rayleigh scattering
	else if (sst < p1 + p2)  evt = GCO_EVT; //Compton scattering
	else evt = EXIT_EVT; //Photoelectric absorption

	return ds;
}
__device__ bool Material::SimpleKnock(size_t it) //return if we should continue simulating this photon
{
	KonckEvent &evt = tv[it].evt;
	if (evt == GCO_EVT) return SimpleCompton(it);
	else if (evt == EXIT_EVT) return false;
	else if (evt == GRA_EVT)
	{
		GRA(it);
		return true;
	}
	return false;
}
__device__ bool Material::SimpleCompton(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	//cache from the slow device memory
	double pE = p.E;

	const double D2 = 1.4142135623731, D1 = 1.0 / D2, D12 = 0.5;

	double Ek = pE / Es;
	double Ek2 = Ek + Ek + 1.0; //2k+1
	double Ek3 = Ek*Ek; //k^2
	double Ek1 = Ek3 - Ek2 - 1.0; //k^2-2k-2
	double tmin = 1 / Ek2; // 1/(2k+1)
	double tmin2 = tmin*tmin; // 1/(2k+1)^2
	double a1 = log(Ek2);
	double a2 = a1 + 2.0*Ek*(1.0 + Ek)*tmin2; //a1+a2 in the manual

	double Ep = 0; //energy of the scattered photon
	double cosTheta = 0; //Direction of the scattered photon
	int iSh = -1; //index of the shell ionized, starting from 0

	double tau = 0;
	double Tcos = 0;
	double S = 0; // S(E,theta)

	if (pE > 5e6) //No Doppler broadening for E greater than 5 MeV
	{
		do
		{
			//sample tau
			do
			{
				if (rng()*a2 < a1) //i=1
					tau = pow(tmin, rng());
				else
					tau = sqrt(1.0 + rng()*(tmin2 - 1.0));
				Tcos = (1.0 + tau*(Ek1 + tau*(Ek2 + tau*Ek3))) / (Ek3*tau*(1.0 + tau*tau));
			} while (rng()>Tcos);
			Ep = tau*pE;
			cosTheta = 1.0 - (1.0 - tau) / (Ek*tau);

			//sample which shell excited
			Tcos = Zt*rng();
			S = 0;
			iSh = -1; //for search judgment
			for (int i = 0; i < NOSCCO; ++i) //for each shell, starting from 0;
			{
				S += FCO[i];
				if (S > Tcos) { iSh = i; break; }
			}
			if (iSh == -1) iSh = NOSCCO - 1; //avoid some critical problem;
		} while (Ep > pE - UICO[iSh]);
	}
	else //should take Doppler broadening into account
	{
		//********sample Direction of the scattered photon
		//calculate S(E,theta=pi), i.e. S0
		double S0 = 0;
		double ni[64];
		double pac[64]; //accumulation funcion of fi*theta(E-Ui)*ni, used to sample which shell being excited
		for (int i = 0; i <NOSCCO; ++i) //for each shell
		{
			if (pE > UICO[i]) //like theta(E-Ui)
			{
				double aux = pE*(pE - UICO[i])*2.0; //2E*(E-Ui), auxiliary variable
				double PzJi0 = FJ0[i] * (aux - Es*UICO[i]) / (Es*sqrt(aux + aux + UICO[i] * UICO[i]));
				if (PzJi0 > 0) ni[i] = 1 - 0.5*exp(D12 - (D1 + D2*PzJi0)*(D1 + D2*PzJi0));
				else ni[i] = 0.5 * exp(D12 - (D1 - D2*PzJi0)*(D1 - D2*PzJi0));
				S0 += FCO[i] * ni[i];
			}
		}

		double cdt1 = 0;
		do
		{
			// Sampling tau
			if (rng()*a2 < a1) tau = pow(tmin, rng());
			else tau = sqrt(1.0 + rng()*(tmin2 - 1.0));
			cdt1 = (1.0 - tau) / (Ek*tau); //1-cos(theta)
			//calculate S(E, theta)
			S = 0;
			for (int i = 0; i <NOSCCO; ++i) //for each shell
			{
				if (pE > UICO[i]) //like theta(E-Ui)
				{
					double aux = pE*(pE - UICO[i])*cdt1; //2E*(E-Ui), auxiliary variable
					double PzJi0 = FJ0[i] * (aux - Es*UICO[i]) / (Es*sqrt(aux + aux + UICO[i] * UICO[i]));
					if (PzJi0 > 0) ni[i] = 1 - 0.5*exp(D12 - (D1 + D2*PzJi0)*(D1 + D2*PzJi0));
					else ni[i] = 0.5 * exp(D12 - (D1 - D2*PzJi0)*(D1 - D2*PzJi0));
					S += FCO[i] * ni[i];
					pac[i] = S;
				}
				else pac[i] = S - 1e-6; //avoid sampling the E<Ui shell
			}
			Tcos = S*(1.0 + tau*(Ek1 + tau*(Ek2 + tau*Ek3))) / (Ek3*tau*(1.0 + tau*tau));
		} while (rng()*S0 > Tcos);
		cosTheta = 1.0 - cdt1; //got the Direction of the scattered photon

		//***************sampling pz;
		double pz = 0;
		double Fmax = 0;
		double Fpz = 0;
		double RS = 0;

		do
		{
			do
			{
				//sample which shell is excited
				iSh = -1; //the index of excited shell
				RS = S*rng(); //random sampled part of S
				for (int i = 0; i < NOSCCO; ++i)
				{
					if (RS < pac[i]) { iSh = i; break; }
				}
				if (iSh == -1) iSh = NOSCCO - 1; //avoid some critical problem;

				//sample pz using formula (2.56)
				double A = rng()*ni[iSh];
				if (A < 0.5) pz = (D1 - sqrt(D12 - log(A + A))) / (D2*FJ0[iSh]);
				else pz = (sqrt(D12 - log(2.0 - A - A)) - D1) / (D2*FJ0[iSh]);
			} while (pz < -1.0);

			//calculate Fmax
			double xqc = 1.0 + tau*(tau - 2.0*cosTheta);
			double AF = sqrt(xqc)*(1.0 + tau*(tau - cosTheta) / xqc);
			if (AF > 0.0) Fmax = 1.0 + AF*0.2;
			else Fmax = 1.0 - AF*0.2;
			Fpz = 1.0 + AF*max(min(pz, 0.2), -0.2);
		} while (Fmax*rng() > Fpz); //rejection selecting pz

		//******************calculate Energy of the scattered photon Ep
		double tz = pz*pz;
		double b1 = 1.0 - tz*tau*tau;
		double b2 = 1.0 - tz*tau*cosTheta;
		int signpz = (pz > 0) ? 1 : -1;
		Ep = pE*(tau / b1)*(b2 + signpz*sqrt(fabs(b2*b2 - b1*(1.0 - tz))));
	}
	//by here, we have got Ep, cosTheta, and iSh(start from 0)

	if (Ep < Eabs[photon]) return false; // can stop the simulation of this photon
	p.changeByCos(cosTheta, 2 * PI*rng());
	p.E = Ep;
	return true;
}
#endif

__device__ void Material::GJump(size_t it)
{
	/*this function move this particle and label what interaction will happen next */
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	double pE = p.E; //cache the variable 
	
	double XE = (log(pE) - ELowL) * invDLE;
	int IE = (int)XE;
	double XIE = XE - IE;
	//write back the energy index
	tv[it].IE = IE;
	tv[it].XIE = XIE;

#if	MULTI_MATERIAL == 1
	double sm = exp((1 - XIE)*WoodcockLogSigma[IE] + XIE*WoodcockLogSigma[IE + 1]);
	double ds = -log(rng()) / sm;//sample the jump distance
#else
	double p1 = exp((1 - XIE)*SGRA[IE] + XIE*SGRA[IE + 1]); //N*sigma for Rayleigh scattering
	double p2 = exp((1 - XIE)*SGCO[IE] + XIE*SGCO[IE + 1]); //N*sigma for Compton scattering
	double p3 = SGPH[IE];
	double p4 = 0;
	if (pE > TEs) p4 = exp((1 - XIE)*SGPP[IE] + XIE*SGPP[IE + 1]);
	double st = (p1 + p2 + p3 + p4); //sum of N*sigma for 4 mechanism
	double ds = -log(rng()) / (st*MaxDensity); //sample the jump distance
#endif

	if (!movePhoton(p, ds))//still in the phantom
	{
#ifdef STEP_DEBUG
		//trackStatus(p);
#endif
		SFloat den = ph[p.iabsv];
		if (den == 0) den = 1.29e-3f;//set the density to be same as air
		//determine what interaction happens next
#if	MULTI_MATERIAL == 1
		//check if it moves into a new material
		char mid = matid[p.iabsv];
		if (mid != tv[it].cMatID)  tv[it].cMatID = matid[p.iabsv]; // update current material id, and cross sections

		double p1 = exp((1 - XIE)*mat[mid].SGRA[IE] + XIE*mat[mid].SGRA[IE + 1]); //N*sigma for Rayleigh scattering
		double p2 = exp((1 - XIE)*mat[mid].SGCO[IE] + XIE*mat[mid].SGCO[IE + 1]); //N*sigma for Compton scattering
		double p3 = mat[mid].SGPH[IE];
		double p4 = 0;
		if (pE > TEs) p4 = exp((1 - XIE)*mat[mid].SGPP[IE] + XIE*mat[mid].SGPP[IE + 1]);
		double st = p1 + p2 + p3 + p4; //sum of N*sigma for 4 mechanism
		sm /= den; // to make it consistent to following code
#else
		double sm = st *(MaxDensity / den);
#endif
		double sst = sm*rng();
		if (sst > st) evt = DELTA_EVT;
		else if (sst < p1) evt = GRA_EVT; //Rayleigh scattering
		else if (sst < p1 + p2)  evt = GCO_EVT; //Compton scattering
		else if (sst < st - p4) evt = GPH_EVT; //Photoelectric absorption
		else evt = GPP_EVT; //Electron-positron pair production
	}
	else evt = EXIT_EVT;
}

__device__ void Material::EJump(size_t it)
{
	/*this function move this particle and label what interaction will happen next */
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;

	GRNG &rng = tv[it].rng;
	int &IE = tv[it].IE; //left index of interval
	double &XIE = tv[it].XIE; //portion in that interval
	KnockMode &Mode = tv[it].Mode;
	double &st = tv[it].st;
	double &DST = tv[it].DST;
	double &DS1 = tv[it].DS1;
	double &W1 = tv[it].W1;
	double &W2 = tv[it].W2;
	double &T1 = tv[it].T1;
	double &T2 = tv[it].T2;
	bool& KDelta = tv[it].KDelta;

	double XE = (log(p.E) - ELowL) * invDLE;
	IE = (int)XE;
	XIE = XE - IE;

	double ds = 0;
	double p2 = exp((1 - XIE)*SEHEL[IE] + XIE*SEHEL[IE + 1]);
	double p3 = exp((1 - XIE)*SEHIN[IE] + XIE*SEHIN[IE + 1]);
	double p4 = exp((1 - XIE)*SEHBR[IE] + XIE*SEHBR[IE + 1]);
	double p5 = exp((1 - XIE)*SEISI[IE] + XIE*SEISI[IE + 1]);
		
	if (Mode == HARD) ds = DS1;
	else //this would be soft interaction
	{
		if (W1E[IE + 1] > -78.3)
		{
			W1 = exp((1 - XIE)*W1E[IE] + XIE*W1E[IE + 1]);
			W2 = exp((1 - XIE)*W2E[IE] + XIE*W2E[IE + 1]);
		}
		else
		{
			W1 = 0;
			W2 = 0;
		}
		if (T1E[IE + 1] > -78.3)
		{
			T1 = exp((1 - XIE)*T1E[IE] + XIE*T1E[IE + 1]);
			T2 = exp((1 - XIE)*T2E[IE] + XIE*T2E[IE + 1]);
		}
		else
		{
			T1 = 0;
			T2 = 0;
		}
		st = p2 + p3 + p4 + p5;
		//get the real soft st
		double DSMaxp = DSMax;
		if (W1 > 1e-20)
		{
			//get the real maximum of DS
			double DSmc = 4 / st;
			if (DSmc < DSMaxp) DSMaxp = DSmc;
			else if (DSMaxp < 1e-8) DSMaxp = DSmc;
			//The value of DSMAXP is randomized to eliminate dose artifacts at the end of the first step
			DSMaxp = (0.5 + rng()*0.5)*DSMaxp;
			//Upper bound for the interaction probability along the step (including artificial energy straggling).
			double EDE0 = W1*DSMaxp;
			double VDE0 = W2*DSMaxp;
			double KW1E = (W1E[IE + 1] - W1E[IE]) / (Etab[IE + 1] - Etab[IE]);
			double KW2E = (W2E[IE + 1] - W2E[IE]) / (Etab[IE + 1] - Etab[IE]);
			double EDEM = EDE0*max(1 - 0.5*KW1E*EDE0, 0.75);
			double VDEM = VDE0*max(1 - (KW1E + 0.5*KW2E)*EDE0, 0.75);
			double W21 = VDEM / EDEM;
			double ELower;
			if (EDEM>9.0*W21) ELower = max(p.E - (EDEM + 3.0*sqrt(VDEM)), Eabs[electron]);
			else if (EDEM > 3.0*W21) ELower = max(p.E - (EDEM + sqrt(3 * VDEM)), Eabs[electron]);
			else ELower = max(p.E - 1.5*(EDEM + W21), Eabs[electron]);
			double XE1 = (log(ELower) - ELowL) * invDLE;
			int IE1 = (int)XE1;
			double XIE1 = XE1 - IE1;
			st = max(st, exp(SETOT[IE1] + (SETOT[IE1 + 1] - SETOT[IE1])*XIE1));
		}

		DST = -log(rng()) / st; //sample the distance

		//modify jump distance with delta interaction
		if (DST < DSMaxp) KDelta = false;
		else
		{
			DST = DSMaxp;
			KDelta = true;
		}

		if (W1 < 1e-20 && T1 < 1e-20) //hard interaction
		{
			Mode = HARD;
			DS1 = 0;
			ds = DST;
		}
		else //soft interaction, random hinge
		{
			ds = DST*rng();
			DS1 = DST - ds;
		}
	}

	if (1==moveEP(p, ds)) evt = EXIT_EVT;
#if	MULTI_MATERIAL == 1
	else if (2 == moveEP(p, ds)) // Moved into a new material so prepare new jump
	{
		Mode = SOFT;
		evt = NEWMAT_EVT; //repeat EJump
	}
#endif
	else // still inside the phantom
	{
#ifdef STEP_DEBUG
		//trackStatus(p);
#endif
		if (Mode == SOFT)
		{
			evt = ESOFT_EVT;
			Mode = HARD;
		}
		else // a hard event will happen
		{
			Mode = SOFT;
			if (KDelta)
			{
				evt = DELTA_EVT;//the maximum allowed step length is exceeded
				return;
			}
			double stnow = p2 + p3 + p4 + p5;
			double sst = max(stnow, st)*rng();
				
			if (sst < p2) evt = EEL_EVT;
			else if (sst < p2 + p3)  evt = EIN_EVT;
			else if (sst < p2 + p3 + p4)  evt = EBR_EVT;
			else if (sst < stnow) evt = ESI_EVT;
			else evt = DELTA_EVT;
		}
	}
}

__device__ void Material::PJump(size_t it)
{
	/*this function move this particle and label what interaction will happen next */
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	GRNG &rng = tv[it].rng;
	int &IE = tv[it].IE; //left index of interval
	double &XIE = tv[it].XIE; //portion in that interval
	KnockMode &Mode = tv[it].Mode;
	double &st = tv[it].st;
	double &DST = tv[it].DST;
	double &DS1 = tv[it].DS1;
	double &W1 = tv[it].W1;
	double &W2 = tv[it].W2;
	double &T1 = tv[it].T1;
	double &T2 = tv[it].T2;
	bool& KDelta = tv[it].KDelta;
		

	double XE = (log(p.E) - ELowL) * invDLE;
	IE = (int)XE;
	XIE = XE - IE;

	double ds = 0;
	double p2 = exp((1 - XIE)*SPHEL[IE] + XIE*SPHEL[IE + 1]);
	double p3 = exp((1 - XIE)*SPHIN[IE] + XIE*SPHIN[IE + 1]);
	double p4 = exp((1 - XIE)*SPHBR[IE] + XIE*SPHBR[IE + 1]);
	double p5 = exp((1 - XIE)*SPISI[IE] + XIE*SPISI[IE + 1]);
	double p6 = exp((1 - XIE)*SPAN[IE] + XIE*SPAN[IE + 1]);
		
	if (Mode == HARD) ds = DS1;
	else //this would be soft interaction
	{
		if (W1P[IE + 1] > -78.3)
		{
			W1 = exp((1 - XIE)*W1P[IE] + XIE*W1P[IE + 1]);
			W2 = exp((1 - XIE)*W2P[IE] + XIE*W2P[IE + 1]);
		}
		else
		{
			W1 = 0;
			W2 = 0;
		}
		if (T1P[IE + 1] > -78.3)
		{
			T1 = exp((1 - XIE)*T1P[IE] + XIE*T1P[IE + 1]);
			T2 = exp((1 - XIE)*T2P[IE] + XIE*T2P[IE + 1]);
		}
		else
		{
			T1 = 0;
			T2 = 0;
		}

		st = p2 + p3 + p4 + p5 + p6;
		//get the real soft st
		double DSMaxp = DSMax;
		if (W1 > 1e-20)
		{
			//get the real maximum of DS
			double DSmc = 4 / st;
			if (DSmc < DSMaxp) DSMaxp = DSmc;
			else if (DSMaxp < 1e-8) DSMaxp = DSmc;
			//The value of DSMAXP is randomized to eliminate dose artifacts at the end of the first step
			DSMaxp = (0.5 + rng()*0.5)*DSMaxp;
			//Upper bound for the interaction probability along the step (including artificial energy straggling).
			double EDE0 = W1*DSMaxp;
			double VDE0 = W2*DSMaxp;
			double KW1P = (W1P[IE + 1] - W1P[IE]) / (Etab[IE + 1] - Etab[IE]);
			double KW2P = (W2P[IE + 1] - W2P[IE]) / (Etab[IE + 1] - Etab[IE]);
			double EDEM = EDE0*max(1 - 0.5*KW1P*EDE0, 0.75);
			double VDEM = VDE0*max(1 - (KW1P + 0.5*KW2P)*EDE0, 0.75);
			double W21 = VDEM / EDEM;
			double ELower;
			if (EDEM>9.0*W21) ELower = max(p.E - (EDEM + 3.0*sqrt(VDEM)), Eabs[positron]);
			else if (EDEM > 3.0*W21) ELower = max(p.E - (EDEM + sqrt(3 * VDEM)), Eabs[positron]);
			else ELower = max(p.E - 1.5*(EDEM + W21), Eabs[positron]);
			double XE1 = (log(ELower) - ELowL) * invDLE;
			int IE1 = (int)XE1;
			double XIE1 = XE1 - IE1;
			st = max(st, exp(SPTOT[IE1] + (SPTOT[IE1 + 1] - SPTOT[IE1])*XIE1));
		}

		DST = -log(rng()) / st; //sample the distance

		//modify jump distance with delta interaction
		if (DST < DSMaxp) KDelta = false;
		else
		{
			DST = DSMaxp;
			KDelta = true;
		}

		if (W1 < 1e-20 && T1 < 1e-20) //hard interaction
		{
			Mode = HARD;
			DS1 = 0;
			ds = DST;
		}
		else //soft interaction, random hinge
		{
			ds = DST*rng();
			DS1 = DST - ds;
		}
	}

	if (1 == moveEP(p, ds)) evt = EXIT_EVT;
#if	MULTI_MATERIAL == 1
	else if (2 == moveEP(p, ds)) // moved into a new material so prepare a new jump
	{
		Mode = SOFT;
		evt = NEWMAT_EVT; //repeat PJump
	}
#endif
	else // still inside the phantom
	{
#ifdef STEP_DEBUG
		//trackStatus(p);
#endif
		if (Mode == SOFT)
		{
			evt = PSOFT_EVT;
			Mode = HARD;
		}
		else // a hard event will happen
		{
			Mode = SOFT;
			if (KDelta)
			{
				evt = DELTA_EVT;//the maximum allowed step length is exceeded
				return;
			}
			double stnow = p2 + p3 + p4 + p5 + p6;
			double sst = max(stnow, st)*rng();
				
			if (sst < p2) evt = PEL_EVT;
			else if (sst < p2 + p3)  evt = PIN_EVT;
			else if (sst < p2 + p3 + p4)  evt = EBR_EVT;
			else if (sst < p2 + p3 + p4 + p5) evt = PSI_EVT;
			else if(sst <stnow) evt = PAN_EVT;
			else evt = DELTA_EVT;
		}
	}
}

__device__ double Material::ESoftKnock(size_t it)
{
	GRNG &rng = tv[it].rng;
	KonckEvent &evt = tv[it].evt;
	Particle &p = tv[it].p;
	int &IE = tv[it].IE; //left index of interval
	double &XIE = tv[it].XIE; //portion in that interval
	double &DST = tv[it].DST;
	double &W1 = tv[it].W1;
	double &W2 = tv[it].W2;
	double &T1 = tv[it].T1;
	double &T2 = tv[it].T2;

	double DE = 0.0;

	if (W1 > 1e-20)
	{
		double EDE0 = W1*DST;
		double VDE0 = W2*DST;
		double K1 = (W1E[IE + 1] - W1E[IE]) / (Etab[IE + 1] - Etab[IE]);
		double K2 = (W2E[IE + 1] - W2E[IE]) / (Etab[IE + 1] - Etab[IE]);
		double EDE = EDE0*max(1 - 0.5*K1*EDE0, 0.75); //average of DE
		double VDE = VDE0*max(1 - (K1 + 0.5*K2)*EDE0, 0.75); //variance of DE
		//Generation of random values DE distributed from 0 to infinity
		//with mean EDE and variance VDE
		double sigma = sqrt(VDE);
		if (sigma < 0.333333333*EDE) //sampling from Gaussian distribution (Box-Muller method) to accelerate the calculation
		{
			double r = 0;
			do
			{
				r = sqrt(-2.0*log(rng()))*sin(2 * PI*rng())* TRUNC;
			} while (fabs(r) > 3);
			DE = EDE + r*sigma;
		}
		else
		{
			double r = rng();
			double EDE2 = EDE*EDE;
			double VDE3 = 3 * VDE;
			if (EDE2 < VDE3)
			{
				double pnull = (VDE3 - EDE2) / (VDE3 + 3.0*EDE2);
				if (r < pnull) DE = 0;
				else DE = 1.5*(EDE + VDE / EDE)*(r - pnull) / (1.0 - pnull);
			}
			else DE = EDE + (2.0*r - 1.0)*sqrt(VDE3);
		}
	}

	p.E -= DE;//lost certain energy
	if (p.E < Eabs[electron])
	{
		DE = p.E + DE; //lost all energy
		evt = EXIT_EVT;
		return DE;
	}
	if (T1 <= 1e-20) return DE;

	// Bielajew's randomly alternate hinge
	if (rng() > 0.5&&DE > 1.0e-3)
	{
		//calculate the energy index
		double XE = (log(p.E) - ELowL) * invDLE;
		IE = (int)XE;
		XIE = XE - IE;
		if (T1E[IE + 1] > -78.3)
		{
			T1 = exp(T1E[IE] + XIE*(T1E[IE + 1] - T1E[IE]));
			T2 = exp(T2E[IE] + XIE*(T2E[IE + 1] - T2E[IE]));
		}
		else T1 = T2 = 0;
		if (T1 < 1e-35) return DE;
	}
	// 1st and 2nd moments of the angular distribution
	double emu1 = 0.5*(1.0 - exp(-DST*T1)); //average mu 1
	double emu2 = emu1 - (1.0 - exp(-DST*T2)) / 6.0; //average mu 2
	//double numerator = 2 * emu1 - 3 * emu2;
	double denominator = 1 - 2 * emu1;
	double b = (2 * emu1 - 3 * emu2) / denominator;
	double a = denominator + b;
	double r = rng();
	double cosTheta;
	if (r < a) cosTheta = 1 - 2 * b*(r / a);
	else cosTheta = 1 - 2 * (b + (1 - b)*(r - a) / (1 - a));
	p.changeByCos(cosTheta, 2 * PI*rng());
	return DE;
}

__device__ double Material::PSoftKnock(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	int &IE = tv[it].IE; //left index of interval
	double &XIE = tv[it].XIE; //portion in that interval
	double &DST = tv[it].DST;
	double &W1 = tv[it].W1;
	double &W2 = tv[it].W2;
	double &T1 = tv[it].T1;
	double &T2 = tv[it].T2;

	double DE = 0;
	if (W1 > 1e-20)
	{
		double EDE0 = W1*DST;
		double VDE0 = W2*DST;
		double K1 = (W1P[IE + 1] - W1P[IE]) / (Etab[IE + 1] - Etab[IE]);
		double K2 = (W2P[IE + 1] - W2P[IE]) / (Etab[IE + 1] - Etab[IE]);
		double EDE = EDE0*max(1 - 0.5*K1*EDE0, 0.75); //average of DE
		double VDE = VDE0*max(1 - (K1 + 0.5*K2)*EDE0, 0.75); //variance of DE
		//Generation of random values DE distributed from 0 to infinity
		//with mean EDE and variance VDE
		double sigma = sqrt(VDE);
		if (sigma < 0.333333333*EDE) //sampling from Gaussian distribution (Box-Muller method) to accelerate the calculation
		{
			double r = 0;
			do
			{
				r = sqrt(-2.0*log(rng()))*sin(2 * PI*rng())* TRUNC;
			} while (fabs(r) > 3);
			DE = EDE + r*sigma;
		}
		else
		{
			double r = rng();
			double EDE2 = EDE*EDE;
			double VDE3 = 3 * VDE;
			if (EDE2 < VDE3)
			{
				double pnull = (VDE3 - EDE2) / (VDE3 + 3.0*EDE2);
				if (r < pnull) DE = 0;
				else DE = 1.5*(EDE + VDE / EDE)*(r - pnull) / (1.0 - pnull);
			}
			else DE = EDE + (2.0*r - 1.0)*sqrt(VDE3);
		}
	}

	p.E -= DE;
	if (p.E < Eabs[positron])
	{
		DE += (p.E+TEs); //lost the kinetic energy of the positron
		PANR(); //the latent energy was taken away by the two photons, so TEs is not added here
		evt = EXIT_EVT;
		return DE;
	}
	if (T1 <= 1e-20) return DE;

	// Bielajew's randomly alternate hinge
	if (rng() > 0.5 && DE > 1.0e-3)
	{
		double XE = (log(p.E) - ELowL) * invDLE;
		IE = (int)XE;
		XIE = XE - IE;
		if (T1E[IE + 1] > -78.3)
		{
			T1 = exp((1 - XIE)*T1P[IE] + XIE*T1P[IE + 1]);
			T2 = exp((1 - XIE)*T2P[IE] + XIE*T2P[IE + 1]);
		}
		else T1 = T2 = 0;
		if (T1 < 1e-35) return DE;
	}
	// 1st and 2nd moments of the angular distribution
	double emu1 = 0.5*(1.0 - exp(-DST*T1)); //average mu 1
	double emu2 = emu1 - (1.0 - exp(-DST*T2)) / 6.0; //average mu 2
	//double numerator = 2 * emu1 - 3 * emu2;
	double denominator = 1 - 2 * emu1;
	double b = (2 * emu1 - 3 * emu2) / denominator;
	double a = denominator + b;
	double r = rng();
	double cosTheta;
	if (r < a) cosTheta = 1 - 2 * b*(r / a);
	else cosTheta = 1 - 2 * (b + (1 - b)*(r - a) / (1 - a));
	p.changeByCos(cosTheta, 2 * PI*rng());
	return DE;
}

__device__ double Material::GRA(size_t it) //Rayleigh scattering, depends on X2COH[241], PDCOH[241], FLCOH, see eq(2.18) in PENELOPE 2005
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;

	//no energy loss, only Direction change
	double cosTheta = 0;
	double x2max = 2.0*log(41.2148*p.E / Es); //define log(xmax^2), xmax=41.2148*E/Es
	int jm = 0, j = 0, jt = 0;
	if (x2max < X2COH[1]) jm = 0;
	else if (x2max > X2COH[239]) jm = 239;
	else jm = (int)((x2max - X2COH[0]) / FLCOH); //j must belong to [0,239]
	//linear interpolation of PImax, actually they're all log values
	double PImax = PDCOH[jm] + (PDCOH[jm + 1] - PDCOH[jm]) / (X2COH[jm + 1] - X2COH[jm]) * (x2max - X2COH[jm]);
	do
	{
		double ru = PImax + log(rng());
		//binary search for the corresponding x^2 of ru
		j = 0;
		int jr = jm + 1;
		do
		{
			jt = (j + jr) / 2;
			if (ru > PDCOH[jt]) j = jt;
			else jr = jt;
		} while (jr - j > 1);

		double dPI = PDCOH[j + 1] - PDCOH[j];
		double x2log = 0;
		if (dPI > 1e-12) x2log = X2COH[j] + (X2COH[j + 1] - X2COH[j]) / dPI *(ru - PDCOH[j]) - x2max;
		else x2log = X2COH[j] - x2max;
		cosTheta = 1.0 - 2.0*exp(x2log);
	} while (rng() > 0.5*(1.0 + cosTheta*cosTheta)); //rejection method

	p.changeByCos(cosTheta, 2 * PI*rng());//only need to change the Direction of the scattered photon

	return 0; //means there's no energy lost
}

__device__ double Material::GCO(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	ParticleStack &pStack = tv[it].pStack;

	//cache from the slow device memory
	double pE = p.E;

	const double D2 = 1.4142135623731, D1 = 1.0 / D2, D12 = 0.5;

	double Ek = pE / Es;
	double Ek2 = Ek + Ek + 1.0; //2k+1
	double Ek3 = Ek*Ek; //k^2
	double Ek1 = Ek3 - Ek2 - 1.0; //k^2-2k-2
	double tmin = 1 / Ek2; // 1/(2k+1)
	double tmin2 = tmin*tmin; // 1/(2k+1)^2
	double a1 = log(Ek2);
	double a2 = a1 + 2.0*Ek*(1.0 + Ek)*tmin2; //a1+a2 in the manual

	double Ep = 0; //energy of the scattered photon
	double cosTheta = 0; //Direction of the scattered photon
	int iSh = -1; //index of the shell ionized, starting from 0

	double tau = 0;
	double Tcos = 0;
	double S = 0; // S(E,theta)

	if (pE > 5e6) //No Doppler broadening for E greater than 5 MeV
	{
		do
		{
			//sample tau
			do
			{
				if (rng()*a2 < a1) //i=1
					tau = pow(tmin, rng());
				else
					tau = sqrt(1.0 + rng()*(tmin2 - 1.0));
				Tcos = (1.0 + tau*(Ek1 + tau*(Ek2 + tau*Ek3))) / (Ek3*tau*(1.0 + tau*tau));
			} while (rng()>Tcos);
			Ep = tau*pE;
			cosTheta = 1.0 - (1.0 - tau) / (Ek*tau);

			//sample which shell excited
			Tcos = Zt*rng();
			S = 0;
			iSh = -1; //for search judgment
			for (int i = 0; i < NOSCCO; ++i) //for each shell, starting from 0;
			{
				S += FCO[i];
				if (S > Tcos) { iSh = i; break; }
			}
			if (iSh == -1) iSh = NOSCCO - 1; //avoid some critical problem;
		} while (Ep > pE - UICO[iSh]);
	}
	else //should take Doppler broadening into account
	{
		//********sample Direction of the scattered photon
		//calculate S(E,theta=pi), i.e. S0
		double S0 = 0;
		double ni[64];
		double pac[64]; //accumulation funcion of fi*theta(E-Ui)*ni, used to sample which shell being excited
		for (int i = 0; i <NOSCCO; ++i) //for each shell
		{
			if (pE > UICO[i]) //like theta(E-Ui)
			{
				double aux = pE*(pE - UICO[i])*2.0; //2E*(E-Ui), auxiliary variable
				double PzJi0 = FJ0[i] * (aux - Es*UICO[i]) / (Es*sqrt(aux + aux + UICO[i] * UICO[i]));
				if (PzJi0 > 0) ni[i] = 1 - 0.5*exp(D12 - (D1 + D2*PzJi0)*(D1 + D2*PzJi0));
				else ni[i] = 0.5 * exp(D12 - (D1 - D2*PzJi0)*(D1 - D2*PzJi0));
				S0 += FCO[i] * ni[i];
			}
		}

		double cdt1 = 0;
		do
		{
			// Sampling tau
			if (rng()*a2 < a1) tau = pow(tmin, rng());
			else tau = sqrt(1.0 + rng()*(tmin2 - 1.0));
			cdt1 = (1.0 - tau) / (Ek*tau); //1-cos(theta)
			//calculate S(E, theta)
			S = 0;
			for (int i = 0; i <NOSCCO; ++i) //for each shell
			{
				if (pE > UICO[i]) //like theta(E-Ui)
				{
					double aux = pE*(pE - UICO[i])*cdt1; //2E*(E-Ui), auxiliary variable
					double PzJi0 = FJ0[i] * (aux - Es*UICO[i]) / (Es*sqrt(aux + aux + UICO[i] * UICO[i]));
					if (PzJi0 > 0) ni[i] = 1 - 0.5*exp(D12 - (D1 + D2*PzJi0)*(D1 + D2*PzJi0));
					else ni[i] = 0.5 * exp(D12 - (D1 - D2*PzJi0)*(D1 - D2*PzJi0));
					S += FCO[i] * ni[i];
					pac[i] = S;
				}
				else pac[i] = S - 1e-6; //avoid sampling the E<Ui shell
			}
			Tcos = S*(1.0 + tau*(Ek1 + tau*(Ek2 + tau*Ek3))) / (Ek3*tau*(1.0 + tau*tau));
		} while (rng()*S0 > Tcos);
		cosTheta = 1.0 - cdt1; //got the Direction of the scattered photon

		//***************sampling pz;
		double pz = 0;
		double Fmax = 0;
		double Fpz = 0;
		double RS = 0;

		do
		{
			do
			{
				//sample which shell is excited
				iSh = -1; //the index of excited shell
				RS = S*rng(); //random sampled part of S
				for (int i = 0; i < NOSCCO; ++i)
				{
					if (RS < pac[i]) { iSh = i; break; }
				}
				if (iSh == -1) iSh = NOSCCO - 1; //avoid some critical problem;

				//sample pz using formula (2.56)
				double A = rng()*ni[iSh];
				if (A < 0.5) pz = (D1 - sqrt(D12 - log(A + A))) / (D2*FJ0[iSh]);
				else pz = (sqrt(D12 - log(2.0 - A - A)) - D1) / (D2*FJ0[iSh]);
			} while (pz < -1.0);

			//calculate Fmax
			double xqc = 1.0 + tau*(tau - 2.0*cosTheta);
			double AF = sqrt(xqc)*(1.0 + tau*(tau - cosTheta) / xqc);
			if (AF > 0.0) Fmax = 1.0 + AF*0.2;
			else Fmax = 1.0 - AF*0.2;
			Fpz = 1.0 + AF*max(min(pz, 0.2), -0.2);
		} while (Fmax*rng() > Fpz); //rejection selecting pz

		//******************calculate Energy of the scattered photon Ep
		double tz = pz*pz;
		double b1 = 1.0 - tz*tau*tau;
		double b2 = 1.0 - tz*tau*cosTheta;
		int signpz = (pz > 0) ? 1 : -1;
		Ep = pE*(tau / b1)*(b2 + signpz*sqrt(fabs(b2*b2 - b1*(1.0 - tz))));
	}
	//by here, we have got Ep, cosTheta, and iSh(start from 0)

	double DE = pE - Ep;//Energy loss of incident photon
	//*************calculate the energy of secondary electron
	double EE = pE - Ep; //energy of secondary electron
	int IZA = KZCO[iSh];
	int ISA = KSCO[iSh];
	if (ISA < 10)
	{
		if (UICO[iSh] > Ecutr) EE = EE - UICO[iSh];
		if (IZA>0) relax(IZA, ISA);
	}

	double phi = 2 * PI*rng();
	if (simuSecondary && EE > Eabs[electron]) // need to simulate the secondary electron
	{
		//*************calculate the Direction of the secondary electron
		double cosThetaE = 0;
		double Q2 = pE*pE + Ep*(Ep - 2 * pE*cosTheta);
		if (Q2 > 1e-12) cosThetaE = (pE - Ep*cosTheta) / sqrt(Q2);
		else cosThetaE = 1.0;

		//***********by here we got Ep, cosTheta, DE, EE, cosThetaE, iSh
		Particle pnew = p;
		pnew.type = electron;
		pnew.changeByCos(cosThetaE, phi + PI);
		pnew.E = EE;
		//store the secondary electron
		pStack.push(pnew);
	}

	if (Ep < Eabs[photon]) // can stop the simulation of this photon
	{
		DE = pE;
		evt = EXIT_EVT;
	}
	else //continue simulating this photon
	{
		p.changeByCos(cosTheta, phi);
		p.E = Ep;
	}

	return DE;
}

__device__ double Material::GPH(size_t it)
{
	double pE = tv[it].p.E;
	if (!simuSecondary)
	{
		tv[it].evt = EXIT_EVT; //means to end the simulation of this photon
		return pE; //lost all energy locally
	}

	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	int IE = tv[it].IE;
	ParticleStack &pStack = tv[it].pStack;

	//to calculate
	int IZZ = 0;
	int ISH = 0;
	double EE = 0;//energy of photoelectric electron
	//double EL = log(p.E);
	double EL = (IE + tv[it].XIE)/invDLE + ELowL; //faster than log operation

	double pac[35]; //accumulation probability for selecting which element is simulated
	int iEel[35]; //index of corresponding energy for each element

	double ptot = 0.0;

	int iL, iR, iM;
	double dEg = 0;
	double pL = 0;
	for (int i = 0; i < elemNum; ++i)
	{
		IZZ = IZ[i];
		//binary search to locate the energy
		iL = IPHF[IZZ - 1] - 1;
		iR = IPHL[IZZ - 1] - 1;
		do
		{
			iM = (iL + iR) / 2;
			if (EL > EPH[iM]) iL = iM;
			else iR = iM;
		} while (iR - iL > 1);
		iEel[i] = iL; //store the location for future use
		dEg = EPH[iL + 1] - EPH[iL];

		if (dEg > 1e-15)	 pL = XPH[iL][0] + (XPH[iL + 1][0] - XPH[iL][0]) / dEg *(EL - EPH[iL]);
		else pL = XPH[iL][0];
		ptot += weight[i] * exp(pL);
		pac[i] = ptot;
	}

	//****  Sample the active element
	int iex = -1; //index of excited element
	double randp = rng()*SGPH[IE] / numDensity;
	for (int i = 0; i < elemNum; ++i)
	{
		if (randp < pac[i])
		{
			iex = i;
			IZZ = IZ[iex];
			break;
		}
	}
	if (iex == -1) //delta interaction
	{
		IZZ = ISH = 0;
		EE = 0;

		return 0.0; // nothing happen
	}
	//****  Sample the active shell of that element
	int iexsh = -1;
	if (NPHS[IZZ - 1] > 1) //need to select which shell is activated
	{
		int iE = iEel[iex];
		dEg = EPH[iE + 1] - EPH[iE];
		//element::gXSect &g=mo[iex].gXS[iE];
		double pt = 0;
		//prepare accumulation table
		if (dEg > 1e-15) //linear interpolation
		{
			ptot = exp(XPH[iE][0] + (XPH[iE + 1][0] - XPH[iE][0]) / dEg*(EL - EPH[iE]));
			for (int j = 1; j <= NPHS[IZZ - 1]; ++j)
			{
				pL = XPH[iE][j] + (XPH[iE + 1][j] - XPH[iE][j]) / dEg*(EL - EPH[iE]);
				pt += exp(pL);
				pac[j - 1] = pt;
			}
		}
		else //don't need linear interpolation
		{
			ptot = exp(XPH[iE][0]);
			for (int j = 0; j < NPHS[IZZ - 1]; ++j)
			{
				pt += exp(XPH[iE][j + 1]);
				pac[j] = pt;
			}
		}

		//do the selection
		randp = rng()*ptot;
		for (int i = 0; i < NPHS[IZZ - 1]; ++i)
		{
			if (randp < pac[i])
			{
				iexsh = i;
				ISH = i + 1;
				break;
			}
		}
		if (iexsh == -1) ISH = 10;
	}
	else ISH = 1;

	//by here, we've got IZZ, ISH, calculate EE next
	if (ISH < 10)
	{
		double EBB = EB[IZZ - 1][ISH - 1];
		if (EBB > Ecutr) EE = pE - EBB;
		else
		{
			EE = pE;
			ISH = 10;
		}
	}
	else EE = pE;

	double DE = pE;
	//by here, we've got IZZ, ISH, EE, calculate scattering angle of the electron
	if (EE > Eabs[electron])
	{
		double cosTheta = 0; //scattered Direction
		if (EE > 1e9) cosTheta = 1.0;
		else
		{
			double gamma = 1 + EE / Es;
			double gamma2 = gamma*gamma;
			double beta = sqrt(1 - 1.0 / gamma2);
			double A = 1 / beta - 1.0;
			double Aux = 0.5*beta*gamma*(gamma - 1)*(gamma - 2); //auxiliary variable
			double g0 = 2.0*(1.0 / A + Aux);
			double nu = 0, xi = 0, gnu = 0;
			do
			{
				xi = rng();
				nu = 2.0 * A / ((A + 2.0)*(A + 2.0) - 4.0 * xi)*(2.0 * xi + (A + 2.0)*sqrt(xi));
				gnu = (2.0 - nu)*(1 / (A + nu) + Aux);
			} while (rng()*g0 > gnu);
			cosTheta = 1 - nu;
		}
		Particle pe = p; //information of the electron
		pe.type = electron;
		pe.E = EE;
		pe.changeByCos(cosTheta, 2 * PI*rng());
		pStack.push(pe);//store the photo-electron in stack
	}

	if (ISH < 9) //the inner shell is excited, and need relax simulation
	{
		relax(IZZ, ISH);
	}

	evt = EXIT_EVT; //means to end the simulation of this photon
	return DE;
}

__device__ double Material::GPP(size_t it)
{
	double pE = tv[it].p.E;
	if (!simuSecondary)
	{
		tv[it].evt = EXIT_EVT; //means to end the simulation of this photon
		return pE; //lost all energy locally
	}
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	ParticleStack &pStack = tv[it].pStack;

	double DE = pE;

	double rki = Es / pE;
	double eps = 0; //energy loss fraction from photon to electron
	if (pE < 1.1e6) //very close to the lower limit
	{
		eps = rki + (1.0 - 2.0*rki)*rng();
	}
	else
	{
		double a = ZEQPP*alpha;
		double a2 = a*a;
		double t = sqrt(2 * rki);
		double t2 = t*t;
		double F0 = (-1.774 - 1.210e1*a + 1.118e1*a2)*t
			+ (8.523 + 7.326e1*a - 4.441e1*a2)*t2
			- (1.352e1 + 1.211e2*a - 9.641e1*a2)*t*t2
			+ (8.946 + 6.205e1*a - 6.341e1*a2)*t2*t2;
		double g0 = Fpp + F0;
		double bm = 4 * rki / Bpp;
		double g1, g2;
		calcG12(bm, g1, g2);
		double phi1m = g1 + g0;
		double phi2m = g2 + g0;
		double rx = 0.5 - rki;
		double P1 = 2.0 / 3.0*rx*rx*phi1m;
		P1 = P1 / (P1 + phi2m);

		double b = 0;
		//sample eps
		while (true)
		{
			if (rng() <= P1) //case 1
			{
				double rnd = 2 * rng() - 1;
				if (rnd < 0)  eps = 0.5 - rx*pow(fabs(rnd), 1.0 / 3.0);
				else eps = 0.5 + rx*pow(rnd, 1.0 / 3.0);
				b = rki / (Bpp*eps*(1 - eps));
				calcG12(b, g1, g2);
				double phi1 = max(g1 + g0, 0.0);
				if (rng()*phi1m <= phi1) break;
			}
			else //case 2
			{
				eps = rki + 2.0*rx*rng();
				b = rki / (Bpp*eps*(1 - eps));
				calcG12(b, g1, g2);
				double phi2 = max(g2 + g0, 0.0);
				if (rng()*phi2m <= phi2) break;
			}
		}
	}

	//here we've got eps
	double cosThetaE = 0, cosThetaP; //Direction of the electron/positron
	double beta = 0;
	//produced electron
	double Ee = eps*pE - Es; //kinetic energy of electron
	if (Ee > Eabs[electron]) //need record new electron
	{
		cosThetaE = 2.0 * rng() - 1.0;
		beta = sqrt(Ee*(Ee + TEs)) / (Ee + Es);
		cosThetaE = (cosThetaE + beta) / (cosThetaE*beta + 1);
	}
	else cosThetaE = 1;

	//produced positron
	double Ep = (1 - eps)*pE - Es; //kinetic energy of positron
	if (Ep > Eabs[positron]) //need record new positron
	{
		cosThetaP = 2.0 * rng() - 1.0;
		beta = sqrt(Ep*(Ep + TEs)) / (Ep + Es);
		cosThetaP = (cosThetaP + beta) / (cosThetaP*beta + 1);
	}
	else cosThetaP = 1;

	Particle pnew;
	if (Ee > Eabs[electron]) //need record new electron
	{
		pnew = p; //copy some information
		pnew.type = electron;
		pnew.changeByCos(cosThetaE, 2 * PI*rng());
		pnew.E = Ee;
		pStack.push(pnew);
	}

	
	if (Ep > Eabs[positron]) //need record new positron
	{
		pnew = p;
		pnew.type = positron;
		pnew.changeByCos(cosThetaP, 2 * PI*rng());
		pnew.E = Ep;
		if (simuPositron)
		{
			pStack.push(pnew);
#ifdef _DEBUG
			p.u = pnew.u;
			p.v = pnew.v;
			p.w = pnew.w;
#endif
			DE -= TEs; //the positron carries away latent energy of 2Es
		}
		else
		{
			PANR();
		}
	}
	else
	{
		PANR(); // the kinetic energy is too low
	}

	evt = EXIT_EVT; //stop tracking the photon
	return DE;
}

__device__ double Material::EEL(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	int IE = tv[it].IE; //left index of interval
	double XIE = tv[it].XIE; //portion in that interval

	double pc = (1 - XIE)*RNDCE[IE] + XIE*RNDCE[IE + 1];
	double cosTheta = 0;
	if (p.E > EELMAX) cosTheta = EELt(exp((1 - XIE)*AE[IE] + XIE*AE[IE + 1]), (1 - XIE)*BE[IE] + XIE*BE[IE + 1], pc, it);
	else cosTheta = EELd(pc, it);
	p.changeByCos(cosTheta, 2 * PI*rng());
	return 0; // no energy lost
}

__device__ double Material::PEL(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	int IE = tv[it].IE; //left index of interval
	double XIE = tv[it].XIE; //portion in that interval

	double pc = (1 - XIE)*RNDCP[IE] + XIE*RNDCP[IE + 1];
	double cosTheta = 0;
	if (p.E > PELMAX) cosTheta = EELt(exp((1 - XIE)*AP[IE] + XIE*AP[IE + 1]), (1 - XIE)*BP[IE] + XIE*BP[IE + 1], pc, it);
	else cosTheta = PELd(pc, it);
	p.changeByCos(cosTheta, 2 * PI*rng());
	return 0; // no energy lost
}

__device__ double Material::EIN(size_t itt)
{
	GRNG &rng = tv[itt].rng;
	Particle &p = tv[itt].p;
	KonckEvent &evt = tv[itt].evt;
	ParticleStack &pStack = tv[itt].pStack;
	int IE = tv[itt].IE; //left index of interval
	double XIE = tv[itt].XIE; //portion in that interval

	//Random sampling of hard inelastic collisions of electrons.
	//Sternheimer-Liljequist GOS model
	//DE energy deposited, EP energy of the scattered electron, cosTheta angle of the scattered electron
	//ES energy of the secondary electron, cosThetaS angle of the secondary electron
	double DE, EP, cosTheta, ES, cosThetaS;//outputs
	double E = p.E; //current energy of electron
	double delta = (1 - XIE)*DEL[IE] + XIE*DEL[IE + 1];
	if (E < Wcc) //delta interaction
	{
		DE = 0;
		return DE;
	}

	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;

	//selecting the active oscillator, i0 is the selected index
	double r = rng();
	int i0 = 0;
	int j0 = NOSC;
	int it = 0;
	do
	{
		it = (i0 + j0) / 2;
		if (r > EINAC[JE][it]) i0 = it;
		else j0 = it;
	} while (j0 - i0 > 1);

	//constants
	double RB = E + TEs;
	double gamma = 1 + E / Es;
	double gamma2 = gamma*gamma;
	double beta2 = 1.0 - 1.0 / gamma2;
	double amol = (1.0 - 1.0 / gamma)*(1.0 - 1.0 / gamma);

	//partial cross sections of the active oscillator
	double Wi = WRI[i0];
	double cps, cp, qm, xhdl, xhdt;
	//distant excitation
	if (Wi > Wcc && Wi < E)
	{
		cps = E*RB;
		cp = sqrt(cps);
		double xhdt0 = max(log(gamma2) - beta2 - delta, 0.0);
		if (Wi>1.0e-6 * E)
		{
			double cpp = sqrt((E - Wi)*(E - Wi + TEs));
			qm = sqrt((cp - cpp)*(cp - cpp) + Es*Es) - Es;
		}
		else
		{
			qm = Wi*Wi / (beta2 * TEs);
			qm = qm*(1.0 - qm / TEs);
		}
		if (qm < Wi)
		{
			xhdl = log(Wi*(qm + TEs) / (qm*(Wi + TEs))) / Wi;
			xhdt = xhdt0 / Wi;
		}
		else
		{
			qm = Wi;
			xhdl = 0.0;
			xhdt = 0.0;
		}
	}
	else
	{
		qm = Wi;
		cps = 0.0;
		cp = 0.0;
		xhdl = 0.0;
		xhdt = 0.0;
	}

	//Close collisions
	double WMaxC = 0.5*E;
	double WCL = max(Wcc, Wi);
	double RCL = WCL / E;
	double xhc; //hard cross section of close interaction
	if (WCL < WMaxC)
	{
		double RL1 = 1.0 - RCL;
		xhc = (amol*(0.5 - RCL) + 1.0 / RCL - 1.0 / RL1 + (1.0 - amol)*log(RCL / RL1)) / E;
	}
	else xhc = 0;

	double xhtot = xhc + xhdl + xhdt;//total cross section
	if (xhtot < 1e-35) return 0.0; //no interaction will happen

	double rxh = rng()*xhtot;
	if (rxh < xhc)// close collision is chosen
	{
		double A = 5.0*amol;
		double ARCL = A*0.5*RCL;
		double FB, RK, RK2, RKF, PHI;
		do
		{
			FB = (1.0 + ARCL)*rng();
			if (FB < 1.0) RK = RCL / (1.0 - FB*(1.0 - (RCL + RCL)));
			else RK = RCL + (FB - 1.0)*(0.5 - RCL) / ARCL;
			RK2 = RK*RK;
			RKF = RK / (1.0 - RK);
			PHI = 1.0 + RKF*RKF - RKF + amol*(RK2 + RKF);
		} while (rng()*(1.0 + A*RK2)>PHI);

		// Energy and scattering angle(primary electron).
		DE = RK*E;
		EP = E - DE;
		cosTheta = sqrt(EP*RB / (E*(RB - DE)));
		// Energy and emission angle of the delta ray.
		if (KS[i0] < 10) //?????
		{
			if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
			else ES = DE;
		}
		else ES = DE;
		cosThetaS = sqrt(DE*RB / (E*(DE + TEs)));
	}
	else if (rxh < xhc + xhdl) //Hard distant longitudinal collision
	{
		DE = Wi;
		EP = E - DE;
		double QS = qm / (1.0 + qm / TEs);
		double Q = QS / (pow((QS / DE)*(1.0 + DE / TEs), rng()) - (QS / TEs));
		double QTREV = Q*(Q + TEs);
		double CPPS = EP*(EP + TEs);
		cosTheta = (CPPS + cps - QTREV) / (2.0*cp*sqrt(CPPS));
		if (cosTheta>1.0) cosTheta = 1.0;
		if (KS[i0] < 10) ES = DE - Ui[i0];
		else ES = DE;
		cosThetaS = 0.5*(DE*(E + RB - DE) + QTREV) / sqrt(cps*QTREV);
		if (cosThetaS>1.0) cosThetaS = 1.0;
	}
	else //Hard distant transverse collision
	{
		DE = Wi;
		EP = E - DE;
		cosTheta = 1.0;
		if (KS[i0] < 10) //?????
		{
			if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
			else ES = DE;
		}
		else ES = DE;
		cosThetaS = 0.5;
	}


	double phi = 2 * PI*rng();
	if (ES > Eabs[electron]) //record secondary electron
	{
		Particle psec = p;
		psec.changeByCos(cosThetaS, phi + PI);
		psec.E = ES;
		psec.type = electron;
		pStack.push(psec);
		
	}
	if (EP > Eabs[electron])
	{
		p.E = EP;
		p.changeByCos(cosTheta, phi);
		
		return DE; //will continue simulating this electron
	}

	evt = EXIT_EVT; //stop simulating this electron
	return E;
}

__device__ double Material::PIN(size_t itt)
{
	GRNG &rng = tv[itt].rng;
	Particle &p = tv[itt].p;
	KonckEvent &evt = tv[itt].evt;
	ParticleStack &pStack = tv[itt].pStack;
	int IE = tv[itt].IE; //left index of interval
	double XIE = tv[itt].XIE; //portion in that interval

	//Random sampling of hard inelastic collisions of electrons.
	//Sternheimer-Liljequist GOS model
	//DE energy deposited, EP energy of the scattered electron, cosTheta angle of the scattered electron
	//ES energy of the secondary electron, cosThetaS angle of the secondary electron
	double DE, EP, cosTheta, ES, cosThetaS;//outputs
	double E = p.E; //current energy of electron
	double delta = (1 - XIE)*DEL[IE] + XIE*DEL[IE + 1];
	if (E < Wcc) //delta interaction
	{
		DE = 0;
		return DE;
	}

	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;

	//selecting the active oscillator, i0 is the selected index
	double r = rng();
	int i0 = 0;
	int j0 = NOSC;
	int it = 0;
	do
	{
		it = (i0 + j0) / 2;
		if (r > PINAC[JE][it]) i0 = it;
		else j0 = it;
	} while (j0 - i0 > 1);

	//constants
	double RB = E + TEs;
	double gamma = 1 + E / Es;
	double gamma2 = gamma*gamma;
	double beta2 = 1.0 - 1.0 / gamma2;
	double G12 = (gamma + 1.0)*(gamma + 1.0);
	double amol = (E / (E + Es))*(E / (E + Es));
	double BHA1 = amol*(2.0*G12 - 1.0) / (gamma2 - 1.0);
	double BHA2 = amol*(3.0 + 1.0 / G12);
	double BHA3 = amol*2.0*gamma*(gamma - 1.0) / G12;
	double BHA4 = amol*(gamma - 1.0)*(gamma - 1.0) / G12;

	//partial cross sections of the active oscillator
	double Wi = WRI[i0];
	double cps, cp, qm, xhdl, xhdt;
	//distant excitation
	if (Wi > Wcc && Wi < E)
	{
		cps = E*RB;
		cp = sqrt(cps);
		double xhdt0 = max(log(gamma2) - beta2 - delta, 0.0);
		if (Wi>1.0e-6 * E)
		{
			double cpp = sqrt((E - Wi)*(E - Wi + TEs));
			qm = sqrt((cp - cpp)*(cp - cpp) + Es*Es) - Es;
		}
		else
		{
			qm = Wi*Wi / (beta2 * TEs);
			qm = qm*(1.0 - qm / TEs);
		}
		if (qm < Wi)
		{
			xhdl = log(Wi*(qm + TEs) / (qm*(Wi + TEs))) / Wi;
			xhdt = xhdt0 / Wi;
		}
		else
		{
			qm = Wi;
			xhdl = 0.0;
			xhdt = 0.0;
		}
	}
	else
	{
		qm = Wi;
		cps = 0.0;
		cp = 0.0;
		xhdl = 0.0;
		xhdt = 0.0;
	}

	//Close collisions
	double WMaxC = E;
	double WCL = max(Wcc, Wi);
	double RCL = WCL / E;
	double xhc; //hard cross section of close interaction
	if (WCL < WMaxC)
	{
		double RL1 = 1.0 - RCL;
		xhc = ((1.0 / RCL - 1.0) + BHA1*log(RCL) + BHA2*RL1
			+ (BHA3 / 2.0)*(RCL*RCL - 1.0) + (BHA4 / 3.0)*(1.0 - RCL*RCL*RCL)) / E;
	}
	else xhc = 0.0;

	double xhtot = xhc + xhdl + xhdt;//total cross section
	if (xhtot < 1e-35) return 0.0; //no interaction will happen

	double rxh = rng()*xhtot;
	if (rxh < xhc)// close collision is chosen
	{
		double RK, PHI;
		do
		{
			RK = RCL / (1.0 - rng()*(1.0 - RCL));
			PHI = 1.0 - RK*(BHA1 - RK*(BHA2 - RK*(BHA3 - BHA4*RK)));
		} while (rng() > PHI);

		// Energy and scattering angle(primary electron).
		DE = RK*E;
		EP = E - DE;
		cosTheta = sqrt(EP*RB / (E*(RB - DE)));
		// Energy and emission angle of the delta ray.
		if (KS[i0] < 10) //?????
		{
			if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
			else ES = DE;
		}
		else ES = DE;
		cosThetaS = sqrt(DE*RB / (E*(DE + TEs)));
	}
	else
	{
		DE = Wi;
		EP = E - DE;
		if (rxh < xhc + xhdl) //Hard distant longitudinal collision
		{
			double QS = qm / (1.0 + qm / TEs);
			double Q = QS / (pow((QS / DE)*(1.0 + DE / TEs), rng()) - (QS / TEs));
			double QTREV = Q*(Q + TEs);
			double CPPS = EP*(EP + TEs);
			cosTheta = (CPPS + cps - QTREV) / (2.0*cp*sqrt(CPPS));
			if (cosTheta > 1.0) cosTheta = 1.0;
			if (KS[i0] < 10) ES = DE - Ui[i0];
			else ES = DE;
			cosThetaS = 0.5*(DE*(E + RB - DE) + QTREV) / sqrt(cps*QTREV);
			if (cosThetaS>1.0) cosThetaS = 1.0;
		}
		else //Hard distant transverse collision
		{
			cosTheta = 1.0;
			if (KS[i0] < 10) //?????
			{
				if (Ui[i0] > Ecutr) ES = DE - Ui[i0];
				else ES = DE;
			}
			else ES = DE;
			cosThetaS = 0.5;
		}
	}
	double phi = 2 * PI*rng();
	if (ES > Eabs[electron]) //record secondary electron
	{
		Particle psec = p;
		psec.changeByCos(cosThetaS, phi + PI);
		psec.E = ES;
		psec.type = electron;
		pStack.push(psec);
	}
	if (EP > Eabs[positron])
	{
		p.E = EP;
		p.changeByCos(cosTheta, phi);
		return DE;
	}
	
	DE = E + TEs;
	evt = EXIT_EVT; //stop simulating current positron
	PANR();
	return DE;
}

__device__ double Material::EBR(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	ParticleStack &pStack = tv[it].pStack;
	int &IE = tv[it].IE; //left index of interval
	double &XIE = tv[it].XIE; //portion in that interval
	double pE = p.E;

	if (pE < Wcr) return 0.0;// no interaction happen
	int i, j, k, m; //index of energy
	if (rng() < XIE) i = IE + 1;
	else i = IE;

	double W = 0;//(reduced) photon energy
	double pt, W1, W2, bj, aj, pmax;
	do
	{
		pt = PBcut[i] + rng()*(PACB[i][NBW - 1] - PBcut[i]); //uniform random number from PACBcut[i] to PACB[i][NWB-1]
		//locate where the kappa is
		j = 0;
		k = NBW - 1;
		do
		{
			m = (j + k) / 2;
			if (pt > PACB[i][m]) j = m;
			else k = m;
		} while (k - j > 1);
		//j = binarySearch(&PACB[i][0], NBW, pt);
		//
		W1 = WB[j];
		W2 = WB[j + 1];
		bj = (PDFB[i][j + 1] - PDFB[i][j]) / (W2 - W1);
		aj = PDFB[i][j] - bj*W1;
		if (W1 < WBcut[i]) W1 = WBcut[i];
		if (W2 < W1)
		{
			printf("Warning: Conflicting end-point values in EBR()");
			W = W1;
			break;
		}
		pmax = max(aj + bj*W1, aj + bj*W2);
		do //rejection sampling
		{
			W = W1*pow(W2 / W1, rng());
		} while (rng()*pmax > aj + bj*W);
		W *= pE; //this is the radiated photon energy
	} while (W < Wcr); //radiated energy should be larger than Wcr, may not end???


	//here we got the radiated energy W. Next, do the angular sampling
	if (W > Eabs[photon]) //sample the direction of the photon
	{
		double beta = sqrt(pE*(pE + TEs)) / (pE + Es);
		double cosTheta = 1.0;
		if (pE > 500e3) //use a simplified dipole distribution for E>500keV
		{
			cosTheta = 2.0*rng() - 1.0;
			if (rng() > 0.75)
			{
				if (cosTheta < 0.0) cosTheta = -pow(-cosTheta, 1.0 / 3.0);
				else cosTheta = pow(cosTheta, 1.0 / 3.0);
			}
			cosTheta = (cosTheta + beta) / (1.0 + beta*cosTheta);
		}
		else //use  fitting expression (3.169)
		{
			int IEE = 0;
			if (beta >BET[5]) IEE = 5;
			else if (beta < BET[0]) IEE = 0;
			else //search for the region where beta located
			{
				IEE = 0;
				int IEE1 = 5;
				int IEEM;
				do
				{
					IEEM = (IEE + IEE1) / 2;
					if (beta>BET[IEEM]) IEE = IEEM;
					else IEE1 = IEEM;
				} while (IEE1 - IEE>1);
			}

			//calculate A and betap first
			double rk = W / pE * 20;
			int ik = min((int)rk, 19);
			double qL = BP1[IEE][ik][0] + beta*(BP1[IEE][ik][1] + beta*(BP1[IEE][ik][2] + beta*BP1[IEE][ik][3]));
			double qR = BP1[IEE][ik + 1][0] + beta*(BP1[IEE][ik + 1][1] + beta*(BP1[IEE][ik + 1][2] + beta*BP1[IEE][ik + 1][3]));
			double q1 = qL + (rk - ik)*(qR - qL);

			qL = BP2[IEE][ik][0] + beta*(BP2[IEE][ik][1] + beta*(BP2[IEE][ik][2] + beta*BP2[IEE][ik][3]));
			qR = BP2[IEE][ik + 1][0] + beta*(BP2[IEE][ik + 1][1] + beta*(BP2[IEE][ik + 1][2] + beta*BP2[IEE][ik + 1][3]));
			double q2 = qL + (rk - ik)*(qR - qL);

			double A = min(exp(q1) / beta, 1.0);
			double betap = min(max(beta*(1.0 + q2 / beta), 0.0), 0.999999999);
			if (rng() < A)
			{
				do
				{
					cosTheta = 2.0*rng() - 1.0;
				} while (2 * rng() > 1 + cosTheta*cosTheta);
			}
			else
			{
				do
				{
					cosTheta = 2.0*rng() - 1.0;
				} while (rng() > 1 - cosTheta*cosTheta);
			}
			cosTheta = (cosTheta + betap) / (1 + betap*cosTheta);
		}

		//here we got the scattering angle of photon
		Particle pnew = p;
		pnew.changeByCos(cosTheta, 2 * PI*rng());
		pnew.E = W;
		pnew.type = photon;
		pStack.push(pnew);
	}

	pE -= W; //current energy of electron;
	
	if (pE > Eabs[p.type])
	{
		p.E = pE;//write back this variable
		return W;
	}
	//else finish simulating this electron/positron
	evt = EXIT_EVT; 
	W += pE;
	if (p.type == positron)
	{
		W += TEs;
		PANR();
	}
	return W;
}

__device__ double Material::ESI(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	int IE = tv[it].IE; //left index of interval
	double XIE = tv[it].XIE; //portion in that interval

	//Ionization of inner shells by impact of electrons.
	//Output arguments, both start from 1 :
	//IZZ ....atomic number of the element where ionization has ocurred.
	//ISH ....atomic electron shell that has been  ionized.
	double PACSI[120];
	int IZSI[120], ISHSI[120];
	int kt = 0;
	PACSI[0] = 0;
	int IZZ, ISH, INDC;
	double SITOT = 0, XSI = 0;

	//prepare accumulation probability
	for (int j = 0; j < elemNum; ++j)
	{
		IZZ = IZ[j] - 1;
		INDC = IESIF[IZZ] - 1;
		for (int ish = 0; ish < NSESI[IZZ]; ++ish)
		{
			double Wcut = EB[IZZ][ish];
			if (Wcut>Ecutr && Wcut < p.E)
			{
				XSI = exp(XESI[INDC + IE][ish] + XIE*(XESI[INDC + IE + 1][ish] - XESI[INDC + IE][ish]));
				if (XSI > 1.1e-35)
				{
					SITOT += XSI*weight[j];
					IZSI[kt] = IZZ + 1;
					ISHSI[kt] = ish + 1;
					++kt;
					PACSI[kt] = SITOT;
				}
			}
		}
	}

	//sampling IZZ and ISH pair
	if (kt == 0) return 0.0;//do nothing

	double rx = PACSI[kt] * rng();
	int is = 0, im;
	do
	{
		im = (is + kt) / 2;
		if (rx > PACSI[im]) is = im;
		else kt = im;
	} while (kt - is>1);

	IZZ = IZSI[is];
	ISH = ISHSI[is];

	relax(IZZ, ISH);

	return 0;
}

__device__ double Material::PSI(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	int IE = tv[it].IE; //left index of interval
	double XIE = tv[it].XIE; //portion in that interval

	//Ionization of inner shells by impact of electrons.
	//Output arguments, both start from 1 :
	//IZZ ....atomic number of the element where ionization has ocurred.
	//ISH ....atomic electron shell that has been  ionized.
	double PACSI[120];
	int IZSI[120], ISHSI[120];
	int kt = 0;
	PACSI[0] = 0;
	int IZZ, ISH, INDC;
	double SITOT = 0, XSI = 0;

	//prepare accumulation probability
	for (int j = 0; j < elemNum; ++j)
	{
		IZZ = IZ[j] - 1;
		INDC = IPSIF[IZZ] - 1;
		for (int ish = 0; ish < NSPSI[IZZ]; ++ish)
		{
			double Wcut = EB[IZZ][ish];
			if (Wcut>Ecutr && Wcut < p.E)
			{
				XSI = exp(XPSI[INDC + IE][ish] + XIE*(XPSI[INDC + IE + 1][ish] - XPSI[INDC + IE][ish]));
				if (XSI > 1.1e-35)
				{
					SITOT += XSI*weight[j];
					IZSI[kt] = IZZ + 1;
					ISHSI[kt] = ish + 1;
					++kt;
					PACSI[kt] = SITOT;
				}
			}
		}
	}

	//sampling IZZ and ISH pair
	if (kt == 0) return 0.0;//do nothing

	double rx = PACSI[kt] * rng();
	int is = 0, im;
	do
	{
		im = (is + kt) / 2;
		if (rx > PACSI[im]) is = im;
		else kt = im;
	} while (kt - is > 1);

	IZZ = IZSI[is];
	ISH = ISHSI[is];

	relax(IZZ, ISH);

	return 0;
}

__device__ double Material::PAN(size_t it)
{
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	KonckEvent &evt = tv[it].evt;
	ParticleStack &pStack = tv[it].pStack;
	double pE = p.E;
	double DE = pE + TEs;

	double E1, cosTheta1, E2, cosTheta2;
	if (pE < Eabs[positron]) //Slow positrons (assumed at rest)
	{
		E1 = 0.5*(pE + TEs);
		E2 = E1;
		cosTheta1 = -1.0 + 2.0*rng();
		cosTheta2 = -cosTheta1;
	}
	else
	{
		//Annihilation in flight(two photons with energy and directions
		//determined from the DCS and energy - momentum conservation).
		double gamma = 1.0 + max(pE, 1.0) / Es;
		double GAM21 = sqrt(gamma*gamma - 1.0);
		double ANI = 1.0 + gamma;
		double CHIMIN = 1.0 / (ANI + GAM21);
		double RCHI = (1.0 - CHIMIN) / CHIMIN;
		double GT0 = ANI*ANI - 2.0;
		double CHI, GREJ;
		do
		{
			CHI = CHIMIN*pow(RCHI, rng());
			GREJ = ANI*ANI*(1.0 - CHI) + gamma + gamma - 1.0 / CHI;
		} while (rng()*GT0 > GREJ);

		double DET = pE + TEs;
		E1 = CHI*DET;
		cosTheta1 = (ANI - 1.0 / CHI) / GAM21;
		double CHIP = 1.0 - CHI;
		E2 = DET - E1;
		cosTheta2 = (ANI - 1.0 / CHIP) / GAM21;
	}

	//need to store the photons
	double phi = 2 * PI*rng();
	if (E1 > Eabs[photon])
	{
		Particle pnew = p;
		pnew.changeByCos(cosTheta1, phi);
		pnew.E = E1;
		pnew.type = photon;
		pStack.push(pnew);
	}
	if (E2 > Eabs[photon])
	{
		Particle pnew = p;
		pnew.changeByCos(cosTheta2, phi + PI);
		pnew.E = E2;
		pnew.type = photon;
		pStack.push(pnew);
	}

	DE = pE + TEs;
	evt = EXIT_EVT;

	return DE;
}

__device__ double Material::PAN(size_t it, Particle& p) //called by the GPP only, it instead uses a customized positron to do annihilation
{
	GRNG &rng = tv[it].rng;
	ParticleStack &pStack = tv[it].pStack;
	double pE = p.E;
	double DE = pE + TEs;

	double E1, cosTheta1, E2, cosTheta2;
	if (pE < Eabs[positron]) //Slow positrons (assumed at rest)
	{
		E1 = 0.5*DE;
		E2 = E1;
		cosTheta1 = -1.0 + 2.0*rng();
		cosTheta2 = -cosTheta1;
	}
	else
	{
		//Annihilation in flight(two photons with energy and directions
		//determined from the DCS and energy - momentum conservation).
		double gamma = 1.0 + max(pE, 1.0) / Es;
		double GAM21 = sqrt(gamma*gamma - 1.0);
		double ANI = 1.0 + gamma;
		double CHIMIN = 1.0 / (ANI + GAM21);
		double RCHI = (1.0 - CHIMIN) / CHIMIN;
		double GT0 = ANI*ANI - 2.0;
		double CHI, GREJ;
		do
		{
			CHI = CHIMIN*pow(RCHI, rng());
			GREJ = ANI*ANI*(1.0 - CHI) + gamma + gamma - 1.0 / CHI;
		} while (rng()*GT0 > GREJ);

		double DET = pE + TEs;
		E1 = CHI*DET;
		cosTheta1 = (ANI - 1.0 / CHI) / GAM21;
		double CHIP = 1.0 - CHI;
		E2 = DET - E1;
		cosTheta2 = (ANI - 1.0 / CHIP) / GAM21;
	}

	//need to store the photons
	double phi = 2 * PI*rng();
	if (E1 > Eabs[photon])
	{
		Particle pnew = p;
		pnew.changeByCos(cosTheta1, phi);
		pnew.E = E1;
		pnew.type = photon;
		pStack.push(pnew);
		DE -= E1;
	}
	if (E2 > Eabs[photon])
	{
		Particle pnew = p;
		pnew.changeByCos(cosTheta2, phi + PI);
		pnew.E = E2;
		pnew.type = photon;
		pStack.push(pnew);
		DE -= E2;
	}
	
	return DE;
}

//other process called by above interactions
__device__ void Material::PANR()
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	ParticleStack &pStack = tv[it].pStack;

	//simulate the static Annihilation of electron and positron
	Particle p1 = p; //copy position
	p1.type = photon;
	double cosTheta = -1 + 2 * rng();
	p1.changeByCos(cosTheta, 2 * PI*rng());
	p1.E = Es;

	Particle p2 = p1; //recoiled particle with opposite velocity
	p2.u = -p1.u;
	p2.v = -p1.v;
	p2.w = -p1.w;

	pStack.push(p1);
	pStack.push(p2);
}

__device__ double Material::relax(int IZ, int IS)
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	GRNG &rng = tv[it].rng;
	Particle &p = tv[it].p;
	ParticleStack &pStack = tv[it].pStack;
	double DE = 0;

	if (IZ < 6 || IS > 9) return 0;
	--IZ;
	--IS;
	//If the shell ionization energy is less than Ecutr, the cascade is not followed.
	if (EB[IZ][IS] < Ecutr) return 0;
	const int stackDepth = 256;
	int isv[stackDepth]; //vacant shell index stack
	int nv = 0; // stack index, pointing to the first element
	isv[0] = IS; //index of initial vacant shell
	int isp = 0;
	int KF, KL, K1, KS;
	double RN, RX;
	int is1, is2;
	ParticleType ptype;// type of generated particle
	do
	{
		//fetch one vacancy from the stack
		isp = isv[nv];
		KF = IFIRST[IZ][isp] - 1;
		KL = ILAST[IZ][isp] - 1;
		--nv;
		//select one transition
		if (KL > KF)
		{
			RN = rng()*(KL - KF + 1); // KL-KF+1 is the point number
			K1 = (int)RN;
			RX = RN - K1; //the decimal part
			if (RX > FR[KF + K1]) KS = IAL[KF + K1] - 1;
			else KS = KF + K1;
		}
		else KS = KF;
		//KS is the index 
		is1 = IS1[KS];
		is2 = IS2[KS];
		if (is2 == 0) // radiative 
		{
			ptype = photon;
			if (is1 < 10 && EB[IZ][is1 - 1] > Ecutr)
			{
				++nv;
				isv[nv] = is1 - 1;
			}
		}
		else // Non-radiative 
		{
			ptype = electron;
			if (is1 < 10 && EB[IZ][is1 - 1] > Ecutr)
			{
				++nv;
				isv[nv] = is1 - 1; // note -1 is necessary since isv starts from 0
			}
			if (is2 < 10 && EB[IZ][is2 - 1] > Ecutr)
			{
				++nv;
				isv[nv] = is2 - 1; // note -1 is necessary since isv starts from 0
			}
		}

		//store the emitted particle
		if (ET[KS]>Eabs[ptype])
		{
			Particle pnew = p; //Position and material index copied
			pnew.type = ptype;
			pnew.w = -1.0 + 2.0*rng();
			double sinTheta = sqrt(1 - pnew.w*pnew.w);
			double phi = 2 * PI* rng();
			pnew.u = sinTheta*cos(phi);
			pnew.v = sinTheta*sin(phi);
			pnew.E = ET[KS];

			pStack.push(pnew); //store second particle

			DE -= pnew.E; //this particle will take away energy
		}

	} while (nv >= 0);

	return DE;
}

__device__ double Material::EELt(double A, double B, double rndc, size_t it)
{
	//used for both electron and positron
	//return the cosTheta
	GRNG &rng = tv[it].rng;
	double A1 = A + 1;
	double B1;
	double mu;
	if (B > 0) //case I
	{
		double muav = A*A1*log(A1 / A) - A;
		B1 = 1 - B;
		double rnd0 = B1*A1*muav / (A + muav);
		double rnd = rndc + rng()*(1 - rndc);
		if (rnd < rnd0) mu = rnd*A / (B1*A1 - rnd);
		else if (rnd > rnd0 + B)
		{
			double rndmB = rnd - B;
			mu = rndmB*A / (B1*A1 - rndmB);
		}
		else mu = muav;
	}
	else //case II
	{
		double BB = -B;
		B1 = 1.0 - BB;
		double muc = rndc*A / (B1*A1 - rndc);
		double PW = B1*A*(1.0 - muc) / (A + muc);
		if (rng()*(BB + PW) < BB) mu = 0.5*(1.0 + sqrt(rng()));
		else
		{
			double rndrc = rng()*(1 - muc);
			mu = (A*rndrc + A1*muc) / (A1 - rndrc);
		}
	}

	return 1 - 2 * mu;
}

__device__ double Material::EELd(double rndc, size_t it)
{
	GRNG &rng = tv[it].rng;
	int &IE = tv[it].IE; //left index of interval
	double &XIE = tv[it].XIE; //portion in that interval

	double mu;
	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;
	double ru = rndc + rng()*(1 - rndc); //range from rndc to 1, hard event
	int itn = (int)(ru*(NP - 1));
	int i = ITLE[itn][JE] - 1;
	int j = ITUE[itn][JE] - 1;
	if (j - i >= 2)
	{
		int k;
		do
		{
			k = (i + j) / 2;
			if (ru > PSE[k][JE]) i = k;
			else j = k;
		} while (j - i > 1);
	}
	double rr = ru - PSE[i][JE];
	if (rr > 1e-16)
	{
		double xx = XSE[i][JE];
		double aa = ASE[i][JE];
		double bb = BSE[i][JE];
		double d = PSE[i + 1][JE] - PSE[i][JE];
		double cd = (1 + aa + bb)*d;
		mu = xx + (cd*rr / (d*d + (aa*d + bb*rr)*rr))*(XSE[i + 1][JE] - xx);
	}
	else mu = XSE[i][JE];

	return 1 - 2 * mu;
}

__device__ double Material::PELd(double rndc, size_t it)
{
	GRNG &rng = tv[it].rng;
	int &IE = tv[it].IE; //left index of interval
	double &XIE = tv[it].XIE; //portion in that interval

	double mu;
	int JE;
	if (rng() < XIE) JE = IE + 1;
	else JE = IE;
	double ru = rndc + rng()*(1 - rndc); //range from rndc to 1, hard event
	int itn = (int)(ru*(NP - 1));
	int i = ITLP[itn][JE] - 1;
	int j = ITUP[itn][JE] - 1;
	if (j - i >= 2)
	{
		int k;
		do
		{
			k = (i + j) / 2;
			if (ru > PSP[k][JE]) i = k;
			else j = k;
		} while (j - i > 1);
	}
	double rr = ru - PSP[i][JE];
	if (rr > 1e-16)
	{
		double xx = XSP[i][JE];
		double aa = ASP[i][JE];
		double bb = BSP[i][JE];
		double d = PSP[i + 1][JE] - PSP[i][JE];
		double cd = (1 + aa + bb)*d;
		mu = xx + (cd*rr / (d*d + (aa*d + bb*rr)*rr))*(XSP[i + 1][JE] - xx);
	}
	else mu = XSP[i][JE];

	return 1 - 2 * mu;
}

/*>>>>>>>>>>>>>>>>>>>>>>>>> end: Material methods definitions >>>>>>>>>>>>>>>>>>>>>>>>>*/



/*<<<<<<<<<<<<<<<<<<<<<<<<< start: tool functions of cuda <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
void cudaErrorCheck(cudaError_t cudaStatus,const char* func=NULL)
{
	if (cudaSuccess != cudaStatus)
	{
		if (func) Log("cuda error: %s in function module-> %s\n", cudaGetErrorString(cudaStatus), func);
		else Log("cuda error: %s\n", cudaGetErrorString(cudaStatus));
		exitApp("cuda function call failed!");
	}
}

void inline cudaKernelCheck(int i, int it = -1)
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
			else if (compute_perf == max_compute_perf)
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


/*<<<<<<<<<<<<<<<<<<<<<<<<< start: PENELOPE method definitions <<<<<<<<<<<<<<<<<<<<<<<<<<*/
void PENELOPE::loadCommonData()
{
	double Etab_h[NEGrid];
	double ELow_h, ELowL_h, EUp_h; //energy up and low bound
	double invDLE_h; //1.0 / delta log(E)
	double FLCOH_h;

// 	double EPH_h[8000], XPH_h[8000][10];
// 	int IPHF_h[99], IPHL_h[99], NPHS_h[99];
// 	double ET_h[15000], FR_h[15000];
// 	int IAL_h[15000], IS1_h[15000], IS2_h[15000], IFIRST_h[99][9], ILAST_h[99][9];
// 	double XESI_h[6000][9], XPSI_h[6000][9];
//  double EB_h[99][30];

	ArrayMgr<double> EPH_h(8000);
	ArrayMgr<double> XPH_h(8000, 10);	
	ArrayMgr<double> ET_h(15000);
	ArrayMgr<double> FR_h(15000);
	ArrayMgr<int> IAL_h(15000);
	ArrayMgr<int> IS1_h(15000);
	ArrayMgr<int> IS2_h(15000);
	ArrayMgr<int> IFIRST_h(99, 9);
	ArrayMgr<int> ILAST_h(99, 9);
	ArrayMgr<double> XESI_h(6000, 9);
	ArrayMgr<double> XPSI_h(6000, 9);
	ArrayMgr<double> EB_h(99, 30);

	int IPHF_h[99], IPHL_h[99], NPHS_h[99];
	int IESIF_h[99], NSESI_h[99], IPSIF_h[99], NSPSI_h[99];
	double BET_h[6];
	//init BET_h[6], used in EBR()
	double Ex[6] = { 1.0e3, 5.0e3, 1.0e4, 5.0e4, 1.0e5, 5.0e5 };
	for (int i = 0; i < 6; ++i) BET_h[i] = sqrt(Ex[i] * (Ex[i] + TEs)) / (Ex[i] + Es);

	FILE *fp = fopen("Materials/common.dat", "rb");
	if (fp == NULL) exitApp("Cannot open the file common.dat!");
	fread(&ELow_h, sizeof(double), 1, fp);
	fread(&EUp_h, sizeof(double), 1, fp);
	fread(&FLCOH_h, sizeof(double), 1, fp);
	fread(EPH_h.getP(), sizeof(double), EPH_h.getInnerLength(), fp);

	fread(XPH_h.getP(), sizeof(double), XPH_h.getInnerLength(), fp);
	XPH_h.toRowFirstLayout();
// 	for (int i = 0; i < 10; ++i)
// 	{
// 		for (int j = 0; j < 8000; ++j) fread(&XPH_h.a(j, i), sizeof(double), 1, fp);
// 	}
	fread(IPHF_h, sizeof(int), 99, fp);
	fread(IPHL_h, sizeof(int), 99, fp);
	fread(NPHS_h, sizeof(int), 99, fp);
	fread(ET_h.getP(), sizeof(double), 15000, fp);
	fread(FR_h.getP(), sizeof(double), 15000, fp);
	fread(IAL_h.getP(), sizeof(int), 15000, fp);
	fread(IS1_h.getP(), sizeof(int), 15000, fp);
	fread(IS2_h.getP(), sizeof(int), 15000, fp);

	fread(IFIRST_h.getP(), sizeof(int), IFIRST_h.getInnerLength(), fp);
	IFIRST_h.toRowFirstLayout();
// 	for (int i = 0; i < 9; ++i)
// 	{
// 		for (int j = 0; j < 99; ++j) fread(&IFIRST_h.a(j,i), sizeof(int), 1, fp);
// 	}

	fread(ILAST_h.getP(), sizeof(int), ILAST_h.getInnerLength(), fp);
	ILAST_h.toRowFirstLayout();
// 	for (int i = 0; i < 9; ++i)
// 	{
// 		for (int j = 0; j < 99; ++j) fread(&ILAST_h[j][i], sizeof(int), 1, fp);
// 	}

	fread(XESI_h.getP(), sizeof(double), XESI_h.getInnerLength(), fp);
	XESI_h.toRowFirstLayout();
// 	for (int i = 0; i < 9; ++i)
// 	{
// 		for (int j = 0; j < 6000; ++j) fread(&XESI_h[j][i], sizeof(double), 1, fp);
// 	}
	fread(IESIF_h, sizeof(int), 99, fp);
	fread(NSESI_h, sizeof(int), 99, fp);

	fread(XPSI_h.getP(), sizeof(double), XPSI_h.getInnerLength(), fp);
	XPSI_h.toRowFirstLayout();
// 	for (int i = 0; i < 9; ++i)
// 	{
// 		for (int j = 0; j < 6000; ++j) fread(&XPSI_h[j][i], sizeof(double), 1, fp);
// 	}
	fread(IPSIF_h, sizeof(int), 99, fp);
	fread(NSPSI_h, sizeof(int), 99, fp);

	fread(EB_h.getP(), sizeof(double), EB_h.getInnerLength(), fp);
	EB_h.toRowFirstLayout();
// 	for (int i = 0; i < 30; ++i)
// 	{
// 		for (int j = 0; j < 99; ++j) fread(&EB_h[j][i], sizeof(double), 1, fp);
// 	}
	fclose(fp);

	//init energy grid
	ELowL_h = log(ELow_h);
	double EUpL = log(EUp_h);
	double DLE_h = (EUpL - ELowL_h) / (NEGrid - 1);
	invDLE_h = 1.0 / DLE_h;
	for (int i = 0; i < NEGrid; ++i)
	{
		Etab_h[i] = exp(ELowL_h + i*DLE_h);
		comConst.Etab_h[i] = Etab_h[i];
	}

	for (unsigned int i = 0; i < gc.size(); ++i)
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));
		//copy data to the GPU
		cudaErrorCheck(cudaMemcpyToSymbol(Etab, &Etab_h, sizeof(Etab_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(ELow, &ELow_h, sizeof(ELow_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(ELowL, &ELowL_h, sizeof(ELowL_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(EUp, &EUp_h, sizeof(EUp_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(invDLE, &invDLE_h, sizeof(invDLE_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(FLCOH, &FLCOH_h, sizeof(FLCOH_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(EPH, EPH_h.getP(), EPH_h.getByteNum()));
		cudaErrorCheck(cudaMemcpyToSymbol(XPH, XPH_h.getP(), XPH_h.getByteNum())); //? must debug to ensure
		cudaErrorCheck(cudaMemcpyToSymbol(IPHF, &IPHF_h, sizeof(IPHF_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(IPHL, &IPHL_h, sizeof(IPHL_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(NPHS, &NPHS_h, sizeof(NPHS_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(ET, ET_h.getP(), ET_h.getByteNum()));
		cudaErrorCheck(cudaMemcpyToSymbol(FR, FR_h.getP(), FR_h.getByteNum()));
		cudaErrorCheck(cudaMemcpyToSymbol(EPH, &EPH_h, sizeof(EPH_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(IAL, IAL_h.getP(), IAL_h.getByteNum()));
		cudaErrorCheck(cudaMemcpyToSymbol(IS1, IS1_h.getP(), IS1_h.getByteNum()));
		cudaErrorCheck(cudaMemcpyToSymbol(IS2, IS2_h.getP(), IS2_h.getByteNum()));
		cudaErrorCheck(cudaMemcpyToSymbol(IFIRST, IFIRST_h.getP(), IFIRST_h.getByteNum())); //?
		cudaErrorCheck(cudaMemcpyToSymbol(ILAST, ILAST_h.getP(), ILAST_h.getByteNum())); //?
		cudaErrorCheck(cudaMemcpyToSymbol(XESI, XESI_h.getP(), XESI_h.getByteNum())); //?
		cudaErrorCheck(cudaMemcpyToSymbol(XPSI, XPSI_h.getP(), XPSI_h.getByteNum())); //?
		cudaErrorCheck(cudaMemcpyToSymbol(IESIF, &IESIF_h, sizeof(IESIF_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(NSESI, &NSESI_h, sizeof(NSESI_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(IPSIF, &IPSIF_h, sizeof(IPSIF_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(NSPSI, &NSPSI_h, sizeof(NSPSI_h)));
		cudaErrorCheck(cudaMemcpyToSymbol(EB, EB_h.getP(), EB_h.getByteNum())); //?
		cudaErrorCheck(cudaMemcpyToSymbol(BET, &BET_h, sizeof(BET_h)));
	}
}

void PENELOPE::phantom2GPU() //copy phantom info to GPU. Make sure the phantom has loaded data
{
	int NVoxel = _phant->getVoxelNum();
	SFloat* dp = new SFloat[NVoxel]; // to initialize the dose score in GPU to zero
	for (int i = 0; i < NVoxel; ++i) dp[i] = 0;
	void* pGPU = NULL;
	for (unsigned int i = 0; i < gc.size(); ++i)
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));

		cudaErrorCheck(cudaMalloc(&pGPU, NVoxel*sizeof(SFloat)), "resize ph");
		gc[i].addGPUPointer(pGPU);
		cudaErrorCheck(cudaMemcpyToSymbol(ph, &pGPU, sizeof(SFloat *)), "copy ph pointer to GPU constant"); //set the array pointer
		cudaErrorCheck(cudaMemcpy(pGPU, _phant->ph.getP(), sizeof(SFloat)*NVoxel, cudaMemcpyHostToDevice), "init the value of ph in GPU");

		cudaErrorCheck(cudaMalloc(&pGPU, NVoxel*sizeof(SFloat)), "resize doseScore");
		gc[i].addGPUPointer(pGPU);
		cudaErrorCheck(cudaMemcpyToSymbol(doseScore, &pGPU, sizeof(SFloat *)), "copy doseScore pointer to GPU constant"); //set the array pointer
		cudaErrorCheck(cudaMemcpy(pGPU, dp, sizeof(SFloat)*NVoxel, cudaMemcpyHostToDevice), "init the value of doseScore in GPU");
		gc[i].d_doseScore = (SFloat*)pGPU; // it will be used to copy back the dose from GPU to CPU

		//copy the rest constant
		double temp;
		cudaErrorCheck(cudaMemcpyToSymbol(NX, &_phant->NX, sizeof(int)), "copy NX to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(NY, &_phant->NY, sizeof(int)), "copy NY to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(NZ, &_phant->NZ, sizeof(int)), "copy NZ to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(DX, &_phant->DX, sizeof(double)), "copy DX to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(DY, &_phant->DY, sizeof(double)), "copy DY to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(DZ, &_phant->DZ, sizeof(double)), "copy DZ to GPU");
		temp = 1.0 / _phant->DX;
		cudaErrorCheck(cudaMemcpyToSymbol(invDX, &temp, sizeof(double)), "copy invDX to GPU");
		temp = 1.0 / _phant->DY;
		cudaErrorCheck(cudaMemcpyToSymbol(invDY, &temp, sizeof(double)), "copy invDY to GPU");
		temp = 1.0 / _phant->DZ;
		cudaErrorCheck(cudaMemcpyToSymbol(invDZ, &temp, sizeof(double)), "copy invDZ to GPU");

		cudaErrorCheck(cudaMemcpyToSymbol(LX, &_phant->LX, sizeof(double)), "copy LX to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(LY, &_phant->LY, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(LZ, &_phant->LZ, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(xo, &_phant->xo, sizeof(double)), "copy xo to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(yo, &_phant->yo, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(zo, &_phant->zo, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(Bx, &_phant->Bx, sizeof(double)), "copy Bx to GPU");
		cudaErrorCheck(cudaMemcpyToSymbol(By, &_phant->By, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(Bz, &_phant->Bz, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(MaxDensity, &_phant->MaxDensity, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(rf, &_phant->rf, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(uniform, &_phant->uniform, sizeof(int)), "copy uniform to GPU");
	}

	delete[] dp;
}

void PENELOPE::initMaterial(ConfigFile *cf)
{
	//determine if we need to prepare the data tables
	string useLast;
	cf->getValue("useLastMats", useLast);
	bool useLastMats = useLast.compare("yes") == 0 ? true : false;
	bool fileExist = existFile("Materials/common.dat");
	ConfigFile *subCF = NULL;
	_MatList.resize(0);
	cf->resetSearchIndex();
	for (_NMAT = 1;; ++_NMAT)
	{
		subCF = cf->getBlock("MAT", true);
		if (NULL == subCF)
		{
			if (1 == _NMAT)
			{
				exitApp("cannot load any material parameters!"); //cannot even find the first material
			}
			else
			{
				--_NMAT;
				break;//at least loaded one material
			}
		}
		char ic[20];
		sprintf(ic, "Materials/mat%d.dat", _NMAT);
		fileExist = existFile(ic);
		int IDNumber = 0;
		subCF->getValue("ID number", IDNumber);
		_MatList.push_back(IDNumber);
	}

	if (useLastMats&&fileExist)
	{
		Log("Penelope tries to load material files generated last time...\n");
		Log.warningColor();
		Log("Make sure that EMax and material parameters weren't changed! \n");
		Log.normalColor();
	}
	else //generate the data file
	{
		if (pID == 0) // Only the main process takes care of data table preparation
		{
			RunTimeCounter rc;
			Log("Penelope will prepare material data tables in advance...\n");
			//<<load the lib functions
#ifdef WIN32
			//SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, FALSE); //try to unload the control-C handler
			HINSTANCE hlib = ::LoadLibrary("PINIT.dll");
#else
			void* hlib = dlopen("./PINIT.so", RTLD_LAZY);
#endif
			if (hlib == NULL) exitApp("Cannot open PINIT.dll or PINIT.so");
			typedef void (CallStyle *SETPARA)(int*, double*, double*, double*, double*, double*, double*, double*);
			typedef void (CallStyle *PINIT)(double*, int*, int*, int*);
			SETPARA setPara = (SETPARA)getLibAddress(hlib, "setpara");
			if (setPara == NULL) exitApp("cannot load function setpara");
			PINIT initMat = (PINIT)getLibAddress(hlib, "pinit");
			if (initMat == NULL) exitApp("cannot load function pinit");
			//>>

			double EMax = SourceHead_Emax();
			cf->getValue("EMax", EMax); //may overwrite EMax


			ConfigFile *subCF = NULL;
			_NMAT = 1;
			double Eabs_e, Eabs_g, Eabs_p, C1, C2, Wcc, Wcr;
			int IDNumber;

			//read parameters of each material
			FILE* fp = NULL;
			fp = fopen("temp_filename.txt", "w");
			if (NULL == fp) exitApp("cannot create temp_filename.txt");
			fprintf(fp, "Materials/\n"); //tell the FORTRAN code in which directory it should export

			if (!fs::exists(fs::path("Materials"))) fs::create_directories(fs::path("Materials"));
			string dbpath;
			cf->getValue("material database", dbpath);
			//search for the definition of beams
			cf->resetSearchIndex();
			for (_NMAT = 1;; ++_NMAT)
			{
				subCF = cf->getBlock("MAT", true);
				if (NULL == subCF)
				{
					if (1 == _NMAT)
					{
						exitApp("cannot load any material parameters!"); //cannot even find the first material
					}
					else
					{
						--_NMAT;
						break;//at least loaded one material
					}
				}
				subCF->getValue("Eabs_e", Eabs_e);
				subCF->getValue("Eabs_g", Eabs_g);
				subCF->getValue("Eabs_p", Eabs_p);
				subCF->getValue("C1", C1);
				subCF->getValue("C2", C2);
				subCF->getValue("Wcc", Wcc);
				subCF->getValue("Wcr", Wcr);
				subCF->getValue("ID number", IDNumber);
				
				if (Eabs_e == 0 || Eabs_g == 0 || Eabs_p == 0 || Wcc == 0 || Wcr == 0) //set the default cutoff energy
				{
					Eabs_e = Eabs_p = 1e-2*EMax;
					Eabs_g = 1e-3*EMax;
					Wcc = Eabs_e;
					Wcr = Eabs_g;
				}

				//search for the material file
				char ic[10];
				sprintf(ic, "%d", IDNumber);
				string filter(ic);
				filter += " ";
				string match = filter;
				matchFile("Materials", match);
				if (match.empty()) //need to exact from the database
				{
					if (1 <= IDNumber && IDNumber <= 280)
					{
						string cmd = "7z e -oMaterials \"";
						cmd += dbpath;
						cmd += "\" \"";
						cmd += ic;
						cmd += " *\"";
						system(cmd.c_str());

						int numTry = 0;
						do
						{
							++numTry;
							// wait some time for the exaction to be finished
							std::this_thread::sleep_for(std::chrono::milliseconds(500));
							match = filter;
							matchFile("Materials", match);
							if (numTry>5) exitApp("Failed to exact the material files from PENELOPE database. Please check the database file's path");
						} while (match.empty());
					}
					else exitApp("Cannot find the designated material file in the database or in the Materials directory");
				}
				fprintf(fp, "%s\n", match.c_str()); //only main process write the file
				setPara(&_NMAT, &Eabs_e, &Eabs_g, &Eabs_p, &C1, &C2, &Wcc, &Wcr);
			}
			fclose(fp); // close the file

			int info = 0; // output information level
			int useBinaryFile = 1;
			initMat(&EMax, &_NMAT, &info, &useBinaryFile); // generate the data files
			if (hlib) freeLib(hlib);
			fs::remove(fs::path("temp_filename.txt"));

			Log("\nIt costs %f seconds to prepare Penelope materials. ", rc.stop());
		}
	}
	
#ifdef USE_MPI
	//let's make sure the data files have been created by synchronizing all processes using MPI_Bcast
	char ic[10];
	MPI_Bcast(ic, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
#endif

	cf->getValue("DSMax", _DSMax);
	cf->getValue("particle stack depth", _NStackDepth);

	//get the loop cycle number for each interaction
	cf->getValue("GRA", simu_GRA);
	cf->getValue("GCO", simu_GCO);
	cf->getValue("GPH", simu_GPH);
	cf->getValue("GPP", simu_GPP);

	cf->getValue("EEL", simu_EEL);
	cf->getValue("EIN", simu_EIN);
	cf->getValue("EBR", simu_EBR);
	cf->getValue("ESI", simu_ESI);

	cf->getValue("SEC", simu_SEC);
	cf->getValue("POS", simu_POS);

	cf->getValue("PEL", simu_PEL);
	cf->getValue("PIN", simu_PIN);
	cf->getValue("PSI", simu_PSI);
	cf->getValue("PAN", simu_PAN);

	//load the material data to CPU memory
	_pMat = new Material[_NMAT];
	for (int i = 0; i < _NMAT; ++i) _pMat[i].load(i + 1);

}

void PENELOPE::initGPU()
{
	RunTimeCounter rc;
	int NGPU = (int)gc.size();

	//load the common data to GPU
	loadCommonData();
	void* pGPU = NULL;
	for (int i = 0; i < NGPU; ++i) //for each GPU
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));

		int NGPUThread = gc[i].NBlock*gc[i].BlockSize;

		//resize initial particle memory in GPU
		cudaErrorCheck(cudaMalloc(&pGPU, NGPUThread * gc[i].NBatch * sizeof(Particle)));
		gc[i].addGPUPointer(pGPU);
		cudaErrorCheck(cudaMemcpyToSymbol(InitPars, &pGPU, sizeof(Particle *)));
		gc[i].d_InitPars = (Particle*)pGPU; //it will be used when copying the initial particles to GPU

		//resize memory and copy Material data to GPU
		cudaErrorCheck(cudaMalloc(&pGPU, sizeof(Material)*_NMAT));
		gc[i].addGPUPointer(pGPU);
		cudaErrorCheck(cudaMemcpyToSymbol(mat, &pGPU, sizeof(Material *)));
		cudaErrorCheck(cudaMemcpy(pGPU, _pMat, sizeof(Material)*_NMAT, cudaMemcpyHostToDevice));

		//resize memory for thread specific variables in GPU
		cudaErrorCheck(cudaMalloc(&pGPU, NGPUThread * sizeof(ThreadVars)));
		gc[i].addGPUPointer(pGPU);
		cudaErrorCheck(cudaMemcpyToSymbol(tv, &pGPU, sizeof(ThreadVars *)));

		//resize memory for the particle stack in GPU
		cudaErrorCheck(cudaMalloc(&pGPU, NGPUThread * _NStackDepth* sizeof(Particle)));
		gc[i].addGPUPointer(pGPU);
		cudaErrorCheck(cudaMemcpyToSymbol(stackBuff, &pGPU, sizeof(Particle *)));

		cudaErrorCheck(cudaMemcpyToSymbol(NStackDepth, &_NStackDepth, sizeof(int)));
		cudaErrorCheck(cudaMemcpyToSymbol(NBatch, &gc[i].NBatch, sizeof(int)));
		cudaErrorCheck(cudaMemcpyToSymbol(DSMax, &_DSMax, sizeof(double)));
		cudaErrorCheck(cudaMemcpyToSymbol(simuPositron, &simu_POS, sizeof(int)));
		cudaErrorCheck(cudaMemcpyToSymbol(simuSecondary, &simu_SEC, sizeof(int)));
	}

	Log("\nIt costs %f seconds to init GPU ", rc.stop());
}


#if(TRANSPORT_OCTREE == 1)
void PENELOPE::initOctree(ConfigFile *octree_cf)
{
	//check if to go through the octree
	string valid;
	octree_cf->getValue("go through", valid);
	if (valid.compare("yes") != 0)
	{
		_UseOctree = false;
		return;
	}
	_UseOctree = true;

	string fname;
	int maxLevel = 9, maxTriangles = 6;
	octree_cf->getValue("3D mesh file", fname);
	octree_cf->getValue("max level", maxLevel);
	octree_cf->getValue("max triangles", maxTriangles);
	createCPUOctree(fname, maxLevel, maxTriangles, _MatList);
	
	int NGPU = gc.size();
#ifdef USE_OPENMP
#pragma omp parallel num_threads(NGPU)
#endif
	{
		int it = omp_get_thread_num(); //thread index starting from 0
		cudaErrorCheck(cudaSetDevice(gc[it].id));
		gc[it].d_gOctreeRoot = creatGPUOctree();
		cudaErrorCheck(cudaMemcpyToSymbol(gOctreeRoot, &gc[it].d_gOctreeRoot, sizeof(goctree_node *)), "copy goctree_node pointer to GPU constant");
		cudaErrorCheck(cudaMemcpyToSymbol(ThroughOctree, &_UseOctree, sizeof(_UseOctree)), "copy ThroughOctree to GPU constant");
	}

	deleteCPUOctree();
}
#endif

void PENELOPE::init(ConfigFile *cf)
{
	_cf = cf;//back up the whole config file
	ConfigFile* pe_cf = cf->getBlock("PENELOPE");
	initMaterial(pe_cf);
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

	int NVox = _phant->getVoxelNum();
	char pid2mid[300];
	for (size_t i = 0; i < _MatList.size(); ++i) pid2mid[_MatList[i]] = char(i);
#if MULTI_MATERIAL == 1
	if (_phant->NMAT == 1) 
	{
		Log.warningColor();
		Log("Warning: Only one material detected in this phantom, please use single material version for better performance");
		Log.normalColor();
	}
	double WoodcockSigma[NEGrid];
	
	char* id = new char[NVox];
	for (int IE = 0; IE < NEGrid; ++IE) // for each energy grid point
	{
		double maxSigma = 0;
		for (int i = 0; i < NVox; ++i) //for each voxel
		{
			id[i] = pid2mid[_phant->matid.a(i)]; //covert from penelope id to material sequence index
			double tSigma = _phant->ph.a(i)*_pMat[id[i]].totalSigmaG(IE);
			if (tSigma > maxSigma) maxSigma = tSigma;
		}
		//now we have the max sigma, but we need to take log of them in order to do interpolation
		WoodcockSigma[IE] = log(maxSigma);
	}
	//copy the woodcock data to GPU
	int NGPU = (int)gc.size();
	void* pGPU = NULL;
	for (int i = 0; i < NGPU; ++i) //for each GPU
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));
		cudaErrorCheck(cudaMemcpyToSymbol(WoodcockLogSigma, &WoodcockSigma, sizeof(WoodcockSigma)));

		cudaErrorCheck(cudaMalloc(&pGPU, NVox*sizeof(char)), "allocate matid in GPU");
		gc[i].addGPUPointer(pGPU);
		cudaErrorCheck(cudaMemcpyToSymbol(matid, &pGPU, sizeof(char *)), "copy matid pointer to GPU constant"); //set the array pointer
		cudaErrorCheck(cudaMemcpy(pGPU, id, sizeof(char)*NVox, cudaMemcpyHostToDevice), "init the value of matid in GPU");
	}
	delete[] id;
#endif
	phantom2GPU();
	// Change the relative density to absolute density in order to calculate the dose in CPU
	for (int i = 0; i < NVox; ++i)
	{
		char mid = pid2mid[_phant->matid.a(i)];
		_phant->ph.a(i) *= SFloat(_pMat[mid].massDensity);
	}
#if(TRANSPORT_OCTREE == 1)
	initOctree(cf->getBlock("OCTREE"));
#endif

}

int PENELOPE::getGPUConfig(ConfigFile* gcf)
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

	int NBlock = 128, BlockSize = 256, NBatch = 100, Source_Reuse_Times = 10;
	string rngStat;
	gcf->getValue("GPU Block Num", NBlock);
	gcf->getValue("GPU Block Dim", BlockSize);
	gcf->getValue("GPU Batch Num", NBatch);
	//gcf->getValue("GPU RNG Statistic", rngStat);
	//gcf->getValue("GRNG Refill Period", GRNG_Refill_Period);
	gcf->getValue("Source Reuse Times", Source_Reuse_Times);

	GPUConfig gpuc;
	gpuc.NBlock = NBlock;
	gpuc.BlockSize = BlockSize;
	gpuc.NBatch = NBatch;
// 	gpuc.refillPeriod = GRNG_Refill_Period;
// 	if (rngStat.compare("yes") == 0) gpuc.rngStat = true;
// 	else gpuc.rngStat = false;
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


	Log("/******************* The following GPU(s) will be used ************************/");
	for (unsigned int i = 0; i < gc.size(); ++i)	printGPUProperties(gc[i].id);
	Log("/************************ End GPU description *********************************/\n\n");

	//create streams of GPU control
	return main_id;
}
/*>>>>>>>>>>>>>>>>>>>>>>>>>> end: PENELOPE method definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


//simple statistic extension of source pool
#ifdef STAT_HEAD
static unsigned long st[63][63] = { 0 };
void SourcePool::statHead()
{
	double x0, y0, z0;
	double u, v, w;
	double xp, yp;
	int i, j;
	Particle* op;
	for (int it = 0; it < _NGenerateThread; ++it)
	{
		op = _PPool[it];
		for (int ip = 0; ip < _PoolCapacity; ++ip)
		{
			x0 = op[ip].x;
			y0 = op[ip].y;
			z0 = op[ip].z;
			u = op[ip].u;
			v = op[ip].v;
			w = op[ip].w;
			xp = -z0 / w*u + x0;
			yp = -z0 / w*v + y0;
			i = int((fabs(xp) + 0.05) / 0.1);
			j = int((fabs(yp) + 0.05) / 0.1);
			if (xp < 0) i = 31 - i;
			else i = 31 + i;
			if (yp < 0) j = 31 - j;
			else j = 31 + j;
			if (i < 0) i = 0;
			if (i > 62) i = 62;
			if (j < 0) j = 0;
			if (j > 62) j = 62;
			++st[i][j];
		}
	}
}
void SourcePool::statOutput()
{
	FILE *fp = fopen("source_stat.txt", "w");
	for (int i = 0; i < 63; ++i)
	{
		for (int j = 0; j < 63; ++j)
		{
			fprintf(fp, "%lu", st[i][j]);
			if (j == 62) fprintf(fp, "\n");
			else fprintf(fp, ";");
		}
	}
	fclose(fp);
}
#endif

RunTimeCounter sourceCounter;

void getSource(SourcePool* sp, volatile int* hist)
{
	sourceCounter.start();
	*hist = sp->prepareCopy();//prepare one batch before any run

	Log("time cost to generate particle = %f s", sourceCounter.stop());
}


int NThread, NProcess, pID;
bool b_abort = false;
bool b_peek = false;
bool b_thread_active = false;
BinaryFile gBF; //to store the raw dose data

ProgressCallBack progressCB = NULL;
PeekDoseCallBack peekDoseCB = NULL;
JobFinishedCallBack jobFinishedCB = NULL;
LogCallBack exitCB = NULL;

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

extern "C" DEXPORT void peekDose() //mark to generate intermediate dose
{
	if(b_thread_active) b_peek = true; //only if the thread is on
}

extern "C" DEXPORT void stopSimulation()
{
	if(b_thread_active) b_abort = true; //only if the thread is on
}

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
	int executionType = 0;
	cf.getValue("execution type", executionType);
	string outname;
	cf.getValue("output file name", outname);

	int NPeek = 0; // if the simulation takes too long time, we can peek the intermediate data.
	cf.getValue("peek output num", NPeek);
	int CPeek = 0; // current peek index
	fs::path cpath = fs::path(outname).parent_path(); //get current directory
	fs::directory_iterator end_itr;//search all files ended with "%.dose" and delete them.
	for (fs::directory_iterator ifile(cpath); ifile != end_itr; ++ifile)
	{
		if (fs::is_regular_file(ifile->status())) //check if it's a file
		{
			string fname = ifile->path().string();
			if (string::npos != fname.find("%.dose"))
			{
#ifdef WIN32
				string cmd = "del /q /f ";
#else //LINUX
				string cmd = "rm -f ";
#endif
				cmd += fname;
				system(cmd.c_str());
			}
		}
	}

	//search and config the GPU part ->> get gc
	PENELOPE PE;
	vector<GPUConfig>& gc = PE.gc;//makes the name shorter
	ConfigFile *gcf = cf.getBlock("GPU");
	int main_id = PE.getGPUConfig(gcf);

	//initialize PENELOPE and SourceHead by configurations
	ConfigFile *scf = cf.getBlock("SOURCEHEAD");
	SourceHead_Init(scf);
	PE.init(&cf);

	fNSIMU -= PE._phant->Hist; //subtract previous histories
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
	SourcePool sp(NGT, NOneFetch, NGStack);

	Log("\nCalculating dose, please wait patiently...\n\n");
	cf.getValue("log short mode", logDescription);
	if (logDescription.compare("yes") == 0) Log.shortMode(true);
	RunTimeCounter rc; //count the calculating time

	//note history number != particle number in the source pool because the particle in source pool isn't the very original one
	const int NGPU = (int)gc.size();
	//history number generated once by the source pool, only modified by one thread,
	//but accessed by multiple threads. It can be read only after the generating thread finished.
	volatile int histNew = 0; //histNew isn't always the same as hisAdd because histNew is modified in an isolated thread
	volatile int histAdd = 0; //this variable is shared by all threads, so add the key word "volatile" for safe
	const int firstSeed = PE._phant->seedBegin(2 * NGPU *gc[0].NBlock*gc[0].BlockSize);
	std::thread sthread; //source generating thread, unattached

	RunTimeCounter kernelCounter;
	RunTimeCounter copyCounter;
	
	sthread = std::thread(&getSource, &sp, &histNew); //start a fetch prepare, return how many new histories generated in histNew

	vector<SFloat*> dose(NGPU);//to store dose from all GPU cards
	vector<SFloat*> uncertainty(NGPU); //to store uncertainty from all GPU cards
	int NVoxel = PE._phant->getVoxelNum();
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

	extern bool b_abort;//it will be set true when the user tries to close the program manually
	bool b_targetErrReached = false;
	//each CPU thread takes care of one GPU
#ifdef USE_OPENMP
#pragma omp parallel num_threads(NGPU)
#endif
	{
		int it = omp_get_thread_num(); //if it==main_id, it's the main thread
		cudaErrorCheck(cudaSetDevice(gc[it].id)); //set which GPU this thread will operate on
		
		const double fHMax = fNSIMU / NGPU; //max history number this GPU will simulate
		double hist = 0; //how many histories has been simulated

		SFloat* gdose = new SFloat[NVoxel]; //to fetch temporary dose from GPU
		memset(gdose, 0, sizeof(SFloat)*NVoxel);
		// resize memory for CPU end storage
		dose[it] = new SFloat[NVoxel]; //to store the final dose of this thread
		uncertainty[it] = new SFloat[NVoxel]; //to store the uncertainty
		// initialize the dose score
		for (int i = 0; i < NVoxel; ++i)
		{
			dose[it][i] = PE._phant->dose[i] / NGPU;
			uncertainty[it][i] = PE._phant->uncertainty[i] / NGPU;
		}

		//generate a random thread for each working thread
		int seed = firstSeed + 2*it*NBlock*BlockSize; //make sure all seeds are unique
		initThreads << <NBlock, BlockSize >> >(seed); // init GRNG and particle stack
		cudaKernelCheck(0);

		int proceed = 1; //indicator about whether to continue
		int source_reuse = Source_Reuse_Times; //count how many times has been reused; force to initially generate incident particles
		int nGRA = 1;
		int nGPH = 1;
		int nGPP = 1;
		int nEEL = 1;
		int nEIN = 1;
		int nEBR = 1;
		int nESI = 1;
		int nPOS = 1;
		

		do //calculating main loop, end when hist >= fHMax
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
				if (it == main_id) copyCounter.start();
				sp.copy2GPU(gc[it].d_InitPars); //copy particles to the corresponding GPU
				if (it == main_id) copyCounter.stop();
				source_reuse = 0; //reset the reuse counter
#pragma omp barrier //wait until all GPUs received the data
				//if (it == main_id) sthread = std::thread(&getSource, &sp, &histNew);
				if (it == main_id) sthread = std::thread(&getSource, &sp, &histNew); //start a fetch prepare
			}
			++source_reuse; //count how many time the source has been used
			hist += histAdd; // histAdd more histories will be simulated

			/****************** Begin a batch run on GPU ********************/
			if (it == main_id) kernelCounter.start(); //only count the kernel time for the main GPU

			int nLoop = 0; //for performance profiling
#if(EXECUTE_TYPE < 0)
			if (executionType == 0) // my "kernel-by-kernel" approach
#endif
#if(EXECUTE_TYPE < 0 || EXECUTE_TYPE == 0)
			{

				start << <NBlock, BlockSize >> >(); //prepare run once in GPU
				cudaKernelCheck(1, it);
				proceed = 1; //set to start one GPU batch
				
				while (proceed)//loop until this batch is finished
				{
					proceed = 0;
					cudaErrorCheck(cudaMemcpyToSymbol(ContinueSignal, &proceed, sizeof(int)));//reset the variable in GPU to zero

					//determine the interaction type of photon
					photonJump << <NBlock, BlockSize >> >();
					cudaKernelCheck(2, it);

					if (PE.simu_GRA > 0)
					{
						if (nGRA >= PE.simu_GRA)
						{
							photonRA << <NBlock, BlockSize >> >();
							cudaKernelCheck(5, it);
							nGRA = 0;
						}
						++nGRA;
					}

					if (PE.simu_GCO) photonCO << <NBlock, BlockSize >> >();
					cudaKernelCheck(6, it);

					if (PE.simu_GPH > 0)
					{
						if (nGPH >= PE.simu_GPH)
						{
							photonPH << <NBlock, BlockSize >> >();
							cudaKernelCheck(7, it);
							nGPH = 0;
						}
						++nGPH;
					}

					if (PE.simu_GPP > 0)
					{
						if (nGPP >= PE.simu_GPP)
						{
							photonPP << <NBlock, BlockSize >> >();
							cudaKernelCheck(8, it);
							nGPP = 0;
						}
						++nGPP;
					}

					if (PE.simu_SEC)
					{
						electronJump << <NBlock, BlockSize >> >();
						cudaKernelCheck(3, it);
						electronSoft << <NBlock, BlockSize >> >();
						cudaKernelCheck(9, it);

						if (PE.simu_EEL > 0)
						{
							if (nEEL >= PE.simu_EEL)
							{
								electronEL << <NBlock, BlockSize >> >();
								cudaKernelCheck(10, it);
								nEEL = 0;
							}
							++nEEL;
						}
						if (PE.simu_EIN > 0)
						{
							if (nEIN >= PE.simu_EIN)
							{
								electronIN << <NBlock, BlockSize >> >();
								cudaKernelCheck(11, it);
								nEIN = 0;
							}
							++nEIN;
						}

						if (PE.simu_EBR > 0)
						{
							if (nEBR >= PE.simu_EBR)
							{
								electronBR << <NBlock, BlockSize >> >();
								cudaKernelCheck(12, it);
								nEBR = 0;
							}
							++nEBR;
						}

						if (PE.simu_ESI > 0)
						{
							if (nESI >= PE.simu_ESI)
							{
								electronSI << <NBlock, BlockSize >> >();
								cudaKernelCheck(13, it);
								nESI = 0;
							}
							++nESI;
						}

						if (PE.simu_POS > 0)
						{
							if (nPOS >= PE.simu_POS)
							{
								positronJump << <NBlock, BlockSize >> >();
								cudaKernelCheck(4, it);
								positronSoft << <NBlock, BlockSize >> >();
								cudaKernelCheck(14, it);
								if (PE.simu_PEL) positronEL << <NBlock, BlockSize >> >();
								cudaKernelCheck(15, it);
								if (PE.simu_PIN) positronIN << <NBlock, BlockSize >> >();
								cudaKernelCheck(16, it);
								if (PE.simu_PSI) positronSI << <NBlock, BlockSize >> >();
								cudaKernelCheck(17, it);
								if (PE.simu_PAN) positronAN << <NBlock, BlockSize >> >();
								cudaKernelCheck(18, it);
								nPOS = 0;
							}
							++nPOS;
						}
					}

					fillParticle << <NBlock, BlockSize >> >(); //refill the "exited" thread with stacked particle or new particle
					cudaKernelCheck(19, it);
					++nLoop; 

					cudaDeviceSynchronize();//make sure all kernels are finished
					//check if we can terminate the simulation in GPU
					cudaErrorCheck(cudaMemcpyFromSymbol(&proceed, ContinueSignal, sizeof(int)));
#ifdef STEP_DEBUG_STATUS
					if (-1 == proceed) //dump the status array
					{
						int len = STEP_NUM;
						double* dumpArray = new double[2*len];
						cudaMemcpyFromSymbol(dumpArray, DebugTack, sizeof(double) *2* len);
						FILE* fp = fopen("gPEN_debug.dat", "w");
						//fwrite(dumpArray, 7 * sizeof(double), len, fp);
						for (int i = 0; i < len; ++i)
							fprintf(fp, "%7d\t%.8g\n", (int)dumpArray[2 * i], dumpArray[2 * i + 1]);
						fclose(fp);
						delete[] dumpArray;
						exitApp("Exported the debug data!");
					}
#endif

#ifdef ACCU_DEBUG
					if (-2 == proceed) break;
#endif
				}

			}
#endif
#if(EXECUTE_TYPE < 0)
			else // traditional "one thread one history" approach
#endif
#if(EXECUTE_TYPE < 0 || EXECUTE_TYPE == 1)
			{

				PENELOPE_Kernel << <NBlock, BlockSize >> >();
				cudaDeviceSynchronize();
#ifdef STEP_DEBUG_STATUS
				cudaErrorCheck(cudaMemcpyFromSymbol(&proceed, ContinueSignal, sizeof(int)));
				if (-1 == proceed) //dump the status array
				{
					int len = STEP_NUM;
					double* dumpArray = new double[2 * len];
					cudaMemcpyFromSymbol(dumpArray, DebugTack, sizeof(double) * 2 * len);
					FILE* fp = fopen("gPEN_debug.dat", "w");
					//fwrite(dumpArray, 7 * sizeof(double), len, fp);
					for (int i = 0; i < len; ++i)
						fprintf(fp, "%7d\t%.8g\n", (int)dumpArray[2 * i], dumpArray[2 * i + 1]);
					fclose(fp);
					delete[] dumpArray;
					exitApp("Exported the debug data!");
				}
#endif
			}
#endif

			if (it == main_id)
			{
#if(EXECUTE_TYPE < 0 || EXECUTE_TYPE == 0)
				if (executionType == 0) Log("Loop number = %d", nLoop);
#endif
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
					Log.warningColor();
					Log("Source time is set to be %fs, kernel time is %fs", sourceCounter.getStoredTime(), kernelCounter.getStoredTime());
					Log("Source_Reuse_Times is set to be %d", Source_Reuse_Times);
					Log.normalColor();
				}
			}
			/****************** End a batch run on GPU *********************/

			cudaErrorCheck(cudaMemcpy(gdose, gc[it].d_doseScore, sizeof(SFloat)*NVoxel, cudaMemcpyDeviceToHost)); //fetch the batch dose from GPU
			for (int i = 0; i < NVoxel; ++i)
			{
				dose[it][i] += gdose[i];
				uncertainty[it][i] += gdose[i] * gdose[i];
				gdose[i] = 0;
			}	
			cudaErrorCheck(cudaMemcpy(gc[it].d_doseScore, gdose, sizeof(SFloat)*NVoxel, cudaMemcpyHostToDevice)); //reset the dose counter in GPU


#ifdef ACCU_DEBUG
			if (-2 == proceed) break;
#endif

#pragma omp barrier //wait until all GPU threads arrive here because Source_Reuse_Times maybe changed
			
			if (it == main_id)
			{
				//calculate if we're going the peek the dose due to NPeek setting
				if (NPeek > 0)
				{
					double fNP = (hist / fHMax)*(NPeek + 1);
					if (CPeek < int(fNP) && int(fNP) <= NPeek)
					{
						++CPeek; //increase current process indicator
						b_peek = true;
					}
				}

				if (b_peek) //need to output the dose to file or in gBF
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
					Phantom tempPhant(*PE._phant);
					double norm = 1.0889e15 * SourceHead_BeamOnTime();
					tempPhant.addDose(d, u, NGPU * nBatch, NGPU * hist, norm, thresholdRegion);
					delete[] d;
					delete[] u;
					tempPhant.getBinaryFile(gBF);

					b_peek = false;
					if (NPeek > 0)
					{
						char chs[50];
						sprintf(chs, "_%.0f%%", hist / fHMax * 100);
						string outnamet = outname;
						outnamet += chs;
						outnamet += ".dose";
						tempPhant.output(outnamet.c_str());
					}
					else
					{
						if(peekDoseCB) peekDoseCB(gBF);
					}
				}
			}
			
			if (targetErr > 0) //calculate the uncertainty
			{
				if (it == main_id)
				{
					double err2 = PE._phant->peekUncertainty(dose[it], uncertainty[it], nBatch + PE._phant->nBatch / NGPU, thresholdRegion);
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

		} while (true);

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
	PE._phant->addDose(d, u, NGPU * nBatch, totHist, norm, thresholdRegion);
	PE._phant->getBinaryFile(gBF);
	delete[] d;
	delete[] u;

	if (b_abort) outname += "_abort";// make sure it wouldn't overwrite the last file
	outname += ".dose";
	string vrFormat;
	cf.getValue("ViewRay format", vrFormat);
	if (vrFormat.compare("yes") == 0) PE._phant->output(outname.c_str(), 1);
	else PE._phant->output(outname.c_str());
	Log("Wait for the source generating thread to finish...");
	sthread.join();
	SourceHead_Delete(); //release the source head resource safely
	Log("Source reuse times = %d", Source_Reuse_Times);
	Log("\nTime statistics for main GPU:");
	Log("Total copy time =%.2f s", copyCounter.getStoredTime());
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