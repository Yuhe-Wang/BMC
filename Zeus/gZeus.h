#pragma once
#ifndef _GZEUS_H_
#define _GZEUS_H_

#include "../Tools/Tools.h"

#include "../SourceHead/SourceHead.h"
#include "../ZeusData/ZeusData.h"

//resolve the compiler's warning
#ifndef __CUDACC__
#define __device__ 
#define __host__
#define __forceinline__
#define __constant__
#endif
//#define USE_MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

//define USE_SINGLE_PRECISION

//notice GlueF will operate on float constant or float function like sin, sqrt
#ifdef USE_SINGLE_PRECISION
typedef float ZFloat;
#define GlueF(a) a##f
#else
typedef double ZFloat;
#define GlueF(a) a
#endif
/************************Start: Constant Definitions **************************/
#define PI GlueF(3.141592653589793238463)
#define Es GlueF(5.10998902e5)
#define TEs GlueF(1.021997804e6)
#define alpha (7.297352536e-3)
#define INV_ELECTRON_MASS GlueF(1.956951337e-6)
#define FLTEPSILON GlueF(1.192092896e-07)

#define NMatMax 5 //the maximum material number
#define NWOODCOCK 4096
#define NLAMDAPHOTON 2048
#define NLAMDACOMPTON 256
#define NLAMDAPAIR 1024
#define NQSURFACE_E 256
#define NQSURFACE_Q 128
#define NSCREENINGPARAMETER 512
#define NRANGE 1024
#define NINVERSERANGE 1024
/************************End: Constant Definitions **************************/

/*************************Start: type/class definition *********************************/
class ParticleR //store status information of a particle
{
public:
	void changeByCos(ZFloat cosTheta, ZFloat phi)
	{
		ZFloat cosPhi = GlueF(cos)(phi);
		ZFloat sinPhi = GlueF(sqrt)(1 - cosPhi*cosPhi);
		if (phi > PI) sinPhi = -sinPhi;
		ZFloat dxy = u*u + v*v;
		ZFloat dxyz = dxy + w*w;
		if (GlueF(fabs)(dxyz - GlueF(1.0)) > GlueF(1e-14))
		{
			ZFloat norm = 1 / GlueF(sqrt)(dxyz);
			u *= norm;
			v *= norm;
			w *= norm;
			dxy = u*u + v*v;
		}

		if (dxy > GlueF(1.0e-28))
		{
			ZFloat sinTheta = GlueF(sqrt)((1 - cosTheta*cosTheta) / dxy);
			ZFloat up = u;
			u = u*cosTheta + sinTheta*(up*w*cosPhi - v*sinPhi);
			v = v*cosTheta + sinTheta*(v*w*cosPhi + up*sinPhi);
			w = w*cosTheta - dxy*sinTheta*cosPhi;
		}
		else
		{
			ZFloat sinTheta = GlueF(sqrt)(1 - cosTheta*cosTheta);
			if (w > 0)
			{
				u = sinTheta*cosPhi;
				v = sinTheta*sinPhi;
				w = cosTheta;
			}
			else
			{
				u = -sinTheta*cosPhi;
				v = -sinTheta*sinPhi;
				w = -cosTheta;
			}
		}
	}
	//data
	ZFloat x, y, z; //position, unit cm
	ZFloat u, v, w; //direction vector
	int ivx, ivy, ivz, iabsv; //voxel index for current particle

	ParticleType type;

	ZFloat E; //energy, unit eV
	ZFloat weight;
};

class ParticleStack
{
	//method
public:
	int maxDepth; //for debug purpose
	ParticleStack(int n)
	{
		ss = NULL;
		ss = new ParticleR[n];
		cur = 0;
		Depth = n;
		maxDepth = 0;
	}
	~ParticleStack()
	{
		if (ss != NULL) delete[] ss;
	}
	__device__ void push(ParticleR &par)
	{
		if (cur < Depth)
		{
			ss[cur] = par; //copy and store the particle
			++cur;
			if (cur > maxDepth) maxDepth = cur;
		}
		else
		{
			printf("exceed the max depth of the GPU stack!!!!!!!!!!!!!!!\n");
		}
	}
	__device__ __forceinline__ bool empty()
	{
		if (cur > 0) return false;
		else return true;
	}
	__device__ __forceinline__ void pop()
	{
		if (cur <= 0) printf("no element to pop in the stack!");
		--cur;
	}
	__device__ __forceinline__ ParticleR& top(){ return ss[cur - 1]; }
	__device__ void init(ParticleR* pp){ cur = 0; ss = pp; }
	//data
private:
	ParticleR* ss;
	int cur; //point to next free space, 0 at initial time
	int Depth;
};



/*************************End: type/class definition *********************************/

void executeJob(const char* configFileName, MPS& configMacro); //execute one job according to the config file

#endif