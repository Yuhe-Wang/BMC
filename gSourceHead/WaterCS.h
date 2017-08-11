#pragma once
#ifndef _WATERCS_H_
#define _WATERCS_H_

#define NWOODCOCK 4096
#define NLAMDAPHOTON 2048
#define NLAMDACOMPTON 256
#define NLAMDAPAIR 1024
#define NQSURFACE_E 256
#define NQSURFACE_Q 128
#define NSCREENINGPARAMETER 512
#define NRANGE 1024
#define NINVERSERANGE 1024

__constant__ ZFloat PhotonCS[NLAMDAPHOTON];
__constant__ ZFloat ComptonCS[NLAMDACOMPTON];
__constant__ ZFloat PairCS[NLAMDAPAIR];
__constant__ ZFloat ScreenCS[NSCREENINGPARAMETER];
__constant__ ZFloat RangeCS[NRANGE];
__constant__ ZFloat InverseRangeCS[NINVERSERANGE];
Texture2DFloat(WaterQS)
__device__ ZFloat WaterQSData[NQSURFACE_E][NQSURFACE_Q];
ZFloat h_WaterQSData[NQSURFACE_E][NQSURFACE_Q];

#define NLOGCACHE 1024
#define NSINECACHE 2048
#ifdef USE_SINGLE_PRECISION
__constant__ ZFloat NegLogX01[NLOGCACHE]; //to make -log(x), 0<x<1 faster
__constant__ ZFloat CosCache[NSINECACHE]; // cache cos(phi)
__constant__ ZFloat SinCache[NSINECACHE]; // cache sin(phi)
#else
__device__ ZFloat NegLogX01[NLOGCACHE]; //to make -log(x), 0<x<1 faster
__device__ ZFloat CosCache[NSINECACHE]; // cache cos(phi)
__device__ ZFloat SinCache[NSINECACHE]; // cache sin(phi)
#endif

bool initWaterCS(const char* filename, vector<GPUConfig>& gc)
{
	FILE* fp = fopen(filename, "rb");
	if (NULL == fp) return false;
	ZFloat* h_PhotonCS = new ZFloat[NLAMDAPHOTON];
	ZFloat* h_ComptonCS = new ZFloat[NLAMDACOMPTON];
	ZFloat* h_PairCS = new ZFloat[NLAMDAPAIR];
	ZFloat* h_ScreenCS = new ZFloat[NSCREENINGPARAMETER];
	ZFloat* h_RangeCS = new ZFloat[NRANGE];
	ZFloat* h_InverseRangeCS = new ZFloat[NINVERSERANGE];
	ZFloat* h_NegLogX01 = new ZFloat[NLOGCACHE];
	ZFloat* h_CosCache = new ZFloat[NSINECACHE];
	ZFloat* h_SinCache = new ZFloat[NSINECACHE];

	double aIndex, bIndex;
	double cs;
	int n;
	//photon
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		h_PhotonCS[i] = ZFloat(cs);
	}

	//Compton
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		h_ComptonCS[i] = ZFloat(cs);
	}

	//pair production
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		h_PairCS[i] = ZFloat(cs);
	}

	//screening
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		h_ScreenCS[i] = ZFloat(cs);
	}

	//range
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		h_RangeCS[i] = ZFloat(cs);
	}

	//inverse range
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		h_InverseRangeCS[i] = ZFloat(cs);
	}

	//FILE* fpw = fopen("QS.txt", "w");
	//qsurface
	double fq, ieq0;
	int ne, nuq;
	fread(&fq, sizeof(double), 1, fp);
	fread(&ieq0, sizeof(double), 1, fp);
	fread(&ne, sizeof(int), 1, fp);
	fread(&nuq, sizeof(int), 1, fp);
	if (ne != NQSURFACE_E) return false;
	if (nuq != NQSURFACE_Q) return false;
	for (int ie = 0; ie < ne; ++ie)
		for (int iq = 0; iq < nuq; ++iq)
		{
			fread(&cs, sizeof(double), 1, fp);
			h_WaterQSData[ie][iq] = ZFloat(cs);
// 			if(iq!=nuq-1) fprintf(fpw, "%g,", cs);
// 			else fprintf(fpw, "%g\n", cs);
		}
	fclose(fp);

	//fclose(fpw);

	double dx = (1.00000001 - 0.01) / (NLOGCACHE - 1);
	for (int i = 0; i < NLOGCACHE; ++i)
	{
		h_NegLogX01[i] = -ZFloat(log(0.01 + i*dx));
	}

	dx = 2 * PI / (NSINECACHE - 1);
	for (int i = 0; i < NSINECACHE; ++i)
	{
		h_SinCache[i] = ZFloat(sin(i * dx));
		h_CosCache[i] = ZFloat(cos(i * dx));
	}

	//Let's copy these data to GPU
	for (unsigned int i = 0; i < gc.size(); ++i)
	{
		cudaErrorCheck(cudaSetDevice(gc[i].id));
		cudaErrorCheck(cudaMemcpyToSymbol(PhotonCS, h_PhotonCS, sizeof(ZFloat)*NLAMDAPHOTON), "copy PhotonCS");
		cudaErrorCheck(cudaMemcpyToSymbol(ComptonCS, h_ComptonCS, sizeof(ZFloat)*NLAMDACOMPTON), "copy ComptonCS");
		cudaErrorCheck(cudaMemcpyToSymbol(PairCS, h_PairCS, sizeof(ZFloat)*NLAMDAPAIR), "copy PairCS");
		cudaErrorCheck(cudaMemcpyToSymbol(ScreenCS, h_ScreenCS, sizeof(ZFloat)*NSCREENINGPARAMETER), "copy ScreenCS");
		cudaErrorCheck(cudaMemcpyToSymbol(RangeCS, h_RangeCS, sizeof(ZFloat)*NRANGE), "copy RangeCS");
		cudaErrorCheck(cudaMemcpyToSymbol(InverseRangeCS, h_InverseRangeCS, sizeof(ZFloat)*NINVERSERANGE), "copy InverseRangeCS");
		cudaErrorCheck(cudaMemcpyToSymbol(WaterQSData, h_WaterQSData, sizeof(ZFloat)*NQSURFACE_E*NQSURFACE_Q), "copy WaterQSData");
		cudaErrorCheck(cudaMemcpyToSymbol(NegLogX01, h_NegLogX01, sizeof(ZFloat)*NLOGCACHE), "copy NegLogX01");
		cudaErrorCheck(cudaMemcpyToSymbol(SinCache, h_SinCache, sizeof(ZFloat)*NSINECACHE), "SinCache");
		cudaErrorCheck(cudaMemcpyToSymbol(CosCache, h_CosCache, sizeof(ZFloat)*NSINECACHE), "CosCache");
	}

	delete[] h_PhotonCS;
	delete[] h_ComptonCS;
	delete[] h_PairCS;
	delete[] h_ScreenCS;
	delete[] h_RangeCS;
	delete[] h_InverseRangeCS;
	//delete[] h_WaterQSData;
	delete[] h_NegLogX01;
	delete[] h_CosCache;
	delete[] h_SinCache;
	return true;
}

__device__ ZFloat WaterPhotonCS(ZFloat E)
{
	ZFloat r = (E - GlueF(4.9e4))*GlueF(0.0014107512060647829);
	int i = int(r);
	r -= i;
	return (1 - r)*PhotonCS[i] + r*PhotonCS[i + 1];
}

__device__ ZFloat WaterComptonCS(ZFloat E)
{
	ZFloat r = (E - GlueF(4.9e4))*GlueF(0.0001757408683666437);
	int i = int(r);
	r -= i;
	return (1 - r)*ComptonCS[i] + r*ComptonCS[i + 1];
}

__device__ ZFloat WaterPairCS(ZFloat E)
{
	ZFloat r = (E - GlueF(1.024e6))*GlueF(0.0021491596638655462);
	int i = int(r);
	r -= i;
	return (1 - r)*PairCS[i] + r*PairCS[i + 1];
}

__device__ ZFloat WaterScreenCS(ZFloat E)
{
	ZFloat r = (E + GlueF(2.0408163265306123e-005)) *GlueF(25884562.370778777);
	int i = int(r);
	r -= i;
	return (1 - r)*ScreenCS[i] + r*ScreenCS[i + 1];
}

__device__ ZFloat WaterRangeCS(ZFloat E)
{
	ZFloat r = (E - GlueF(10e3))*GlueF(0.00068657718120805373);
	int i = int(r);
	r -= i;
	return (1 - r)*RangeCS[i] + r*RangeCS[i + 1];
}

__device__ ZFloat WaterInverseRangeCS(ZFloat E)
{
	ZFloat r = E *GlueF(1447.3007247284104);
	int i = int(r);
	r -= i;
	return (1 - r)*InverseRangeCS[i] + r*InverseRangeCS[i + 1];
}

__device__ void WaterQSGetEnergyIndex1(ZFloat ie, int &energyIndex, ZFloat &probability)
{
	ZFloat rle = GlueF(12916953.824948313) * (ie + GlueF(2.0408163265306123e-5));
	if (rle > 0)
	{
		energyIndex = (int)rle;
		probability = rle - energyIndex;
		if (energyIndex >= NQSURFACE_E - 1)
		{
			energyIndex = NQSURFACE_E - 1;
			probability = -1;
		}
	}
	else
	{
		energyIndex = 0;
		probability = -1;
	}
}

__device__ ZFloat WaterQSurface(int anEnergyIndex, ZFloat u)
{
	ZFloat ru = u *(NQSURFACE_Q - 1);
#ifdef USE_SINGLE_PRECISION
	int ju = (int)ru;
	ru -= ju;
	ZFloat a = getWaterQS(ju, anEnergyIndex);
	ZFloat b = getWaterQS(ju + 1, anEnergyIndex);
	return a*(1 - ru) + b*ru;
	//return getWaterQS(ru, anEnergyIndex);
#else
	if (u == 1) return WaterQSData[anEnergyIndex][NQSURFACE_Q - 1];
	int ju = (int)ru;
	ru -= ju;
	ZFloat a = WaterQSData[anEnergyIndex][ju];
	ZFloat b = WaterQSData[anEnergyIndex][ju + 1];
	return a * (1 - ru) + b * ru;
#endif
};

__device__ ZFloat lambdaCalculator(ZFloat x)
{
	//x = linspace(0.01, 1.00000001, 1024);
	if (x >= GlueF(0.01))
	{
		ZFloat r = (x - GlueF(0.01)) * ZFloat((NLOGCACHE - 1) / (1.00000001 - 0.01));
		int i = int(r);
		r -= i;
		return  NegLogX01[i] * (1 - r) + NegLogX01[i + 1] * r;
	}
	else return -GlueF(log)(x);
};

__device__ void randomAzimuth(ZFloat rndno, ZFloat &cphi, ZFloat &sphi)
{
// 	int i = rndno*(NSINECACHE - 1);
// 	cphi = CosCache[i];
// 	sphi = SinCache[i];

	ZFloat phi = 2 * PI * rndno;
	cphi = GlueF(cos)(phi);
	//sphi = (phi < PI ? 1 : -1)*GlueF(sqrt)(1 - cphi*cphi);
	sphi = GlueF(sin)(phi);
}

#endif