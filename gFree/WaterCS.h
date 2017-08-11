#pragma once
#ifndef _WATERCS_H_
#define _WATERCS_H_

__constant__ ZFloat PhotonCS[NLAMDAPHOTON];
__constant__ ZFloat ComptonCS[NLAMDACOMPTON];
__constant__ ZFloat PairCS[NLAMDAPAIR];
__constant__ ZFloat ScreenCS[NSCREENINGPARAMETER];
__constant__ ZFloat RangeCS[NRANGE];
__constant__ ZFloat InverseRangeCS[NINVERSERANGE];
Texture2DFloat(WaterQS)
ZFloat WaterQSData[NQSURFACE_E][NQSURFACE_Q];
__constant__ ZFloat NegLogX01[1024]; //to make -log(x), 0<x<1 faster

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
	//ZFloat* h_WaterQSData = new ZFloat[NQSURFACE_E*NQSURFACE_Q];
	ZFloat* h_NegLogX01 = new ZFloat[1024];

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
			fwrite(&cs, sizeof(double), 1, fp);
			WaterQSData[ie][iq] = ZFloat(cs);
		}
	fclose(fp);

	//x = linspace(0.01, 1.00000001, 1024);
	double dx = (1.00000001 - 0.01) / 1023;
	for (int i = 0; i < 1024; ++i)
	{
		h_NegLogX01[i] = -ZFloat(log(0.01 + i*dx));
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
		//cudaErrorCheck(cudaMemcpyToSymbol(WaterQSData, h_WaterQSData, sizeof(ZFloat)*NQSURFACE_E*NQSURFACE_Q), "copy WaterQSData");
		cudaErrorCheck(cudaMemcpyToSymbol(NegLogX01, h_NegLogX01, sizeof(ZFloat)*1024), "copy NegLogX01");
	}

	delete[] h_PhotonCS;
	delete[] h_ComptonCS;
	delete[] h_PairCS;
	delete[] h_ScreenCS;
	delete[] h_RangeCS;
	delete[] h_InverseRangeCS;
	//delete[] h_WaterQSData;
	delete[] h_NegLogX01;
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
	return (1 - r)*InverseRangeCS[i] + r*InverseRangeCS[i + 1];
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

__device__ __forceinline__ void WaterQSGetEnergyIndex1(ZFloat ie, int &energyIndex, ZFloat &probability)
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

__device__ __forceinline__ ZFloat WaterQSurface(int anEnergyIndex, ZFloat u)
{
	int nuqp = NQSURFACE_Q - 1;
	ZFloat ru = u * nuqp;
// 	int ju = (int)ru;
// 	if (ju > nuqp - 1) ju = nuqp - 1;
// 	ru -= ju;
	//return WaterQSData[anEnergyIndex][ju] * (1 - ru) + WaterQSData[anEnergyIndex][ju + 1] * ru;
	return getWaterQS(anEnergyIndex, ru);
};

__device__ __forceinline__ ZFloat lambdaCalculator(ZFloat x)
{
	//x = linspace(0.01, 1.00000001, 1024);
	if (x >= GlueF(0.01))
	{
		ZFloat r = (x - GlueF(0.01)) *GlueF(1033.333322895623);
		int i = int(r);
		r -= i;
		return  NegLogX01[i] * (1 - r) + NegLogX01[i + 1] * r;
	}
	else return -GlueF(log)(x);
};

#endif