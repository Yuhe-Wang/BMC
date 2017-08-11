#pragma once
#ifndef _WATERCS_H_
#define _WATERCS_H_

ZFloat PhotonCS[NLAMDAPHOTON];
ZFloat ComptonCS[NLAMDACOMPTON];
ZFloat PairCS[NLAMDAPAIR];
ZFloat ScreenCS[NSCREENINGPARAMETER];
ZFloat RangeCS[NRANGE];
ZFloat InverseRangeCS[NINVERSERANGE];
ZFloat WaterQSData[NQSURFACE_E][NQSURFACE_Q];

bool initWaterCS(const char* filename)
{
	FILE* fp = fopen(filename, "rb");
	if (NULL == fp) return false;

	
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
		PhotonCS[i] = ZFloat(cs);
	}

	//Compton
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		ComptonCS[i] = ZFloat(cs);
	}

	//pair production
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		PairCS[i] = ZFloat(cs);
	}

	//screening
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		ScreenCS[i] = ZFloat(cs);
	}

	//range
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		RangeCS[i] = ZFloat(cs);
	}

	//inverse range
	fread(&bIndex, sizeof(double), 1, fp);
	fread(&aIndex, sizeof(double), 1, fp);
	fread(&n, sizeof(int), 1, fp);
	for (int i = 0; i < n; ++i)
	{
		fread(&cs, sizeof(double), 1, fp);
		InverseRangeCS[i] = ZFloat(cs);
	}

	//qsurface
	double fq, ieq0;
	int ne, nuq;
	fread(&fq, sizeof(double), 1, fp);
	fread(&ieq0, sizeof(double), 1, fp);
	fread(&ne, sizeof(int), 1, fp);
	fread(&nuq, sizeof(int), 1, fp);
	for (int ie = 0; ie < ne; ++ie)
		for (int iq = 0; iq < nuq; ++iq)
		{
			fread(&cs, sizeof(double), 1, fp);
			WaterQSData[ie][iq] = ZFloat(cs);
		}


	fclose(fp);
	return true;
}

ZFloat WaterPhotonCS(ZFloat E)
{
	ZFloat r = (E - GlueF(4.9e4))*GlueF(0.0014107512060647829);
	int i = int(r);
	r -= i;
	return (1 - r)*PhotonCS[i] + r*PhotonCS[i + 1];
}

ZFloat WaterComptonCS(ZFloat E)
{
	ZFloat r = (E - GlueF(4.9e4))*GlueF(0.0001757408683666437);
	int i = int(r);
	r -= i;
	return (1 - r)*ComptonCS[i] + r*ComptonCS[i + 1];
}

ZFloat WaterPairCS(ZFloat E)
{
	ZFloat r = (E - GlueF(1.024e6))*GlueF(0.0021491596638655462);
	int i = int(r);
	r -= i;
	return (1 - r)*PairCS[i] + r*PairCS[i + 1];
}

ZFloat WaterScreenCS(ZFloat E)
{
	ZFloat r = (E + GlueF(2.0408163265306123e-005)) *GlueF(25884562.370778777);
	int i = int(r);
	r -= i;
	return (1 - r)*ScreenCS[i] + r*ScreenCS[i + 1];
}

ZFloat WaterRangeCS(ZFloat E)
{
	ZFloat r = (E - GlueF(10e3))*GlueF(0.00068657718120805373);
	int i = int(r);
	r -= i;
	return (1 - r)*RangeCS[i] + r*RangeCS[i + 1];
}

ZFloat WaterInverseRangeCS(ZFloat E)
{
	ZFloat r = E *GlueF(1447.3007247284104);
	int i = int(r);
	r -= i;
	return (1 - r)*InverseRangeCS[i] + r*InverseRangeCS[i + 1];
}

void WaterQSGetEnergyIndex1(ZFloat ie, int &energyIndex, ZFloat &probability)
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

ZFloat WaterQSurface(int anEnergyIndex, ZFloat u)
{
	int nuqp = NQSURFACE_Q - 1;
	ZFloat ru = u * nuqp;
	int ju = (int)ru;
	if (ju > nuqp - 1) ju = nuqp - 1;
	ru -= ju;
	return WaterQSData[anEnergyIndex][ju] * (1 - ru) + WaterQSData[anEnergyIndex][ju + 1] * ru;
};

#endif