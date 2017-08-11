
#include "../Tools/Tools.h"

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		printf("usage: %s zeusDoseFilename\n\npress enter to exit...", argv[0]);
		getchar();
		return 0;
	}
	FILE *fp = fopen(argv[1], "rb");
	if (NULL == fp)
	{
		printf("cannot open zeus dose file!");
		getchar();
		return -1;
	}
	int version = 0;
	fread(&version, sizeof(int), 1, fp);
	fread(&version, sizeof(int), 1, fp);
	int NX, NY, NZ;
	double DX, DY, DZ;
	double OX, OY, OZ;
	fread(&NX, sizeof(int), 1, fp);
	fread(&DX, sizeof(double), 1, fp);
	fread(&OX, sizeof(double), 1, fp);

	fread(&NY, sizeof(int), 1, fp);
	fread(&DY, sizeof(double), 1, fp);
	fread(&OY, sizeof(double), 1, fp);

	fread(&NZ, sizeof(int), 1, fp);
	fread(&DZ, sizeof(double), 1, fp);
	fread(&OZ, sizeof(double), 1, fp);

	int ncount = NX*NY*NZ;
	float *dose = new float[ncount];
	float *uncertainty = new float[ncount];

	ArrayMgr<SFloat> phantom(NX, NY, NZ);
	phantom.setValue(1.0);
	ArrayMgr<SFloat> pdose(NX, NY, NZ);
	ArrayMgr<SFloat> puncertainty(NX, NY, NZ);
	
	fread(dose, sizeof(float), ncount, fp);
	int uncertainRead = (int)fread(uncertainty, sizeof(float), ncount, fp);
	
	for (int ix = 0; ix < NX; ++ix)
	{
		for (int iy = 0; iy < NY; ++iy)
		{
			for (int iz = 0; iz < NZ; ++iz)
			{
				pdose.a(ix, iy, iz) = dose[ix + NX*iy + NX*NY*iz];
			}
		}
	}
	if (uncertainRead)
	{
		for (int ix = 0; ix < NX; ++ix)
		{
			for (int iy = 0; iy < NY; ++iy)
			{
				for (int iz = 0; iz < NZ; ++iz)
				{
					puncertainty.a(ix, iy, iz) = uncertainty[ix + NX*iy + NX*NY*iz] / dose[ix + NX*iy + NX*NY*iz];
				}
			}
		}
	}

	delete[] dose;
	delete[] uncertainty;
	fclose(fp);


	//let out put the dose file in my format
	BinaryFile BF;
	double CT_kVp = 110;
	BF.push("CT_kVp", &CT_kVp, sizeof(double));
	int uniform = 1;
	BF.push("uniform", &uniform, sizeof(int));
	double origin = -DX*NX / 2;
	BF.push("xo", &origin, sizeof(double));
	origin = -DY*NY / 2;
	BF.push("yo", &origin, sizeof(double));
	origin = -DZ*NZ / 2;
	BF.push("zo", &origin, sizeof(double));
	BF.push("NX", &NX, sizeof(int));
	BF.push("NY", &NY, sizeof(int));
	BF.push("NZ", &NZ, sizeof(int));
	BF.push("DX", &DX, sizeof(double));
	BF.push("DY", &DY, sizeof(double));
	BF.push("DZ", &DZ, sizeof(double));
	
	BF.push("phantom", (void*)phantom.getP(), sizeof(SFloat)*NX*NY*NZ);
	BF.push("dose", (void*)pdose.getP(), sizeof(SFloat)*NX*NY*NZ);
	if (uncertainRead)
	{
		BF.push("uncertainty", (void*)puncertainty.getP(), sizeof(SFloat)*NX*NY*NZ);
	}
	if (argc == 3) BF.write(argv[2]);
	else BF.write("myDose.dose");

	return 0;

}