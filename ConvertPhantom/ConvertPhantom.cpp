#include "../Tools/Tools.h"
#ifdef WIN32
//windows specific, hide the cmd window
#pragma comment(linker, "/subsystem:windows /entry:mainCRTStartup")
#include "Shellapi.h" // to call ShellExecute()
#endif

void exitApp(const char *inf)
{
	printf("fatal error: %s", inf);
	getchar();
	exit(-1);
}

int main(int argc, char *argv[])
{
	//open the config file and get phantom config
	if (argc < 2)
	{
		printf("usage: %s newConfigFileName\n\npress enter to exit...", argv[0]);
		getchar();
		return 0;
	}

	ConfigFile cff;
	if (!cff.parse(argv[1])) exitApp("cannot parse the config file");
	ConfigFile* cf = cff.getBlock("PHANTOM");

	double isoPX = 0, isoPY = 0, isoPZ = 0;//iso center position in patient coordinate
	double DX, DY, DZ;
	double COPX, COPY, COPZ;
	double xo, yo, zo;
	int NX, NY, NZ;
	ArrayMgr<SFloat> ph;

	vector<double> dc;
	cf->getValue("DX DY DZ", dc);
	DX = dc[0]; DY = dc[1]; DZ = dc[2];

	cf->resetSearchIndex();
	cf->getValue("matrix", dc, true);
	NX = int(dc[0]); NY = int(dc[1]); NZ = int(dc[2]);
	if (NX <= 0 || NY <= 0 || NZ <= 0 || dc.size() != 7) exitApp("Incorrect custom phantom definition");
	//the origin of the patient coordinates lies at the center of the cuboid,
	COPX = -NX*DX / 2;
	COPY = -NY*DY / 2;
	COPZ = -NZ*DZ / 2;
	isoPX = -dc[3];
	isoPY = -dc[4];
	isoPZ = -dc[5];
	xo = COPX - isoPX;
	yo = COPY - isoPY;
	zo = COPZ - isoPZ;
	
	ph.resize(NX, NY, NZ);
	ph.setValue(SFloat(dc[6]));

	while (cf->getValue("matrix", dc, true)) //fill inside with different density
	{
		int snx = int(dc[0]), sny = int(dc[1]), snz = int(dc[2]);
		if (snx > NX || sny > NY || snz > NZ || dc.size() != 7) exitApp("Incorrect custom phantom definition");
		//need to find the origin index of the new cuboid
		int cix = lround((dc[3] - snx*DX*0.5 - xo) / DX);
		int ciy = lround((dc[4] - sny*DY*0.5 - yo) / DY);
		int ciz = lround((dc[5] - snz*DZ*0.5 - zo) / DZ);

		int pix, piy, piz;
		for (int i = 0; i < snx; ++i)
			for (int j = 0; j < sny; ++j)
				for (int k = 0; k < snz; ++k)
				{
					pix = i + cix;
					piy = j + ciy;
					piz = k + ciz;
					if (0 <= pix && pix < NX && 0 <= piy && piy < NY && 0 <= piz && piz < NZ)
						ph.a(pix, piy, piz) = SFloat(dc[6]);
				}
	}

	//may adjust the isoPX, isoPY, isoPZ
	double SSD = 0;
	if (cf->getValue("SSD", SSD)) isoPY = SAD - SSD - NY*DY / 2;
	vector<double> centerXZ;
	if (cf->getValue("centerXZ", centerXZ))
	{
		isoPX = -centerXZ[0];
		isoPZ = -centerXZ[1];
	}

	//update xo,yo,zo
	xo = COPX - isoPX;
	yo = COPY - isoPY;
	zo = COPZ - isoPZ;

	//overwrite couch density into above phantom
	string couchFile;
	if (cf->getValue("couchConfig", couchFile)) //has the density information
	{
		FILE* fp = fopen(couchFile.c_str(), "rb");
		if (NULL == fp) exitApp("Cannot open the couch density file you assigned");
		int CNX, CNY, CNZ;
		float CDX, CDY, CDZ;
		float offset_x, offset_y, offset_z;
		fread(&CNX, sizeof(int), 1, fp);
		fread(&CNY, sizeof(int), 1, fp);
		fread(&CNZ, sizeof(int), 1, fp);
		fread(&CDX, sizeof(float), 1, fp);
		fread(&CDY, sizeof(float), 1, fp);
		fread(&CDZ, sizeof(float), 1, fp);
		if (fabs(CDX - DX) > 1e-5 || fabs(CDY - DY) > 1e-5 || fabs(CDZ - DZ) > 1e-5) exitApp("Inconsistent couch voxel size to the phantom");
		fread(&offset_x, sizeof(float), 1, fp);
		fread(&offset_y, sizeof(float), 1, fp);
		fread(&offset_z, sizeof(float), 1, fp);
		ArrayMgr<SFloat> cMat(CNX, CNY, CNZ);
		float elem = 0;
		for (int k = 0; k < CNZ; ++k)
			for (int j = 0; j < CNY; ++j)
				for (int i = 0; i < CNX; ++i)
				{
					fread(&elem, sizeof(float), 1, fp);
					cMat.a(i, j, k) = SFloat(elem);
				}
		fclose(fp);

		//find the couch's location
		vector<double> pos;
		cf->getValue("couchXYZ", pos);
		if (pos.size() != 3) exitApp("Incorrect couch position!");
		//fit the couch density matrix to above phantom

		//first locate the index of couch's origin
		int cix = lround((pos[0] - CNX*CDX / 2 - xo) / DX);
		int ciy = lround((pos[1] - yo) / DY);
		int ciz = lround((pos[2] - CNZ*CDZ / 2 - zo) / DZ);
		//now the indx(couch)+cix == indx(phantom), and so forth
		//we may need to extend z direction of couch to fit above phantom
		int pix, piy, piz;
		if (ciz > 0) //need a face copy
		{
			for (int piz = 0; piz < ciz; ++piz)
			{
				for (int ix = 0; ix < CNX; ++ix)
					for (int iy = 0; iy < CNY; ++iy)
					{
						pix = ix + cix;
						piy = iy + ciy;
						if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY)
							ph.a(pix, piy, piz) = cMat.a(ix, iy, 0);
					}
			}
		}
		if (ciz + CNZ < NZ) //need another face copy
		{
			for (int piz = ciz + CNZ; piz < NZ; ++piz)
			{
				for (int ix = 0; ix < CNX; ++ix)
					for (int iy = 0; iy < CNY; ++iy)
					{
						pix = ix + cix;
						piy = iy + ciy;
						if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY)
							ph.a(pix, piy, piz) = cMat.a(ix, iy, CNZ - 1);
					}
			}
		}
		//do the major replacement
		for (int ix = 0; ix < CNX; ++ix)
			for (int iy = 0; iy < CNY; ++iy)
				for (int iz = 0; iz < CNZ; ++iz)
				{
					pix = ix + cix;
					piy = iy + ciy;
					piz = iz + ciz;
					if (pix >= 0 && pix < NX &&piy >= 0 && piy < NY && piz >= 0 && piz < NZ)
						ph.a(pix, piy, piz) = cMat.a(ix, iy, iz);
				}
	}

	//may change the couch density
	double couchDensityFactor = 1.0;
	if (cf->getValue("couch density factor", couchDensityFactor))
	{
		int ty = NY - 45;//scan for the position of ViewRa's nominal couch
		for (; ty >= 0; --ty) //the couch's height takes about 48 voxels
		{
			if (fabs(ph.a(NX / 2, ty, NZ / 2) - 1.819) < 1e-4 && fabs(ph.a(NX / 2, ty + 1, NZ / 2) - 0.0778) < 1e-4) break;
		}
		if (ty != -1) //we got the couch top position
		{
			for (int iy = ty; iy <= ty + 48 && iy < NY; ++iy) //the couch's height takes about 48 voxels
			{
				for (int ix = 0; ix < NX; ++ix)
					for (int iz = 0; iz < NZ; ++iz)
					{
						ph.a(ix, iy, iz) *= SFloat(couchDensityFactor);
					}
			}
		}
		else printf("Warning: failed to modify the couch's density. Make sure this couch exist.");
	}
	else printf("Warning: cannot get couch density factor in config file; will not modify the couch density!");

	//
	VIEWRAY_FORMAT vf;
	vf.nx = NX;
	vf.ny = NY;
	vf.nz = NZ;
	vf.dx = DX;
	vf.dy = DY;
	vf.dz = DZ;
	vf.offset_x = COPX + DX / 2;
	vf.offset_y = COPY + DY / 2;
	vf.offset_y = COPZ + DZ / 2;
	vf.m = ph;
	fs::path cp(fs::current_path());
	fs::current_path(fs::path(argv[1]).parent_path());//change current path
	vf.write("mfQAData.red");
	//windows specific
#ifdef WIN32
	ShellExecute(NULL, NULL, "mfQAData.red", NULL, NULL, SW_SHOW);
#endif
	if (argc >= 3)
	{
		pauseSeconds(3);
#ifdef WIN32
		system("del /q mfQAData.red");
#else
        system("rm /f mfQAData.red");
#endif
	}
	fs::current_path(cp);//change current path
	return 0;
}