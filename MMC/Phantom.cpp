
#include "Phantom.h"
#include "../SourceHead/SourceHead.h"
// note that the simulation should be done in the cuboid coordinates system
// so p.x, p.y, p.z may need transfer from the source samples

void Phantom::loadPhantom(ConfigFile *cf)
{
	RunTimeCounter rc;
	
	//read NX,NY,NZ; DX,DY,DZ(LX, LY, LZ); density matrix
	bool b_fromFile = true;
	string fname;
	if (!cf->getValue("phantomName", fname)) b_fromFile = false; //no phantom file specified

	cf->getValue("trim dose threshold", trimDensityThreshold);
	if (trimDensityThreshold < 0) trimDensityThreshold = 0;

	double isoPX = 0, isoPY = 0, isoPZ = 0;//iso center position in patient coordinate

	//first try to load from the file
	if (b_fromFile)
	{
		BinaryFile BF;
		if (BF.read(fname.c_str()))//try to load infomation from *.phtm file
		{
			BF.get("CT_kVp", &CT_kVp); 
			BF.get("NX", &NX);
			BF.get("NY", &NY);
			BF.get("NZ", &NZ);
			BF.get("DX", &DX);
			BF.get("DY", &DY);
			BF.get("DZ", &DZ);
			BF.get("COPX", &COPX);
			BF.get("COPY", &COPY);
			BF.get("COPZ", &COPZ);
			ph.resize(NX, NY, NZ);
			BF.get("phantom", ph.getP());
		}
		else //try to load RED file then
		{
			VIEWRAY_FORMAT VF;
			if (VF.read(fname.c_str()))
			{
				NX = VF.nx;
				NY = VF.ny;
				NZ = VF.nz;
				DX = VF.dx;
				DY = VF.dy;
				DZ = VF.dz;
				COPX = VF.offset_x - DX / 2;
				COPY = VF.offset_y - DY / 2;
				COPZ = VF.offset_z - DZ / 2;
				ph = VF.m;
			}
			else
			{
				Log("Warning: Cannot load phantom from file. Will try to use customized phantom in config file");
				b_fromFile = false;
			}
		}
	}
		
	if(!b_fromFile)//need to get regular phantom from user's definition. Patient's coordinate center is the cuboid center
	{
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
							ph.a(pix, piy, piz)= SFloat(dc[6]);
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
	}
	else
	{
		SourceHead_GetIsoCenter(&isoPX, &isoPY, &isoPZ);	
	}
	//it's correct regardless b_fromFile
	xo = COPX - isoPX;
	yo = COPY - isoPY;
	zo = COPZ - isoPZ;

	LX = NX*DX; LY = NY*DY; LZ = NZ*DZ;

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
						ph.a(ix, iy, iz)*= SFloat(couchDensityFactor);
					}
			}
		}
		else Log("Warning: failed to modify the couch's density. Make sure this couch exist.");
	}
	else Log("Warning: cannot get couch density factor in config file; will not modify the couch density!");
	
	//read configuration of the magnetic field
	cf->getValue("magnetic field strength", Bs); 
	if (Bs <= 0) rf = 0;
	else rf = 100 / (lightSpeed*Bs);

	vector<double> dir;
	cf->getValue("magnetic field direction", dir);
	if (dir.size() == 2) //give two polar angles: theta and phi
	{
		dir[0] *= PI / 180;
		dir[1] *= PI / 180;
		Bx = sin(dir[0])*cos(dir[1]);
		if (fabs(Bx) < 1e-8) Bx = 0;
		By = sin(dir[0])*sin(dir[1]);
		if (fabs(By) < 1e-8) By = 0;
		Bz = cos(dir[0]);
		if (fabs(Bz) < 1e-8) Bz = 0;
	}
	else if (dir.size() == 3) //give unit vector
	{
		if (fabs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] - 1.0) > 1e-8) exitApp("Wrong unit vector!");
		Bx = dir[0];
		By = dir[1];
		Bz = dir[2];
	}
	else exitApp("Wrong magnetic field direction parameters!");

	maxDensity = (SFloat)getMaxDensity(); //get max density and check if it's uniform
	dose.resize(NX, NY, NZ);
	dose.setValue(0);
	uncertainty = dose;

	density = ph.getP();
	Log("\nIt costs %f seconds to init Penelope. ", rc.stop());
}
void Phantom::output(const char* fname) //output the dose counter
{
	int LEN = getNVoxel();
	BinaryFile BF;

	BF.push("CT_kVp", &CT_kVp, sizeof(double));
	BF.push("Hist", &Hist, sizeof(double));
	BF.push("nBatch", &nBatch, sizeof(int));
	BF.push("uniform", &uniform, sizeof(int));
	BF.push("NX", &NX, sizeof(int));
	BF.push("NY", &NY, sizeof(int));
	BF.push("NZ", &NZ, sizeof(int));
	BF.push("DX", &DX, sizeof(double));
	BF.push("DY", &DY, sizeof(double));
	BF.push("DZ", &DZ, sizeof(double));
	BF.push("xo", &xo, sizeof(double)); //tell DoseViewer where is the ISO center
	BF.push("yo", &yo, sizeof(double));
	BF.push("zo", &zo, sizeof(double));
	BF.push("COPX", &COPX, sizeof(double));//keep original DICOM coordinate(may load extra DICOM data later)
	BF.push("COPY", &COPY, sizeof(double));
	BF.push("COPZ", &COPZ, sizeof(double));
	BF.push("Bs", &Bs, sizeof(double));
	BF.push("Bx", &Bx, sizeof(double));
	BF.push("By", &By, sizeof(double));
	BF.push("Bz", &Bz, sizeof(double));
	BF.push("prescriptionDose", &prescriptionDose, sizeof(double));
	BF.push("treatmentFraction", &treatmentFraction, sizeof(int));

	//may trim the dose just based on the density
	int len = NX*NY*NZ;
	for (int i = 0; i < len; ++i)
	{
		if (ph.a(i) < trimDensityThreshold) dose.a(i) = 0;
	}

	BF.push("phantom", ph.getP(), sizeof(float)*LEN);
	BF.push("dose", dose.getP(), sizeof(float)*LEN);
	BF.push("uncertainty", uncertainty.getP(), sizeof(float)*LEN);
	if(!BF.write(fname)) exitApp("Cannot create the dose file. Please check the file name or writing access.");

	Log("Accumulated history number in dose file = %d", Hist);
	Log("The dose file named <%s> has been generated ^_^\n", fname);
}
bool Phantom::readDose(const char* fname)
{
	BinaryFile BF;
	if (!BF.read(fname)) return false;
	BF.get("CT_kVp", &CT_kVp);

	BF.get("Hist", &Hist);
	BF.get("nBatch", &nBatch);
	BF.get("uniform", &uniform);
	BF.get("NX", &NX);
	BF.get("NY", &NY);
	BF.get("NZ", &NZ);
	BF.get("DX", &DX);
	BF.get("DY", &DY);
	BF.get("DZ", &DZ);
	BF.get("xo", &xo); //this can tell DoseViewer where is the ISO center
	BF.get("yo", &yo);
	BF.get("zo", &zo);
	LX = NX*DX; LY = NY*DY; LZ = NZ*DZ;
	BF.get("Bs", &Bs);
	BF.get("Bx", &Bx);
	BF.get("By", &By);
	BF.get("Bz", &Bz);
	if (0 == Bs) rf = 0;
	else rf = 100 / (lightSpeed*Bs);
	ph.resize(NX, NY, NZ);
	BF.get("phantom", ph.getP());

	getMaxDensity(); //get max density and check if it's uniform
	dose.resize(NX, NY, NZ);
	uncertainty.resize(NX, NY, NZ);
	BF.get("dose", dose.getP());
	BF.get("uncertainty", uncertainty.getP());

	return true;
}
