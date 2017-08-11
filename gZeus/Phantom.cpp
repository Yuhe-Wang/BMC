
#include "Phantom.h"
#include "../SourceHead/SourceHead.h"
// note that the simulation should be done in the cuboid coordinates system
// so p.x, p.y, p.z may need transfer from the source samples
LogProvider Log;

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
			matid.resize(NX, NY, NZ);
			if (!BF.get("matid", matid.getP())) matid = 278; //default water
			else //need to convert to relative density first
			{
				PENELOPE_MAT penMat;
				if (!penMat.load()) exitApp("Cannot open PENELOPE_densities");
				int LEN = NX*NY*NZ;
				for (int i = 0; i < LEN; ++i) ph.a(i) /= penMat.density(matid.a(i));
			}
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
				matid.resize(NX, NY, NZ);
				matid = 278; //default water
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
		if (NX <= 0 || NY <= 0 || NZ <= 0) exitApp("Incorrect voxel dimensions!");
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
		ph = (dc.size() >= 7 ? SFloat(dc[6]) : SFloat(1)); //1 means use the default density of the corresponding material
		matid.resize(NX, NY, NZ);
		matid = (dc.size() >= 8 ? short(dc[7]) : short(278));

		while (cf->getValue("matrix", dc, true)) //fill inside with different density
		{
			int snx = int(dc[0]), sny = int(dc[1]), snz = int(dc[2]);
			if (snx > NX || sny > NY || snz > NZ) exitApp("Incorrect voxel dimensions!");
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
						{
							ph.a(pix, piy, piz) = (dc.size() >= 7 ? SFloat(dc[6]) : SFloat(1)); //1 means use the default density of the corresponding material
							matid.a(pix, piy, piz) = (dc.size() >= 8 ? short(dc[7]) : short(278)); //default water
						}
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
	
	if (NX*NY*NZ> PhantomBufferSize) exitApp("the voxel number is too big. Please modify PhantomBufferSize in source code");

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

	NMAT = getMatNum(); //get the number of materials
	
	//Note here we get relative density to each corresponding material
	if (NMAT==1) getMaxDensity(); //only one material =>  get max relative density and check if it's uniform
	//else need extra calculation for Woodcock tracking
	dose.resize(NX, NY, NZ);
	dose.setValue(0);
	uncertainty = dose;
	
	Log("\nIt costs %f seconds to load phantom ", rc.stop());
}

void Phantom::output(const char* fname, int FormatType) //output the dose counter
{
	int LEN = getVoxelNum();
	//may trim the dose just based on the density
	for (int i = 0; i < LEN; ++i)
	{
		if (ph.a(i) < trimDensityThreshold) dose.a(i) = 0;
	}

	if (0 == FormatType)
	{
		BinaryFile BF;
		// usual output items
		BF.push("NX", &NX, sizeof(int));
		BF.push("NY", &NY, sizeof(int));
		BF.push("NZ", &NZ, sizeof(int));
		BF.push("DX", &DX, sizeof(double));
		BF.push("DY", &DY, sizeof(double));
		BF.push("DZ", &DZ, sizeof(double));
		BF.push("phantom", ph.getP(), sizeof(SFloat)*LEN);
		BF.push("matid", matid.getP(), sizeof(short)*LEN);
		BF.push("dose", dose.getP(), sizeof(SFloat)*LEN);
		BF.push("uncertainty", uncertainty.getP(), sizeof(SFloat)*LEN);

		BF.push("Hist", &Hist, sizeof(double));
		BF.push("nBatch", &nBatch, sizeof(int));
		BF.push("Bs", &Bs, sizeof(double));
		BF.push("Bx", &Bx, sizeof(double));
		BF.push("By", &By, sizeof(double));
		BF.push("Bz", &Bz, sizeof(double));
		BF.push("xo", &xo, sizeof(double)); //tell DoseViewer where is the ISO center
		BF.push("yo", &yo, sizeof(double));
		BF.push("zo", &zo, sizeof(double));
		BF.push("COPX", &COPX, sizeof(double));//keep original DICOM coordinate(may load extra DICOM data later)
		BF.push("COPY", &COPY, sizeof(double));
		BF.push("COPZ", &COPZ, sizeof(double));
		
		BF.push("prescriptionDose", &prescriptionDose, sizeof(double));
		BF.push("treatmentFraction", &treatmentFraction, sizeof(int));
		BF.push("trimDensityThreshold", &trimDensityThreshold, sizeof(double));
		BF.push("rngSeed", &rngSeed, sizeof(int)); // used for continuous simulation

		//optional output; just for viewer, don't read it in function readDose()
		BF.push("uniform", &uniform, sizeof(int));
		BF.push("NMAT", &NMAT, sizeof(int));

		if (!BF.write(fname)) exitApp("Cannot create the dose file. Please check the file name or writing access.");
	}
	else //ViewRay's format
	{
		VIEWRAY_FORMAT VF;
		VF.nx = NX;
		VF.ny = NY;
		VF.nz = NZ;
		VF.dx = DX;
		VF.dy = DY;
		VF.dz = DZ;
		VF.offset_x = COPX + DX / 2;
		VF.offset_y = COPY + DY / 2;
		VF.offset_z = COPZ + DZ / 2;
		VF.m = dose;
		VF.err = uncertainty;
		VF.Marjor_Header = 1;
		VF.Minor_Header = 0;
		VF.write(fname);
	}
	

	Log("Accumulated history number in dose file = %.0f", Hist);
	Log("The dose file named <%s> has been generated ^_^\n\n", fname);
}

void Phantom::getBinaryFile(BinaryFile& BF)
{
	int LEN = getVoxelNum();
	BF.push("NX", &NX, sizeof(int));
	BF.push("NY", &NY, sizeof(int));
	BF.push("NZ", &NZ, sizeof(int));
	BF.push("DX", &DX, sizeof(double));
	BF.push("DY", &DY, sizeof(double));
	BF.push("DZ", &DZ, sizeof(double));
	BF.push("phantom", ph.getP(), sizeof(SFloat)*LEN);
	BF.push("matid", matid.getP(), sizeof(short)*LEN);
	BF.push("dose", dose.getP(), sizeof(SFloat)*LEN);
	BF.push("uncertainty", uncertainty.getP(), sizeof(SFloat)*LEN);

	BF.push("Hist", &Hist, sizeof(double));
	BF.push("nBatch", &nBatch, sizeof(int));
	BF.push("Bs", &Bs, sizeof(double));
	BF.push("Bx", &Bx, sizeof(double));
	BF.push("By", &By, sizeof(double));
	BF.push("Bz", &Bz, sizeof(double));
	BF.push("xo", &xo, sizeof(double)); //tell DoseViewer where is the ISO center
	BF.push("yo", &yo, sizeof(double));
	BF.push("zo", &zo, sizeof(double));
	BF.push("COPX", &COPX, sizeof(double));//keep original DICOM coordinate(may load extra DICOM data later)
	BF.push("COPY", &COPY, sizeof(double));
	BF.push("COPZ", &COPZ, sizeof(double));

	BF.push("prescriptionDose", &prescriptionDose, sizeof(double));
	BF.push("treatmentFraction", &treatmentFraction, sizeof(int));
	BF.push("trimDensityThreshold", &trimDensityThreshold, sizeof(double));
	BF.push("rngSeed", &rngSeed, sizeof(int)); // used for continuous simulation

	//optional output; just for viewer, don't read it in function readDose()
	BF.push("uniform", &uniform, sizeof(int));
	BF.push("NMAT", &NMAT, sizeof(int));
}

bool Phantom::previousDose(const char* fname) // usually call it after loadPhantom() to continue previous simulation
{
	BinaryFile BF;
	if (!BF.read(fname)) return false; //may be not the right format
	int nx, ny, nz;
	double dx, dy, dz;
	BF.get("NX", &nx);
	BF.get("NY", &ny);
	BF.get("NZ", &nz);
	BF.get("DX", &dx);
	BF.get("DY", &dy);
	BF.get("DZ", &dz);
	if (nx != NX || ny != NY || nz != NZ || dx != DX || dy != DY || dz != DZ) return false;

	
	if (!BF.get("dose", dose.getP())) return false;
	if (!BF.get("uncertainty", uncertainty.getP())) return false;
	if (!BF.get("Hist", &Hist)) return false;
	if (!BF.get("nBatch", &nBatch)) return false;
	if (!BF.get("rngSeed", &rngSeed)) return false;

	double invNorm = (DX*DY*DZ)/(1.602e-16*1.0889e15 * SourceHead_BeamOnTime());
	int NVoxel = getVoxelNum();
	for (int i = 0; i < NVoxel; ++i)
	{
		double D = Hist*dose[i];//previous total dose
		double U = (nBatch*uncertainty[i] * uncertainty[i] + 1) * D * D / nBatch; //previous sum dose^2
		double k = ph[i] * invNorm;
		dose.a(i) = SFloat(D*k); //DE
		uncertainty.a(i) = SFloat(U*k*k); //DE^2
	}
	return true;
}

// double Phantom::addDose(SFloat* d, SFloat* u, int NB, double fH, double aNorm) //return average uncertainty of 50% dose region
// {
// 	SFloat norm = SFloat(1.602e-16*aNorm / (DX*DY*DZ));
// 	int NVoxel = getVoxelNum();
// 	int CB = nBatch + NB; //current total batch number
// 	double CH = Hist + fH; //current history number
// 	SFloat D = 0, U = 0;
// 	SFloat maxDose = 0;
// 	for (int i = 0; i < NVoxel; ++i)
// 	{
// 		D = SFloat(Hist*dose[i]);//previous total dose
// 		if (nBatch) U = (nBatch*uncertainty[i] * uncertainty[i] + 1) * D * D / nBatch; //previous sum dose^2
// 		else U = 0;
// 
// 		//normalize the additional dose and uncertainty
// 		d[i] *= norm;
// 		u[i] *= norm*norm;
// 		if (ph[i] > 0)
// 		{
// 			SFloat dNorm = 1.0f / ph[i];
// 			d[i] *= dNorm;
// 			u[i] *= dNorm*dNorm;
// 			D += d[i]; //current total dose
// 			U += u[i]; //current sum dose^2
// 		}
// 		else
// 		{
// 			D = 0;
// 			U = 0;
// 		}
// 
// 
// 		dose.a(i) = SFloat(D / CH); //current average dose per history
// 		if (D) uncertainty.a(i) = SFloat(sqrt(U / (D*D) - 1.0 / CB)); //current standard deviation % based on batches
// 		else uncertainty.a(i) = -1;
// 
// 		if (dose.a(i) > maxDose) maxDose = dose.a(i);
// 	}
// 	//update batch number and history number
// 	nBatch = CB;
// 	Hist = CH;
// 
// 	//find out the uncertainty of 50% dose region
// 	int n50 = 0;
// 	double errN50 = 0;
// 	for (int i = 0; i < NVoxel; ++i)
// 	{
// 		if (dose[i] > 0.5*maxDose)
// 		{
// 			++n50;
// 			errN50 += uncertainty[i] * uncertainty[i];
// 		}
// 	}
// 	double ave_err = sqrt(errN50 / n50);
// 	Log("Total batch number = %d", CB);
// 	Log("D>50%%*Dmax region occupies %.1f%% all voxels, with average uncertainty = %.1f%%", 100.0*double(n50) / NVoxel, ave_err*100.0);
// 	return ave_err;
// }

double Phantom::addDose(SFloat* d, SFloat* u, int NB, double fH, double aNorm,double threshold) //return average uncertainty of 50% dose region
{
	SFloat norm = SFloat(1.602e-16*aNorm / (DX*DY*DZ));
	int NVoxel = getVoxelNum();
	SFloat maxDose = 0;

	nBatch += NB;
	Hist += fH;
	for (int i = 0; i < NVoxel; ++i)
	{
		SFloat D = d[i];
		SFloat U = u[i];
		if (ph[i] > 0)
		{
			SFloat dNorm = norm / ph[i];
			D *= dNorm;//current total dose
			U *= dNorm*dNorm; //current sum dose^2
		}
		else
		{
			D = 0;
			U = 0;
		}

		dose.a(i) = SFloat(D / Hist); //current average dose per history
		if (D > 0) uncertainty.a(i) = SFloat(sqrt(U / (D*D) - 1.0 / nBatch)); //current standard deviation % based on batches
		else uncertainty.a(i) = -1;

		if (dose.a(i) > maxDose) maxDose = dose.a(i);
	}
	
	//find out the uncertainty of 50% dose region
	int n50 = 0;
	double errN50 = 0;
	for (int i = 0; i < NVoxel; ++i)
	{
		if (dose[i] > threshold*maxDose)
		{
			++n50;
			errN50 += uncertainty[i] * uncertainty[i];
		}
	}
	double ave_err = sqrt(errN50 / n50);
	Log("\nTotal batch number = %d", nBatch);
	Log("D>%d%%*Dmax region occupies %.1f%% all voxels, with average uncertainty = %.1f%%", int(100*threshold), 100.0*double(n50) / NVoxel, ave_err*100.0);
	return ave_err;
}

double Phantom::peekUncertainty(SFloat* d, SFloat* u, int NB, double threshold)
{
	if (NB == 1) NB = 2;
	//assume NB is the total batch for one GPU work thread
	int NVoxel = getVoxelNum();
	SFloat maxDose = 0;
	for (int i = 0; i < NVoxel; ++i)
	{
		if (u[i] > 0) uncertainty.a(i) = SFloat(sqrt(u[i] / (d[i] * d[i]) - 1.0 / NB));
		else uncertainty.a(i) = -1;

		if (d[i] > maxDose) maxDose = d[i];
	}

	//find out the uncertainty of the above threshold dose region
	int n = 0;
	double err = 0;
	threshold *= maxDose; //relative to absolute value
	for (int i = 0; i < NVoxel; ++i)
	{
		if (d[i] > threshold)
		{
			++n;
			err += uncertainty[i] * uncertainty[i];
		}
	}
	//return sqrt(err / n);
	return err / n;
}

bool Phantom::lineInPhantom(Particle& p)
{
	//first translate the coordinate system
	p.x -= xo;
	p.y -= yo;
	p.z -= zo;
	//assuming the incident particle go straight line, judge if it will enter the phantom.
	//if true, give the interaction position
	const double Delta = 1e-5;
	if (p.x < 0 || p.x >= LX || p.y < 0 || p.y >= LY || p.z < 0 || p.z >= LZ) //initial position lays outside the phantom
	{
		if (p.x < 0)
		{
			//now it's outside the phantom	
			if (p.u > 0)
			{
				double t = p.x / p.u;
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
				double t = (LX - Delta - p.x) / p.u;
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
				double t = p.y / p.v;
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
				double t = (LY - Delta - p.y) / p.v;
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
				double t = p.z / p.w;
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
				double t = (LZ - Delta - p.z) / p.w;
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