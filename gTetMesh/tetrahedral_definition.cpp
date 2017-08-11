#include "tetrahedral_definition.h"
#undef at
#include "matio.h"
#include <unordered_map>
#include <ctime>
#include "../SourceHead/SourceHead.h"

int Tetrahedral::load(Tetrahedral*& ts, octree_node*& root, const char* fname, ArrayMgr<SFloat>& node, ArrayMgr<int>& elem, ArrayMgr<int>& f2t, vector<double> shift) //return the number of tetrahedral
{
	//RunTimeCounter timeCounter;
	mat_t *matfp;
	matvar_t *matvar;
	matfp = Mat_Open(fname, MAT_ACC_RDWR);
	if (NULL == matfp) {
		fprintf(stderr, "Error opening MAT file \"%s\"!\n", fname);
		return 0;
	}

	matvar = Mat_VarRead(matfp,"node");
	if (matvar == NULL || matvar->rank != 2 || matvar->dims[1] != 3 || matvar->data_type != MAT_T_DOUBLE || matvar->class_type != MAT_C_DOUBLE)
	{
		fprintf(stderr, "Cannot load variable node!\n");
		return 0;
	};
	ArrayMgr<double> temp;
	
	temp.resize((int32_t)matvar->dims[0], (int32_t)matvar->dims[1]);
	temp.copy((double*)matvar->data);
	temp *= 0.1;//change the unit from mm to cm
	Mat_VarFree(matvar);
	temp.cast<SFloat>(node);
	
	//calculate the center coordinates
	SFloat xmin, ymin, zmin;
	SFloat xmax, ymax, zmax;
	xmin = xmax = node.a(0, 0);
	ymin = ymax = node.a(0, 1);
	zmin = zmax = node.a(0, 2);
	int nNode = node.getWidth(1);
	for (int i = 1; i < nNode; ++i)
	{
		xmin = min(xmin, node.a(i, 0));
		xmax = max(xmax, node.a(i, 0));
		ymin = min(ymin, node.a(i, 1));
		ymax = max(ymax, node.a(i, 1));
		zmin = min(zmin, node.a(i, 2));
		zmax = max(zmax, node.a(i, 2));
	}
	SFloat cx = (xmin + xmax) / 2;
	SFloat cy = (ymin + ymax) / 2;
	SFloat cz = (zmin + zmax) / 2;
	//scale and shift the mesh
	SFloat scaleFactor = 1;
	SFloat shiftx = 0, shifty = 0, shiftz = 0;
	if (shift.size() >= 4) scaleFactor = (SFloat)shift[3];
	if (shift.size() >= 3)
	{
		shiftx = (SFloat)shift[0];
		shifty = (SFloat)shift[1];
		shiftz = (SFloat)shift[2];
	}
	for (int i = 1; i < nNode; ++i)
	{
		node.a(i, 0) = (node.a(i, 0) - cx)*scaleFactor + cx + shiftx;
		node.a(i, 1) = (node.a(i, 1) - cy)*scaleFactor + cy + shifty;
		node.a(i, 2) = (node.a(i, 2) - cz)*scaleFactor + cz + shiftz;
	}

	bool b_f2t = true;
	matvar = Mat_VarRead(matfp, "f2t");
	//ArrayMgr<int> f2t;
	
	if (matvar == NULL || matvar->rank != 2 || matvar->dims[1] != 4 || matvar->data_type != MAT_T_DOUBLE || matvar->class_type != MAT_C_DOUBLE)
	{
		//cannot load the variable f2t
		b_f2t = false;
	}
	else
	{
		temp.resize((int32_t)matvar->dims[0], (int32_t)matvar->dims[1]);
		temp.copy((double*)matvar->data);
		temp -= 1; // from MATLAB index to c index
		temp.cast<int>(f2t);
	}

	matvar = Mat_VarRead(matfp, "elem");
	if (matvar == NULL || matvar->rank != 2 || matvar->dims[1] != 5 || matvar->data_type != MAT_T_DOUBLE || matvar->class_type != MAT_C_DOUBLE)
	{
		fprintf(stderr, "Cannot load the variable elem!\n");
		return 0;
	};
	temp.resize((int32_t)matvar->dims[0], (int32_t)matvar->dims[1]);
	temp.copy((double*)matvar->data);
	temp -= 1; // from MATLAB index to c index
	temp.cast<int>(elem);

	int NT = elem.getWidth(1);
	ts = new Tetrahedral[NT];

	vector<triangle_octree> triangles; //surface triangles

	if (b_f2t)
	{
		for (int i = 0; i < NT; ++i)
		{
			ts[i].mid = short(elem.a(i, 4));
			int n0 = elem.a(i, 0);
			int n1 = elem.a(i, 1);
			int n2 = elem.a(i, 2);
			int n3 = elem.a(i, 3);
			Vector v0((ZFloat)node.a(n0, 0), (ZFloat)node.a(n0, 1), (ZFloat)node.a(n0, 2));
			Vector v1((ZFloat)node.a(n1, 0), (ZFloat)node.a(n1, 1), (ZFloat)node.a(n1, 2));
			Vector v2((ZFloat)node.a(n2, 0), (ZFloat)node.a(n2, 1), (ZFloat)node.a(n2, 2));
			Vector v3((ZFloat)node.a(n3, 0), (ZFloat)node.a(n3, 1), (ZFloat)node.a(n3, 2));
			ts[i].face[0].v[0] = v1;
			ts[i].face[0].v[1] = v2;
			ts[i].face[0].v[2] = v3;

			ts[i].face[1].v[0] = v0;
			ts[i].face[1].v[1] = v3;
			ts[i].face[1].v[2] = v2;

			ts[i].face[2].v[0] = v0;
			ts[i].face[2].v[1] = v1;
			ts[i].face[2].v[2] = v3;

			ts[i].face[3].v[0] = v0;
			ts[i].face[3].v[1] = v2;
			ts[i].face[3].v[2] = v1;
			for (char j = 0; j < 4; ++j) ts[i].face[j].calcPluckerCoordinates();
		}
		int NF = f2t.getWidth(1);
		for (int i = 0; i < NF; ++i)
		{
			if (f2t.a(i, 0) == -1)
			{
				printf("Error: a face doesn't belong to any tetrahedral.\n");
				getchar();
				return -1;
			}
			ts[f2t.a(i, 0)].face[f2t.a(i, 1)].iadj = f2t.a(i, 2);
			if (f2t.a(i, 2) == -1) //fv.idt0 & fv.idf0 points to a boundary triangle, we need to extract it
			{
				TriangleFace& f = ts[f2t.a(i, 0)].face[f2t.a(i, 1)];
				triangles.push_back(triangle_octree(f.v[0].x, f.v[0].y, f.v[0].z, f.v[1].x, f.v[1].y, f.v[1].z, f.v[2].x, f.v[2].y, f.v[2].z, f2t.a(i, 0)));
			}
			else ts[f2t.a(i, 2)].face[f2t.a(i, 3)].iadj = f2t.a(i, 0);
		}
	}
	else // no variable f2t available
	{
		//std::unordered_map<string, FaceValue> fmap; //face map
		auto pmap = new std::unordered_map<string, FaceValue>(); //face map
		std::unordered_map<string, FaceValue>& fmap = *pmap;

		//scan the tetrahedral to correct the vertex sequence
		for (int i = 0; i < NT; ++i)
		{
			int n0 = elem.a(i, 0);
			int n1 = elem.a(i, 1);
			int n2 = elem.a(i, 2);
			int n3 = elem.a(i, 3);
			Vector v0((ZFloat)node.a(n0, 0), (ZFloat)node.a(n0, 1), (ZFloat)node.a(n0, 2));
			Vector v1((ZFloat)node.a(n1, 0), (ZFloat)node.a(n1, 1), (ZFloat)node.a(n1, 2));
			Vector v2((ZFloat)node.a(n2, 0), (ZFloat)node.a(n2, 1), (ZFloat)node.a(n2, 2));
			Vector v3((ZFloat)node.a(n3, 0), (ZFloat)node.a(n3, 1), (ZFloat)node.a(n3, 2));

			Vector v01 = v1 - v0;
			Vector v02 = v2 - v0;
			Vector vn = v01&v02;
			double up = vn*(v3 - v0);
			if (up < 0) //need to adjust the vertex sequence of the tetrahedral in the original data
			{
				swapData<int>(elem.a(i, 0), elem.a(i, 1));
				//printf("up<0\n");
			}
		}
	
		for (int i = 0; i < NT; ++i)
		{
			ts[i].mid = short(elem.a(i, 4));
			int n0 = elem.a(i, 0);
			int n1 = elem.a(i, 1);
			int n2 = elem.a(i, 2);
			int n3 = elem.a(i, 3);
			Vector v0((ZFloat)node.a(n0, 0), (ZFloat)node.a(n0, 1), (ZFloat)node.a(n0, 2));
			Vector v1((ZFloat)node.a(n1, 0), (ZFloat)node.a(n1, 1), (ZFloat)node.a(n1, 2));
			Vector v2((ZFloat)node.a(n2, 0), (ZFloat)node.a(n2, 1), (ZFloat)node.a(n2, 2));
			Vector v3((ZFloat)node.a(n3, 0), (ZFloat)node.a(n3, 1), (ZFloat)node.a(n3, 2));

			FaceKey fkey;
			ts[i].face[0].v[0] = v1;
			ts[i].face[0].v[1] = v2;
			ts[i].face[0].v[2] = v3;
			fkey.set(n1, n2, n3);
			fmap[fkey()].set(i, 0);

			ts[i].face[1].v[0] = v0;
			ts[i].face[1].v[1] = v3;
			ts[i].face[1].v[2] = v2;
			fkey.set(n0, n2, n3);
			fmap[fkey()].set(i, 1);

			ts[i].face[2].v[0] = v0;
			ts[i].face[2].v[1] = v1;
			ts[i].face[2].v[2] = v3;
			fkey.set(n1, n0, n3);
			fmap[fkey()].set(i, 2);

			ts[i].face[3].v[0] = v0;
			ts[i].face[3].v[1] = v2;
			ts[i].face[3].v[2] = v1;
			fkey.set(n1, n2, n0);
			fmap[fkey()].set(i, 3);

			for (char j = 0; j < 4; ++j) ts[i].face[j].calcPluckerCoordinates();
		}

		//iterate the face map to set the adjacent tetrahedral index for each face
		//collect those surface triangles at the same time
		int NF = (int)fmap.size();
		f2t.resize(NF, 4);
		int idf = 0;
		for (auto it = fmap.begin(); it != fmap.end(); ++it)
		{
			FaceValue& fv = it->second;
			if (fv.idt0 == -1)
			{
				printf("Error: a face doesn't belong to any tetrahedral.\n");
				getchar();
				return -1;
			}
			f2t.a(idf, 0) = fv.idt0;
			f2t.a(idf, 1) = fv.idf0;
			f2t.a(idf, 2) = fv.idt1;
			f2t.a(idf, 3) = fv.idf1;
			++idf;

			ts[fv.idt0].face[fv.idf0].iadj = fv.idt1;
			if (fv.idt1 == -1) //fv.idt0 & fv.idf0 points to a boundary triangle, we need to extract it
			{
				TriangleFace& f = ts[fv.idt0].face[fv.idf0];
				triangles.push_back(triangle_octree(f.v[0].x, f.v[0].y, f.v[0].z, f.v[1].x, f.v[1].y, f.v[1].z, f.v[2].x, f.v[2].y, f.v[2].z, fv.idt0));
			}
			else ts[fv.idt1].face[fv.idf1].iadj = fv.idt0;
		}

		elem.cast<double>(temp);
		temp += 1; // c index to FORTRAN index
		memcpy(matvar->data, temp.getP(), temp.getInnerLength()*sizeof(double));
		Mat_VarDelete(matfp, "elem");
		Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
		Mat_VarFree(matvar);

		f2t.cast<double>(temp);
		temp += 1; // c index to FORTRAN index
		size_t dims[2] = { (size_t)NF, 4 };
		matvar = Mat_VarCreate("f2t", MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, temp.getP(), 0);
		if (NULL == matvar) fprintf(stderr, "Error creating variable for 'f2t'\n");
		else
		{
			Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
		}

#ifdef NDEBUG
		//fmap.clear(); //it's very slow monitored by debugger
		delete pmap;
#endif
	}

	Mat_VarFree(matvar);
	Mat_Close(matfp);
	
	vector<double> bound(6);
	bound[0] = xmin; bound[1] = xmax;
	bound[2] = ymin; bound[3] = ymax;
	bound[4] = zmin; bound[5] = zmax;
	root = octree_node::create_octree(triangles, 9, 6, bound);

	return NT;
}

SFloat Tetrahedral::getMass()
{
	Vector v0 = face[0].v[0] - face[0].v[1];
	Vector v1 = face[0].v[0] - face[0].v[2];
	Vector v2 = face[0].v[0] - face[1].v[0];
	return SFloat(density*fabs((v0&v1)*v2)/6);
}

Tetrahedral* TetMesh::tet = NULL;
octree_node* TetMesh::root = NULL;

bool TetMesh::load(ConfigFile* cf)
{
	release();
	string fname;
	if (!cf->getValue("mesh file", fname)) return false;
	vector<double> adjust;
	cf->getValue("mesh adjustment", adjust);
	RunTimeCounter timeCounter;
	nTet = Tetrahedral::load(tet, root, fname.c_str(), node, elem, f2t, adjust);
	dose.resize(nTet);
	dose = 0;
	uncertainty = dose;
	Log("It takes %.1f seconds to init tetrahedral mesh",timeCounter.stop());
	if (nTet <= 0) return false;
	cf->getValue("trim dose threshold", trimDensityThreshold);
	if (trimDensityThreshold < 0) trimDensityThreshold = 0;

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
		if (fabs(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] - 1.0) > 1e-8)
		{
			printf("Wrong unit vector!\n");
			return false;
		}
		Bx = dir[0];
		By = dir[1];
		Bz = dir[2];
	}
	else
	{
		printf("Wrong magnetic field direction parameters!\n");
		return false;
	}
	return true;
}

bool TetMesh::previousDose(const char* fname)
{
	BinaryFile BF;
	if (!BF.read(fname)) return false; //may be not the right format
	int nt;
	BF.get("nTet", &nt);

	if (nt != nTet) return false;

	if (!BF.get("dose", dose.getP())) return false;
	if (!BF.get("uncertainty", uncertainty.getP())) return false;
	if (!BF.get("Hist", &Hist)) return false;
	if (!BF.get("nBatch", &nBatch)) return false;
	if (!BF.get("rngSeed", &rngSeed)) return false;

	double invNorm = 1.0 / (1.602e-16*1.0889e15 * SourceHead_BeamOnTime());
	int NVoxel = nTet;
	for (int i = 0; i < NVoxel; ++i)
	{
		double D = Hist*dose[i];//previous total dose
		double U = (nBatch*uncertainty[i] * uncertainty[i] + 1) * D * D / nBatch; //previous sum dose^2
		double k = tet[i].getMass()* invNorm;
		dose.a(i) = SFloat(D*k); //DE
		uncertainty.a(i) = SFloat(U*k*k); //DE^2
	}
	return true;
}

void TetMesh::output(const char* fname, int FormatType) //output the dose counter
{
	int LEN = nTet;
	//may trim the dose just based on the density
	for (int i = 0; i < LEN; ++i)
	{
		if (tet[i].density < trimDensityThreshold) dose.a(i) = 0;
	}

	if (0 == FormatType)
	{
		BinaryFile BF;
		BF.push("nTet", &LEN, sizeof(int));
		int nNode = node.getWidth(1);
		BF.push("nNode", &nNode, sizeof(int));
		BF.push("node", node.getP(), sizeof(SFloat)*nNode * 3); //(x,y,z)
		BF.push("elem", elem.getP(), sizeof(int)*LEN*elem.getWidth(2));
		BF.push("dose", dose.getP(), sizeof(SFloat)*LEN);
		BF.push("uncertainty", uncertainty.getP(), sizeof(SFloat)*LEN);
		int nFace = f2t.getWidth(1);
		BF.push("nFace", &nFace, sizeof(int));
		BF.push("f2t", f2t.getP(), sizeof(int)*f2t.getInnerLength());

		BF.push("Hist", &Hist, sizeof(double));
		BF.push("nBatch", &nBatch, sizeof(int));
		BF.push("Bs", &Bs, sizeof(double));
		BF.push("Bx", &Bx, sizeof(double));
		BF.push("By", &By, sizeof(double));
		BF.push("Bz", &Bz, sizeof(double));
		
		BF.push("prescriptionDose", &prescriptionDose, sizeof(double));
		BF.push("treatmentFraction", &treatmentFraction, sizeof(int));
		BF.push("trimDensityThreshold", &trimDensityThreshold, sizeof(double));
		BF.push("rngSeed", &rngSeed, sizeof(int)); // used for continuous simulation

		if (!BF.write(fname)) exitApp("Cannot create the dose file. Please check the file name or writing access.");
	}
	Log("Accumulated history number in dose file = %.0f", Hist);
	Log("The dose file named <%s> has been generated ^_^\n\n", fname);
}

void TetMesh::getBinaryFile(BinaryFile& BF)
{
	int LEN = nTet;
	BF.push("nTet", &LEN, sizeof(int));
	int nNode = node.getWidth(1);
	BF.push("nNode", &nNode, sizeof(int));
	BF.push("node", node.getP(), sizeof(SFloat)*nNode * 3); //(x,y,z)
	BF.push("elem", elem.getP(), sizeof(int)*LEN * elem.getWidth(2));
	BF.push("dose", dose.getP(), sizeof(SFloat)*LEN);
	BF.push("uncertainty", uncertainty.getP(), sizeof(SFloat)*LEN);
	int nFace = f2t.getWidth(1);
	BF.push("nFace", &nFace, sizeof(int));
	BF.push("f2t", f2t.getP(), sizeof(int)*f2t.getInnerLength());

	BF.push("Hist", &Hist, sizeof(double));
	BF.push("nBatch", &nBatch, sizeof(int));
	BF.push("Bs", &Bs, sizeof(double));
	BF.push("Bx", &Bx, sizeof(double));
	BF.push("By", &By, sizeof(double));
	BF.push("Bz", &Bz, sizeof(double));

	BF.push("prescriptionDose", &prescriptionDose, sizeof(double));
	BF.push("treatmentFraction", &treatmentFraction, sizeof(int));
	BF.push("trimDensityThreshold", &trimDensityThreshold, sizeof(double));
	BF.push("rngSeed", &rngSeed, sizeof(int)); // used for continuous simulation
}

double TetMesh::addDose(SFloat* d, SFloat* u, int NB, double fH, double aNorm, double threshold) //return average uncertainty of 50% dose region
{
	SFloat norm = SFloat(1.602e-16*aNorm);
	int NVoxel = nTet;
	SFloat maxDose = 0;

	nBatch += NB;
	Hist += fH;
	for (int i = 0; i < NVoxel; ++i)
	{
		SFloat D = d[i];
		SFloat U = u[i];

		SFloat dNorm = norm / tet[i].getMass();
		D *= dNorm;//current total dose
		U *= dNorm*dNorm; //current sum dose^2

		dose.a(i) = SFloat(D / Hist); //current average dose per history
		if (D > 0) uncertainty.a(i) = SFloat(sqrt(U / (D*D) - 1.0 / nBatch)); //current standard deviation % based on batches
		else uncertainty.a(i) = -1;

		if (dose.a(i) > maxDose) maxDose = dose.a(i);
	}

	//find out the uncertainty of threshold dose region
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
	Log("D>%d%%*Dmax region occupies %.1f%% all voxels, with average uncertainty = %.1f%%", int(100 * threshold), 100.0*double(n50) / NVoxel, ave_err*100.0);
	return ave_err;
}

double TetMesh::peekUncertainty(SFloat* d, SFloat* u, int NB, double threshold)
{
	if (NB == 1) NB = 2;
	//assume NB is the total batch for one GPU work thread
	int NVoxel = nTet;
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
