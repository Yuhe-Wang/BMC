#ifndef  _TETRAHEDRAL_DEFINITION_H
#define  _TETRAHEDRAL_DEFINITION_H

#include "../Tools/Tools.h"
#include "octree_definition.h"
#include "common.h"



struct TriangleFace
{
	int iadj; //adjacent tetrahedral. if it's -1, means this face is the boundary
	Vector pv[3][2]; // plucker coordinates
	Vector v[3];
	Vector vn;
	void calcPluckerCoordinates()
	{
		pv[0][0] = v[2] - v[1];
		pv[0][1] = v[2] & v[1];
		pv[1][0] = v[0] - v[2];
		pv[1][1] = v[0] & v[2];
		pv[2][0] = v[1] - v[0];
		pv[2][1] = v[1] & v[0];

		vn = (pv[0][1] + pv[1][1] + pv[2][1])*(-1);
		vn.normalizeToUnitLength();
	}
	__device__ __host__ ZFloat dist(ZFloat x, ZFloat y, ZFloat z)
	{
		return (v[0].x - x)*vn.x + (v[0].y - y)*vn.y + (v[0].z - z)*vn.z;
	}
};

struct Tetrahedral
{
	Tetrahedral() :density(1){};
	//Vector v[4]; //used to calculate intersection point.
	TriangleFace face[4];
	short mid;
	SFloat density;
	SFloat getMass();
	static int load(Tetrahedral*& tet, octree_node*& root, const char* fname, ArrayMgr<SFloat>& node, ArrayMgr<int>& elem, ArrayMgr<int>& f2t, vector<double> shift = vector<double>()); //remember to delete[] tet and octree_node::delete_octree(root) after use
	
	void reverseClockwise()
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				for (int k = 0; k < 2; ++k)
				{
					face[i].pv[j][k] *= -1;
				}
			}
		}
	}

	void copyToStoreTet(StoreTet& tet)
	{
		tet.mid = mid;
		tet.density = density;
		for (short i = 0; i < 4; ++i)
		{
			for (short j = 0; j < 3; ++j) tet.face[i].v[j] = face[i].v[j];
			tet.face[i].iadj = face[i].iadj;
		}
	}
	void copyFromStoreTet(StoreTet& tet)
	{
		mid = tet.mid;
		density = tet.density;
		
		for (short i = 0; i < 4; ++i)
		{
			for (short j = 0; j < 3; ++j) face[i].v[j] = tet.face[i].v[j];
			face[i].iadj = tet.face[i].iadj;
		}
	}
};

struct FaceKey // the key value for hash table
{
	FaceKey(){};
	FaceKey(int i0, int i1, int i2){ iv0 = i0; iv1 = i1; iv2 = i2; sort(); }
	void set(int i0, int i1, int i2){ iv0 = i0; iv1 = i1; iv2 = i2; sort(); }
	int iv0, iv1, iv2;
	void sort()
	{
		if (iv0 > iv1) swapData<int>(iv0, iv1);
		if (iv1 > iv2) swapData<int>(iv1, iv2);
		if (iv0 > iv1) swapData<int>(iv0, iv1);
	}
	string operator () () const
	{
		char ch[48];
		sprintf(ch, "%d,%d,%d", iv0, iv1, iv2);
		return string(ch);
	}
	bool operator < (const FaceKey& rk) const
	{
		if (iv0 < rk.iv0) return true;
		else if (iv0 == rk.iv0)
		{
			if (iv1 < rk.iv1) return true;
			else if (iv1 == rk.iv1)
			{
				if (iv2 < rk.iv2) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	}
};

struct FaceValue
{
	FaceValue() :idt0(-1), idt1(-1){};
	void set(int idt, int idf)
	{
		if (idt0 == -1)
		{
			idt0 = idt;
			idf0 = idf;
			return;
		}
		else if (idt1 == -1)
		{
			idt1 = idt;
			idf1 = idf;
			return;
		}
		else
		{
			printf("error: this face belongs to more than two tetrahedrals.\n");
			getchar();
			exit(-1);
		}
	}
	int idt0;
	int idf0; // which face of this tetrahedral

	int idt1;
	int idf1; // which face of this tetrahedral
};

class TetMesh // this class handles mesh loading and dose generating
{
public:
	TetMesh() : trimDensityThreshold(0), Bs(0), Bx(0), By(0), Bz(1), rf(0),
		nBatch(0), Hist(0)
	{
	}
	bool load(ConfigFile* cf);
	static void release()
	{
		if (tet) delete[] tet;
		tet = NULL;
		if (root) octree_node::delete_octree(root);
		root = NULL;
	}
	bool previousDose(const char* fname);
	double addDose(SFloat* d, SFloat* u, int NB, double fH, double aNorm, double threhold = 0.5); //return the square of average uncertainty of threshold% dose region
	double peekUncertainty(SFloat* d, SFloat* u, int NB, double threshold);
	void getBinaryFile(BinaryFile& BF);
	void output(const char* fname, int FormatType = 0);
	int seedBegin(int shift) //it returns previous seed and update rngSeed
	{
		int ret = rngSeed;
		rngSeed += shift;
		return ret;
	}
	
	static Tetrahedral* tet;
	static octree_node* root;

	ArrayMgr<SFloat> dose;
	ArrayMgr<SFloat> uncertainty;
	ArrayMgr<SFloat> node;
	ArrayMgr<int> elem;
	ArrayMgr<int> f2t;//nFace* 4 matrix

	int nTet;
	int nBatch; // default 0
	double Hist; // default 0
	double trimDensityThreshold;
	double prescriptionDose;
	int treatmentFraction;
	double Bs, Bx, By, Bz;
	double rf;
	int rngSeed;
};
#endif // _TETRAHEDRAL_DEFINITION_H