#include "Phantom.h"
#include "gZeus.h"

//#define PRECISE

class RandomAzimuthHelper {

public:

	RandomAzimuthHelper(int nbin)
	{
		_table.resize(nbin);
		double dphi = 1. / (nbin - 1); _invBin = 1 / dphi;
		double cold = 1, sold = 0;
		for (int i = 1; i < nbin; ++i) {
			double phi = dphi*i;
			double c = cos(2 * PI*phi), s = sin(2 * PI*phi);
			_table[i - 1]._bc = (c - cold)*_invBin;
			_table[i - 1]._ac = c - _table[i - 1]._bc*phi;
			_table[i - 1]._bs = (s - sold)*_invBin;
			_table[i - 1]._as = s - _table[i - 1]._bs*phi;
			cold = c; sold = s;
		}
		_table[nbin - 1] = _table[nbin - 2];
	}

	inline void compute(double phi, double &cphi, double &sphi) const {
		int bin = (int)(phi*_invBin);
		cphi = _table[bin]._ac + _table[bin]._bc*phi;
		sphi = _table[bin]._as + _table[bin]._bs*phi;
	};

private:

	struct Data {
		double _ac, _bc, _as, _bs;
		Data() {};
	};

	vector<Data> _table;
	double       _invBin;

};

struct LambdaWoodcock
{
	double aIndex, bIndex;
	double a[NWOODCOCK];
	double b[NWOODCOCK];
	__device__ double lambdawck(double E)
	{
		int i = int(aIndex* (E - bIndex));
		return a[i] + b[i] * E;
	}
	bool init(const SFloat *dens, const char *matids, int nvoxels)
	{
		vector<double> av, bv;
		ZeusData_lambdaWoodcock(NWOODCOCK, dens, matids, nvoxels, aIndex, bIndex, av, bv);
		int n = av.size();
		if (n > NWOODCOCK) return false;
		for (int i = 0; i < n; ++i)
		{
			a[i] = av[i];
			b[i] = bv[i];
		}
		return true;
	}
};

struct LambdaPhoton
{
	double aIndex, bIndex;
	int n;
	double a[NLAMDAPHOTON];
	double b[NLAMDAPHOTON];
	double c[NLAMDAPHOTON];
	double d[NLAMDAPHOTON];
	__device__ double lambda(double E)
	{
		int i = int(aIndex* (E - bIndex));
		double E2 = E*E;
		return a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
	}
	bool init(int matid)
	{
		vector<double> av, bv, cv, dv;
		ZeusData_lambdaPhoton(matid, aIndex, bIndex, av, bv, cv, dv);
		n = av.size();
		if (n > NLAMDAPHOTON) return false;
		for (int i = 0; i < n; ++i)
		{
			a[i] = av[i];
			b[i] = bv[i];
			c[i] = cv[i];
			d[i] = dv[i];
		}
		double em = (n - 1) / aIndex + bIndex;
		return true;
	}
	void dumpData(FILE* fp)
	{
		fwrite(&bIndex, sizeof(double), 1, fp);
		fwrite(&aIndex, sizeof(double), 1, fp);
		fwrite(&n, sizeof(int), 1, fp);
		for (int i = 0; i < n; ++i)
		{
			double E = i / aIndex + bIndex;
			double E2 = E*E;
			double cs = a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
			fwrite(&cs, sizeof(double), 1, fp);
		}
	}
};

struct LambdaCompton
{
	double aIndex, bIndex;
	int n;
	double a[NLAMDACOMPTON];
	double b[NLAMDACOMPTON];
	double c[NLAMDACOMPTON];
	double d[NLAMDACOMPTON];
	__device__ double lambda(double E)
	{
		int i = int(aIndex* (E - bIndex));
		double E2 = E*E;
		return a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
	}
	bool init(int matid)
	{
		vector<double> av, bv, cv, dv;
		ZeusData_lambdaCo(matid, aIndex, bIndex, av, bv, cv, dv);
		n = av.size();
		if (n > NLAMDACOMPTON) return false;
		for (int i = 0; i < n; ++i)
		{
			a[i] = av[i];
			b[i] = bv[i];
			c[i] = cv[i];
			d[i] = dv[i];
		}
		return true;
	}
	void dumpData(FILE* fp)
	{
		fwrite(&bIndex, sizeof(double), 1, fp);
		fwrite(&aIndex, sizeof(double), 1, fp);
		fwrite(&n, sizeof(int), 1, fp);
		for (int i = 0; i < n; ++i)
		{
			double E = i / aIndex + bIndex;
			double E2 = E*E;
			double cs = a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
			fwrite(&cs, sizeof(double), 1, fp);
		}
	}
};

struct LambdaPair
{
	double aIndex, bIndex;
	int n;
	double a[NLAMDAPAIR];
	double b[NLAMDAPAIR];
	double c[NLAMDAPAIR];
	double d[NLAMDAPAIR];
	__device__ double lambda(double E)
	{
		int i = int(aIndex* (E - bIndex));
		double E2 = E*E;
		return a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
	}
	bool init(int matid)
	{
		vector<double> av, bv, cv, dv;
		ZeusData_lmbdapp(matid, aIndex, bIndex, av, bv, cv, dv);
		n = av.size();
		if (n > NLAMDAPAIR) return false;
		for (int i = 0; i < n; ++i)
		{
			a[i] = av[i];
			b[i] = bv[i];
			c[i] = cv[i];
			d[i] = dv[i];
		}
		return true;
	}
	void dumpData(FILE* fp)
	{
		fwrite(&bIndex, sizeof(double), 1, fp);
		fwrite(&aIndex, sizeof(double), 1, fp);
		fwrite(&n, sizeof(int), 1, fp);
		for (int i = 0; i < n; ++i)
		{
			double E = i / aIndex + bIndex;
			double E2 = E*E;
			double cs = a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
			fwrite(&cs, sizeof(double), 1, fp);
		}
	}
};

struct QSurface
{
	double fq, ieq0;
	double qs[NQSURFACE_E][NQSURFACE_Q];
	int ne, nuq;
	bool init(int matid)
	{
		vector<vector<double>> qsv;
		ZeusData_QSurface(matid, fq, ieq0, qsv);
		ne = qsv.size();
		if (ne > NQSURFACE_E) return false;
		nuq = qsv[0].size();
		if (nuq > NQSURFACE_Q) return false;
		for (int ie = 0; ie < ne; ++ie)
			for (int iq = 0; iq < nuq; ++iq)
				qs[ie][iq] = qsv[ie][iq];
		return true;
	}

	__device__ void getEnergyIndex1(double ie, int &energyIndex, double &probability)
	{
		double rle = fq * (ie - ieq0);
		if (rle > 0)
		{
			energyIndex = (int)rle;
			probability = rle - energyIndex;
			if (energyIndex >= ne - 1)
			{
				energyIndex = ne - 1;
				probability = -1;
			}
		}
		else
		{
			energyIndex = 0;
			probability = -1;
		}
	}

	__device__ double qsurf(int anEnergyIndex, double u)
	{
		int nuqp = nuq - 1;
		double ru = u * nuqp;
		int ju = (int)ru;
		if (ju > nuqp - 1) ju = nuqp - 1;
		ru -= ju;
		return qs[anEnergyIndex][ju] * (1 - ru) + qs[anEnergyIndex][ju + 1] * ru;
	};

	void dumpData(FILE* fp)
	{
		fwrite(&fq, sizeof(double), 1, fp);
		fwrite(&ieq0, sizeof(double), 1, fp);
		fwrite(&ne, sizeof(int), 1, fp);
		fwrite(&nuq, sizeof(int), 1, fp);
		for (int ie = 0; ie < ne; ++ie)
			for (int iq = 0; iq < nuq; ++iq)
			{
				fwrite(&(qs[ie][iq]), sizeof(double), 1, fp);
			}
	}
};

struct ScreeningParameter
{
	double aIndex, bIndex;
	int n;
	double a[NSCREENINGPARAMETER];
	double b[NSCREENINGPARAMETER];
	double c[NSCREENINGPARAMETER];
	double d[NSCREENINGPARAMETER];
	__device__ double bw1(double E)
	{
		int i = int(aIndex* (E - bIndex));
		double E2 = E*E;
		return a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
	}
	bool init(int matid)
	{
		vector<double> av, bv, cv, dv;
		ZeusData_ScreeningParameter(matid, aIndex, bIndex, av, bv, cv, dv);
		n = av.size();
		if (n > NSCREENINGPARAMETER) return false;
		for (int i = 0; i < n; ++i)
		{
			a[i] = av[i];
			b[i] = bv[i];
			c[i] = cv[i];
			d[i] = dv[i];
		}
		return true;
	}
	void dumpData(FILE* fp)
	{
		fwrite(&bIndex, sizeof(double), 1, fp);
		fwrite(&aIndex, sizeof(double), 1, fp);
		fwrite(&n, sizeof(int), 1, fp);
		for (int i = 0; i < n; ++i)
		{
			double E = i / aIndex + bIndex;
			double E2 = E*E;
			double cs = a[i] + b[i] * E + c[i] * E2 + d[i] * E*E2;
			fwrite(&cs, sizeof(double), 1, fp);
		}
	}
};

struct Range
{
	double aIndex, bIndex;
	int n;
	double a[NRANGE];
	double b[NRANGE];
	__device__ double lambda(double E)
	{
		int i = int(aIndex* (E - bIndex));
		return a[i] + b[i] * E;
	}
	bool init(int matid)
	{
		vector<double> av, bv;
		ZeusData_Range(matid, aIndex, bIndex, av, bv);
		n = av.size();
		if (n > NRANGE) return false;
		for (int i = 0; i < n; ++i)
		{
			a[i] = av[i];
			b[i] = bv[i];
		}
		return true;
	}
	void dumpData(FILE* fp)
	{
		fwrite(&bIndex, sizeof(double), 1, fp);
		fwrite(&aIndex, sizeof(double), 1, fp);
		fwrite(&n, sizeof(int), 1, fp);
		for (int i = 0; i < n; ++i)
		{
			double E = i / aIndex + bIndex;
			double cs = a[i] + b[i] * E;
			fwrite(&cs, sizeof(double), 1, fp);
		}
	}
};

struct InverseRange
{
	double aIndex, bIndex;
	int n;
	double a[NINVERSERANGE];
	double b[NINVERSERANGE];
	__device__ __host__ double lambda(double E)
	{
		int i = int(aIndex* (E - bIndex));
		return a[i] + b[i] * E;
	}
	bool init(int matid)
	{
		vector<double> av, bv;
		ZeusData_InverseRange(matid, aIndex, bIndex, av, bv);
		n = av.size();
		if (n > NINVERSERANGE) return false;
		for (int i = 0; i < n; ++i)
		{
			a[i] = av[i];
			b[i] = bv[i];
		}
		return true;
	}

	void dumpData(FILE* fp)
	{
		fwrite(&bIndex, sizeof(double), 1, fp);
		fwrite(&aIndex, sizeof(double), 1, fp);
		fwrite(&n, sizeof(int), 1, fp);
		for (int i = 0; i < n; ++i)
		{
			double E = i / aIndex + bIndex;
			double cs = a[i] + b[i] * E;
			fwrite(&cs, sizeof(double), 1, fp);
		}
	}
};

//this class is also accessed by the GPU so there's no construction function
struct Material
{
	int _matid;
	LambdaPhoton lmdPhoton;
	LambdaCompton lmdCompton;
	LambdaPair lmdPair;
	QSurface qsurface;
	ScreeningParameter scnpara;
	Range range;
	InverseRange invRange;
	bool init(int matid)
	{
		_matid = matid;
		if (!lmdPhoton.init(_matid)) return false;
		if (!lmdCompton.init(_matid)) return false;
		if (!lmdPair.init(_matid)) return false;
		if (!qsurface.init(_matid)) return false;
		if (!scnpara.init(_matid)) return false;
		if (!range.init(_matid)) return false;
		if (!invRange.init(_matid)) return false;
		return true;
	}
	bool dumpData(const char* filename)
	{
		FILE* fp = fopen(filename, "wb");
		if (NULL == fp) return false;
		lmdPhoton.dumpData(fp);
		lmdCompton.dumpData(fp);
		lmdPair.dumpData(fp);
		scnpara.dumpData(fp);
		range.dumpData(fp);
		invRange.dumpData(fp);
		qsurface.dumpData(fp);
		fclose(fp);
		return true;
	}
};

/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: variables in device memory <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
int NStackDepth;
int simuPositron;
int simuSecondary;
int NMaxSplit = 50;

ZFloat EAbsPhoton = 50e3;
ZFloat EAbsElectron = 50e3;
ZFloat ERangeCut = 10e3;
ZFloat EMaxCSDA = 200e3;
//}}

int NMAT = 1;
LambdaWoodcock lambdaWoodcock;
Material mat;
Phantom pht;
RandomAzimuthHelper randAzimuth(2048);
//}}
/*>>>>>>>>>>>>>>>>>>>>>>>>> end: variables in device memory >>>>>>>>>>>>>>>>>>>>>>>>>>>*/

#include "WaterCS.h"

// class PRNG //simple version, need test to verify it's uniform if force calling fillbuffer
// {
// public:
// 	__device__ void init(int seed) //init the generator by its index in an array or threads
// 	{
// 		long In = 1234 + seed;
// 
// 		In = 16807 * (In % 127773) - 2836 * (In / 127773);
// 		if (In < 0) In += 2147483647;
// 		ISEED1 = In;
// 		In = 16807 * (In % 127773) - 2836 * (In / 127773);
// 		if (In < 0) In += 2147483647;
// 		ISEED2 = In;
// 	}
// 	__device__ __forceinline__ float operator () () //default return random number in (0,1)
// 	{
// 		const double USCALE = 1.0 / 2.147483563e9;
// 		long I1 = ISEED1 / 53668;
// 		ISEED1 = 40014 * (ISEED1 - I1 * 53668) - I1 * 12211;
// 		if (ISEED1 < 0) ISEED1 += 2147483563;
// 
// 		long I2 = ISEED2 / 52774;
// 		ISEED2 = 40692 * (ISEED2 - I2 * 52774) - I2 * 3791;
// 		if (ISEED2 < 0) ISEED2 += 2147483399;
// 
// 		long IZ = ISEED1 - ISEED2;
// 		if (IZ < 1) IZ += 2147483562;
// 		return float(IZ*USCALE);
// 	}
// 	__device__ void printStatus()
// 	{
// 		printf("IS1 = %d, IS2 = %d\n", ISEED1, ISEED2);
// 	}
// private:
// 	int ISEED1;
// 	int ISEED2;
// };



/*<<<<<<<<<<<<<<<<<<<<<<<<<<<<< start: phantom method definitions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
/*


__device__ bool lineIn(ParticleR & p)
{
	double px = p.x;
	double py = p.y;
	double pz = p.z;
	//first translate the coordinate system
	px -= xo;
	py -= yo;
	pz -= zo;
	//assuming the incident particle go straight line, judge if it will enter the phantom.
	//if true, give the interaction position
	const double Delta = 1e-5;
	double sx, sy, sz;
	double x = -1.0, y = -1.0, z = -1.0;
	if (px < 0 || px >= LX || py < 0 || py >= LY || pz<0 || pz>LZ) //initial position lays outside the phantom
	{
		double pu = p.u;
		double pv = p.v;
		double pw = p.w;
		if (px < 0)
		{
			//now it's outside the phantom	
			if (pu > 0)
			{
				sx = 0;
				sy = py - px*pv / pu;
				sz = pz - px*pw / pu;
				if (0 <= sy &&sy < LY && 0 <= sz&&sz < LZ)
				{
					x = sx;
					y = sy;
					z = sz;
				}
			}
			else return false;
		}
		else if (px >= LX)
		{
			if (pu < 0)
			{
				sx = LX - Delta;
				sy = py + (sx - px)*pv / pu;
				sz = pz + (sx - px)*pw / pu;
				if (0 <= sy &&sy < LY && 0 <= sz&&sz < LZ)
				{
					x = sx;
					y = sy;
					z = sz;
				}
			}
			else return false;
		}

		if (py < 0)
		{
			//now it's outside the phantom
			if (pv > 0)
			{
				sy = 0;
				sx = px - py*pu / pv;
				sz = pz - py*pw / pv;
				if (0 <= sx &&sx < LX && 0 <= sz&&sz < LZ)
				{
					x = sx;
					y = sy;
					z = sz;
				}
			}
			else return false;
		}
		else if (py >= LY)
		{
			if (pv < 0)
			{
				sy = LY - Delta;
				sx = px + (sy - py)*pu / pv;
				sz = pz + (sy - py)*pw / pv;
				if (0 <= sx &&sx < LX && 0 <= sz&&sz < LZ)
				{
					x = sx;
					y = sy;
					z = sz;
				}
			}
			else return false;
		}

		if (pz < 0)
		{
			//now it's outside the phantom
			if (pw > 0)
			{
				sz = 0;
				sy = py - pz*pv / pw;
				sx = px - pz*pu / pw;
				if (0 <= sy &&sy < LY && 0 <= sx&&sx < LX)
				{
					x = sx;
					y = sy;
					z = sz;
				}
			}
			else return false;
		}
		else if (pz >= LZ)
		{
			if (pw < 0)
			{
				sz = LZ - Delta;
				sy = py + (sz - pz)*pv / pw;
				sx = px + (sz - pz)*pu / pw;
				if (0 <= sy &&sy < LY && 0 <= sx&&sx < LX)
				{
					x = sx;
					y = sy;
					z = sz;
				}
			}
			else return false;
		}

		if (x != -1.0) //got effective intersection point (x,y,z)
		{
			const double step = min(DX, min(DY, DZ));
			const int preNum = 1; //how many step it will forwardly detect
			const double dmx = step*pu;
			const double dmy = step*pv;
			const double dmz = step*pw;
			while (true)
			{

				int ix = int((x + preNum*dmx) / DX);
				int iy = int((y + preNum*dmy) / DY);
				int iz = int((z + preNum*dmz) / DZ);
				if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) return false;//it will leave the phantom without scattering
				if (ph[at(ix, iy, iz)] > 0.04) break; //stop when it get close to the target
				x += dmx;
				y += dmy;
				z += dmz;
			}
			px = x; //store the position
			py = y;
			pz = z;
		}
		else return false;
	}

	//prepare the voxel index
	int ix = int(px / DX);
	int iy = int(py / DY);
	int iz = int(pz / DZ);
	px -= ix*DX;
	py -= iy*DY;
	pz -= iz*DZ;

	//write back to global memory
	p.x = px;
	p.y = py;
	p.z = pz;
	p.ivx = ix;
	p.ivy = iy;
	p.ivz = iz;
	p.iabsv = at(ix, iy, iz);

	return true; //now it's ready to transport
}

__device__ bool movePhoton(ParticleR & p, double ds) //return whether particle leaves the phantom
{
	double px = p.x + ds*p.u + p.ivx*DX;
	double py = p.y + ds*p.v + p.ivy*DY;
	double pz = p.z + ds*p.w + p.ivz*DZ;
	if (px < 0 || px >= LX || py < 0 || py >= LY || pz < 0 || pz >= LZ) return true;
	//calculate the voxel index
	int ix = int(px / DX);
	int iy = int(py / DY);
	int iz = int(pz / DZ);

	p.x = px - DX*ix;
	p.y = py - DY*iy;
	p.z = pz - DZ*iz;
	p.ivx = ix;
	p.ivy = iy;
	p.ivz = iz;
	p.iabsv = at(ix, iy, iz);
	return false;
}

__device__ void ParticleStack::push(ParticleR &par)
{
	if (cur == NStackDepth)
	{
		printf("!!!!!!!!!!!!!!!exceed the max depth of the GPU stack!!!!!!!!!!!!!!!\n\
			   		!!!!!!!!!!!!!!!will abandon this particle now!!!!!!!!!!!!!!!!!!!!!!\n\
							!!!!!!!!!!!!!!!current maximum stack depth is %d !!!!!!!!!!!!!!!!\n\
									!!!!!!!!!!!!!!!you should increase NStackDepth later!!!!!!!!!!!!!!!\n", NStackDepth);
		//exitKernel("exceed the max depth of the GPU stack!");
	}
	else
	{

		ss[cur] = par; //copy and store the particle
		++cur;
	}
}


__device__ bool refill(size_t it)
{
	if (tv[it].pStack.empty())
	{
		int cur = InitCur[it];
		while (cur < NBatch)
		{
			ParticleR p = InitPars[it*NBatch + cur];
			if (!lineIn(p))
			{
				++cur;
				continue;
			}
			else
			{
				tv[it].p = p;
				cur++;
				InitCur[it] = cur;
				return true;
			}
		}
		return false;
	}
	else
	{
		tv[it].p = tv[it].pStack.top(); //it must be already in the phantom
		tv[it].pStack.pop();
		return true;
	}	
}

__device__ void samcom(size_t it, double energy, double &efrac, double &costh) 
{
	PRNG &rng = tv[it].rng;

	double ko = energy*INV_ELECTRON_MASS;
	double broi = 1 + 2 * ko; double bro = 1 / broi;
	double br, temp;
	if (ko < 10) {
		// "low" energy case: uniformly between bro and bro1.
		double bro1 = 1 - bro;
		double ko2 = ko*ko;
		double rejmax = ko2*(broi + bro);
		double br2;
		do {
			br = bro + bro1*rng(); br2 = br*br;
		} while (rng()*br2*rejmax > ko2*br*(br2 + 1) - (1 - br)*(br*broi - 1));
		temp = (1 - br) / (ko*br);
	}
	else {
		// "high" energy case: the usual way 
		double broi2 = broi*broi;
		double alpha1 = log(broi);
		double alpha2 = ko*(broi + 1)*bro*bro;
		double alphaS = alpha1 + alpha2;
		double sint;
		do {
			br = rng()*alphaS < alpha1 ? exp(alpha1*rng())*bro : sqrt(rng()*(broi2 - 1) + 1)*bro;
			temp = (1 - br) / (ko*br); sint = temp*(2 - temp);
		} while (rng()*(1 + br*br) < br*sint);
	}
	efrac = br;
	costh = 1 - temp;
}


__device__ 	int flightInConstantMagneticField(ParticleR& p, double Eend)
{
	//move for electron/photon. Note that coordinates are relative to each voxel
	double pE = p.E;
	//fetch data from the slow global memory
	double px = p.x;
	double py = p.y;
	double pz = p.z;
	double pu = p.u;
	double pv = p.v;
	double pw = p.w;
	int ivx = p.ivx;
	int ivy = p.ivy;
	int ivz = p.ivz;
	int iabsv = p.iabsv;
	//////////////////////////////////////////////
	const double deltax = 0.01;
	double voxden = ph[iabsv];
	double invDen = 1 / voxden;
	int matid = mids[iabsv] - 1;
	double range = (pE - ERangeCut)*mat[matid].range.lambda(pE);
	double finalRange = Eend > ERangeCut ? (Eend - ERangeCut)*mat[matid].range.lambda(Eend):0;
	double e = 0.5*(pE + Eend);
	if (e < EAbsElectron) e = EAbsElectron; // limit to Eabs

	double uwx = Bx;
	double uwy = By;
	double uwz = Bz;
	double Rb = rf*sqrt(pE*(e + TEs));
	double maxStep = sqrt(2 * Rb*deltax); //max allowed distance to move to ensure accuracy

	int loopInd = 0;
	while (true)
	{
		loopInd++;
		double step = (range - finalRange) *invDen;

		bool finalStep = true;
		if (step > maxStep)
		{
			step = maxStep;
			finalStep = false;
		}

		//check if it intersect with the boundary of current voxel
		short idex, dvox; //idex = 1,2,3 means x,y,z direction; dvox = +1, -1 means moving in positive or negative direction
		bool intersect = false;

		if (pv > 0)
		{
			double next = DY - py;
			if (pv*step > next)
			{
				step = next / pv; idex = 2; dvox = 1; intersect = true;
			}
		}
		else if (pv < 0)
		{
			double next = -py;
			if (pv*step < next)
			{
				step = next / pv; idex = 2; dvox = -1; intersect = true;
			}
		}

		if (pw > 0)
		{
			double next = DZ - pz;
			if (pw*step > next)
			{
				step = next / pw; idex = 3; dvox = 1; intersect = true;
			}
		}
		else if (pw < 0)
		{
			double next = -pz;
			if (pw*step < next)
			{
				step = next / pw; idex = 3; dvox = -1; intersect = true;
			}
		}

		if (pu > 0)
		{
			double next = DX - px;
			if (pu*step > next)
			{
				step = next / pu; idex = 1; dvox = 1; intersect = true;
			}
		}
		else if (pu < 0)
		{
			double next = -px;
			if (pu*step < next)
			{
				step = next / pu; idex = 1; dvox = -1; intersect = true;
			}
		}

		if (intersect) finalStep = false;

		double newEnergy = Eend;
		if (!finalStep) 
		{
			range -= step*voxden;
			newEnergy = mat[matid].invRange.lambda(range);
			if (newEnergy < ERangeCut) newEnergy = 0;
		}
// 		if (pE < newEnergy)
// 		{
// 			printf("Negative Energy\n");
// 		}
		deposit((pE - newEnergy)*p.weight,iabsv);
		//up the energy to the new one
		pE = newEnergy;
		if (pE < ERangeCut) return 1; // if this is the final step to local absorption, we don't need to update the position and the direction

		//move the electron/positron
		px += pu*step;
		py += pv*step;
		pz += pw*step;

		double vuw = pu*uwx + pv*uwy + pw*uwz;
		double vperpx = pu - vuw * uwx,
			vperpy = pu - vuw * uwy,
			vperpz = pu - vuw * uwz;
		double vxwx = vperpy*uwz - vperpz*uwy,
			vxwy = vperpz*uwx - vperpx*uwz,
			vxwz = vperpx*uwy - vperpy*uwx;
		// The step-length dependent variables f1 & f2
		double f1, f2;
		double arg = step / Rb;
		if (arg < 0.2)
		{
			// arg is small, so use power series expansion of sine and cosine
			double arg2 = arg*arg;
			f1 = -0.5*arg2 + 0.0416666667*arg2*arg2;  // for 0.2, relative error is 2.2e-6
			f2 = arg - 0.16666667*arg*arg2;           // for 0.2, relative error is 1.3e-5, absolute error is 2.6e-6
		}
		else
		{
			f1 = cos(arg) - 1;
			f2 = sin(arg);
		}
		// Direction change 
		double dvx = f1*vperpx - f2*vxwx;  // would simplify to f1*_v.x - f2*_v.y;
		double dvy = f1*vperpy - f2*vxwy;  // would simplify to f1*_v.y + f2*_v.x;
		double dvz = f1*vperpz - f2*vxwz;  // would simplify to 0 (i.e., component along the magnetic field remains constant).

		if (finalStep)
		{
			//update direction and break the loop
			pu += dvx;
			pv += dvy;
			pw += dvz;

			break;
		}

		//not the final step, we may need to check the direction
		if (intersect)
		{
			//
			// We are entering a new voxel. But because we are also changing the direction, we need to verify that the direction has not changed 
			// in a way that we are actually coming back into the voxel from which we came. The condition for coming back into the same voxel 
			// is that v*(v + dv) < 0, where v and dv are the direction and the direction change of the component crossing the voxel boundary
			//
			switch (idex)
			{
			case 3:
				if (pw*pw + pw*dvz >= 0) //enter a new voxel
				{
					ivz += dvox;
					if (dvox > 0) {
						if (ivz >= NZ) return 1;
						pz = 0;
					}
					else {
						if (ivz < 0) return 1;
						pz = DZ;
					}
				}
				else intersect = false;
				break;
			case 2:
				if (pv*pv + pv*dvy >= 0)
				{
					ivy += dvox;
					if (dvox > 0) {
						if (ivy >= NY) return 1;
						py = 0;
					}
					else {
						if (ivy < 0) return 1;
						py = DY;
					}
				}
				else intersect = false;
				break;
			default:
				if (pu*pu + pu*dvx >= 0)
				{
					ivx += dvox;
					if (dvox > 0) {
						if (ivx >= NX) return 1;
						px = 0;
					}
					else {
						if (ivx < 0) return 1;
						px = DX;
					}
				}
				else intersect = false;
			}
		}// end if (intersect)

		//still intersect after the direction check, need to update the voxel density and index
		if (intersect)
		{
			iabsv = at(ivx, ivy, ivz);
			voxden = ph[iabsv];
			invDen = 1 / voxden;

			int newMat = mids[iabsv] - 1;
			if (newMat != matid) 
			{
				// New material, we simply recompute the ranges in the new material and keep going
				matid = newMat;
				range = (pE - ERangeCut)*mat[matid].range.lambda(pE);
				finalRange = Eend > ERangeCut ? (Eend - ERangeCut)*mat[matid].range.lambda(Eend) : 0;
			}
		}

		//update direction and keep going
		pu += dvx;
		pv += dvy;
		pw += dvz;
	}

	//write back energy, direction and position
	p.E = pE;
	p.u = pu;
	p.v = pv;
	p.w = pw;
	p.x = px;
	p.y = py;
	p.z = pz;

	p.ivx = ivx;
	p.ivy = ivy;
	p.ivz = ivz;
	p.iabsv = iabsv;
	return 0;
}

__device__ 	int freeFlight(ParticleR& p, double Eend)
{
	return 1;
}

__device__ __forceinline__ double tperp(ParticleR& p)
{
	double px = p.x;
	double py = p.y;
	double pz = p.z;

	double tx = DX - px;
	if (px < tx) tx = px;
	double ty = DY - py; 
	if (py < ty) ty = py;
	double tz = DZ - pz;
	if (pz < tz) tz = pz;
	return tx < ty && tx < tz ? tx : ty < tz ? ty : tz; //min of (tx, ty, tz)
}

__device__ double samsca(size_t it, double e, int matid)
{
	PRNG &rng = tv[it].rng;
	// The screening parameter at this energy.
	double ie = -1 / e;
	double b = mat[matid].scnpara.bw1(ie);
	double oneb = 1.0 + b;

	// Determine energy bin of the pre-computed q surface.
	// Instead of linearly interpolating between energies in each qsurf() evaluation below, we 
	// use the lower or higher energy of the bin in which the energy falls with the correponding probability
	int je; double pe;
	mat[matid].qsurface.getEnergyIndex1(ie, je, pe);
	if (pe > 0) 
	{
		if (rng() < pe) ++je;
	}

	double u;
	while (1) 
	{
		u = rng();//u is in [0, 1)
		if (rng() < mat[matid].qsurface.qsurf(je, u)) break;
	}
	double mu = (oneb - u*(oneb + b)) / (oneb - u); // oneb>=1.0, so (oneb-u) is always positive
	return mu;
}

__device__ void electronRun(size_t it)
{
	PRNG &rng = tv[it].rng;
	ParticleR &p = tv[it].p;
	double E = p.E;
	if (E <= EAbsElectron)
	{
		if (E < ERangeCut) return;
		if (rf) flightInConstantMagneticField(p, 0);
		else freeFlight(p, 0);
		return;
	}
	double fuelel = E > EMaxCSDA ? EMaxCSDA : E;
	// Random hinge energy
	double fuelxt = fuelel * rng(); //the left fuel to exhaust for next flight
	fuelel -= fuelxt; //the fuel to be exhaust during the flight
	double ebefor = E;
	while (true)
	{
		double nextEnergy = E - fuelel;
		if (nextEnergy <= ERangeCut) nextEnergy = 0;
		int status = rf ? flightInConstantMagneticField(p, nextEnergy) : freeFlight(p, nextEnergy);
		if (status == 1) return; //below the cut-off energy or exit the phantom

		// Check if we can discontinue transport because the electron cannot escape the current voxel
		double tp = tperp(p);
		int iabsv = p.iabsv;
		int matid = mids[iabsv] - 1;
		E = p.E;
		double range = (E - ERangeCut)*mat[matid].range.lambda(E);
		if (range < tp*ph[iabsv])
		{
			deposit(E*p.weight, iabsv);
			return; //end the simulation of this electron
		}

		//do elastic multi-scattering
		double costhe = samsca(it, ebefor, matid);
		if (costhe < 1) p.changeByCos(costhe, 2 * PI*rng());


		// energy at the end of this elastic scattering step becomes initial energy of next step
		ebefor = E - fuelxt;
		if (ebefor > EAbsElectron)
		{
			// still above threshold => determine random hinge of next step
			double nextEloss = ebefor > EMaxCSDA ? EMaxCSDA : ebefor;
			double nextFuel = rng()*nextEloss;
			// now combine the remaining elastic fuel of this step (fuelxt) with the pre-hinge fuel of the next step
			// in this way there is no need to come back to this function and the electron can be tracked directly until the 
			// next elastic event.
			fuelel = fuelxt + nextFuel;
			// and also remember the post-hinge fuel of the next step
			fuelxt = nextEloss - nextFuel;
		}
		else 
		{
			if (E < ERangeCut) return;
			if (rf) flightInConstantMagneticField(p, 0);
			else freeFlight(p, 0);
			return;
		}
	}
	
	

}
__device__ void photonRun(size_t it)
{
	PRNG &rng = tv[it].rng;
	ParticleR &p = tv[it].p;
	double E = p.E;
	
	double e = E*1e-6;
	double pInitial = e > 0.3 ? -0.053139 + 1.0695 * e - 0.24783 * e*e + 0.028566 * e*e*e : 0.25;
	double nsplit = int(pInitial*NMaxSplit);
	double lammin = lambdaWoodcock.lambdawck(E);
	double rnno1 = rng();
	double rnno2 = rng();
	double delta = 1.0 / nsplit;

	int keepCompton = -1;
	int keepRejected = -1;
	int keepAnnih = -1;

	double ECompt, uCompt, vCompt, wCompt;
	double ESec, uSec, vSec, wSec;

	double Epair1, Epair2;

	double _lamph[1], _lamco[1], _lampair[1]; 
	for (int i = 0; i < NMAT; ++i)
	{
		_lamph[i] = -1; _lamco[i] = -1; _lampair[i] = -1;
	}

	
	ParticleR sp = p;//remember the initial particle status
	
	for (int isplit = 0; isplit < nsplit; ++isplit)
	{
		ParticleR pi = sp;
		double rnnoi = 1 - delta*(rnno1 + isplit);
		if (rnnoi <= 0) break;
		double s = -lammin*log(rnnoi);
		if (movePhoton(pi, s)) break; //may move out of the phantom
		
		int absvox = pi.iabsv;
		double voxden = ph[absvox];
		double lamden = lammin * voxden;
		int matid = mids[absvox] - 1;

		double lamph = _lamph[matid];
		if (lamph < 0)
		{
			lamph = mat[matid].lmdPhoton.lambda(E);
			_lamph[matid] = lamph;
		}

		if (rnno2 < lamden*lamph) //accept the knock 
		{
			//calculate cross-section of Compton scattering
			double lamco = _lamco[matid];
			if (lamco < 0) 
			{
				lamco = mat[matid].lmdCompton.lambda(E); 
				_lamco[matid] = lamco;
			}

			if (rnno2 < lamden * lamco)// It's a Compton interaction
			{
				if (keepCompton < 0) 
				{
					// Haven't sampled a Compton interaction yet, so do it now.
					// The index of the randomly selected surviving scattered photon
					keepCompton = (int)(rng()*nsplit);
					// Sample the interaction
					double efracCompt, costheCompt;
					samcom(it, E, efracCompt, costheCompt);
					// Azimuthal scattering angle
					double phiCompt = 2 * PI*rng();
					pi.newDirection(costheCompt, phiCompt, uCompt, vCompt, wCompt);
					ECompt = E*efracCompt;
					
					
					// Compute Compton electron direction
					double e0 = E * INV_ELECTRON_MASS;
					double efrac1 = 1 - efracCompt;
					double cost;
					if (efrac1 > 1.192092896e-07)
					{
						cost = (1.0 + e0) * sqrt(efrac1 / (e0*(2.0 + e0*efrac1)));
						if (cost > 1) cost = 1;
					}
					else cost = 0;
					pi.newDirection(cost, PI + phiCompt, uSec, vSec, wSec);
					ESec = E*(1 - efracCompt);
				}

				if (isplit == keepCompton && ECompt > EAbsPhoton) 
				{
					pi.E = ECompt;
					pi.u = uCompt;
					pi.v = vCompt;
					pi.w = wCompt;
					tv[it].pStack.push(pi);
				}
				
				pi.type = electron;
				pi.E = ESec;
				pi.weight *= delta;
				pi.u = uSec;
				pi.v = vSec;
				pi.w = wSec;
				p = pi; //make it to the current electron to simulate
				electronRun(it);
			}
			else // Not a Compton, so check if pair or photo.
			{
				bool doPair = E > TEs;
				if (doPair) 
				{
					double lampair = _lampair[matid];
					if (lampair < 0) 
					{
						lampair = mat[matid].lmdPair.lambda(E); 
						_lampair[matid] = lampair;
					}
					if (rnno2 > lamden*(lamco + lampair)) doPair = false;
				}
				if (doPair) // It's a pair production -> the photon disappears (but we add annihilation photons below as needed).
				{
					if (keepAnnih < 0) // Haven't sampled a pair event yet, so do it now. 
					{
						Epair1 = rng() * (E - TEs);
						Epair2 = E - TEs - Epair1;
						// The index of the randomly selected surviving annihilation photons
						keepAnnih = (int)(rng()*nsplit);
					}
					if (isplit == keepAnnih)
					{
						pi.type = photon;
						double cosTheta = -1 + 2 * rng();
						double u, v, w;
						pi.newDirection(cosTheta, 2 * PI*rng(), u, v, w);
						pi.u = u;
						pi.v = v;
						pi.w = w;
						pi.E = Es;
						tv[it].pStack.push(pi);

						pi.u = -u;
						pi.v = -v;
						pi.w = -w;
						tv[it].pStack.push(pi);
					}
					// Put an e+/e- pair on the stack. We do not distinguish between electrons and positrons at this point.
					pi.E = Epair1;
					pi.type = electron;
					pi.weight *= delta;
					tv[it].pStack.push(pi);

					pi.E = Epair2;
					tv[it].pStack.push(pi);
				}
				else // It's a photo absorption -> the photon disappears
				{
					// Put resulting electron on the stack.
					//_stack.setScndry(-1,energy,eweight,_r.x, _r.y, _r.z, dir.x, dir.y, dir.z,_xvox,_yvox,_zvox,_absvox);
					// Instead of putting the resulting electron on the stack, transport it here.
					pi.type = electron;
					pi.weight *= delta;
					p = pi; //make it as the current electron to simulate
					electronRun(it);
				}
			}
		}
		else // The interaction was rejected. 
		{
			if (keepRejected < 0) // The index of the randomly selected surviving rejected photon interaction
			{
				keepRejected = (int)(rng()*nsplit);
			}
			if (isplit == keepRejected) tv[it].pStack.push(pi); //we need to put the original photon at its new location on the stack
		}
	}
}
__global__ void gZeusRun() //let's implement in a simple way first
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	InitCur[it] = 0;//reset the current particle index
	while (refill(it))
	{
		//now the particle is in phantom now
		if (tv[it].p.type == photon) photonRun(it);
		else electronRun(it);
	}
}

__global__ void fillParticle() //handle threads which are labeled with EXIT_EVT
{
	size_t it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	bool hasSource = true;
	++tv[it].loop;

	if (tv[it].evt == EXIT_EVT)
	{
		if (!simuSecondary || tv[it].pStack.empty())
		{
			//try to load a photon from the source pool
			if (InitCur[it] < NBatch)
			{
				//fetch the primary particles
				tv[it].p = InitPars[it*NBatch + InitCur[it]];
				++InitCur[it];
				tv[it].evt = INIT_EVT;
				//move the primary particle to the phantom
				if (!lineIn(tv[it].p)) tv[it].evt = EXIT_EVT;
			}
			else hasSource = false;//cannot get any incident photon for this thread
		}
		else //need to simulate the secondary particles stored in the stack
		{
			tv[it].p = tv[it].pStack.top();
			tv[it].pStack.pop();
			tv[it].evt = INIT_EVT;
		}
	}

	if (hasSource) ++ContinueSignal; //if ContinueSignal!=0, then the loop continues
}
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end: Kernel definitions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
int Phantom::electronFlight(ParticleR& p, ZFloat Eend)
{
	//move for electron/photon. Note that coordinates are relative to each voxel
	const ZFloat deltax = GlueF(0.01);
	ZFloat voxden = ph.a(p.iabsv);
	ZFloat invDen = 1 / voxden;
#ifdef PRECISE
	ZFloat range = (p.E - ERangeCut)*mat.range.lambda(p.E);
	ZFloat finalRange = Eend > ERangeCut ? (Eend - ERangeCut)*mat.range.lambda(Eend) : 0;
#else
	ZFloat range = (p.E - ERangeCut)*WaterRangeCS(p.E);
	ZFloat finalRange = Eend > ERangeCut ? (Eend - ERangeCut)*WaterRangeCS(Eend) : 0;
#endif
	ZFloat e = GlueF(0.5)*(p.E + Eend);
	if (e < EAbsElectron) e = EAbsElectron; // limit to Eabs

	ZFloat uwx = Bx;
	ZFloat uwy = By;
	ZFloat uwz = Bz;
	ZFloat Rb = rf*GlueF(sqrt)(e*(e + TEs));
	ZFloat Rbi = 1 / Rb; //in case this value being used many times
	ZFloat maxStep = GlueF(sqrt)(2 * Rb*deltax); //max allowed distance to move to ensure accuracy

	int loopInd = 0;
	while (true)
	{
		loopInd++;
		ZFloat step = (range - finalRange) *invDen;

		bool finalStep = true;
		if (step > maxStep)
		{
			step = maxStep;
			finalStep = false;
		}

		//check if it intersect with the boundary of current voxel
		short idex, dvox; //idex = 1,2,3 means x,y,z direction; dvox = +1, -1 means moving in positive or negative direction
		bool intersect = false;

		if (p.v > 0)
		{
			ZFloat next = DY - p.y;
			if (p.v*step > next)
			{
				step = next / p.v; idex = 2; dvox = 1; intersect = true;
			}
		}
		else if (p.v < 0)
		{
			ZFloat next = -p.y;
			if (p.v*step < next)
			{
				step = next / p.v; idex = 2; dvox = -1; intersect = true;
			}
		}

		if (p.w > 0)
		{
			ZFloat next = DZ - p.z;
			if (p.w*step > next)
			{
				step = next / p.w; idex = 3; dvox = 1; intersect = true;
			}
		}
		else if (p.w < 0)
		{
			ZFloat next = -p.z;
			if (p.w*step < next)
			{
				step = next / p.w; idex = 3; dvox = -1; intersect = true;
			}
		}

		if (p.u > 0)
		{
			ZFloat next = DX - p.x;
			if (p.u*step > next)
			{
				step = next / p.u; idex = 1; dvox = 1; intersect = true;
			}
		}
		else if (p.u < 0)
		{
			ZFloat next = -p.x;
			if (p.u*step < next)
			{
				step = next / p.u; idex = 1; dvox = -1; intersect = true;
			}
		}

		if (intersect) finalStep = false;

		ZFloat newEnergy = Eend;
		if (!finalStep)
		{
			range -= step*voxden;
#ifdef PRECISE
			newEnergy = mat.invRange.lambda(range);
#else
			newEnergy =WaterInverseRangeCS(range);
#endif
			if (newEnergy < ERangeCut) newEnergy = 0;
		}
		deposit(p.iabsv, (p.E - newEnergy)*p.weight);
		//up the energy to the new one
		p.E = newEnergy;
		if (p.E < ERangeCut) return 1; // if this is the final step to local absorption, we don't need to update the position and the direction

		//move the electron/positron
		p.x += p.u*step;
		p.y += p.v*step;
		p.z += p.w*step;

		ZFloat vuw = p.u*uwx + p.v*uwy + p.w*uwz;
		ZFloat vperpx = p.u - vuw * uwx,
			vperpy = p.v - vuw * uwy,
			vperpz = p.w - vuw * uwz;
		ZFloat vxwx = vperpy*uwz - vperpz*uwy,
			vxwy = vperpz*uwx - vperpx*uwz,
			vxwz = vperpx*uwy - vperpy*uwx;
		// The step-length dependent variables f1 & f2
		ZFloat f1, f2;
		ZFloat arg = step * Rbi;
		if (arg < GlueF(0.2))
		{
			// arg is small, so use power series expansion of sine and cosine
			ZFloat arg2 = arg*arg;
			f1 = -GlueF(0.5)*arg2 + GlueF(0.0416666667)*arg2*arg2;  // for 0.2, relative error is 2.2e-6
			f2 = arg - GlueF(0.16666667)*arg*arg2;           // for 0.2, relative error is 1.3e-5, absolute error is 2.6e-6
		}
		else
		{
			f1 = GlueF(cos)(arg) - 1;
			f2 = GlueF(sin)(arg);
		}
		// Direction change 
		ZFloat dvx = f1*vperpx - f2*vxwx;  // would simplify to f1*_v.x - f2*_v.y;
		ZFloat dvy = f1*vperpy - f2*vxwy;  // would simplify to f1*_v.y + f2*_v.x;
		ZFloat dvz = f1*vperpz - f2*vxwz;  // would simplify to 0 (i.e., component along the magnetic field remains constant).

		if (finalStep)
		{
			//update direction and break the loop
			p.u += dvx;
			p.v += dvy;
			p.w += dvz;

			break;
		}

		//not the final step, we may need to check the direction
		if (intersect)
		{
			//
			// We are entering a new voxel. But because we are also changing the direction, we need to verify that the direction has not changed 
			// in a way that we are actually coming back into the voxel from which we came. The condition for coming back into the same voxel 
			// is that v*(v + dv) < 0, where v and dv are the direction and the direction change of the component crossing the voxel boundary
			//
			switch (idex)
			{
			case 3:
				if (p.w*p.w + p.w*dvz >= 0) //enter a new voxel
				{
					p.ivz += dvox;
					if (dvox > 0) {
						if (p.ivz >= NZ) return 1;
						p.z = 0;
					}
					else {
						if (p.ivz < 0) return 1;
						p.z = DZ;
					}
				}
				else intersect = false;
				break;
			case 2:
				if (p.v*p.v + p.v*dvy >= 0)
				{
					p.ivy += dvox;
					if (dvox > 0) {
						if (p.ivy >= NY) return 1;
						p.y = 0;
					}
					else {
						if (p.ivy < 0) return 1;
						p.y = DY;
					}
				}
				else intersect = false;
				break;
			default:
				if (p.u*p.u + p.u*dvx >= 0)
				{
					p.ivx += dvox;
					if (dvox > 0) {
						if (p.ivx >= NX) return 1;
						p.x = 0;
					}
					else {
						if (p.ivx < 0) return 1;
						p.x = DX;
					}
				}
				else intersect = false;
			}
		}// end if (intersect)

		//still intersect after the direction check, need to update the voxel density and index
		if (intersect)
		{
			p.iabsv = at(p.ivx, p.ivy, p.ivz);
			voxden = ph[p.iabsv];
			invDen = 1 / voxden;
		}

		//update direction and keep going
		p.u += dvx;
		p.v += dvy;
		p.w += dvz;
	}
	return 0;
}

class LambdaCalculator {

public:

	LambdaCalculator(int nbin){
		_table.resize(nbin);
		_xmin = GlueF(0.01); ZFloat xmax = GlueF(1.00000001);
		ZFloat dx = (xmax - _xmin) / (nbin - 1);
		ZFloat fold = GlueF(log)(_xmin);
		_bbin = 1 / dx; _abin = -_xmin*_bbin;
		for (int i = 1; i < nbin; ++i) {
			ZFloat x = _xmin + dx*i;
			ZFloat f = GlueF(log)(x);
			_table[i - 1].first = (f - fold)*_bbin;
			_table[i - 1].second = f - _table[i - 1].first*x;
			fold = f;
		}
		_table[nbin - 1] = _table[nbin - 2];
	};

	inline ZFloat compute(ZFloat rnno) const {
		ZFloat result;
		if (rnno > _xmin) {
			int bin = (int)(rnno*_bbin + _abin);
			result = _table[bin].first*rnno + _table[bin].second;
		}
		else result = GlueF(log)(rnno);
		return result;
	};

private:

	vector<pair<ZFloat, ZFloat> >  _table;
	ZFloat                        _xmin;
	ZFloat                        _abin, _bbin;

};

LambdaCalculator _lambdaCalculator(1024);


void Rotate(ZFloat&x, ZFloat& y, ZFloat& z, ZFloat costh, ZFloat cosph, ZFloat sinph) {
	ZFloat costh2 = costh*costh;
	ZFloat rho2 = x * x + y * y;
	if (rho2 > 0 && costh2 < 1) {
		ZFloat a = GlueF(sqrt)((1 - costh2) / rho2);
		ZFloat xrho = x * a;
		ZFloat yrho = y * a;
		ZFloat ztmp = z * cosph;
		x = x * costh - yrho * sinph + ztmp * xrho;
		y = y * costh + xrho * sinph + ztmp * yrho;
		z = z * costh - rho2 * a * cosph;
	}
	else {
		if (costh2 >= 1) {
			if (costh < 0) {
				x = -x; y = -y; z = -z;
			}
			return;
		}
		ZFloat b = GlueF(sqrt)(1 - costh2);
		y = b * sinph;
		if (z > 0) {
			x = b * cosph;
			z = costh;
		}
		else {
			x = -b * cosph;
			z = -costh;
		}
	}
};
void samcom(ZFloat energy, ZFloat &efrac, ZFloat &costh, PRNG& rng)
{
	ZFloat ko = energy*INV_ELECTRON_MASS;
	ZFloat broi = 1 + 2 * ko; ZFloat bro = 1 / broi;
	ZFloat br, temp;
	if (ko < 10) {
		// "low" energy case: uniformly between bro and bro1.
		ZFloat bro1 = 1 - bro;
		ZFloat ko2 = ko*ko;
		ZFloat rejmax = ko2*(broi + bro);
		ZFloat br2;
		do {
			br = bro + bro1*rng(); br2 = br*br;
		} while (rng()*br2*rejmax > ko2*br*(br2 + 1) - (1 - br)*(br*broi - 1));
		temp = (1 - br) / (ko*br);
	}
	else {
		// "high" energy case: the usual way 
		ZFloat broi2 = broi*broi;
		ZFloat alpha1 = GlueF(log)(broi);
		ZFloat alpha2 = ko*(broi + 1)*bro*bro;
		ZFloat alphaS = alpha1 + alpha2;
		ZFloat sint;
		do {
			br = rng()*alphaS < alpha1 ? GlueF(exp)(alpha1*rng())*bro : GlueF(sqrt)(rng()*(broi2 - 1) + 1)*bro;
			temp = (1 - br) / (ko*br); sint = temp*(2 - temp);
		} while (rng()*(1 + br*br) < br*sint);
	}
	efrac = br;
	costh = 1 - temp;
}

ZFloat tperp(ParticleR& p)
{
	ZFloat tx = pht.DX - p.x;
	if (p.x < tx) tx = p.x;
	ZFloat ty = pht.DY - p.y;
	if (p.y < ty) ty = p.y;
	ZFloat tz = pht.DZ - p.z;
	if (p.z < tz) tz = p.z;
	return tx < ty && tx < tz ? tx : ty < tz ? ty : tz; //min of (tx, ty, tz)
}

ZFloat samsca(ZFloat e, PRNG& rng)
{

	// The screening parameter at this energy.
	ZFloat ie = -1 / e;
#ifdef PRECISE
	ZFloat b = mat.scnpara.bw1(ie);
#else
	ZFloat b = WaterScreenCS(ie);
#endif
	ZFloat oneb = GlueF(1.0) + b;

	// Determine energy bin of the pre-computed q surface.
	// Instead of linearly interpolating between energies in each qsurf() evaluation below, we 
	// use the lower or higher energy of the bin in which the energy falls with the correponding probability
	int je; ZFloat pe;
#ifdef PRECISE
	mat.qsurface.getEnergyIndex1(ie, je, pe);
#else
	WaterQSGetEnergyIndex1(ie, je, pe);
#endif
	if (pe > 0)
	{
		if (rng() < pe) ++je;
	}

	ZFloat u;
	while (1)
	{
		u = rng();//u is in [0, 1)
#ifdef PRECISE
		if (rng() < mat.qsurface.qsurf(je, u)) break;
#else
		if (rng() < WaterQSurface(je, u)) break;
#endif
	}
	ZFloat mu = (oneb - u*(oneb + b)) / (oneb - u); // oneb>=1.0, so (oneb-u) is always positive
	return mu;
}

double sE = 0;
double aE = 0;

void electronRun(ParticleR& p, PRNG& rng)
{
	if (p.E < 743e3) sE += 1;
	aE += 1;

	if (p.E <= EAbsElectron)
	{
		if (p.E < ERangeCut) return;
		pht.electronFlight(p, 0);
		return;
	}
	ZFloat fuelel = p.E > EMaxCSDA ? EMaxCSDA : p.E;
	// Random hinge energy
	ZFloat fuelxt = fuelel * rng(); //the left fuel to exhaust for next flight
	fuelel -= fuelxt; //the fuel to be exhaust during the flight
	ZFloat ebefor = p.E;
	while (true)
	{
		ZFloat nextEnergy = p.E - fuelel;
		if (nextEnergy <= ERangeCut) nextEnergy = 0;
		int status = pht.electronFlight(p, nextEnergy);
		if (status == 1) return; //below the cut-off energy or exit the phantom

		// Check if we can discontinue transport because the electron cannot escape the current voxel
		ZFloat tp = tperp(p);

#ifdef PRECISE
		ZFloat range = (p.E - ERangeCut)*mat.range.lambda(p.E);
#else
		ZFloat range = (p.E - ERangeCut)*WaterRangeCS(p.E);
#endif
		if (range < tp*pht.getDensity(p))
		{
			pht.deposit(p.iabsv, p.E*p.weight);
			return; //end the simulation of this electron
		}

		//do elastic multi-scattering
		ZFloat costhe = samsca(ebefor, rng);

		if (costhe < 1) //p.changeByCos(costhe, 2 * PI*rng());
		{
			ZFloat cphi, sphi; 
			randAzimuth.compute(rng(), cphi, sphi);
			Rotate(p.u, p.v, p.w, costhe, cphi, sphi);
		}


		// energy at the end of this elastic scattering step becomes initial energy of next step
		ebefor = p.E - fuelxt;
		if (ebefor > EAbsElectron)
		{
			// still above threshold => determine random hinge of next step
			ZFloat nextEloss = ebefor > EMaxCSDA ? EMaxCSDA : ebefor;
			ZFloat nextFuel = rng()*nextEloss;
			// now combine the remaining elastic fuel of this step (fuelxt) with the pre-hinge fuel of the next step
			// in this way there is no need to come back to this function and the electron can be tracked directly until the 
			// next elastic event.
			fuelel = fuelxt + nextFuel;
			// and also remember the post-hinge fuel of the next step
			fuelxt = nextEloss - nextFuel;
		}
		else
		{
			if (p.E < ERangeCut) return;
			pht.electronFlight(p, 0);
			return;
		}
	}
}

void smartPhoton(ParticleR& p, ParticleStack& pStack, PRNG& rng)
{
	while (true)
	{
		//splitting number
		ZFloat e = p.E*GlueF(1e-6);
		ZFloat pInitial = e > GlueF(0.3) ? GlueF(-0.053139) + GlueF(1.0695) * e - GlueF(0.24783) * e*e + GlueF(0.028566) * e*e*e : GlueF(0.25);
		ZFloat asplit = pInitial*NMaxSplit;
		int nsplit = (int)asplit;
		asplit -= nsplit;
		if (rng() < asplit) ++nsplit;
		if (nsplit < 1) return;

		//ZFloat lammin = lambdaWoodcock.lambdawck(p.E);
		ZFloat rnno1 = rng();
		ZFloat rnno2 = rng();
		ZFloat delta = 1 / ZFloat(nsplit);
		ZFloat eweight = p.weight*delta;

#ifdef PRECISE
		ZFloat lammin = lambdaWoodcock.lambdawck(p.E);
		ZFloat lamph = -1;
#else
		ZFloat lamph = WaterPhotonCS(p.E);
		ZFloat lammin = 1 / (lamph*pht.MaxDensity);
#endif
		ZFloat lamco = -1;
		ZFloat lampair = -1;

		// 	int keepCompton = -1;
		// 	int keepRejected = -1;
		// 	int keepAnnih = -1;
		int keepID = (int)(rng()*nsplit);
		// These will be used to remember the relevant Compton scattering variables
		ZFloat eCompElec = -1;
		ZFloat eu = p.u, ev = p.v, ew = p.w;
		ZFloat costheCompt, cphiCompt, sphiCompt; //used to calculate the direction of the scattered photon
		
		ZFloat Epair1 = -1; // , Epair2 = -1;

		ParticleR pOld = p; //remember the initial status
		p.E = 0; //default to exit
		for (int isplit = 0; isplit < nsplit; ++isplit)
		{
			ParticleR pi = pOld; //start with the initial status
			ZFloat rnnoi = 1 - delta*(rnno1 + isplit);
			if (rnnoi <= 0) break;
			ZFloat s = -lammin*_lambdaCalculator.compute(rnnoi);

			if (pht.photonFlight(pi, s)) break;

			// update lammin with density in the voxel and get material id.
			ZFloat voxden = pht.getDensity(pi);
			ZFloat lamden = lammin * voxden;
#ifdef PRECISE
			if(lamph<0) lamph = mat.lmdPhoton.lambda(pi.E);
#else
			if (lamph<0) lamph = WaterPhotonCS(pi.E);
#endif

			// Check if a real interaction
			if (rnno2 < lamden*lamph) // yes. 
			{
#ifdef PRECISE
				if (lamco < 0) lamco = mat.lmdCompton.lambda(pi.E);
#else
				if (lamco < 0) lamco = WaterComptonCS(pi.E);
#endif

				if (rnno2 < lamden * lamco) // It's a Compton interaction
				{
					if (eCompElec < 0) // Haven't sampled a Compton interaction yet, so do it now.
					{
						//keepCompton = (int)(rng()*nsplit);
						// Sample the interaction
						ZFloat efracCompt = 1;
						samcom(pi.E, efracCompt, costheCompt, rng);
// 						ZFloat phi = 2 * PI*rng();
// 						cphiCompt = GlueF(cos)(phi);
// 						if (phi < PI) sphiCompt = GlueF(sqrt)(1 - cphiCompt*cphiCompt);
// 						else sphiCompt = -GlueF(sqrt)(1 - cphiCompt*cphiCompt);
						randAzimuth.compute(rng(), cphiCompt, sphiCompt);

						// Electron energy
						efracCompt = 1 - efracCompt;
						eCompElec = pi.E*efracCompt;
						// Compute Compton electron direction
						ZFloat e0 = pi.E * INV_ELECTRON_MASS;
						//ZFloat efrac1 = 1 - efracCompt;
						ZFloat cost;
						if (efracCompt > FLTEPSILON)
						{
							cost = (1 + e0) * GlueF(sqrt)(efracCompt / (e0*(2.0 + e0*efracCompt)));
							if (cost > 1) cost = 1;
						}
						else cost = 0;
						Rotate(eu, ev, ew, cost, -cphiCompt, -sphiCompt);
					}

					//if (isplit == keepCompton && eCompGamma > EAbsPhoton)
					if (isplit == keepID)
					{
						ZFloat eCompGamma = pi.E - eCompElec;
						if (eCompGamma > EAbsPhoton)
						{
							Rotate(pi.u, pi.v, pi.w, costheCompt, cphiCompt, sphiCompt);
							pi.E = eCompGamma;
							p = pi;
						}
						else p.E = 0; //indicate that the energy is too low, so we can abandon simulating this photon
					}

					// Now, instead of first pushing the Compton electron onto the stack and later getting it back, we simply transport it here.
					pi.type = electron;
					pi.E = eCompElec;
					pi.weight = eweight;
					pi.u = eu;
					pi.v = ev;
					pi.w = ew;
					electronRun(pi, rng);
				}
				else // Not a Compton, so check if pair or photo.
				{
					bool doPair = pi.E > TEs;
					if (doPair)
					{
#ifdef PRECISE
						if (lampair < 0) lampair = mat.lmdPair.lambda(pi.E);
#else
						if (lampair < 0) lampair = WaterPairCS(pi.E);
#endif

						if (rnno2 > lamden*(lamco + lampair)) doPair = false;
					}
					if (doPair) // It's a pair production -> the photon disappears (but we add annihilation photons below as needed).
					{
						if (Epair1 < 0) // Haven't sampled a pair event yet, so do it now. 
						{
							Epair1 = rng() * (pi.E - 2.0*510998.910); //Es = 510998.910 //debug purpose
						}

						if (isplit == keepID)
						{
							ZFloat vz = 2 * rng() - 1;
							ZFloat sinthe = GlueF(sqrt)(1.0 - vz*vz);
							ZFloat phi = 2 * PI * rng();
							p.u = sinthe * GlueF(cos)(phi);
							p.v = sinthe * GlueF(sin)(phi);
							p.w = vz;
							p.E = 510998.910; //for debug, it should be Es
							//copy the position
							p.x = pi.x; p.y = pi.y; p.z = pi.z;
							p.ivx = pi.ivx; p.ivy = pi.ivy; p.ivz = pi.ivz; p.iabsv = pi.iabsv;
							pStack.push(p);

							//The other photon has opposite direction
							p.u = -p.u;
							p.v = -p.v;
							p.w = -p.w;
						}

						// Put an e+/e- pair on the stack. We do not distinguish between electrons and positrons at this point.
						ZFloat bx = pi.x; ZFloat by = pi.y; ZFloat bz = pi.z; //backup position
						int bix = pi.ivx; int biy = pi.ivy; int biz = pi.ivz; int biabs = pi.iabsv;

						pi.type = electron;
						pi.E = Epair1;
						pi.weight = eweight;
						electronRun(pi, rng);
						//restore the position, direction
						pi.x = bx; pi.y = by; pi.z = bz;
						pi.ivx = bix; pi.ivy = biy; pi.ivz = biz; pi.iabsv = biabs;
						pi.u = pOld.u; pi.v = pOld.v; pi.w = pOld.w;
						pi.E = pOld.E - 2.0*510998.910 - Epair1;
						electronRun(pi, rng);
					}
					else
					{
						// It's a photo absorption -> the photon disappears
						pi.type = electron;
						pi.weight = eweight;
						electronRun(pi, rng);
// 						if (isplit == keepID)
// 						{
// 							p.E = 0; // no energy now. It's useless since we have set p.E = 0 before this loop
// 						}

					}
				}
			}
			else // The interaction was rejected. 
			{
				if (isplit == keepID) //copy the position, and energy
				{
					p.x = pi.x; p.y = pi.y; p.z = pi.z;
					p.ivx = pi.ivx; p.ivy = pi.ivy; p.ivz = pi.ivz; p.iabsv = pi.iabsv;
					p.E = pi.E;
				}
			}
		}
		if (p.E < EAbsPhoton) return;
	}
	
}

void executeJob(const char* configFileName, MPS& configMacro) //execute one job according to the config file
{
	RunTimeCounter totalTime;
	ConfigFile cf(configFileName); //parse the total config file
	//find out where is the config file located
	string cdir(configFileName);
	size_t pos = cdir.find_last_of("\\/");
	cdir.erase(++pos, string::npos);
	cf.macroReplace(string("$cdir$"), cdir);
	cf.macroReplace(configMacro);

	string logDir, logAppend, logDescription;
	cf.getValue("log file directory", logDir);
	cf.getValue("log file append", logAppend);
	cf.getValue("log description", logDescription);
	if (0 == logAppend.compare("no")) logAppend = "w";
	else logAppend = "a";
	Log.setLogName(0, logAppend, logDir);//start log recording for this job
	Log("The job config file name = %s", configFileName);
	if (logDescription.compare("NA") != 0) Log("Short description: %s", logDescription.c_str());
	Log("Start log time = %s\n\n", Log.timeNow());

	double fNSIMU = 1e7;
	cf.getValue("NSIMU", fNSIMU);

	//initialize GZEUS and SourceHead by configurations
	ConfigFile *scf = cf.getBlock("SOURCEHEAD");
	SourceHead_Init(scf);

	//config the source particle pool
	int NGT = 1;
	scf->getValue("NThread", NGT);
	
	string dataDir;
	scf->getValue("DataDir",dataDir);
	if (!ZeusData_load(dataDir.c_str())) exitApp("Cannot load Zeus cross-sections correctly!");
	
	ConfigFile* ph_cf = cf.getBlock("PHANTOM");
	if (ph_cf == NULL) exitApp("Cannot find the phantom configuration");
	pht.loadPhantom(ph_cf);
	SourceHead_GetPrescrition(&(pht.prescriptionDose), &(pht.treatmentFraction));
	int NVoxel = pht.getNVoxel();

	Log("\nCalculating dose, please wait patiently...\n\n");
	RunTimeCounter rc; //count the calculating time
	vector<SFloat*> dose(NGT);//to store dose from all GPU cards
	vector<double> histArray(NGT);
	lambdaWoodcock.init(pht.ph.getP(), pht.mids.getP(), NVoxel);
	mat.init(0);
	//mat.dumpData("ZeusCrossSections.binary");
	initWaterCS("ZeusCrossSections.binary");
	//each thread takes care of one GPU
#ifdef USE_OPENMP
#pragma omp parallel num_threads(NGT)
#endif
	{
		int it = omp_get_thread_num();
		double fHMax = fNSIMU/NGT;
		int NPerBath = int(fHMax / 48);
		fHMax = NPerBath * 48;
		double percentHist = fHMax / 100;
		int NPercent = 1;

		double hist = 0;
		PRNG rng;
		//PRNG rng_src;
		//rng_src.init(1234 + it); //init the RNG
		rng.init(1234 + it);

		Particle* pool = new Particle[300];
		ParticleStack pStack(100);
		do //calculating main loop, end when hist >= fHMax
		{
			int np = SourceHead_Sample(&rng, pool);
			hist++;
			
			for (int ip = 0; ip < np; ++ip)
			{
				ParticleR p;
				p.E = ZFloat(pool[ip].E);
				p.x = ZFloat(pool[ip].x);
				p.y = ZFloat(pool[ip].y);
				p.z = ZFloat(pool[ip].z);
				p.u = ZFloat(pool[ip].u);
				p.v = ZFloat(pool[ip].v);
				p.w = ZFloat(pool[ip].w);
				p.weight = ZFloat(pool[ip].weight);
				p.type = pool[ip].type;
				if (!pht.lineIn(p)) continue;
				while (true)
				{
					if (p.type == photon)
					{
						smartPhoton(p, pStack, rng);
					}
					else
					{
						electronRun(p, rng);
					}
					if (pStack.empty()) break;
					else
					{
						p = pStack.top();
						pStack.pop();
					}
				}
			}


			//print the speed information
			if (it == 0 && hist > NPercent*percentHist)
			{
				NPercent++;
				double time = rc.stop(true);
				double speed = hist / time;
				double rest = 0;
				if (fHMax > hist) rest = (fHMax - hist) / speed;
				else rest = 0;
				Log("CPU processed ------------------------ %3.1f%%,   speed = %d h/s\n", hist*100.0 / fHMax, int(speed));
				Log("Time escaped = %.1f min, left time expected = %.1f min", time / 60.0, rest / 60.0);
			}

		
		} while (hist < fHMax);

		//finish in this CPU thread
		histArray[it] = hist;
		if (it == 0) Log("GPU processed ------------------------ 100%%,   speed = %d h/s\n", int(hist / rc.stop(true)));
		if (it == 0) Log("\nWait all GPUs to finish their job...\n");
		delete[] pool;
#ifdef USE_OPENMP
#pragma omp critical
#endif
		{
			Log("In thread %d, Max stack depth = %d\n", it, pStack.maxDepth);
		}
	} //end openMP

	Log("All GPUs have finished their simulation job! Collecting dose...\n\n");

	//merge dose and dose^2 from all GPU cards
	double totHist = histArray[0];
	for (int i = 1; i < NGT; ++i)
	{
		totHist += histArray[i];
	}

	// Get the dose. The dose comes back as Gray per photon times normalization (the constant passed to the 
	// getDose() function below). The purpose of normalization is to converted to real dose 
	// by multiplying with twice the current source activity (because there are two photons per decay)
	// and the total beam on time. For a fresh 15,000 Ci source and 1 minute beam on time, the factor is 
	// 2 * 1.5e4 * 3.7e10 * 60 = 6.66e16. 
	// Here we simply use a normalization that will produce numerical values typically in the range 1-10.

	double norm = 1.0889e15 * SourceHead_BeamOnTime();
	pht.addDose(totHist, norm);

	string outname;
	cf.getValue("output file name", outname);
	outname += ".dose";
	pht.output(outname.c_str());

	SourceHead_Delete(); //release the source head resource safely
	ZeusData_unload();

	Log("Time statistics for main GPU:");
	Log("Mixed running time = %.1f minutes, total history number = %d", rc.stop(true) / 60.0, totHist);
	Log("The overall simulating speed = %d hist/sec\n\n", int(totHist / rc.stop(true)));
	Log("End log time = %s\n\n", Log.timeNow());
	Log("/##############################################################################/\n\n");
	Log("The ratio of secondary electron < 743 KeV is %f", sE / aE);
	Log.closeFile();
} //end executeJob