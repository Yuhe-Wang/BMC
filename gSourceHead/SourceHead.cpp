
#include "../Tools/Tools.h"
#include "SourceHead.h"
#include "gSourceHead.h"

int NSampleStack = 400;

namespace Constants
{
	/**
	* Electron mass energy equivalent [eV].
	* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?mec2mev|search_for=electron+mass">
	* Click here for the NIST reference, but note the units are MeV.</a>
	*/
	const double ELECTRON_MASS = 510998.910;

	/*! 1/ELECTRON_MASS in eV^-1*/
	const double INV_ELECTRON_MASS = 1.956951337e-6;

	/**
	* Classical electron radius [centimeters].
	* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?re|search_for=electron+radius">
	* Click here for the NIST reference, but note the units are in meters.</a>
	*/
	const double ELECTRON_RADIUS = 2.8179402894e-13;
	/**
	* Avogadro constant [1/mol].
	* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=Avogadro">
	* Click here for the NIST reference.</a>
	*/
	const double AVOGADRO = 6.02214179e23;
	/**
	* Speed of light in a vacuum [centimeters/sec].
	* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?c|search_for=speed+of+light">
	* Click here for the NIST reference, but note the units are [meters/sec].</a>
	*/
	const double SPEED_OF_LIGHT = 2.99792458e10;
	/**
	* Energy of 1 eV [J].
	* <a href="http://physics.nist.gov/cgi-bin/cuu/Value?tevj|search_for=electron+volt">
	* Click here for the NIST reference.</a>
	*/
	const double ELECTRON_VOLT = 1.602176487e-19;
	/**
	* The mathematical constant PI.
	* N[Pi,30] from Mathematica 7.
	*/
	const double PI = 3.14159265358979323846264338328;
	/**
	* The square root of two.
	* N[Sqrt[2],30] from Mathematica 7.
	*/
	const double SQRT2 = 1.41421356237309;
	/**
	* This constant, the max value of a double precision floating point number,
	* is returned by method LambdaMoller::lmbdamo when the input energy is zero or less,
	* and is used to initialize a variable in method LambdaPhoton::SetupWoodcock, before
	* using that variable to hold a minimum value.
	*/
	const double DBLMAX = 1.79769313486231e308;
	/**
	* Zero divides are avoided in DoseWorker::inters by setting quotient to this value.
	*/
	const double FLTMAX = 3.402823466e38;
	/**
	* From microsoft, smallest number such that 1.0+FLTEPSILON != 1.0.
	*/
	const double FLTEPSILON = 1.192092896e-07;
	/**
	* A very small number, used in LambdaPhoton::SetupWoodcock.
	*/
	const double EPS = 1.0e-10;
	/**
	* The length of the table in LambdaPhoton::SetupWoodcock.
	*/
	const int WOODCOCKARRAYSIZE = 4096;
	/**
	* Number of leaf pairs in the MLC. Only the middle 30 are real.
	*/
	const int NLEAF_PAIRS_REAL = 30;
	/**
	* Number of leaf pairs, on each side, which are not real.
	*/
	const int N_LEAF_PAIRS_FAKE = 4;
	/**
	* Number of leaf pairs in the MLC. Only the middle 30 are real.
	*/
	const int NLEAF_PAIRS_TOTAL = 38;
	/**
	* Distance in CM from source to isocenter.
	*/
	const double SOURCE_TO_ISO_DISTANCE = 105;
	/**
	* Leaf thickness in CM.
	*/
	const double LEAF_THICKNESS = 1.05;
	/**
	* Leaf length in CM.
	*/
	const double LEAF_LENGTH = 100.0;
	/**
	* MLC model, the radius of the cylinder representing the inner leaves.
	*/
	const double MLC_INNER_RADIUS = 41;
	/**
	* MLC model, the radius of the cylinder representing the inner leaves.
	*/
	const double MLC_OUTER_RADIUS = 50;
	/**
	* The middle of the MLC, in the Y dimension.
	*/
	const double MLC_MIDDLE = (MLC_INNER_RADIUS + MLC_OUTER_RADIUS) / 2.0;
	/**
	* Minimum X location of MLC.
	*/
	const double MLC_MINX = -15;
	/**
	* Maximum X location of MLC.
	*/
	const double MLC_MAXX = 15;
	/**
	* MLC design constant
	*/
	const double INNER_RIGHT_FOCUS = 0.03;
	/**
	* MLC design constant
	*/
	const double INNER_GAP_FOCUS = -0.035;
	/**
	* MLC design constant
	*/
	const double INNER_Z_SHIFT = -0.015;
	/**
	* MLC design constant
	*/
	const double OUTER_LEFT_FOCUS = -0.03;
	/**
	* MLC design constant
	*/
	const double OUTER_GAP_FOCUS = 0.015;
	/**
	* MLC design constant
	*/
	const double OUTER_Z_SHIFT = 0.035;
}

typedef enum {
	/*! An enumeration used to identify the various materials in the treatment head model of the ViewRay delivery system. */
	Tungsten = 0,
	Air = 1,
	GradientCoil = 2,
	RFShield = 3,
	TxPCB = 4,
	TxFR4 = 5,
	UnknownMaterial = 6
} Medium;


SourceHead* gSH = NULL;
bool SourceHead_Init(const void* cf)
{
	gSH = new SourceHead;
	return gSH->init((ConfigFile*)cf);
}
void SourceHead_Delete()
{
	if (gSH) delete gSH;
}
double SourceHead_BeamOnTime() //get the total beam on time, which is used in dose normalization
{
	if (gSH) return gSH->beamOnTime();
	else return 0;
}
void SourceHead_GetIsoCenter(double* px, double* py, double* pz)
{
	if (gSH) gSH->isoCenter(px, py, pz);
}
void SourceHead_GetPrescrition(double* pDose, int* fraction)
{
	if (gSH) gSH->prescription(pDose, fraction);
}
double SourceHead_Emax()
{
	if (gSH) return gSH->maxEnergy();
	else return 0;
}
bool SourceHead_ConvertBeams(const char* txtbeams, const char* beamsPath)
{
	//this is an independent function that we don't need to call SourceHead_Init
	SourceHead sh;
	SourceHead::LeafAdjustment = -0.07;//it's safe since this interface is only designed for converting beams file
	if (!sh.getViewRayBeams(string(txtbeams))) return false;
	if (!sh.exportViewRayBeams(string(beamsPath))) return false;
	return true;
}
/******************************************* Class Implementation ******************************************/

/***********************************************************************************************************/


/****************************************** SEGMENT & BEAM *************************************************/
const int     SEGMENT::NLeaf = Constants::NLEAF_PAIRS_REAL;
const double  SEGMENT::LeafWidth = Constants::LEAF_THICKNESS; //unit cm
const double  SEGMENT::DX = NLeaf*LeafWidth; //unit cm

SEGMENT::SEGMENT() :primAT(NULL), scatAT(NULL), primMask(NULL), scatMask(NULL)
{
	leafPositions.resize(NLeaf);
	for (int i = 0; i < NLeaf; ++i) //set all the leaves in close state
	{
		leafPositions[i] = make_pair(0, 0);
	}
}
SEGMENT::SEGMENT(ConfigFile *cf) :primAT(NULL), scatAT(NULL), primMask(NULL), scatMask(NULL)
{
	leafPositions.resize(NLeaf);
	for (int i = 0; i < NLeaf; ++i) //set all the leaves in close state
	{
		leafPositions[i] = make_pair(0, 0);
	}
	load(cf);
}
bool SEGMENT::load(ConfigFile *cf) //true means loaded successfully
{
	vector<double> values;
	if (!cf->getValue("on time", onTime)) return false;
	if (!cf->getValue("leaf position", values)) return false;
	bool status = setPos(values);
	if (!status) return status;
	while (cf->getValue("leaf position", values, true)) //if we can find more overlap definition
	{
		status = setPos(values);
	}
	return status;
}
bool SEGMENT::loadViewRaySegment(const string& vrs) //load segment from ViewRay's format
{
	onTime = -1;
	size_t start = 0;
	for (int i = 0; i < 3; ++i)	start = vrs.find("\n", start + 1);
	sscanf(vrs.c_str() + start, "     Total Beam-On Time (s): %lf", &onTime);
	if (onTime < 0)
	{
		printf("Error: cannot read the beam on time for the segment.\n");
		return false;
	}
	for (int i = 0; i < 3; ++i)	start = vrs.find("\n", start + 1);
	double left = -1, right = -1;
	int leafIndex = 0;
	int sr = -100;
	for (int i = 0; i < NLeaf; ++i)
	{
		sscanf(vrs.c_str() + start, "      %d      %lf                %lf\n", &leafIndex, &left, &right);
		start = vrs.find("\n", start + 1); //update the search header index
		if (right - left + SourceHead::LeafAdjustment >= 0) //may apply an overall adjustment on leaf distance
		{
			left -= SourceHead::LeafAdjustment*0.5;
			right += SourceHead::LeafAdjustment*0.5;
		}
		leafPositions[leafIndex - 1] = make_pair(left, right);
	}
	return true;
}
bool SEGMENT::setPos(vector<double> &values)
{
	if (values.size() != 4) exitApp("The number of leaf position parameters isn't 4!");
	int start = int(values[0]);
	int end = int(values[1]);
	if (start<0 || end>NLeaf - 1) return false;
	if (values[3] - values[2] + SourceHead::LeafAdjustment >= 0) //may apply an overall adjustment on leaf distance
	{
		values[2] -= SourceHead::LeafAdjustment*0.5;
		values[3] += SourceHead::LeafAdjustment*0.5;
	}
	pair<double, double> pos = make_pair(values[2], values[3]);
	for (int i = start; i <= end; ++i) leafPositions[i] = pos;
	return true;
}

bool BEAM::load(ConfigFile *cf)
{
	segments.resize(0);
	vector<double> values;
	if (!cf->getValue("gantry angle", gantryAngle)) return false;
	cosPhi = cos(gantryAngle*PI / 180);
	sinPhi = sin(gantryAngle*PI / 180);
	if (!cf->getValue("head index", headIndex)) return false;

	ConfigFile *seg_cf;
	cf->resetSearchIndex();
	for (int NSeg = 1;; ++NSeg)
	{
		seg_cf = cf->getBlock("segment", true);
		if (NULL == seg_cf)
		{
			if (1 == NSeg) return false; //cannot even find the first material
			else break; //at least loaded one material
		}
		segments.push_back(SEGMENT(seg_cf));
	}

	return true;
}
bool BEAM::loadViewRayBeam(const string& vrb)//load beam configuration from ViewRay's plan
{
	//find beam index
	int angle = 0;
	sscanf(vrb.c_str(), "Beam %d:\n    Angle: %d", &headIndex, &angle);
	if (angle < 0 || angle>360)
	{
		printf("Error: incorrect beam angle!\n");
		return false;
	}
	gantryAngle = angle;
	cosPhi = cos(gantryAngle*PI / 180);
	sinPhi = sin(gantryAngle*PI / 180);

	segments.resize(0);//make sure previous segments are released
	//search the definition of each segment
	vector<string> vrSeg;
	size_t spos = vrb.find("Segment Id");
	size_t epos;
	do
	{
		epos = vrb.find("Segment Id", spos + 1);
		vrSeg.push_back(vrb.substr(spos, epos - spos));
		spos = epos;
	} while (epos != string::npos);

	segments.resize(vrSeg.size());
	for (size_t i = 0; i < vrSeg.size(); ++i) segments[i].loadViewRaySegment(vrSeg[i]);

	return true;
}
/***********************************************************************************************************/


/********************************************** CoPhaseSpace ***********************************************/
CoPhaseSpace::CoPhaseSpace(const char *fname) {
	const static char *func = "CoPhaseSpace::CoPhaseSpace";
	_dataFileName = fname;
	_isValid = false;
	_fluenceP = 0;
	_prob1 = 0;
	_fluenceS = 0;
	_pScoarse = 0;
	_pScoarse1 = 0;

	//
	// *** Open phase space file for reading
	//
	ifstream in(fname, std::ios::binary);
	if (!in) {
		exitApp("failed to open the phase space file");
	}

	//
	// *** Read header info
	//
	in.read((char*)&_Ncase, sizeof(double));
	in.read((char*)&_weight, sizeof(double));
	//
	// Tony has used the exact gradient coil geometry, while in the treatment head simulation source 
	// I used an approximation that ended up attenuating too much by 0.55%. I could of course fix the 
	// treatment head simulation source, but that would require also changing the calibration factor in the system.
	// To avoid such changes, we reduce the weight of the particles from the phase space source by 0.9945.
	//
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	_weight *= 0.9945;
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	in.read((char*)&_Nprim1, sizeof(double));
	in.read((char*)&_Nprim2, sizeof(double));
	in.read((char*)&_Nscat, sizeof(double));
	in.read((char*)&_xmin, sizeof(double));
	in.read((char*)&_xmax, sizeof(double));
	in.read((char*)&_ymin, sizeof(double));
	in.read((char*)&_ymax, sizeof(double));
	in.read((char*)&_Emin, sizeof(double));
	in.read((char*)&_Emax, sizeof(double));
	in.read((char*)&_delxP, sizeof(double));
	in.read((char*)&_delyP, sizeof(double));
	in.read((char*)&_delxS, sizeof(double));
	in.read((char*)&_delyS, sizeof(double));
	in.read((char*)&_delxS1, sizeof(double));
	in.read((char*)&_delyS1, sizeof(double));
	in.read((char*)&_delxSc, sizeof(double));
	in.read((char*)&_delySc, sizeof(double));
	in.read((char*)&_sourcePos, sizeof(double));
	in.read((char*)&_phspPos, sizeof(double));
	_phspPos -= 0.0001;  // this is to avoid taking an extra step of zero length during the head simulation
	in.read((char*)&_Nx, sizeof(int));
	in.read((char*)&_Ny, sizeof(int));
	in.read((char*)&_Nx1, sizeof(int));
	in.read((char*)&_Ny1, sizeof(int));
	in.read((char*)&_Nxb, sizeof(int));
	in.read((char*)&_Nyb, sizeof(int));
	in.read((char*)&_Nxb1, sizeof(int));
	in.read((char*)&_Nyb1, sizeof(int));
	in.read((char*)&_Nxbc, sizeof(int));
	in.read((char*)&_Nybc, sizeof(int));
	in.read((char*)&_nE, sizeof(int));
	if (in.fail()) {
		ZeusLog("%s: Failed to read header info from %s\n", func, fname); return;
	}
	if (false) //may output these for debugging
	{
		ZeusLog("Ncase=%g w=%g Np1=%g Np2=%g Nscat=%g xmin=%g xmax=%g ymin=%g ymax=%g Emin=%g Emax=%g spos=%g phsppos=%g\n",
			_Ncase, _weight, _Nprim1, _Nprim2, _Nscat, _xmin, _xmax, _ymin, _ymax, _Emin, _Emax, _sourcePos, _phspPos);
		ZeusLog("Nx=%d Ny=%d Nx1=%d Ny1=%d nE=%d\n", _Nx, _Ny, _Nx1, _Ny1, _nE);
		ZeusLog("Particles per case: %g\n", _particlesPerCase);
	}
	_particlesPerCase = (_Nprim1 + _Nprim2 + _Nscat) / _Ncase;


	//
	// *** Allocate memory and read primary fluence
	//
	int N = _Nx*_Ny;
	_fluenceP = new double[N]; in.read((char*)_fluenceP, N*sizeof(double));
	if (in.fail()) {
		ZeusLog("%s: Failed to read primiary fluence 1 from %s\n", func, fname); return;
	}
	_prob1 = new double[N]; in.read((char*)_prob1, N*sizeof(double));
	if (in.fail()) {
		ZeusLog("%s: Failed to read primiary fluence 2 from %s\n", func, fname); return;
	}
	for (int j = 0; j < N; ++j) {
		_fluenceP[j] += _prob1[j];
		if (_fluenceP[j] > 0) _prob1[j] /= _fluenceP[j];
	}

	//
	// *** Allocate memory and read scatter fluence
	//
	_fluenceS = new double[N]; in.read((char*)_fluenceS, N*sizeof(double));
	if (in.fail()) {
		ZeusLog("%s: Failed to read scatter fluence from %s\n", func, fname); return;
	}

	//
	// *** Now comes the data for the _Nx1 x _Ny1 tiles 
	//
	int N1 = _Nx1*_Ny1; int M = _Nxb*_Nyb; int M1 = _Nxb1*_Nyb1; int Mc = _Nxbc*_Nybc;

	//
	// *** Resize vectors to the number of tiles
	//
	_scatSpectrum.resize(N1, 0);
	_thetaPhiPflu.resize(N1, 0);
	_thetaPhiSflu.resize(N1, 0);
	_thetaPhiSflu1.resize(N1, 0);
	_thetaPhiSfluC.resize(N1, 0);
	_pScoarse = new double[N1];
	_pScoarse1 = new double[N1];

	//
	// *** Loop over tiles
	//
	for (int j = 0; j < N1; ++j) {
		//
		// *** Allocate and read scatter spectrum for tile j
		//
		_scatSpectrum[j] = new unsigned short[2 * _nE];
		in.read((char*)_scatSpectrum[j], 2 * _nE*sizeof(unsigned short));
		if (in.fail()) {
			ZeusLog("%s: Failed reading spectrum for tile %d from %s\n", func, j, fname); return;
		}

		//
		// *** Allocate and read primary angular distribution for tile j
		//
		_thetaPhiPflu[j] = new unsigned short[2 * M];
		in.read((char*)_thetaPhiPflu[j], 2 * M*sizeof(unsigned short));
		if (in.fail()) {
			ZeusLog("%s: Failed reading theta/phi primary fluence for tile %d from %s\n", func, j, fname); return;
		}

		//
		// *** Allocate and read scatter angular distribution for tile j on the _Nxb x _Nyb grid
		//
		_thetaPhiSflu[j] = new unsigned short[2 * M];
		in.read((char*)_thetaPhiSflu[j], 2 * M*sizeof(unsigned short));
		if (in.fail()) {
			ZeusLog("%s: Failed reading theta/phi scatter fluence for tile %d from %s\n", func, j, fname); return;
		}

		//
		// *** Allocate and read scatter 1 angular distribution for tile j on the _Nxb1 x _Nyb1 grid
		//
		_thetaPhiSflu1[j] = new unsigned short[2 * M1];
		in.read((char*)_thetaPhiSflu1[j], 2 * M1*sizeof(unsigned short));
		if (in.fail()) {
			ZeusLog("%s: Failed reading theta/phi scatter 1 fluence for tile %d from %s\n", func, j, fname); return;
		}

		//
		// *** Allocate and read scatter angular distribution for tile j on the coarser _Nxbc x _Nybc grid
		//
		_thetaPhiSfluC[j] = new unsigned short[2 * Mc];
		in.read((char*)_thetaPhiSfluC[j], 2 * Mc*sizeof(unsigned short));
		if (in.fail()) {
			ZeusLog("%s: Failed reading coarse theta/phi scatter fluence for tile %d from %s\n", func, j, fname); return;
		}

		//
		// *** Read probability for the scattered photon falling within the fine grid
		//
		in.read((char *)&_pScoarse[j], sizeof(double));
		if (in.fail()) {
			ZeusLog("%s: Failed reading coarse scatter probability for tile %d from %s\n", func, j, fname); return;
		}

		//
		// *** Read probability for the scattered photon falling within the fine 1 grid
		//
		in.read((char *)&_pScoarse1[j], sizeof(double));
		if (in.fail()) {
			ZeusLog("%s: Failed reading coarse scatter 1 probability for tile %d from %s\n", func, j, fname); return;
		}

	}

	//
	// *** Initialize variables used at run time for sampling particles
	//
	_pPrim1 = _Nprim1 / (_Nprim1 + _Nprim2);
	_dE = (_Emax - _Emin) / _nE;
	_dx = (_xmax - _xmin) / _Nx;
	_dy = (_ymax - _ymin) / _Ny;
	_dx1 = (_xmax - _xmin) / _Nx1;
	_dy1 = (_ymax - _ymin) / _Ny1;
	//_dxti0 = ((double)_Nx1)/(_xmax - _xmin);
	//_dxti1 = -(_xmin*_Nx1)/(_xmax - _xmin);
	//_dyti0 = ((double)_Ny1)/(_ymax - _ymin);
	//_dyti1 = -(_ymin*_Ny1)/(_ymax - _ymin);
	_dxti0 = ((double)_Nx1) / _xmax;
	_dxti1 = 0;
	_dyti0 = ((double)_Ny1) / _ymax;
	_dyti1 = 0;

	_minxP = -_delxP; _delxP *= 2. / _Nxb;
	_minyP = -_delyP; _delyP *= 2. / _Nyb;
	_minxS = -_delxS; _delxS *= 2. / _Nxb;
	_minyS = -_delyS; _delyS *= 2. / _Nyb;
	_minxS1 = -_delxS1; _delxS1 *= 2. / _Nxb1;
	_minyS1 = -_delyS1; _delyS1 *= 2. / _Nyb1;
	_minxSc = -_delxSc; _delxSc *= 2. / _Nxbc;
	_minySc = -_delySc; _delySc *= 2. / _Nybc;

	//
	// *** Mark the phase space source as valid
	//
	_isValid = true;
}

CoPhaseSpace::~CoPhaseSpace() {
	if (_fluenceP) delete[] _fluenceP;
	if (_prob1) delete[] _prob1;
	if (_fluenceS) delete[] _fluenceS;
	if (_pScoarse) delete[] _pScoarse;
	if (_pScoarse1) delete[] _pScoarse1;
	for (int j = 0; j < (int)_scatSpectrum.size(); ++j) if (_scatSpectrum[j]) delete[] _scatSpectrum[j];
	for (int j = 0; j < (int)_thetaPhiPflu.size(); ++j) if (_thetaPhiPflu[j]) delete[] _thetaPhiPflu[j];
	for (int j = 0; j < (int)_thetaPhiSflu.size(); ++j) if (_thetaPhiSflu[j]) delete[] _thetaPhiSflu[j];
	for (int j = 0; j < (int)_thetaPhiSflu1.size(); ++j) if (_thetaPhiSflu1[j]) delete[] _thetaPhiSflu1[j];
	for (int j = 0; j < (int)_thetaPhiSfluC.size(); ++j) if (_thetaPhiSfluC[j]) delete[] _thetaPhiSfluC[j];
}

#define IROUND(x) int( (x > 0.0f) ? floor(x + 0.5f) : ceil(x - 0.5f) );
void CoPhaseSpace::setMask(double x1, double x2, double y1, double y2, char *mask) {
	double dxi = 1. / _dx; double dyi = 1. / _dy;
	int ix1 = IROUND(dxi*(x1 - _xmin));
	if (ix1 < 0) ix1 = 0; else if (ix1 > _Nx - 1) ix1 = _Nx - 1;
	int ix2 = IROUND(dxi*(x2 - _xmin));
	if (ix2 < 0) ix2 = 0; else if (ix2 > _Nx - 1) ix2 = _Nx - 1;
	int iy1 = IROUND(dyi*(y1 - _ymin));
	if (iy1 < 0) iy1 = 0; else if (iy1 > _Ny - 1) iy1 = _Ny - 1;
	int iy2 = IROUND(dyi*(y2 - _ymin));
	if (iy2 < 0) iy2 = 0; else if (iy2 > _Ny - 1) iy2 = _Ny - 1;
	for (int iy = iy1; iy <= iy2; ++iy) for (int ix = ix1; ix <= ix2; ++ix) mask[ix + iy*_Nx] = 1;
}
/***********************************************************************************************************/


/******************************************** MLCattenuation ***********************************************/
MLCattenuation::MLCattenuation(const char* fname) : _data(0), _isValid(false) {
	const static char *func = "MLCattenuation::MLCattenuation";
	//
	// *** Open file for reading and check for error
	//
	ifstream in(fname, std::ios::binary);
	if (!in) {
		ZeusLog("%s: failed to open %s\n", func, fname);
		return;
	}
	//
	// *** Read header info (ymin, ymax, number of bins) and check for errors
	//
	float ymin, ymax; int Ny;
	in.read((char *)&ymin, sizeof(float));
	in.read((char *)&ymax, sizeof(float));
	in.read((char *)&Ny, sizeof(int));
	if (in.fail()) {
		ZeusLog(": failed reading header from %s\n", func);
		return;
	}
	if (ymin >= ymax || Ny < 2) {
		ZeusLog("%s: found bad header ymin=%g ymax=%g Ny=%d in %s\n", func, ymin, ymax, Ny, fname);
		return;
	}
	//
	// *** Allocate and read attenuation data
	//
	_data = new unsigned short[Ny];
	in.read((char *)_data, Ny*sizeof(unsigned short));
	if (in.fail()) {
		ZeusLog("%s: failed reading attenuation data from %s\n", func, fname);
		delete[] _data; _data = 0; return;
	}
	//
	// *** Data was successfully loaded, so initialize interpolation variables and mark as valid.
	//
	_ymin = ymin; _ymax = ymax; _Ny = Ny;
	_attLeft = 6; _attRight = 6;
	_dyi0 = ((double)Ny) / (_ymax - _ymin);
	_dyi1 = -_ymin*_dyi0;
	_isValid = true;
	ZeusLog("%s: successfully initialized from data in %s\n", func, fname);
}

MLCattenuation::~MLCattenuation() {
	if (_data) delete[] _data;
}
/***********************************************************************************************************/


/********************************************** SourceHead *************************************************/
CoPhaseSpace* SourceHead::phsp = NULL;
int SourceHead::NReject[3] = { 1, 10, 34 };
double SourceHead::LeafAdjustment = 0.0;

bool SourceHead::init(ConfigFile *cf)
{
	RunTimeCounter rc;

	cf->getValue("Sample Stack Depth", NSampleStack);
	cf->getValue("Leaf spacing adjustment", LeafAdjustment);
	cf->getValue("Energy scale factor", kE);

	string dir;
	cf->getValue("DataDir", dir);
	int lastInd = (int)dir.size() - 1;
	if (dir[lastInd] != '/' && dir[lastInd] != '\\') dir += '/';

	string File;
	cf->getValue("CoPhaseSpace", File);
	phsp = new CoPhaseSpace(File.c_str());

	if (!cf->getValue("Beams file", File) || !getViewRayBeams(File)) //try to load from ViewRay's beams file
	{
		//use customized definition instead
		vector<double> values;
		if (!cf->getValue("isocenter", values) || values.size() != 3) return false;
		isoX = values[0];
		isoY = values[1];
		isoZ = values[2];
		//search for the definition of beams
		ConfigFile *beamCF;
		cf->resetSearchIndex();
		for (int NBeam = 1;; ++NBeam)
		{
			beamCF = cf->getBlock("Beam", true);
			if (NULL == beamCF)
			{
				if (1 == NBeam) exitApp("cannot load any beam parameter!"); //cannot even find the first material
				else  break;//at least loaded one beam
			}
			beams.resize(NBeam);
			beams[NBeam - 1].load(beamCF);
		}
	}

	//force overwriting the isocenter coordinates relative to patient center
	string overwrite("no");
	cf->getValue("overwrite isocenter", overwrite);
	if (0 == overwrite.compare("yes"))
	{
		vector<double> values;
		if (!cf->getValue("isocenter", values) || values.size() != 3) return false;
		isoX = values[0];
		isoY = values[1];
		isoZ = values[2];
	}

	//this should come from configure file
	NReject[0] = 1;
	NReject[1] = 10;
	NReject[2] = 34;
	//customize the sampling alias table for each segment
	int Nx = phsp->getNx();
	int Ny = phsp->getNy();
	// The margins around the MLC opening. These should come from a config file!!!
	double margin1 = 0.35; //primary photon small margin
	double margin2 = 0.90; //scatter photon small margin
	double margin3 = 1.5;  //scatter photon large margin
	// A gap margin. This one is used for the closed leaves
	double gapMargin = 0.15;
	// Allocate memory for masks for MLC opening plus margins
	char *mask1 = new char[Nx*Ny];
	char *mask1a = new char[Nx*Ny];
	char *mask2 = new char[Nx*Ny];
	char *mask3 = new char[Nx*Ny];
	// The probabilities for the fluence bins depending on mask
	double probIn = 1.0 / NReject[1];
	double probOut = 1.0 / NReject[2];

	// Get the virtual source position from the phase space data.
	double sourcePos = fabs(phsp->getSourcePosition());
	// Scaling factor from the isocenter plane to the phase space plane.
	double scale1 = 1 - phsp->getPhspPosition() / sourcePos;
	// The minimum y-position across the leaves
	double ymin = -0.5*SEGMENT::LeafWidth*SEGMENT::NLeaf;
	// Get the primary and scatter fluence data from the phase space.
	const double *fluP = phsp->getPrimaryFluence();
	const double *fluS = phsp->getScatterFluence();
	// Storage for the fluences modified by the MLC
	vector<double> sampP(Nx*Ny);
	vector<double> sampS(Nx*Ny);
	// The binSampler will be used for preparing alias tables.
	BinSampler binSampler;

	vector<double> segFluence;
	int NBeam = (int)beams.size();
	for (int ib = 0; ib < NBeam; ++ib)
	{
		int NSeg = (int)beams[ib].segments.size();
		for (int is = 0; is < NSeg; ++is)
		{
			SEGMENT &s = beams[ib].segments[is];
			//set mask array zero first
			for (int j = 0; j<Nx*Ny; ++j) { mask1[j] = 0; mask1a[j] = 0; mask2[j] = 0; mask3[j] = 0; }
			//set mask based on the leaves position
			for (int i = 0; i < s.NLeaf; ++i)
			{
				// *** The min/max y-location for this leaf pair in the phase-space plane
				double y1 = scale1*(ymin + s.LeafWidth*i);
				double y2 = scale1*(ymin + s.LeafWidth*(i + 1));
				// *** Is this leaf pair open?
				if (s.leafPositions[i].second - s.leafPositions[i].first > 0.01) { //unit cm
					// *** Yes, it is. Compute min/max x-position of the opening in the phase-space plane and then set the primary and scatter masks
					//     for the opening plus margin(s)
					double x1 = scale1*s.leafPositions[i].first; double x2 = scale1*s.leafPositions[i].second;      // Min/max x.
					phsp->setMask(x1 - margin1, x2 + margin1, y1 - margin1, y2 + margin1, mask1);    // Set primary mask plus margin
					phsp->setMask(x1 - margin2, x2 + margin2, y1 - margin2, y2 + margin2, mask2);    // Set scatter mask plus small margin
					phsp->setMask(x1 - margin3, x2 + margin3, y1 - margin3, y2 + margin3, mask3);    // Set scatter mask plus large margin
				}
				else {
					// *** No, its closed.
					//     Compute position of closing in the phase-space plane and set masks
					double x = 0.5*scale1*(s.leafPositions[i].first + s.leafPositions[i].second);                    // Position x
					phsp->setMask(x - gapMargin, x + gapMargin, y1 - gapMargin, y2 + gapMargin, mask1);          // Set primary mask plus gap margin
					phsp->setMask(x - gapMargin, x + gapMargin, y1 - gapMargin, y2 + gapMargin, mask2);          // Set scatter mask plus gap margin
					phsp->setMask(x - 2 * gapMargin, x + 2 * gapMargin, y1 - 2 * gapMargin, y2 + 2 * gapMargin, mask1a); // Set primary mask plus larger gap margin
					phsp->setMask(x - 2 * gapMargin, x + 2 * gapMargin, y1 - 2 * gapMargin, y2 + 2 * gapMargin, mask3);  // Set scatter mask plus larger gap margin
				}
			}

			// *** Allocate memory for alias tables and masks for this segment
			s.primAT = new unsigned short[2 * Nx*Ny];
			s.scatAT = new unsigned short[2 * Nx*Ny];
			s.primMask = new char[Nx*Ny];
			s.scatMask = new char[Nx*Ny];

			double sumP = 0, sumS = 0, sum0 = 0;
			//set the mask of each bin and get the fluence modified by MLC
			for (int j = 0; j < Nx*Ny; ++j) {
				// Accumulate total fluence
				sum0 += fluP[j] + fluS[j];
				// Set primary sampling probability and primary mask for this bin
				sampP[j] = fluP[j];
				if (mask1[j]) s.primMask[j] = 0;   // i.e., fluence bin is within opening plus margin
				else {
					// fluence bin is outside of opening. 
					if (mask1a[j]) {
						// fluence bin is within the larger margin 
						s.primMask[j] = 1;
						sampP[j] *= probIn;
					}
					else {
						// fluence bin is outside of the larger margin
						s.primMask[j] = 2;
						sampP[j] *= probOut;
					}
				}
				// Set scatter sampling probability and scatter mask for this bin
				sampS[j] = fluS[j];
				if (mask2[j]) s.scatMask[j] = 0;  // i.e., fluence bin is within opening plus small margin
				else {
					// fluence bin is outside of opening plus small margin
					if (mask3[j]) {
						// fluence bin is within the larger margin
						s.scatMask[j] = 1;
						sampS[j] *= probIn;
					}
					else {
						// fluence bin is outside of the larger margin
						s.scatMask[j] = 2;
						sampS[j] *= probOut;
					}
				}
				// Accumulate sampling probabilities
				sumP += sampP[j]; sumS += sampS[j];
			}

			binSampler.initialize(sampP);
			// Copy primary position sampling table to unsigned shorts (to save memory). This still gives us more than enough resolution in sampling probability
			const vector<pair<double, int> > &table = binSampler.getTable();
			for (int j = 0; j < Nx*Ny; ++j) {
				s.primAT[2 * j] = (unsigned short)(65535.*table[j].first);
				s.primAT[2 * j + 1] = table[j].second;
			}
			// Initialize scatter position sampling table
			binSampler.initialize(sampS);
			// Copy scatter position sampling table to unsigned shorts (to save memory). This still gives us more than enough resolution in sampling probability
			const vector<pair<double, int> > &table1 = binSampler.getTable();
			for (int j = 0; j < Nx*Ny; ++j) {
				s.scatAT[2 * j] = (unsigned short)(65535.*table1[j].first);
				s.scatAT[2 * j + 1] = table1[j].second;
			}

			// Number of photons to sample for this segment per call to getParticles().
			s.nSample = (sumP + sumS)*phsp->getParticlesPerCase() / sum0;
			// The probability to sample a primary photon
			s.pPrim = sumP / (sumP + sumS);
			s.beamID = ib;
			//s.fluence = s.onTime*(sumP + sumS); //this is what I thought
			s.fluence = s.onTime; //but Zeus use this
			segFluence.push_back(s.fluence);
			segIndex.push_back(make_pair(ib, is));
		}
	}

	segSampler.initialize(segFluence);

	delete[] mask1;
	delete[] mask1a;
	delete[] mask2;
	delete[] mask3;

	ZeusLog("It takes %f seconds to initiate the SourceHead\n\n", rc.stop());

	return true;
}

bool SourceHead::getViewRayBeams(const string& file) // get beam definition from ViewRay's plan
{
	//let's find the file extension
	string ext = file.substr(file.find_last_of("."), string::npos);
	if (ext == ".txt") // use ViewRay plan overview, which may change due to ViewRay's update
	{
		FILE *fp = fopen(file.c_str(), "r");
		if (NULL == fp)
		{
			printf("Cannot open %s. Will try customized beam configuration.\nPress enter to continue...\n", file.c_str());
			getchar();
			return false;
		}
		string vrPlan;
		char ch = 0;
		do
		{
			ch = fgetc(fp);
			vrPlan += ch;
		} while (EOF != ch);
		fclose(fp);

		size_t spos = 0, epos = 0;

		//find the prescription dose and fraction number
		spos = vrPlan.find("Total Prescription Dose");
		spos = vrPlan.find(":", spos) + 1;
		sscanf(vrPlan.c_str() + spos, "%lf", &prescriptionDose);

		spos = vrPlan.find("Prescription Dose Per Fraction");
		spos = vrPlan.find(":", spos) + 1;
		double fractionDose;
		sscanf(vrPlan.c_str() + spos, "%lf", &fractionDose);
		treatmentFraction = (int)round(prescriptionDose / fractionDose);

		//find the ISO center
		spos = vrPlan.find("Isocenter Coordinate (cm)");
		sscanf(vrPlan.c_str() + spos, "Isocenter Coordinate (cm): (%lf,%lf,%lf)", &isoX, &isoY, &isoZ);
		//divide the setting to groups
		vector<string> vrGroup;
		spos = vrPlan.find("Planning Beam Angle Group");
		do
		{
			epos = vrPlan.find("Planning Beam Angle Group", spos + 1);
			vrGroup.push_back(vrPlan.substr(spos, epos - spos));
			spos = epos;
		} while (epos != string::npos);

		vector<string> vrBeam;
		//divide each group to beams
		for (size_t i = 0; i < vrGroup.size(); ++i)
		{
			size_t b1 = vrGroup[i].find("Beam 1");
			size_t b2 = vrGroup[i].find("Beam 2");
			size_t b3 = vrGroup[i].find("Beam 3");
			string sb1 = vrGroup[i].substr(b1, b2 - b1);
			string sb2 = vrGroup[i].substr(b2, b3 - b2);
			string sb3 = vrGroup[i].substr(b3, string::npos - b3);
			if (string::npos == sb1.find("Angle: N/A")) vrBeam.push_back(sb1);
			if (string::npos == sb2.find("Angle: N/A")) vrBeam.push_back(sb2);
			if (string::npos == sb3.find("Angle: N/A")) vrBeam.push_back(sb3);
		}

		beams.resize(vrBeam.size());
		for (size_t i = 0; i < beams.size(); ++i)
		{
			beams[i].loadViewRayBeam(vrBeam[i]);
		}
		//need to change the coordinates of ISO center. This may vary depending on the verison of ViewRay software
		double isoT = isoY;
		isoY = -isoZ;
		isoZ = isoT;
		return true;
	}
	else if (ext == ".beams")
	{
		FILE *fp = fopen(file.c_str(), "rb");
		if (NULL == fp)
		{
			printf("Cannot open %s. Will try customized beam configuration.\nPress enter to continue...\n", file.c_str());
			getchar();
			return false;
		}

		unsigned int Major_Header, Minor_Header, Num_Of_Beams, Head, Num_Of_Segments;
		double GantryAngle, Isocenter_X, Isocenter_Y, Isocenter_Z, Beam_On_Time, Leaf_Position_Left, Leaf_Position_Right;
		fread(&Major_Header, sizeof(unsigned int), 1, fp);
		fread(&Minor_Header, sizeof(unsigned int), 1, fp);
		fread(&Num_Of_Beams, sizeof(unsigned int), 1, fp);
		beams.resize(Num_Of_Beams);
		for (unsigned int i = 0; i < Num_Of_Beams; ++i) //for each beam
		{
			fread(&GantryAngle, sizeof(double), 1, fp);
			fread(&Isocenter_X, sizeof(double), 1, fp);
			fread(&Isocenter_Y, sizeof(double), 1, fp);
			fread(&Isocenter_Z, sizeof(double), 1, fp);
			fread(&Head, sizeof(unsigned int), 1, fp);
			fread(&Num_Of_Segments, sizeof(unsigned int), 1, fp);
			beams[i].gantryAngle = GantryAngle;
			beams[i].sinPhi = sin(GantryAngle*PI / 180);
			beams[i].cosPhi = cos(GantryAngle*PI / 180);
			beams[i].headIndex = Head;
			beams[i].segments.resize(Num_Of_Segments);
			for (unsigned int j = 0; j < Num_Of_Segments; ++j)
			{
				fread(&Beam_On_Time, sizeof(double), 1, fp);
				beams[i].segments[j].onTime = Beam_On_Time;
				for (int k = 0; k < SEGMENT::NLeaf; ++k)
				{
					fread(&Leaf_Position_Left, sizeof(double), 1, fp);
					fread(&Leaf_Position_Right, sizeof(double), 1, fp);
					if (Leaf_Position_Right - Leaf_Position_Left + SourceHead::LeafAdjustment >= 0) //may apply an overall adjustment on leaf distance
					{
						Leaf_Position_Left -= SourceHead::LeafAdjustment*0.5;
						Leaf_Position_Right += SourceHead::LeafAdjustment*0.5;
					}
					beams[i].segments[j].leafPositions[k] = make_pair(Leaf_Position_Left, Leaf_Position_Right);
				}
			}
		}
		isoX = Isocenter_X;
		isoY = Isocenter_Y;
		isoZ = Isocenter_Z;
		fclose(fp);
		return true;
	}

	return false;
}

bool SourceHead::exportViewRayBeams(const string& file) //only export in .beams format
{
	FILE *fp = fopen(file.c_str(), "wb");
	if (NULL == fp)
	{
		printf("Cannot create file named %s.", file.c_str());
		return false;
	}

	unsigned int Major_Header = 1, Minor_Header = 0, Num_Of_Beams = (unsigned int)beams.size(), Head, Num_Of_Segments;
	double GantryAngle, Isocenter_X, Isocenter_Y, Isocenter_Z, Beam_On_Time, Leaf_Position_Left, Leaf_Position_Right;
	fwrite(&Major_Header, sizeof(unsigned int), 1, fp);
	fwrite(&Minor_Header, sizeof(unsigned int), 1, fp);
	fwrite(&Num_Of_Beams, sizeof(unsigned int), 1, fp);
	Isocenter_X = isoX;
	Isocenter_Y = isoY;
	Isocenter_Z = isoZ;
	for (unsigned int i = 0; i < Num_Of_Beams; ++i) //for each beam
	{
		GantryAngle = beams[i].gantryAngle;
		Head = beams[i].headIndex;
		Num_Of_Segments = (unsigned int)beams[i].segments.size();
		fwrite(&GantryAngle, sizeof(double), 1, fp);
		fwrite(&Isocenter_X, sizeof(double), 1, fp);
		fwrite(&Isocenter_Y, sizeof(double), 1, fp);
		fwrite(&Isocenter_Z, sizeof(double), 1, fp);
		fwrite(&Head, sizeof(unsigned int), 1, fp);
		fwrite(&Num_Of_Segments, sizeof(unsigned int), 1, fp);

		for (unsigned int j = 0; j < Num_Of_Segments; ++j)
		{
			Beam_On_Time = beams[i].segments[j].onTime;
			fwrite(&Beam_On_Time, sizeof(double), 1, fp);

			for (int k = 0; k < SEGMENT::NLeaf; ++k)
			{
				Leaf_Position_Left = beams[i].segments[j].leafPositions[k].first;
				Leaf_Position_Right = beams[i].segments[j].leafPositions[k].second;
				fwrite(&Leaf_Position_Left, sizeof(double), 1, fp);
				fwrite(&Leaf_Position_Right, sizeof(double), 1, fp);
			}
		}
	}

	fclose(fp);
	return true;
}
SourceHead::~SourceHead()
{
	if (phsp) delete phsp;
}
/***********************************************************************************************************/

/************************************************GPU interface**********************************************/
double g_beamOnTime, g_isox, g_isoy, g_isoz, g_pDose;
int g_fraction;

void CoPhaseSpace4GPU::init(void* vps)
{
	CoPhaseSpace* ps = (CoPhaseSpace*)vps;
	for (int i = 0; i < 57600; ++i) _prob1[i] = (ZFloat)ps->_prob1[i];
	for (int i = 0; i < 225; ++i)
	{
		_pScoarse[i] = (ZFloat)ps->_pScoarse[i];
		_pScoarse1[i] = (ZFloat)ps->_pScoarse1[i];
		for (int j = 0; j < 2 * 254; ++j) _scatSpectrum[i][j] = ps->_scatSpectrum[i][j];
		for (int j = 0; j < 2 * 65536; ++j) _thetaPhiPflu[i][j] = ps->_thetaPhiPflu[i][j];
		for (int j = 0; j < 2 * 65536; ++j) _thetaPhiSflu[i][j] = ps->_thetaPhiSflu[i][j];
		for (int j = 0; j < 2 * 16384; ++j) _thetaPhiSflu1[i][j] = ps->_thetaPhiSflu1[i][j];
		for (int j = 0; j < 2 * 4096; ++j) _thetaPhiSfluC[i][j] = ps->_thetaPhiSfluC[i][j];
	}
}

void SEGMENT4GPU::init(void* vbeam, void* vsg)
{
	BEAM* beam = (BEAM*)vbeam;
	SEGMENT* sg = (SEGMENT*)vsg;
	nSample = (ZFloat)sg->nSample;
	pPrim = (ZFloat)sg->pPrim;
	const int nd = 2 * 240 * 240;
	for (int i = 0; i < nd; ++i)
	{
		primAT[i] = sg->primAT[i];
		scatAT[i] = sg->scatAT[i];
	}
	for (int i = 0; i < 240*240; ++i)
	{
		primMask[i] = sg->primMask[i];
		scatMask[i] = sg->scatMask[i];
	}

	gantryAngle = (ZFloat)beam->gantryAngle;
	cosPhi = (ZFloat)beam->cosPhi;
	sinPhi = (ZFloat)beam->sinPhi;
}

void getHostData(ConfigFile* cf, CoPhaseSpace4GPU*&hps, SEGMENT4GPU*&hseg, int& hnSeg, int& nMaxPhoton, vector<vector<pair<double, double> > >& aListofSegments)
{
	SourceHead_Init(cf);
	g_beamOnTime = SourceHead_BeamOnTime();
	SourceHead_GetIsoCenter(&g_isox, &g_isoy, &g_isoz);
	SourceHead_GetPrescrition(&g_pDose, &g_fraction);

	//init the phase space data
	hps = new CoPhaseSpace4GPU;
	hps->init(gSH->getPhaseSpace());

	//init the segment data
	int nb = (int)gSH->beams.size();
	for (int i = 0; i < nb; ++i)
	{
		hnSeg += (int)gSH->beams[i].segments.size();
	}
	hseg = new SEGMENT4GPU[hnSeg];
	hnSeg = 0;
	ZFloat dMaxPhoton = 0;
	const vector<pair<double, int> >& table = gSH->segSampler.getTable();
	for (int i = 0; i < nb; ++i)
	{
		int ns = (int)gSH->beams[i].segments.size();
		for (int j = 0; j < ns; ++j)
		{
			hseg[hnSeg].init(&(gSH->beams[i]), &(gSH->beams[i].segments[j]));
			hseg[hnSeg].wfirst = (ZFloat)table[hnSeg].first;
			hseg[hnSeg].wsecond = table[hnSeg].second;
			dMaxPhoton = max(dMaxPhoton, hseg[hnSeg].nSample);
			++hnSeg;
		}
	}
	nMaxPhoton = int(dMaxPhoton) + 1; // +1 for safe memory allocation

	aListofSegments = gSH->getListofSegments();

	SourceHead_Delete();
}

double SourceHead4GPU_BeamOnTime()
{
	return g_beamOnTime;
}

void SourceHead4GPU_GetIsoCenter(double* px, double* py, double* pz)
{
	*px = g_isox;
	*py = g_isoy;
	*pz = g_isoz;
}

void SourceHeadGPU_GetPrescription(double* pDose, int* fraction)
{
	*pDose = g_pDose;
	*fraction = g_fraction;
}