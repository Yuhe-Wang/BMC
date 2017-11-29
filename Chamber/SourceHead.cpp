
#include "SourceHead.h"

#define ZeusLog printf
int NSampleStack = 400;
void exitApp(const char *inf)
{
	ZeusLog("fatal error: %s", inf);
	getchar();
	exit(-1);
}

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

class Vector {
public:

	double _x;  //!< x-coordinate
	double _y;  //!< y-coordinate
	double _z;  //!< z-coordinate

	/*! \brief Default constructor */
	Vector() : _x(0), _y(0), _z(0) {}

	/*! \brief Construct from input values \a aX, \a aY and \a aZ */
	Vector(double aX, double aY, double aZ) : _x(aX), _y(aY), _z(aZ) {}

	/*! \brief Assignment operator */
	Vector &operator=(const Vector &v) { _x = v._x; _y = v._y; _z = v._z; return *this; };

	/*! \brief Add \a v to this vector and return the result */
	Vector operator+(const Vector &v) const { return Vector(_x + v._x, _y + v._y, _z + v._z); };

	/*! \brief Add vector \a v and assign the result to this vector */
	Vector &operator+=(const Vector &v) { _x += v._x; _y += v._y; _z += v._z; return *this; };

	/*! \brief Returns the scalar product of vector \a v and this vector */
	double operator*(const Vector &v) const { return _x*v._x + _y*v._y + _z*v._z; };

	/*! \brief Multiply all co-ordinates with \a f and assign the result to this vector */
	Vector &operator*=(double f) { _x *= f; _y *= f; _z *= f; return *this; };

	/*! \brief Multiply \a v with \a f and return the result */
	friend Vector operator*(double f, Vector &v) { return v*f; };

	/*! \brief Substract the vector \a v from this vector and return the result */
	Vector operator-(const Vector &v) const { return Vector(_x - v._x, _y - v._y, _z - v._z); };

	/*! \brief Substract \a v and assign the result to this vector */
	Vector &operator-=(const Vector &v) { _x -= v._x; _y -= v._y; _z -= v._z; return *this; };

	/*! \brief Return the multiplication of this vector with \a f */
	Vector operator*(double f) const { return Vector(_x*f, _y*f, _z*f); };

	/*! \brief Get the squared length of the vector */
	double getLengthSquared() const { return _x*_x + _y*_y + _z*_z; };

	/*! \brief Get the length of the vector */
	double getLength() const { return sqrt(getLengthSquared()); };

	/*! \brief Normalize the vector to unit length */
	void normalizeToUnitLength() {
		double length2 = getLengthSquared();
		if (length2) {
			double invLength = 1 / sqrt(length2);
			_x *= invLength; _y *= invLength; _z *= invLength;
		}
	};

	/*! \brief Rotate the vector by the polar angle Theta and the azimutha angle Phi */
	void rotate(double cosTheta, double sinTheta, double cosPhi, double sinPhi) {
		double sinz = _x*_x + _y*_y;
		if (sinz > 1e-20) {
			sinz = sqrt(sinz);
			double c1 = sinTheta / sinz;
			double cphi = _z*cosPhi;
			double cx = _x*cosTheta, cy = _y*cosTheta;
			double cx1 = cphi*_x - _y*sinPhi;
			double cy1 = cphi*_y + _x*sinPhi;
			_x = c1*cx1 + cx; _y = c1*cy1 + cy;
			_z = _z*cosTheta - sinz*sinTheta*cosPhi;
		}
		else { _x = sinTheta*cosPhi; _y = sinTheta*sinPhi; _z *= cosTheta; }
	};

};

class BinSampler{

public:

	/*! Construct an uninitialized instance, which can be later initialized using the initialize() function. */
	BinSampler() {};

	/*! Construct an instance from the probabilities \a probabilities.
	*
	*  The constructor calls the initialize() function to set up the alias table.
	*/
	BinSampler(const vector<double> &probabilities) { initialize(probabilities); };

	/*! Initialize the instance from the bin probabilities \a probabilities.
	*
	* Transforms the discrete probabilities into an alias table suitable for fast random sampling.
	* Pre-existing data (if any) is destroyed.
	* Return value is zero on success. No-zero return codes are
	* -1, if there are less than two bins (it does not make sense to randomly sample from a single probability)
	* -2, if some of the probabilities are negative (i.e., these are not really probabilities)
	* -3, if all probabilities are zero.
	*/
	int initialize(const vector<double> &probabilities);

	/*! Returns a randomly selected bin using the random number \a r.
	*
	* Note: to not waste CPU cycles for checking if the random number really is in the 0...1 interval,
	*       no check is being made that this is satisfied => result will not make sense if this condition
	*       is not satidfied and you may even crash your program if \a r is outside of 0...1.
	*/
	int sampleBin(double r) const {
		int nBin = (int)_table.size();
		if (nBin != 1)
		{
			double aBin = r*nBin;
			int bin = (int)aBin; aBin -= bin;
			return aBin < _table[bin].first ? bin : _table[bin].second;
		}
		else return 0;

	};

	/*! Get a reference to the internal alias table representation of the input data */
	const vector<pair<double, int> >& getTable() const { return _table; };

protected:

	/*! The alias table. */
	vector<pair<double, int> >  _table;

};
int BinSampler::initialize(const vector<double> &probabilities) {
	const static char* func = "BinSampler::initialize";
	_table.clear();
	int nbin = (int)probabilities.size();
	if (0 == nbin) {
		exitApp("no input in BinSampler::initialize");
	}
	if (1 == nbin){
		_table.resize(nbin);
		return 0;
	}
	_table.resize(nbin);
	double sum = 0;
	for (int j = 0; j < nbin; ++j) {
		if (probabilities[j] < 0) {
			ZeusLog("%s: found negative probability %g. All probabilities must be non-negative\n", func, probabilities[j]);
			return -2;
		}
		sum += probabilities[j];
		_table[j] = pair<double, int>(1, -1);
	}
	if (!sum) return -3;
	sum /= nbin;
	vector<double> p(probabilities);
	vector<int> low(nbin), high(nbin); int nl = 0, nh = 0;
	for (int i = 0; i<nbin; ++i) {
		if (p[i] > sum) high[nh++] = i;
		else if (p[i] < sum) low[nl++] = i;
		else { _table[i].first = 1; _table[i].second = i; }
	}
	double sumi = 1 / sum;
	while (nl > 0 && nh > 0) {
		int jl = low[--nl];
		int jh = high[--nh];
		double diff = sum - p[jl];
		p[jh] -= diff;
		_table[jl].first = p[jl] * sumi;
		_table[jl].second = jh;
		if (p[jh] < sum) low[nl++] = jh;
		else if (p[jh] > sum) ++nh;
		else { _table[jh].first = 1; _table[jh].second = jh; }
	}
	for (int i = 0; i < nbin; ++i) if (_table[i].second < 0) { _table[i].first = 1; _table[i].second = i; }
	return 0;
}

class CoPhaseSpace {

public:

	/*! Constructor. \a fname is the name of the phase space file in the appropriate format. */
	CoPhaseSpace(const char *fname);

	/*! Destructor */
	~CoPhaseSpace();

	/*! Is the object valid? */
	bool isValid() const { return _isValid; };

	/*! Sample one particle from the phase space.

	Inputs: \a rndm is the random number generator, \a x the position along the leaves, \a y the position across the leaves,
	\a isPrimary is true if the particle sampled is to be primary and false if scatter.
	Output: The energy is stored in \a E, the azimuthal and polar angles about the direction defined by the particle position and the
	virtual source position are stored in \a sint (the sine of the azimuthal angle) and \a phi.
	Return value: 0 on success, some non-zero error code on failure.
	*/
	int sample(PRNG &rng, double x, double y, bool isPrimary, double &E, double &sint, double &phi);

	/*! Sample one particle from the phase space.

	Similar to the above function, except that it also samples the partcle position, returning it in \a x, and sets \a u to the particle directon.
	\a type is set to zero for E=1.17 primary photons, to 1 for E=1.33 primary photons, and to 2 for scattered photons
	*/
	int sample(PRNG &rng, bool isPrimary, int bin, double &E, Vector &x, Vector &u, int &type);

	/*! Returns the virtual source position (z in the original KMC co-ordinate system, y in the ViewRay system). */
	double getSourcePosition() const { return _sourcePos; };

	/*! Returns the phsp plane position (z in the original KMC co-ordinate system, y in the ViewRay system). */
	double getPhspPosition() const { return _phspPos; };

	/*! Returns the primary particle fluence above the MLC (_Nx x _Ny bins) */
	const double* getPrimaryFluence() const { return _fluenceP; }

	/*! Returns the probability for the 1.17 MeV line of the Co-60 spectrum (_Nx x _Ny bins) */
	const double* getProbability1() const { return _prob1; }
	/*! Returns the scattered particle fluence above the MLC (_Nx x _Ny bins) */
	const double* getScatterFluence() const { return _fluenceS; }

	/*! Returns the number of x-bins used for the particle fluence */
	int getNx() const { return _Nx; };

	/*! Returns the number of y-bins used for the particle fluence */
	int getNy() const { return _Ny; };

	/*! Returns the total number of fluence bins */
	int getFluenceBins() const { return _Nx*_Ny; };

	/*! Given a recatngle x1...x2, z1...z2, this function sets \a mask to 1 for the fluence bins inside the rectangle, and to zero outside.
	Obviously \a mask must be of size of at least _Nx x _Ny bins.
	*/
	void setMask(double x1, double x2, double y1, double y2, char *mask);

	/*! Returns the number of particles scored in the phase space data per initial history */
	inline double getParticlesPerCase() const { return _particlesPerCase; };

	/*! Returns the statistical weight of the particles in the phase space data file */
	inline double getParticleWeight() const { return _weight; };

protected:

	string                  _dataFileName;    //!< The file name of the data file used to intialize this instance.

	double*                 _fluenceP;        //!< Primary particle fluence as a function of x,y position in a plane just above the MLC
	double*                 _prob1;           //!< Probability for a 1.17 MeV photon as a function of x,y position in a plane just above the MLC
	double*                 _fluenceS;        //!< Scatter particle fluence as a function of x,y position in a plane just above the MLC
	double*                 _pScoarse;        //!< Probability for a scattered photon to fall into the first angular table
	double*                 _pScoarse1;       //!< Probability for a scattered photon to fall into the second angular table
	vector<unsigned short*> _scatSpectrum;    //!< Scatter particle energy spectra in _Nx1*_Ny1 tiles prepared for convenient sampling
	vector<unsigned short*> _thetaPhiPflu;    //!< Primary particle fluence as a function of sin(theta),phi in _Nx1*_Ny1 tiles
	vector<unsigned short*> _thetaPhiSflu;    //!< Scatter particle fluence as a function of sin(theta),phi in _Nx1*_Ny1 tiles
	vector<unsigned short*> _thetaPhiSflu1;   //!< Scatter particle fluence as a function of sin(theta),phi in _Nx1*_Ny1 tiles
	vector<unsigned short*> _thetaPhiSfluC;   //!< Used for sampling sin(theta),phi of primary particles

	double                  _xmin;            //!< Minimum x-position of the phase-space data (Note: this is at the MLC position and not in the isocenter)
	double                  _xmax;            //!< Maximum x-position of the phase-space data (Note: this is at the MLC position and not in the isocenter)
	double                  _ymin;            //!< Minimum y-position of the phase-space data (Note: this is at the MLC position and not in the isocenter)
	double                  _ymax;            //!< Maximum y-position of the phase-space data (Note: this is at the MLC position and not in the isocenter)
	double                  _Emin;            //!< Minimum energy for the scattered particle energy spectra
	double                  _Emax;            //!< Maximum energy for the scattered particle energy spectra
	double                  _sourcePos;       //!< The virtual source position (virtual source is at (0,0,_sourcePos)
	double                  _phspPos;         //!< The phsp plane z-position 
	double                  _dE;              //!< Energy bin width
	double                  _dx;              //!< x-bin width for the _fluenceP and _fluenceS histograms
	double                  _dy;              //!< y-bin width for the _fluenceP and _fluenceS histograms
	double                  _dx1;             //!< x-width of the _Nx1 x _Ny1 tiles for spectrum and angular distribution histograms
	double                  _dy1;             //!< y-width of the _Nx1 x _Ny1 tiles for spectrum and angular distribution histograms
	double                  _delxP;           //!< When read from the phase-space file, the x-range for primary angular deflections. Then transformed to bin width
	double                  _delyP;           //!< When read from the phase-space file, the y-range for primary angular deflections. Then transformed to bin width
	double                  _delxS;           //!< When read from the phase-space file, the x-range for scatter angular deflections. Then transformed to bin width
	double                  _delyS;           //!< When read from the phase-space file, the y-range for scatter angular deflections. Then transformed to bin width
	double                  _delxS1;          //!< When read from the phase-space file, the x-range for scatter angular deflections. Then transformed to bin width
	double                  _delyS1;          //!< When read from the phase-space file, the y-range for scatter angular deflections. Then transformed to bin width
	double                  _delxSc;          //!< When read from the phase-space file, the x-range for coarse scatter angular deflections. Then transformed to bin width
	double                  _delySc;          //!< When read from the phase-space file, the y-range for coarse scatter angular deflections. Then transformed to bin width
	double                  _minxP;           //!< Minimum of angular deflection of primary photons in x
	double                  _minyP;           //!< Minimum of angular deflection of primary photons in y
	double                  _minxS;           //!< Minimum of angular deflection of scattered photons in x
	double                  _minyS;           //!< Minimum of angular deflection of scattered photons in y
	double                  _minxS1;          //!< Minimum of angular deflection of scattered photons in x
	double                  _minyS1;          //!< Minimum of angular deflection of scattered photons in y
	double                  _minxSc;          //!< Minimum of angular deflection of scattered photons in x for the coarse grid
	double                  _minySc;          //!< Minimum of angular deflection of scattered photons in y for the coarse grid
	double                  _dxti0;           //!< Used to convert x-positions into a bin value
	double                  _dxti1;           //!< Used to convert x-positions into a bin value
	double                  _dyti0;           //!< Used to convert y-positions into a bin value
	double                  _dyti1;           //!< Used to convert y-positions into a bin value

	double                  _weight;          //!< The statistical weight of the particles in the phase-space data (they all have the same weight)
	double                  _Ncase;           //!< The number of histories that were run to generate the phase space data
	double                  _Nprim1;          //!< The number of 1.17 MeV primary photons that were scored
	double                  _Nprim2;          //!< The number of 1.33 MeV primary photons that were scored
	double                  _Nscat;           //!< The number of scattered photons that were scored
	double                  _pPrim1;          //!< The probability for picking a 1.17 MeV primary photon (vs the 1.33 MeV line).
	double                  _particlesPerCase;//!< Particles per initial case = (_Nprim1 + _Nprim2 + _Nscat)/_Ncase;
	int                     _Nx;              //!< Number of x-bins for the primary and scatter fluence
	int                     _Ny;              //!< Number of y-bins for the primary and scatter fluence
	int                     _nE;              //!< Number of energy bins for the scattered spectrums
	int                     _Nx1;             //!< Number of x-bins for the fluence scored differential in energy and angle
	int                     _Ny1;             //!< Number of y-bins for the fluence scored differential in energy and angle
	int                     _Nxb;             //!< Number of x- angular bins for primary photons and part of the scattered photons
	int                     _Nyb;             //!< Number of y- angular bins for primary photons and part of the scattered photons
	int                     _Nxb1;            //!< Number of x- angular bins for part of the scattered photons
	int                     _Nyb1;            //!< Number of y- angular bins for part of the scattered photons
	int                     _Nxbc;            //!< Number of x- angular bins for part of the scattered photons
	int                     _Nybc;            //!< Number of y- angular bins for part of the scattered photons

	bool                    _isValid;         //!< True, if all intializations went OK

};

class HeadTransport;

class MLCattenuation {

public:

	/*! Constructor. Loads the data from the file \a fname */
	MLCattenuation(const char* fname);

	/*! Destructor */
	~MLCattenuation();

	/*! Returns one over the minimum attenuation for location across the leaves given  \a y. (within integer precision) */
	inline int getInverseRejection(double y) const {
		int iy = (int)(_dyi0*y + _dyi1);
		int result = iy >= 0 && iy < _Ny ? _data[iy] : iy < 0 ? _attLeft : _attRight;
		return result - 5;
	};

	/*! Is this MLCattenuation instance valid? */
	bool isValid() const { return _isValid; };

private:

	double         _ymin;      //!< The minimum initialized y-position (typically the left edge of the first leaf pair)
	double         _ymax;      //!< The maximum initialized y-position (typically the right edge of the first leaf pair)
	double         _dyi0;      //!< Used to determine the bin from a y-position
	double         _dyi1;      //!< Used to determine the bin from a y-position
	unsigned short*  _data;      //!< The actual pre-computed attenuation data
	int              _Ny;        //!< The number of bins in _data
	int              _attLeft;   //!< Value returned for y < _ymin
	int              _attRight;  //!< Value returned for y > _ymax
	bool             _isValid;   //!< Flaf indicating if the instance is valid (i.e., the data has been successfully loaded)
};

class SEGMENT
{
public:
	SEGMENT();
	SEGMENT(ConfigFile *cf);
	~SEGMENT()
	{
		if (primAT) delete[] primAT;
		if (scatAT) delete[] scatAT;
		if (primMask) delete[] primMask;
		if (scatMask) delete[] scatMask;
	}
	bool load(ConfigFile *cf);
	bool loadViewRaySegment(const string& vrs); //load segment from ViewRay's format
	bool setPos(vector<double> &values);
	int sample(PRNG &rng, Particle pars[], int is, CoPhaseSpace *ps, HeadTransport* headTransport, MLCattenuation* mlcAttenuation = NULL);

	//static constant about the MLC
	static const int NLeaf;
	static const double LeafWidth;
	static const double DX;

	vector<pair<double, double> >  leafPositions;  //!< List of left and right leaf positions in cm in the isocenter plane for all leaf pairs
	double                         onTime;         //!< Time (or relative fluence) for this segment

	unsigned short*    primAT;        //!< Alias sampling table for primary photons (used to sample x- and y-position in the plane above the MLC)
	unsigned short*    scatAT;        //!< Same as _primAT but for scattered particles
	char*              primMask;      //!< The primary particles mask 
	char*              scatMask;      //!< The scattered particles mask
	double             nSample;       //!< Number of times to sample
	double             pPrim;         //!< Probability for primary.
	int                beamID;        //!< which beam it belongs to. The ID starts from zero
	double             fluence;       //!< relative fluence during the treatment, used to sample the segment
	// Mask = 0 means within open+small margin
	// Mask = 1 means between small and large margin
	// Mask = 2 means outside large margin
};

class BEAM //describing the beam configuration
{
public:
	BEAM(){};
	BEAM(ConfigFile *cf){ load(cf); }
	bool load(ConfigFile *cf);
	bool loadViewRayBeam(const string& vrb);//load beam configuration from ViewRay's plan
	double beamOnTime()
	{
		double sum = 0;
		for (size_t i = 0; i < segments.size(); ++i)
		{
			sum += segments[i].onTime;
		}
		return sum;
	}
	vector<SEGMENT> segments;     //!< The list of segments that belong to this beam
	double  gantryAngle;          //!< The gantry angle in degree
	double  sinPhi, cosPhi;       //!< sin and cos of the gantry angle to accelerate calculation
// 	double  Xiso;                 //!< The x-position of the isocenter in cm
// 	double  Yiso;                 //!< The y-position of the isocenter in cm
// 	double  Ziso;                 //!< The z-position of the isocenter in cm
	int     headIndex;
};

class SourceHead
{
public:
	//method
	SourceHead():isoX(0),isoY(0),isoZ(0),prescriptionDose(0),treatmentFraction(1), kE(1){}
	~SourceHead();
	//SourceHead(ConfigFile *cf){ init(cf); }

	bool init(ConfigFile *cf);
	vector<BEAM> beams;
	double beamOnTime()
	{
		double sum = 0;
		for (size_t i = 0; i < beams.size(); ++i)
		{
			sum += beams[i].beamOnTime();
		}
		return sum;
	}
	void isoCenter(double* px, double* py, double* pz)
	{
		*px = isoX;
		*py = isoY;
		*pz = isoZ;
	}
	void prescription(double* pDose, int* fration)
	{
		*pDose = prescriptionDose;
		*fration = treatmentFraction;
	}
	int sample(PRNG *rng, Particle pars[]);//one call to his method means getting photons for one history
	//int sample(PRNG *rng, Particle *pars, int N); //generate given number of particles, return how many histories called

	bool getViewRayBeams(const string& file);//get beam config from ViewRay's plan file
	bool exportViewRayBeams(const string& file);
	double maxEnergy() { return kE*1.35e6; } //unit eV
	static CoPhaseSpace* getPhaseSpace(){ return phsp; }
	static HeadTransport* getHeadTransport(){ return headTransport; }

#ifdef STAT_HEAD
	void statHead();
	void statOutput();
#endif
	//data
	static int NReject[3];
	static double LeafAdjustment;
private:
	//method
	
	//data
	static HeadTransport *headTransport;
	static CoPhaseSpace *phsp;//phase space, only need one instance now
	static MLCattenuation* mlcAttenuation;
	BinSampler segSampler;
	vector<pair<int, int> > segIndex; //first = beam index; second = segment index
	double isoX, isoY, isoZ;
	double prescriptionDose;
	int treatmentFraction;

	//experiment variable
	double kE; // energy scaling factor
};

SourceHead* gSH = NULL;
PRNG* gRNG = NULL;
Particle* gArray = NULL;
int* gNP = NULL;

bool SourceHead_Init(ConfigFile* cf, int NT)
{
	gSH = new SourceHead;
	gRNG = new PRNG[NT];
	for (int i = 0; i < NT; ++i) gRNG[i].init(1234 + i);
	gArray = new Particle[NT*NSampleStack];
	gNP = new int[NT];
	for (int i = 0; i < NT; ++i) gNP[i] = 0;
	return gSH->init(cf);
}
void SourceHead_Delete()
{
	if (gSH) delete gSH;
	if (gRNG) delete[] gRNG;
	if (gArray) delete[] gArray;
	if (gNP) delete[] gNP;
}

int SourceHead_Sample(int it, double& E, double& weight, MonteCarlo::Vector& x, MonteCarlo::Vector& v, ParticleType& type)
{
	if (gNP[it] == 0)
	{
		gNP[it] = gSH->sample(gRNG + it, gArray + it*NSampleStack);

		if (gNP[it] == 0) return -1; //fail to get one history
	}
	//here at least we have one particle ready
	--gNP[it];
	Particle& p = gArray[it*NSampleStack + gNP[it]];
	E = p.E;
	weight = p.weight;
	x.x = p.x;
	x.y = p.y;
	x.z = p.z;
	v.x = p.u;
	v.y = p.v;
	v.z = p.w;
	type = p.type;

	return gNP[it]; //return how many particles left for one sampling
}
double SourceHead_Emax()
{
	if (gSH) return gSH->maxEnergy();
	else return 0;
}
double SourceHead_NormFactor()
{
	if (gSH) return 1.0889e15 * gSH->beamOnTime();
	else return 0;
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
class IGeometry {
	/*! \brief An abstract interface class for interrogating the geometry during a Monte Carlo simulation */
public:
	virtual ~IGeometry() {};

	virtual int isWhere(unsigned int aTimeIndex, const Vector &aX, Medium & material) const = 0;

	virtual int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const = 0;

	virtual int getNumRegions() = 0;
};

class Plane {
public:

	/*!
	* \brief Query x component of plane point.
	*
	* Return x component of point used to define this plane.
	*/
	double X() const { return _point[0]; }
	/*!
	* \brief Query y component of plane point.
	*
	* Return y component of point used to define this plane.
	*/
	double Y() const { return _point[1]; }
	/*!
	* \brief Query z component of plane point.
	*
	* Return z component of point used to define this plane.
	*/
	double Z() const { return _point[2]; }

	/*!
	* \brief Constructor for the plane class.
	*
	* Constructor.  Define a plane by specifying a point and a normal.
	* \param point A point on the plane.
	* \param norm Normal to the plane.
	*/
	Plane(Vector point, Vector norm) {
		_point[0] = point._x; _point[1] = point._y; _point[2] = point._z;
		_norm[0] = norm._x; _norm[1] = norm._y; _norm[2] = norm._z;
		_regionCode = 0;
	}

	/*!
	* \brief What is the region number of this plane.
	*/
	void setRegionCode(int aRegionCode) { _regionCode = aRegionCode; }

	/*!
	* \brief Default constructor.
	*/
	Plane() { }

	/**
	* \brief What is the next region the particle will encounter.
	*
	* If the particle will hit this plane, then this plane is the
	* next region along the path.
	* \param origin - position of particle.
	* \param direction - direction of particle.
	* \param distance - intended distance.  This will be reset if we intersect another region in a shorter distance.
	* \param aNewRegion - Index of next region along path of particle.
	* \return Index of next region along path of particle.
	*/
	int nextRegionAlongPath(const Vector &origin, const Vector &direction, double &distance, int& aNewRegion) const {
		double d = distance;
		if (hit(origin, direction, d)) {
			if (d < distance) {
				distance = d;
				aNewRegion = _regionCode;
				return aNewRegion;
			}
		}
		return -1;
	}



	/**
	* \brief Determine if a particle will hit this plane.
	*
	* Equations in How to code geometry by Walter R Helson.  Writing subroutine howfar. Page 9.
	*
	* \param origin - position of particle.
	* \param direction - direction of particle.
	* \param distance - intended distance.  This will be reset if we intersect another region in a shorter distance.
	* \ return true if particle will hit this plane.
	*/
	bool hit(const Vector &origin, const Vector &direction, double &distance) const {

		Vector C(_point[0], _point[1], _point[2]);
		Vector N(_norm[0], _norm[1], _norm[2]);

		Vector X(origin._x, origin._y, origin._z);
		Vector U(direction._x, direction._y, direction._z);

		double numer = (C._x - X._x) * N._x + (C._y - X._y) * N._y + (C._z - X._z) * N._z;
		double denom = U._x * N._x + U._y * N._y + U._z * N._z;

		if (fabs(denom)>1e-8) {
			if (fabs(numer)>1e-8) {
				double t = numer / denom;
				if (t < 0) return false; // travelling away from plane
				distance = t;
				return true;
			}
			else {
				// particle in the plane
				return true;
			}
		}
		else {
			// No component in direction of the plane.
			return false;
		}

	}

	/*!
	* \brief Store the point used to define this plane.
	*/
	double _point[3];

	/*!
	* \brief Store the normal to this plane.
	*/
	double _norm[3];

	/*!
	* \brief Return this code if the plane will be hit.
	*/
	int _regionCode;

};

class Cylinder {
public:

	/*! \brief Construct a cylinder with radius \a aRadius and axis passing through \a aCenter */
	Cylinder(const double aRadius, const Vector & aCenter) : _radius(aRadius) {
		_radiusSquared = _radius * _radius;
		_center._x = aCenter._x;
		_center._y = aCenter._y;
		_center._z = aCenter._z;
		_regionCode = 0;
	}

	/*! \brief Set the region code for this cylinder.
	*
	* The region code is used in the region numbering scheme of geometry classes making use of a z-cylinder
	*/
	void setRegionCode(const int aRegionCode) { _regionCode = aRegionCode; }

	/*! \brief Get the region code of this cylinder */
	int getRegionCode() const { return _regionCode; }

	/*! \brief Construct a dummy cylinder */
	Cylinder() {
		_radius = 0;
		_radiusSquared = 0;
		_center._x = 0;
		_center._y = 0;
		_center._z = 0;
		_regionCode = 0;
	}

	/*! \brief Compute the next intersection distance with this cylinder for the ray defined by \a aX and \a aU.
	*
	* If \a aIsInside is true, the position is known to be inside, else it is known to be outside.
	* Returns true if the ray intersects the cylinder before \a aDistance and false otherwise.
	* In the former case \a aDistance is truncated to the distance to the intersection.
	*/
	bool nextRegionAlongPath(const Vector & aX, const Vector & aU, double& aDistance, bool aIsInside) const;

	/*! \brief Returns true if the psoition \a aX is inside the cylinder, false otherwise. */
	bool isInside(const Vector &aX) const;

	/*! \brief Returns the radius of this cylinder. */
	double getRadius() const { return _radius; }

	/*! \brief Returns the radius squared of this cylinder. */
	double getRadiusSquared() const { return _radiusSquared; }

private:
	Vector _center;             //!< A point on the axis of this cylinder
	double  _radius;             //!< The radius of this cylinder
	double  _radiusSquared;      //!< The radius squared of this cylinder
	int       _regionCode;         //!< The region code of this cylinder
};

class CylinderPair : public IGeometry {

public:

	/*! \brief Construct the pair of cylinders.
	*
	* @param[in]  aRadii  The two radii of the cylinder pair
	* @param[in]  aCenter A point on the axis of the cylinder pair
	* @param[in]  aInnerRegion The inner region index (i.e., the region inside the smaller cylinder)
	* @param[in]  aMiddleRegion The middle region index (i.e., the region between the two cylinders)
	* @param[in]  aOuterRegion The outer region index (i.e., the region outside the larger cylinder)
	*/
	CylinderPair(const pair<double, double> aRadii, const Vector aCenter, const int aInnerRegion, const int aMiddleRegion, const int aOuterRegion) {

		_innerRegion = aInnerRegion;
		_middleRegion = aMiddleRegion;
		_outerRegion = aOuterRegion;

		_radii.first = aRadii.first;
		_radii.second = aRadii.second;

		_radiiSquared.first = _radii.first * _radii.first;
		_radiiSquared.second = _radii.second * _radii.second;

		_center._x = aCenter._x;
		_center._y = aCenter._y;
		_center._z = aCenter._z;

		_cyls.first = new Cylinder(_radii.first, _center);
		_cyls.second = new Cylinder(_radii.second, _center);

	}

	/*! \brief Destructor */
	~CylinderPair() {
		delete _cyls.first; delete _cyls.second;
	}

	/*! \brief Implements the IGeometry::getNumRegions() interface */
	int getNumRegions() { return 3; }

	/*! \brief Implements the IGeometry::isWhere() interface */
	int isWhere(unsigned int, const Vector &aX, Medium&) const {
		double distance = (aX._x - _center._x) * (aX._x - _center._x);
		distance += (aX._z - _center._z) * (aX._z - _center._z);
		if (distance <= _radiiSquared.first) return _innerRegion;
		if (distance <= _radiiSquared.second) return _middleRegion;
		return _outerRegion;
	}

	/*! \brief Implements the IGeometry::nextRegionAlongPath() interface */
	int nextRegionAlongPath(unsigned int, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& aMedium) const {

		int nextRegion = -1;

		Vector shifted(aX);
		shifted -= _center;

		if (aCurrentRegion == _innerRegion) {

			if (_cyls.first->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) { aMedium = Air; return _middleRegion; }

		}
		else if (aCurrentRegion == _middleRegion) {

			double B = shifted._x*aU._x + shifted._z*aU._z;
			if (B >= 0) { // travelling outward
				if (_cyls.second->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) { aMedium = Air; return _outerRegion; }
			}
			else { // travelling inward
				if (_cyls.first->nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) { aMedium = Air; return _innerRegion; }
				if (_cyls.second->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) { aMedium = Air; return _outerRegion; }
			}


		}
		else if (aCurrentRegion == _outerRegion) {

			double B = shifted._x*aU._x + shifted._z*aU._z;
			if (B < 0) { // travelling inward
				if (_cyls.second->nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) { aMedium = Air; return _middleRegion; }
			}

		}

		return nextRegion;

	}

private:

	pair<double, double> _radii;             //!< The two radii
	pair<double, double> _radiiSquared;      //!< The two radii squared
	Vector _center;                          //!< A point on the cylinder pair axis
	int _innerRegion;                           //!< The inner region index
	int _middleRegion;                          //!< The middle region index
	int _outerRegion;                           //!< The outer region index
	pair<Cylinder*, Cylinder*> _cyls;      //!< The actual implementation of the geometry methods.

};
bool Cylinder::nextRegionAlongPath(const Vector &aX, const Vector &aU, double& aDistance, bool aIsInside) const {

	double xs = aX._x - _center._x, zs = aX._z - _center._z;
	double ux = aU._x, uz = aU._z;

	//
	// The quadratic equation to solve is A*t^2 + 2*B*t - C, with A, B, C defined below. The discriminant is D = B^2 + C and the solutions, 
	// if they exist, are (-B +/- sqrt(D))/A.
	// There are no solutions when
	//   a) A = 0 (trajectory is parallel to cylinder axis
	//   b) D < 0 Strictly speaking, this can only happen if aIsInside is false (C < 0 in this case)
	//   c) aIsInside is false and B >= 0
	//
	double A = ux*ux + uz*uz;
	if (A < 1e-10) return false; // travelling parallel to cylinder axis
	double B = ux * xs + uz * zs;
	if (!aIsInside && B >= 0) return false;
	double C = _radiusSquared - xs*xs - zs*zs;
	double D = B*B + A*C;
	if (D < 0) return false;

	if (aIsInside) {
		// when we are inside, C >= 0 (except for numerical roundoff), so the discriminant is always > |B|. 
		// Thus, the positive solution is (-B + sqrt(D))/A. To improve precision when B > 0, we multiply and divide the solution by (sqrt(D) + B), which 
		// then simplifies to C/(sqrt(D) + B).
		// To avoid the relatively expensive sqrt and division evaluations, we check if A*aDistance^2 + 2*B*aDistance - C < 0, and if so, simply return false.
		//
		if (aDistance*(A*aDistance + 2 * B) < C) return false;
		double t = B > 0 ? C / (sqrt(D) + B) : (sqrt(D) - B) / A;
		aDistance = t; return true;
	}
	// If here, we are outside and also B < 0. In that case, the solution (first intersection) is (-B - sqrt(D))/A = -(sqrt(D) + B)/A. To improve precision, 
	// we multiply and divide with (sqrt(D) - B), and the solution simplifies to -C/(sqrt(D) - B)
	double t = -C / (sqrt(D) - B);
	if (t <= aDistance) { aDistance = t; return true; }
	return false;

	return false;

}

class BaseRegion : public IGeometry {
	/**
	* \brief  Base region for MLC and gradient
	*
	* Class describes the geometry of the base region, a large parallelpiped marking the boundaries
	* of the universe, with cylinders for the MLC and gradient.
	*/
public:

	BaseRegion() :
		_top(Vector(0, 0, 70), Vector(0, 0, -1)),
		_bot(Vector(0, 0, -5), Vector(0, 0, +1)),
		_lef(Vector(-42, 0, 0), Vector(+1, 0, 0)),
		_rig(Vector(+42, 0, 0), Vector(-1, 0, 0)),
		_bac(Vector(0, +42, 0), Vector(0, -1, 0)),
		_fro(Vector(0, -42, 0), Vector(0, +1, 0)),
		_mlc(pair<double, double>(41.0, 50.0), Vector(0, 0, Constants::SOURCE_TO_ISO_DISTANCE), 1, 2, 0),
		_gradient(pair<double, double>(35.0, 40.5), Vector(0, 0, 0), 3, 4, 0),
		_baseRegion(0),
		_mlcRegionInner(1), _mlcRegionMiddle(2), _mlcRegionOuter(0),
		_gradientRegionInner(3), _gradientRegionMiddle(4), _gradientRegionOuter(0),
		_maxRegion(5)
	{

		_bot.setRegionCode(-1);
		_top.setRegionCode(-1);
		_lef.setRegionCode(-1);
		_rig.setRegionCode(-1);
		_fro.setRegionCode(-1);
		_bac.setRegionCode(-1);

		_xplanes.push_back(_lef);
		_xplanes.push_back(_rig);
		_yplanes.push_back(_fro);
		_yplanes.push_back(_bac);
		_zplanes.push_back(_bot);
		_zplanes.push_back(_top);

	}

	~BaseRegion() {
		_xplanes.clear(); _yplanes.clear(); _zplanes.clear();
	}

	/**
	*  Maximum z coordinate in the universe.
	*/
	double getZmax() { return _top.Z(); }

	/**
	* Minimum z coordinate in the universe.
	*/
	double getZmin() { return _bot.Z(); }

	/**
	* Maximum number of regions in this class.
	*/
	int getNumRegions() { return _maxRegion; }

	/**
	* Return the base region that the particle is in.
	* \param aX position of the particle.
	* \return Which region particle is in.
	*/
	int isWhere(unsigned int, const Vector &aX, Medium&) const {

		Medium medium;

		if (aX._x > _xplanes[0].X() && aX._x < _xplanes[1].X()) {
			if (aX._y > _yplanes[0].Y() && aX._y < _yplanes[1].Y()) {
				if (aX._z > _zplanes[0].Z() && aX._z < _zplanes[1].Z()) {

					// Check if inside MLC
					int mlcReg = _mlc.isWhere(0, aX, medium);
					if (mlcReg != _mlcRegionOuter) return mlcReg;

					// Check if inside gantry.
					int gradReg = _gradient.isWhere(0, aX, medium);
					if (gradReg != _gradientRegionOuter) return gradReg;

					// Inside universe.
					return _baseRegion;

				}
			}
		}

		return -1;

	}

	/**
	* Return the next region along the path
	*
	* \param aCurrentRegion - which region is the particle in?
	* \param aX - position of particle.
	* \param aU - direction of particle.
	* \param aIntendedDistance - distance particle is to travel.  Will be reset if another region boundary is encountered in a shorter distance.
	* \param aMedium - Returned, the material of the next region.
	* \return Index of next region along path of particle.
	*/
	int nextRegionAlongPath(unsigned int, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& aMedium) const {

		int newRegion = -1;

		if (aCurrentRegion == _baseRegion) {

			for (int i = 0; i < 2; i++) _xplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			double distance = aIntendedDistance;
			int tempMlcRegion = _mlc.nextRegionAlongPath(0, _mlcRegionOuter, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempMlcRegion;
				aIntendedDistance = distance;
				if (newRegion == _mlcRegionOuter) newRegion = _baseRegion;
			}

			distance = aIntendedDistance;
			int tempGradientRegion = _gradient.nextRegionAlongPath(0, _gradientRegionOuter, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempGradientRegion;
				aIntendedDistance = distance;
				if (newRegion == _gradientRegionOuter) newRegion = _baseRegion;
			}

		}
		else if (aCurrentRegion == _mlcRegionInner || aCurrentRegion == _mlcRegionMiddle) {

			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			double distance = aIntendedDistance;
			int tempMlcRegion = _mlc.nextRegionAlongPath(0, aCurrentRegion, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempMlcRegion;
				aIntendedDistance = distance;
			}

		}
		else if (aCurrentRegion == _gradientRegionInner || aCurrentRegion == _gradientRegionMiddle) {

			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			double distance = aIntendedDistance;
			int tempGradientRegion = _gradient.nextRegionAlongPath(0, aCurrentRegion, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempGradientRegion;
				aIntendedDistance = distance;
			}

		}

		return newRegion;

	}

private:
	/**
	* The plane bounding the top of the universe.
	*/
	Plane _top;
	/**
	* The plane bounding the bottom of the universe.
	*/
	Plane _bot;
	/**
	* The plane bounding the left side of the universe.
	*/
	Plane _lef;
	/**
	* The plane bounding the rightt side of the universe.
	*/
	Plane _rig;
	/**
	* The plane bounding the back of the universe.
	*/
	Plane _bac;
	/**
	* The plane bounding the from of the universe.
	*/
	Plane _fro;

	/**
	* The left and right planes are stored in this vector.
	*/
	vector<Plane> _xplanes;
	/**
	* The top and bottom planes are stored in this vector.
	*/
	vector<Plane> _yplanes;
	/**
	* The front and back planes are stored in this vector.
	*/
	vector<Plane> _zplanes;

	/**
	* Two cylinders marking inner and outer MLC boundaries.
	*/
	CylinderPair _mlc;

	/**
	* Two cylinders marking inner and outer gradientn boudnaries.
	*/
	CylinderPair _gradient;

	/**
	* Region number of the base region.  This region is neither
	* inside an MLC cylinder or a Gradient cylinder.
	*/
	int _baseRegion;

	/**
	* Region number of the region inside the innermost MLC cylinder.
	*/
	int _mlcRegionInner;
	/**
	* Region number of the region inside the MLC cylinder pair.
	*/
	int _mlcRegionMiddle;
	/**
	* Region number outside the outermost MLC cylinder.
	*/
	int _mlcRegionOuter;

	/**
	* Region number of the region inside the innermost gradient cylinder.
	*/
	int _gradientRegionInner;
	/**
	* Region number of the region inside the gradient cylinder pair.
	*/
	int _gradientRegionMiddle;
	/**
	* Region number outside the outermost gradient cylinder.
	*/
	int _gradientRegionOuter;

	/**
	* Number of regions represented in this class.
	*/
	int _maxRegion;

};

class BaseRegionMlc : public IGeometry {
	/**
	* \brief  Base region for MLC only.
	*
	* Class describes the geometry of the base region, a large parallelpiped marking the boundaries
	* of the universe, with cylinders for the MLC.
	*/
public:

	/**
	* \brief  Constructor.
	*
	* Here is defined the extents of the universe, and the cylinders marking
	* the boundaries of the MLC.
	*
	*/
	BaseRegionMlc() :
		_top(Vector(0, 0, 70), Vector(0, 0, -1)),
		_bot(Vector(0, 0, -5), Vector(0, 0, +1)),
		_lef(Vector(-42, 0, 0), Vector(+1, 0, 0)),
		_rig(Vector(+42, 0, 0), Vector(-1, 0, 0)),
		_bac(Vector(0, +42, 0), Vector(0, -1, 0)),
		_fro(Vector(0, -42, 0), Vector(0, +1, 0)),
		_mlc(pair<double, double>(41.0, 50.0), Vector(0, 0, Constants::SOURCE_TO_ISO_DISTANCE), 1, 2, 0),
		_baseRegion(0),
		_mlcRegionInner(1), _mlcRegionMiddle(2), _mlcRegionOuter(0),
		_maxRegion(3)
	{

		_bot.setRegionCode(-1);
		_top.setRegionCode(-1);
		_lef.setRegionCode(-1);
		_rig.setRegionCode(-1);
		_fro.setRegionCode(-1);
		_bac.setRegionCode(-1);

		_xplanes.push_back(_lef);
		_xplanes.push_back(_rig);
		_yplanes.push_back(_fro);
		_yplanes.push_back(_bac);
		_zplanes.push_back(_bot);
		_zplanes.push_back(_top);

	}

	/**
	* \brief Destructor.
	*
	* Empty class vectors.
	*/
	~BaseRegionMlc() {
		_xplanes.clear(); _yplanes.clear(); _zplanes.clear();
	}
	/**
	*  Maximum z coordinate in the universe.
	*/
	double getZmax() { return _top.Z(); }

	/**
	* Minimum z coordinate in the universe.
	*/
	double getZmin() { return _bot.Z(); }

	/**
	* Maximum number of regions in this class.
	*/
	int getNumRegions() { return _maxRegion; }

	/**
	* Return the base region that the particle is in.
	* \param aX position of the particle.
	* \return Which region particle is in.
	*/
	int isWhere(unsigned int, const Vector &aX, Medium&) const {

		Medium medium;

		if (aX._x > _xplanes[0].X() && aX._x < _xplanes[1].X()) {
			if (aX._y > _yplanes[0].Y() && aX._y < _yplanes[1].Y()) {
				if (aX._z > _zplanes[0].Z() && aX._z < _zplanes[1].Z()) {

					// Check if inside MLC
					int mlcReg = _mlc.isWhere(0, aX, medium);
					if (mlcReg != _mlcRegionOuter) return mlcReg;

					// Inside universe.
					return _baseRegion;

				}
			}
		}

		return -1;

	}

	/**
	* Return the next region along the path
	*
	* \param aCurrentRegion - which region is the particle in?
	* \param aX - position of particle.
	* \param aU - direction of particle.
	* \param aIntendedDistance - distance particle is to travel.  Will be reset if another region boundary is encountered in a shorter distance.
	* \param aMedium - Returned, the material of the next region.
	* \return Index of next region along path of particle.
	*/
	int nextRegionAlongPath(unsigned int, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& aMedium) const {

		int newRegion = -1;

		if (aCurrentRegion == _baseRegion) {

			for (int i = 0; i < 2; i++) _xplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			int tempMlcRegion = _mlc.nextRegionAlongPath(0, _mlcRegionOuter, aX, aU, aIntendedDistance, aMedium);
			if (tempMlcRegion != -1) newRegion = tempMlcRegion;
			if (newRegion == _mlcRegionOuter) newRegion = _baseRegion;

		}
		else if (aCurrentRegion == _mlcRegionInner || aCurrentRegion == _mlcRegionMiddle) {

			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			double distance = aIntendedDistance;
			int tempMlcRegion = _mlc.nextRegionAlongPath(0, aCurrentRegion, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempMlcRegion;
				aIntendedDistance = distance;
			}

		}

		return newRegion;

	}

private:


	/**
	* The plane bounding the top of the universe.
	*/
	Plane _top;
	/**
	* The plane bounding the bottom of the universe.
	*/
	Plane _bot;
	/**
	* The plane bounding the left side of the universe.
	*/
	Plane _lef;
	/**
	* The plane bounding the rightt side of the universe.
	*/
	Plane _rig;
	/**
	* The plane bounding the back of the universe.
	*/
	Plane _bac;
	/**
	* The plane bounding the from of the universe.
	*/
	Plane _fro;

	/**
	* The left and right planes are stored in this vector.
	*/
	vector<Plane> _xplanes;
	/**
	* The top and bottom planes are stored in this vector.
	*/
	vector<Plane> _yplanes;
	/**
	* The front and back planes are stored in this vector.
	*/
	vector<Plane> _zplanes;

	/**
	* Two cylinders marking inner and outer MLC boundaries.
	*/
	CylinderPair _mlc;

	/**
	* Region number of the base region.  This region is neither
	* inside an MLC cylinder or a Gradient cylinder.
	*/
	int _baseRegion;

	/**
	* Region number of the region inside the innermost MLC cylinder.
	*/
	int _mlcRegionInner;
	/**
	* Region number of the region inside the MLC cylinder pair.
	*/
	int _mlcRegionMiddle;
	/**
	* Region number outside the outermost MLC cylinder.
	*/
	int _mlcRegionOuter;

	/**
	* Number of regions represented in this class.
	*/
	int _maxRegion;


};

class Cylinders : public IGeometry {
	/*! \brief Provides geometry methods for a set of concentric cylinders with axis along the y-direction.
	*
	* The cylinders are infinitely long (they are cut by other geometry components).
	* The set of cylinders is used to model the gradient coils in the ViewRay delivery system.
	*/
public:

	Cylinders(const vector<double> aRadii, const Vector &aCenter, const vector<int> &aRegionCodes, int aOuterRegion) : _outerRegion(aOuterRegion) {

		_ncyl = (int)aRadii.size();

		for (int i = 0; i < _ncyl; i++) {
			_radii.push_back(aRadii[i]);
			_radiiSquared.push_back(aRadii[i] * aRadii[i]);
			_cyls.push_back(new Cylinder(aRadii[i], aCenter));
			_regionCodes.push_back(aRegionCodes[i]);
		}

		_media.push_back(Air);
		_media.push_back(TxFR4);
		_media.push_back(Air);
		_media.push_back(TxPCB);
		_media.push_back(Air);
		_media.push_back(RFShield);
		_media.push_back(Air);
		_media.push_back(GradientCoil);

		_center._x = aCenter._x;
		_center._y = aCenter._y;
		_center._z = aCenter._z;

	}

	~Cylinders() {
		for (int i = 0; i < _ncyl; i++) delete _cyls[i];
		_radii.clear(); _radiiSquared.clear(); _regionCodes.clear();
	}

	int getNumRegions() { return _ncyl + 1; }

	// \return Which region particle is in.
	int isWhere(unsigned int, const Vector &aX, Medium& medium) const {
		double distance = (aX._x - _center._x) * (aX._x - _center._x);
		distance += (aX._z - _center._z) * (aX._z - _center._z);

		if (distance >= _radiiSquared[_ncyl - 1]) { medium = Air; return _outerRegion; }

		for (int i = 0; i < _ncyl; i++) {
			if (distance < _radiiSquared[i]) { medium = _media[i]; return _regionCodes[i]; }
		}

		return _outerRegion;

	}

	// \return Index of next region along path of particle.
	int nextRegionAlongPath(unsigned int, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& medium) const {

		int nextRegion = -1;

		Vector shifted(aX);
		shifted -= _center;

		// 
		if (aCurrentRegion == -1) aCurrentRegion = isWhere(0, aX, medium);

		if (aCurrentRegion == _regionCodes[0]) {

			if (_cyls[0]->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
				medium = _media[1];
				return _regionCodes[1];
			}
			else
			{
				return _regionCodes[0];
			}

		}
		else if (aCurrentRegion == _outerRegion) {

			double B = shifted._x*aU._x + shifted._z*aU._z;
			if (B < 0) { // travelling inward
				if (_cyls[_ncyl - 1]->nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) {
					medium = _media[_ncyl - 1];
					return _regionCodes[_ncyl - 1];
				}
			}

		}
		else if (aCurrentRegion < _regionCodes[_ncyl - 1]) {

			double B = shifted._x*aU._x + shifted._z*aU._z;

			if (B >= 0) { // travelling outward

				if (_cyls[aCurrentRegion]->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					medium = _media[aCurrentRegion + 1];
					return _regionCodes[aCurrentRegion + 1];
				}
				else
				{
					return _regionCodes[aCurrentRegion];
				}

			}
			else { // travelling inward

				if (_cyls[aCurrentRegion - 1]->nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) {
					medium = _media[aCurrentRegion - 1];
					return _regionCodes[aCurrentRegion - 1];
				}
				if (_cyls[aCurrentRegion]->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					medium = _media[aCurrentRegion];
					return _regionCodes[aCurrentRegion];
				}

				return _regionCodes[aCurrentRegion];

			}

		}
		else if (aCurrentRegion == _regionCodes[_ncyl - 1]) {

			double B = shifted._x*aU._x + shifted._z*aU._z;

			if (B >= 0) { // travelling outward

				if (_cyls[aCurrentRegion]->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					return -1;
				}
				else
				{
					return _regionCodes[_ncyl - 1];
				}

			}
			else { // travelling inward

				if (_cyls[aCurrentRegion - 1]->nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) {
					medium = _media[aCurrentRegion - 1];
					return _regionCodes[aCurrentRegion - 1];
				}
				if (_cyls[aCurrentRegion]->nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					medium = _media[aCurrentRegion];
					return _regionCodes[aCurrentRegion];
				}

				return _regionCodes[aCurrentRegion];

			}

		}

		return nextRegion;

	}

private:

	vector<double> _radii;
	vector<double> _radiiSquared;
	vector<Medium> _media;
	Vector _center;
	vector<int> _regionCodes;
	vector<Cylinder*> _cyls;
	int _ncyl;
	int _outerRegion;

};

class ViewRayGantryCoils : public IGeometry {
	/*! \brief Implements the IGeometry interface for the gradient coils of the ViewRay delivery system.
	*
	* This class contains the coils around the patient.
	* Gradient, RF, TxPCB and TxFr.
	*/
public:

	/*! \brief Constructor.
	*
	* Construct one model for the inner and outer cylinders of the four coils
	*/
	ViewRayGantryCoils();

	/*! \brief Destructor */
	~ViewRayGantryCoils();

	/*! \brief Determine region that particle is in.
	*
	* Implements the IGeometry::isWhere() method.
	* \param[in] aTimeIndex Unused
	* \param[in] aX Location of particle.
	* \param[out] material On exit contains the material the particle is in.
	* \return Which region is particle in
	*/
	int isWhere(unsigned int aTimeIndex, const Vector &aX, Medium & material) const;

	/*! \brief Implements the IGeometry::nextRegionAlongPath() interface. */
	int nextRegionAlongPath(unsigned int, int aCurrentRegion, const Vector &aX, const Vector &aU,
		double &aIntendedDistance, Medium & material) const;

	/*! \brief Implements the IGeometry::getNumRegions() interface */
	int getNumRegions() { return int(_ncyl); }

private:

	vector<double> _radii;   //!< The radii of the coild

	vector<int> _regionCodes;  //!< The region codes

	Vector _isocenter;      //!< The isocenter location (which is the center of the coil cylinders)

	size_t _ncyl;              //!< Number of concentric cylinders in the model 

	Cylinders* _cyls;       //!< Actual implementation of the geometry methods

};
int ViewRayGantryCoils::isWhere(unsigned int, const Vector &aX, Medium & material) const {

	return _cyls->isWhere(0, aX, material);

}
int ViewRayGantryCoils::nextRegionAlongPath(unsigned int, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& aNextMedium) const {

	return _cyls->nextRegionAlongPath(0, aCurrentRegion, aX, aU, aIntendedDistance, aNextMedium);

}
ViewRayGantryCoils::ViewRayGantryCoils() : _isocenter(0, 0, 0) {

	_radii.push_back(35.000000); // TxCoilFr4InnerRadius
	_radii.push_back(35.200000); // TxCoilFr4OuterRadius
	_radii.push_back(35.501854); // TxCoilPCBInnerRadius
	_radii.push_back(35.550000); // TxCoilPCBOuterRadius
	_radii.push_back(39.600000); // RfShieldInnerRadius
	_radii.push_back(39.634540); // RfShieldOuterRadius
	_radii.push_back(40.000000); // GradientCoilInnerRadius
	_radii.push_back(40.500000); // GradientCoilOuterRadius

	_ncyl = _radii.size();

	int regionCode = 0;
	for (size_t i = 0; i < _ncyl; i++) { _regionCodes.push_back(regionCode); regionCode++; }

	_cyls = new Cylinders(_radii, _isocenter, _regionCodes, int(_ncyl + 1));

}
ViewRayGantryCoils::~ViewRayGantryCoils() { delete _cyls; }

class ViewRayLeaf {
public:
	/**
	* \brief Construct a ViewRay leaf class.
	*
	*/
	ViewRayLeaf(int aRegionOffset, double aZ, Vector & aLeftFocus, Vector & aRightFocus, Vector & aGapFocus, double aZfocus);
	/**
	* \brief Add MLC segments to the leaf class
	*
	*/
	void setOpenings(const vector<pair<double, double> > &leafPositions);
	/**
	* \brief Which region is particle in?
	*
	* \param[in] aTimeIndex - Which MLC segment to use.  Determines leaf position.
	* \param[in] aX - Position of particle.
	* \param[out] material - Which material is particle in.
	* \return Region the particle is in.
	*/
	int isWhere(unsigned int aTimeIndex, const Vector &aX, Medium & material) const;
	/**
	* \brief Which region is the particle in, but dont consider the air gap.
	*
	* Method is called when particle has left the air gap an intersected the leaf pair.
	*
	* \param[in] aTimeIndex - Which MLC segment to use.  Determines leaf position.
	* \param[in] aX - Position of particle.
	* \param[out] material - Which material is particle in.
	* \return Region the particle is in.
	*/
	int isWhereNoAirGap(unsigned int aTimeIndex, const Vector &aX, Medium & material) const;
	/**
	* \brief Determine next region particle will intersect.
	*
	* \param[in] aTimeIndex Which MLC segment to use.
	* \param[in] aCurrentRegion Which region is particle in.
	* \param[in] aX Particle position.
	* \param[in] aU Particle direction.
	* \param[in,out] aIntendedDistance Max distance to search for next region.  This value is reset to distance to next region
	* if it is less than aIntendedDistance.
	* \param[out] material
	*/
	int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const;
	/**
	* \brief Return the  minimum z coord of the leaf.
	*/
	double getMinZ() const { return _minZ; }
	/**
	* \brief Return the medium this region consists of.
	*/
	Medium getMedium(int aRegion) const;
private:
	/**
	* \brief Minimum z location.
	*/
	double _minZ;
	/**
	* \brief Position of the z focus of the leaf.
	*/
	double _zFocus;
	/**
	* \brief Normal vector to the air gap plane.
	*/
	Vector _airGapNormal;
	/**
	* \brief Position of the air gap plane.
	*/
	double _airGapPosition;
	/**
	* The left position of the leaf.  One element for each MLC segment.
	*/
	vector< double > _leftPos;
	/**
	* The right position of the leaf.  One element for each MLC segment.
	*/
	vector< double > _rightPos;
	/**
	* The minimum x coordinate of the leaf.  One element for each MLC segment.
	*/
	vector< double > _minX;
	/**
	* The maximum x coordinate of the leaf.  One element for each MLC segment.
	*/
	vector< double > _maxX;
	/**
	* The normals of the four planes that represent the leaf.
	*/
	vector< vector< Vector > > _leafNormals;
	/**
	* The positions of the four planes that represent the leaf.
	*/
	vector< vector< double > > _leafPositions;
	/**
	* The left focus of the leaf.
	*/
	Vector _leftFocus;
	/**
	* The right focus of the leaf.
	*/
	Vector _rightFocus;
	/**
	* The first region index of this leaf.
	*/
	int _regionOffset;
	/**
	* For the leaf class, there are three regions.
	* Region 0 is the left leaf, 1 is the air in between the left and right leaf, and 2 is in the right leaf.
	*/
	int _numberOfRegions;
};
ViewRayLeaf::ViewRayLeaf(int aRegionOffset, double aZ, Vector & aLeftFocus, Vector & aRightFocus, Vector & aGapFocus, double aZfocus) : _minZ(aZ), _zFocus(aZfocus), _regionOffset(aRegionOffset) {
	using namespace Constants;
	double z = SOURCE_TO_ISO_DISTANCE;
	double y = _minZ + LEAF_THICKNESS;
	double aux = 1.0 / sqrt(y*y + z*z);
	_airGapNormal = Vector(0, z*aux, y*aux);
	_airGapPosition = aGapFocus * _airGapNormal;
	_leftFocus = aLeftFocus;
	_rightFocus = aRightFocus;
	_numberOfRegions = 3;
}
int ViewRayLeaf::isWhere(unsigned int aTimeIndex, const Vector &aLocation, Medium & aMedium) const {
	// Air gap?
	double normDotLoc = _airGapNormal._x*aLocation._x + _airGapNormal._y*aLocation._y + _airGapNormal._z*aLocation._z;
	if (normDotLoc > _airGapPosition) {
		aMedium = Air;
		return 3;
	}

	double x1 = _leafNormals[aTimeIndex][1] * aLocation;
	double x2 = _leafNormals[aTimeIndex][2] * aLocation;

	if (x1 > _leafPositions[aTimeIndex][1] && x2 < _leafPositions[aTimeIndex][2]) {
		aMedium = Air;
		return (1);
	}
	else if (x1 <= _leafPositions[aTimeIndex][1]) {
		aMedium = Tungsten;
		return 0;
	}
	else {
		aMedium = Tungsten;
		return (2);
	}

}
Medium ViewRayLeaf::getMedium(int iregion) const {
	if (iregion == 0) return Tungsten;
	if (iregion == 1) return Air;
	if (iregion == 2) return Tungsten;
	return Air;
}
// \return Index of next region along path of particle.
int ViewRayLeaf::nextRegionAlongPath(unsigned int timeIndex, int ireg, const Vector &x, const Vector &u, double &t, Medium& newmed)  const {
	if (ireg >= 0) {
		double up = _leafNormals[timeIndex][ireg] * u;
		int inew = ireg;
		if (up < -Constants::FLTEPSILON) {
			double xp = _leafNormals[timeIndex][ireg] * x;
			if (xp < _leafPositions[timeIndex][ireg] - up*t) {
				double tt = (_leafPositions[timeIndex][ireg] - xp) / up; inew = ireg - 1;
				if (tt > t) ZeusLog("Error: Huh1: %g %g", t, tt);
				t = tt;
				if (inew >= 0 && newmed) newmed = getMedium(inew);
			}
		}
		up = _leafNormals[timeIndex][ireg + 1] * u;
		if (up > Constants::FLTEPSILON) {
			double xp = _leafNormals[timeIndex][ireg + 1] * x;
			if (xp > _leafPositions[timeIndex][ireg + 1] - up*t) {
				double tt = (_leafPositions[timeIndex][ireg + 1] - xp) / up; inew = ireg + 1;
				if (tt > t) ZeusLog("Error: Huh2: %g %g", t, tt);
				t = tt;
				if (inew >= _numberOfRegions) inew = -1;
				if (inew >= 0 && newmed) newmed = getMedium(inew);
			}
		}
		return inew;
	}
	double dz = _zFocus - x._z;
	if (dz <= 0) {
		ZeusLog("Error: XplanesVRMLC::nextRegionAlongPath(): does not worh for points above focal point!\n");
	}
	double xp, up; int inew = ireg;
	xp = _leafNormals[timeIndex][0] * x;
	if (xp < _leafPositions[timeIndex][0]) {
		up = _leafNormals[timeIndex][0] * u;
		if (up > Constants::FLTEPSILON) {
			double tt = (_leafPositions[timeIndex][0] - xp) / up;
			if (tt < t) {
				t = tt; inew = 0;
				if (newmed) newmed = getMedium(0);
			}
		}
		return inew;
	}
	xp = _leafNormals[timeIndex][_numberOfRegions] * x;
	if (xp > _leafPositions[timeIndex][_numberOfRegions]) {
		up = _leafNormals[timeIndex][_numberOfRegions] * u;
		if (up < -Constants::FLTEPSILON) {
			double tt = (_leafPositions[timeIndex][_numberOfRegions] - xp) / up;
			if (tt < t) {
				t = tt; inew = _numberOfRegions - 1;
				if (newmed) newmed = getMedium(_numberOfRegions - 1);
			}
		}
		return inew;
	}
	return inew;
}
int ViewRayLeaf::isWhereNoAirGap(unsigned int aTimeIndex, const Vector &aLocation, Medium & aMedium) const {

	double x1 = _leafNormals[aTimeIndex][1] * aLocation;
	double x2 = _leafNormals[aTimeIndex][2] * aLocation;

	if (x1 > _leafPositions[aTimeIndex][1] && x2 < _leafPositions[aTimeIndex][2]) {
		aMedium = Air;
		return (1);
	}
	else if (x1 <= _leafPositions[aTimeIndex][1]) {
		aMedium = Tungsten;
		return 0;
	}
	else {
		aMedium = Tungsten;
		return (2);
	}

}
void ViewRayLeaf::setOpenings(const vector<pair<double, double> > &leafPositions) {
	std::vector< pair< double, double > >::size_type  nsegment = leafPositions.size();
	if (_leftPos.size() != nsegment) {
		_leftPos.resize(nsegment);
		_rightPos.resize(nsegment);
		_minX.resize(nsegment);
		_maxX.resize(nsegment);
		_leafNormals.resize(nsegment);
		_leafPositions.resize(nsegment);
	}
	for (unsigned int j = 0; j < nsegment; ++j) {
		_leftPos[j] = leafPositions[j].first;
		_rightPos[j] = leafPositions[j].second;
		_minX[j] = _leftPos[j] - Constants::LEAF_LENGTH;
		_maxX[j] = _rightPos[j] + Constants::LEAF_LENGTH;
		if (_leafNormals[j].size() < 4) _leafNormals[j].resize(4);
		if (_leafPositions[j].size() < 4) _leafPositions[j].resize(4);
		Vector a;
		double zFocus = Constants::SOURCE_TO_ISO_DISTANCE;
		a = Vector(zFocus, 0, _minX[j]);     a.normalizeToUnitLength(); _leafNormals[j][0] = a; _leafPositions[j][0] = a*_leftFocus;
		a = Vector(zFocus, 0, _leftPos[j]);  a.normalizeToUnitLength(); _leafNormals[j][1] = a; _leafPositions[j][1] = a*_leftFocus;
		a = Vector(zFocus, 0, _rightPos[j]); a.normalizeToUnitLength(); _leafNormals[j][2] = a; _leafPositions[j][2] = a*_rightFocus;
		a = Vector(zFocus, 0, _maxX[j]);     a.normalizeToUnitLength(); _leafNormals[j][3] = a; _leafPositions[j][3] = a*_rightFocus;
	}
}

class ViewRayGaps {
public:
	/**
	* \brief Constructor of the gaps class.
	*
	* \param aGapFocus The focal point of the gap.
	* \param azFocusY The y focus of the gap.
	* \param aZshift The shift of the gap.
	*
	*/
	ViewRayGaps(Vector & aGapFocus, double azFocusY, double aZshift);
	/**
	* \brief Which region is the particle in.
	*
	* \param aTimeIndex Which MLC segmet is being delivered.
	* \param aX Position of particle.
	* \return Which region particle is in.
	*/
	int isWhere(unsigned int aTimeIndex, const Vector &aX) const;
	/**
	* \brief
	*
	* \return Index of next region along path of particle.
	*/
	int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const;
	/**
	* \brief Destructor.
	*
	*/
	~ViewRayGaps();
private:
	/**
	* \brief Number of regions modelled in this class.
	*/
	int _numberOfRegions;
	/**
	* \brief Normal of the plane representing the air gap.
	*/
	vector<Vector> _airGapNormal;
	/**
	* \brief Position of the plane representing the air gap.
	*/
	vector<double> _airGapPosition;
	/**
	* \brief What medium does the current region consist of.
	*/
	Medium getMedium(int aRegion) const;
	/**
	* \brief z-Position of the focal spot.
	*/
	double _zFocus;
	/**
	* \brief True, if there is a gap between the leaves
	*/
	bool _hasGap;
	/**
	* \brief Min/max z-position
	*/
	double _minZ, _maxZ;
	/**
	* \brief The position of the focal point
	*/
	Vector _focus;
	/**
	* \brief _zFocus - _zTop
	*/
	double _focusToTop;
	/**
	* \brief Leaf thickness
	*/
	double _deltaZ;

};
ViewRayGaps::ViewRayGaps(Vector & aGapFocus, double azFocusY, double aZshift) : _zFocus(azFocusY), _hasGap(true), _focus(Vector(0, aZshift, 0)) {
	using namespace Constants;
	_numberOfRegions = 2 * NLEAF_PAIRS_TOTAL;
	_deltaZ = LEAF_THICKNESS;

	_focusToTop = 105;
	Vector focus(0, aZshift, 105);

	double yloc = 0, zloc = 0, aux = 0;
	int n = NLEAF_PAIRS_TOTAL;
	double farLeft = -(n / 2)*LEAF_THICKNESS;
	double farRight = (n / 2)*LEAF_THICKNESS;
	_minZ = farLeft; _maxZ = farRight;

	_airGapNormal.resize(_numberOfRegions + 1);
	_airGapPosition.resize(_numberOfRegions + 1);

	int ii = 0;
	for (int y = 0; y <= NLEAF_PAIRS_TOTAL; y++) {
		zloc = SOURCE_TO_ISO_DISTANCE;
		yloc = farLeft + double(y)*LEAF_THICKNESS;
		aux = 1.0 / sqrt(yloc*yloc + zloc*zloc);
		_airGapNormal[ii] = Vector(0, zloc*aux, yloc*aux);
		_airGapPosition[ii] = (focus * _airGapNormal[ii]);
		ii++;
		if (y != NLEAF_PAIRS_TOTAL) {
			yloc += LEAF_THICKNESS;
			aux = 1.0 / sqrt(yloc*yloc + zloc*zloc);
			_airGapNormal[ii] = Vector(0, zloc*aux, yloc*aux);
			_airGapPosition[ii] = aGapFocus * _airGapNormal[ii];
			ii++;
		}
	}
}
ViewRayGaps::~ViewRayGaps() {
	_airGapNormal.clear();
	_airGapPosition.clear();
}
int ViewRayGaps::isWhere(unsigned int, const Vector &x) const {
	double dz = Constants::SOURCE_TO_ISO_DISTANCE - x._z; double y = x._y - _focus._y;
	int region = (dz <= 0 || _focusToTop*y < dz*_minZ || _focusToTop*y >= dz*_maxZ) ? -1 :
		(int)((_focusToTop*y - dz*_minZ) / (dz*_deltaZ));
	if (_hasGap && region >= 0) {
		region *= 2;
		if (region >= 76) { return -1; }
		if (_airGapNormal[region + 1] * x > _airGapPosition[region + 1]) ++region;
	}
	return region;
};
Medium ViewRayGaps::getMedium(int aRegion) const { int rmod4 = aRegion % 4; return rmod4 == 1 || rmod4 == 3 ? Air : Tungsten; }
int ViewRayGaps::nextRegionAlongPath(unsigned int, int ireg, const Vector &x, const Vector &u, double &t, Medium& newmed) const {

	if (ireg >= 0) {
		double up = _airGapNormal[ireg] * u;
		int inew = ireg;
		if (up < -Constants::FLTEPSILON) {
			double xp = _airGapNormal[ireg] * x;
			if (xp < _airGapPosition[ireg] - up*t) {
				double tt = (_airGapPosition[ireg] - xp) / up; inew = ireg - 1;
				if (tt > t)  ZeusLog("Error: Huh1: t=%g tt=%g", t, tt);
				t = tt;
				if (inew >= 0 && newmed) newmed = getMedium(inew);
			}
		}
		up = _airGapNormal[ireg + 1] * u;
		if (up > Constants::FLTEPSILON) {
			double xp = _airGapNormal[ireg + 1] * x;
			if (xp > _airGapPosition[ireg + 1] - up*t) {
				double tt = (_airGapPosition[ireg + 1] - xp) / up; inew = ireg + 1;
				if (tt > t)  ZeusLog("Error: Huh2: t=%g tt=%g", t, tt);
				t = tt;
				if (inew >= _numberOfRegions) inew = -1;
				if (inew >= 0 && newmed) newmed = getMedium(inew);
			}
		}
		return inew;
	}
	double dz = _zFocus - x._z;
	if (dz <= 0) {
		ZeusLog("Error: YplanesVRMLC::nextRegionAlongPath(): does not work for points above focal point! x=(%g,%g,%g)\n", x._x, x._y, x._z);
	}
	double xp, up; int inew = ireg;
	xp = _airGapNormal[0] * x;
	if (xp < _airGapPosition[0]) {
		up = _airGapNormal[0] * u;
		if (up > Constants::FLTEPSILON) {
			double tt = (_airGapPosition[0] - xp) / up;
			if (tt < t) {
				t = tt; inew = 0;
				if (newmed) newmed = getMedium(0);
			}
		}
		return inew;
	}
	xp = _airGapNormal[_numberOfRegions] * x;
	if (xp > _airGapPosition[_numberOfRegions]) {
		up = _airGapNormal[_numberOfRegions] * u;
		if (up < -Constants::FLTEPSILON) {
			double tt = (_airGapPosition[_numberOfRegions] - xp) / up;
			if (tt < t) {
				t = tt; inew = _numberOfRegions - 1;
				if (newmed) newmed = getMedium(_numberOfRegions - 1);
			}
		}
		return inew;
	}
	return inew;
}


class ViewRayMlcCylinders : public IGeometry {
	/**
	* \brief A class to model cylinders surrounding the leaves.
	*
	* There are two cylinders, each which has a full set of leaf pairs.
	*
	*/
public:
	/**
	* \brief Constructor.
	*
	* Sets the radii of the cyliinders modelling the MLC.
	* There are three cylinders, separating the upper and lower
	* leaf halves.
	*
	*/
	ViewRayMlcCylinders() : _center(0, 0, 105), _numberOfRegions(3) {
		_radiiSquared.push_back(41 * 41);
		_radiiSquared.push_back(45.5*45.5);
		_radiiSquared.push_back(50 * 50);
		_radii.push_back(41);
		_radii.push_back(45.5);
		_radii.push_back(50);
	}
	/**
	* \brief Destructor.
	*
	*/
	~ViewRayMlcCylinders() { _radiiSquared.clear(); _radii.clear(); }
	/**
	* \brief Which region is the particle in.
	*
	* \param[in] aTimeIndex Which MLC segment is currently used.
	* \param[in] aX Position of the particle.
	* \param[out] material Material the region is in.
	* \return The region the particle is in.
	*
	*/
	int isWhere(unsigned int aTimeIndex, const Vector &aX, Medium & material) const;
	/**
	* \brief What is the next region along the particle path.
	*
	* \param[in] aTimeIndex Which MLC segment is currently used.
	* \param[in] aCurrentRegion Which region is particle in.
	* \param[in] aX Position of particle.
	* \param[in] aU Direction of particle.
	* \param[in,out] aIntendedDistance Max distance to search for next region.  This value is reset to distance to next region
	* if it is less than aIntendedDistance.
	* \param[out] material Material of next region along path.
	*/
	int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const;
	/**
	* \brief Query medium this region consists of.
	*
	*/
	Medium getMedium(int aRegion) const;
	/**
	* \brief Query number of regions modelled in this class.
	*
	*/
	int getNumRegions() { return 3; }
private:
	/**
	* \brief The center of all the cylinders.
	*
	*/
	Vector _center;
	/**
	* \brief The square of the cylinder radii.
	*
	*/
	vector<double> _radiiSquared;
	/**
	* \brief The radii of the cylinders.
	*
	*/
	vector<double> _radii;
	/**
	* \brief Number of regions modelled in this class.
	*
	*/
	int _numberOfRegions;
};
int ViewRayMlcCylinders::isWhere(unsigned int, const Vector &aLocation, Medium & aMedium) const {

	double z = aLocation._z - _center._z;
	double x = aLocation._x - _center._x;
	double d = x*x + z*z;

	aMedium = Air;
	if (d < _radiiSquared[0]) return 0;

	aMedium = Tungsten;
	if (d < _radiiSquared[1]) return 1;
	if (d < _radiiSquared[2]) return 2;

	aMedium = UnknownMaterial;
	return -1;

}
// \return Index of next region along path of particle.
int ViewRayMlcCylinders::nextRegionAlongPath(unsigned int, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& nextMedium) const {

	Vector xp(aX - _center);
	double A, B, C;
	A = aU._x*aU._x + aU._z*aU._z;
	B = xp._x*aU._x + xp._z*aU._z;
	C = xp._x*xp._x + xp._z*xp._z;
	if (A < 1e-10) return aCurrentRegion;  // travelling parallel to axis.

	int nextRegion = aCurrentRegion;
	if (aCurrentRegion >= 0) {
		if (B >= 0 || !aCurrentRegion) {  // in the innermost cylinder or travelling outwards
			double C1 = _radiiSquared[aCurrentRegion] - C;
			if (A*aIntendedDistance*aIntendedDistance + 2 * B*aIntendedDistance >= C1) {
				double D = B*B + A*C1;
				// if ( D < 0 ) D = 0;
				aIntendedDistance = B > 0 ? C1 / (sqrt(D) + B) : (sqrt(D) - B) / A;
				++nextRegion;
				if (nextRegion >= _numberOfRegions) nextRegion = -1;
			}
		}
		else {
			// travelling inwards
			// -> first check inner cylinder
			double C1 = _radiiSquared[aCurrentRegion - 1] - C;
			double D = B*B + A*C1;
			if (D > 0) {
				double distanceToInner = -C1 / (sqrt(D) - B);
				if (distanceToInner < aIntendedDistance) {
					aIntendedDistance = distanceToInner; --nextRegion;
				}
			}
			else {
				// missing inner cylinder
				C1 = _radiiSquared[aCurrentRegion] - C;
				if (A*aIntendedDistance*aIntendedDistance + 2 * B*aIntendedDistance >= C1) {
					aIntendedDistance = (sqrt(fabs(D)) - B) / A;
					++nextRegion;
					if (nextRegion >= _numberOfRegions) nextRegion = -1;
				}
			}
		}
	}
	else if (B < 0) {
		// if outside, only need to check for intersection if B<0
		double C1 = _radiiSquared[_numberOfRegions - 1] - C;
		double D = B*B + A*C1;
		if (D > 0) {
			double distanceToInner = -C1 / (sqrt(D) - B);
			//double distanceToInner = (sqrt(D) - B)/A;
			if (distanceToInner < aIntendedDistance) {
				aIntendedDistance = distanceToInner; nextRegion = _numberOfRegions - 1;
			}
		}
	}

	if (nextRegion != aCurrentRegion) {
		if (nextMedium && nextRegion >= 0) nextMedium = getMedium(nextRegion);
	}

	return nextRegion;

}
Medium ViewRayMlcCylinders::getMedium(int aRegion) const {
	if (aRegion == 0) return UnknownMaterial;
	if (aRegion == 1) return Tungsten;
	if (aRegion == 2) return Tungsten;
	return UnknownMaterial;
}

class ViewRayLeaves {
public:
	/**
	* \brief Constructor
	*
	* \param[in] aRegionOffset index of first region.
	* \param[in] aLeafWidth Leaf width CM.
	* \param[in] aLeftFocus Left focus of leaf.
	* \param[in] aRightFocus Right focus of leaf.
	* \param[in] aGapFocus Focus of gap.
	* \param[in] aZshift Zshift
	*
	*/
	ViewRayLeaves(int aRegionOffset, double aLeafWidth, Vector & aLeftFocus, Vector & aRightFocus, Vector & aGapFocus, double aZshift);
	/**
	* \brief Set the MLC segments.
	*
	* \param aListOfSegments
	*/
	void setSegments(const vector<vector<pair<double, double> > > &aListOfSegments);
	/**
	* \brief Which region particle is in.
	*
	* \return Which region particle is in.
	*/
	int isWhere(unsigned int aTimeIndex, const Vector &aX, Medium & material) const;
	/**
	* \brief What is the next region along the particle path.
	*
	* \param[in] aTimeIndex Which MLC segment is currently used.
	* \param[in] aCurrentRegion Which region is particle in.
	* \param[in] aX Position of particle.
	* \param[in] aU Direction of particle.
	* \param[in,out] aIntendedDistance Max distance to search for next region.  This value is reset to distance to next region
	* if it is less than aIntendedDistance.
	* \param[out] material Material of next region along path.
	*/
	int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const;
	/**
	* \brief Min z coord inside MLC model.
	*
	*/
	double getMinZ() { return _leaf[0]->getMinZ(); }
	/**
	* \brief Max z coord inside MLC model.
	*
	*/
	double getMaxZ() { std::vector<ViewRayLeaf*>::size_type n = _leaf.size(); return _leaf[n - 1]->getMinZ() + Constants::LEAF_THICKNESS; }
	/**
	* \brief One class for each leaf.
	*
	*/
	vector<ViewRayLeaf*> _leaf;
	/**
	* \brief Destructor.
	*
	*/
	~ViewRayLeaves();
private:
	/**
	* \brief The z shift of the leaves.
	*
	*/
	double _zshift;
};
ViewRayLeaves::ViewRayLeaves(int aRegionOffset, double aLeafWidth, Vector & aLeftFocus, Vector & aRightFocus, Vector & aGapFocus, double aZshift) : _zshift(aZshift) {
	int n = Constants::NLEAF_PAIRS_TOTAL;
	_leaf.resize(n);
	int nd2 = n / 2;
	double leftMost = -(nd2)*aLeafWidth;
	for (int ii = 0; ii < n; ++ii) _leaf[ii] = new ViewRayLeaf(aRegionOffset + ii * 4, leftMost + ii*aLeafWidth, aLeftFocus, aRightFocus, aGapFocus, aZshift);
}
ViewRayLeaves::~ViewRayLeaves() {
	int n = Constants::NLEAF_PAIRS_TOTAL;
	for (int ii = 0; ii < n; ++ii) delete _leaf[ii];
}
int ViewRayLeaves::isWhere(unsigned int aTimeIndex, const Vector &aLocation, Medium & material) const {
	using namespace Constants;
	double srcToIso = SOURCE_TO_ISO_DISTANCE;

	// Which leaf are we in
	double y = aLocation._y - _zshift;
	double dz = srcToIso - aLocation._z;

	int n = NLEAF_PAIRS_TOTAL;
	double farLeft = -(n / 2)*LEAF_THICKNESS;
	double farRight = (n / 2)*LEAF_THICKNESS;

	int region = (dz <= 0 || srcToIso*y < dz*(farLeft) || srcToIso*y >= dz*farRight) ? -1 :
		int((srcToIso * y - dz * (farLeft)) / (dz * LEAF_THICKNESS));

	if (region >= 0) {
		int localRegion = _leaf[region]->isWhere(aTimeIndex, aLocation, material);
		return region * 4 + localRegion;
	}

	return region;

}
void ViewRayLeaves::setSegments(const vector<vector<pair<double, double> > > &aListOfSegments) {
	using namespace Constants;
	// Validation.
	std::vector<vector<pair<double, double> > >::size_type nsegment = aListOfSegments.size();
	if (nsegment < 1) return;
	for (unsigned int j = 0; j < nsegment; ++j) {
		if ((int)aListOfSegments[j].size() != NLEAF_PAIRS_REAL) {
			ZeusLog("ViewRayLeaves::setSegments: segment %d has %d instead of %d leaf positions\n", j + 1, int(aListOfSegments[j].size()), NLEAF_PAIRS_REAL);
		}
	}

	// Set the positions of the real leaves.
	vector<pair<double, double> > leafPositions(nsegment);
	for (int leafPair = 0; leafPair < NLEAF_PAIRS_REAL; ++leafPair) {
		for (unsigned int j = 0; j < nsegment; ++j) leafPositions[j] = aListOfSegments[j][leafPair];
		_leaf[leafPair + 4]->setOpenings(leafPositions);
	}

	// Set the fake leaves to closed.
	vector<pair<double, double> > outPositions(nsegment);
	for (unsigned int j = 0; j < nsegment; ++j) outPositions[j] = pair<double, double>(-20, -20);
	for (int leafPair = 0; leafPair < Constants::N_LEAF_PAIRS_FAKE; ++leafPair) {
		_leaf[leafPair]->setOpenings(outPositions);
		_leaf[leafPair + 34]->setOpenings(outPositions);
	}

}

class ViewRayMLC : public IGeometry {
	/**
	* \brief A class to model the ViewRay MLC.
	*
	* A leaf is modelled with an innner and outer half, so there is one model for each leaf for each half.
	*/
public:
	/**
	* \brief Constructor.  Initialize the MLC class.
	*
	* Particles that travel above this value are no longer tracked.
	*
	* @param[in] aMaxZCoord The maximum z value in the universe.
	*/
	ViewRayMLC(double aMaxZCoord = 105);
	/**
	* \brief Set the MLC segments inside this class.
	*
	* \param[in] aListOfSegments  A vector of leaef segments.
	*/
	void setSegments(const vector<vector<pair<double, double> > > &aListOfSegments);
	/**
	* \brief Which region is the particle in.
	*
	* \param[in] aTimeIndex Which MLC segment is currently used.
	* \param[in] aX Position of the particle.
	* \param[out] material Material the region is in.
	* \return The region the particle is in.
	*
	*/
	int isWhere(unsigned int aTimeIndex, const Vector &aX, Medium & material) const;
	/**
	* \brief What is the next region along the particle path.
	*
	* \param[in] aTimeIndex Which MLC segment is currently used.
	* \param[in] aCurrentRegion Which region is particle in.
	* \param[in] aX Position of particle.
	* \param[in] aU Direction of particle.
	* \param[in,out] aIntendedDistance Max distance to search for next region.  This value is reset to distance to next region
	* if it is less than aIntendedDistance.
	* \param[out] material Material of next region along path.
	*/
	int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const;
	/**
	* \brief Return the number of regions in this class.
	*
	*/
	int getNumRegions() { return Constants::NLEAF_PAIRS_TOTAL * 4 * 2; }
	/**
	* \brief Set the maximum z coordinate in the universe.
	*
	* @param[in] aMaxZCoord The maximum z value in the universe.
	*/
	void setMaxZCoord(double aMaxZCoord) { _maxZ = aMaxZCoord; }
	/**
	* \brief Destructor.
	*
	*/
	~ViewRayMLC();
private:
	/**
	* \brief Store a model of the inner portion of the leaf.
	* One class for each leaf.
	*
	*/
	ViewRayLeaves *_innerLeaves;
	/**
	* \brief Store a model of the outer portion of the leaf.
	* One class for each leaf.
	*
	*/
	ViewRayLeaves *_outerLeaves;
	/**
	* \brief
	*
	* Leaves are modelled with an upper and lower half.
	* Each half fits inside a cylinder.
	*/
	vector<ViewRayLeaves*> _cylinderOfLeaves;
	/**
	* \brief A vector of air gaps between the leaves.
	*
	* Each air gap is one class inside this vector.
	*
	*/
	vector<ViewRayGaps*> _gaps;
	/**
	* \brief The radius of the inner cylinder.
	*
	* The source is inside this radius.
	*
	*/
	double _innerRadiusSquared;
	/**
	* \brief The radius of the cylinder separating the upper and lower
	* halves of the leaf.
	*
	*/
	double _middleRadiusSquared;
	/**
	* \brief The radius of the cylinder bounding the outside of the MLC.
	*
	* Outside this radius is the base region, the region between the MLC
	* and the gantry.
	*/
	double _outerRadiusSquared;
	/**
	* \brief Two cylinders lining the inside and the outside of the MLC.
	*
	*/
	ViewRayMlcCylinders _boundingCylinders;
	/**
	* \brief The maximum z value in the universe.
	*
	*/
	double _maxZ;

};
int ViewRayMLC::isWhere(unsigned int aTimeIndex, const Vector &aLocation, Medium & material) const {

	material = UnknownMaterial;
	int regionIndex = -1;

	double zc = aLocation._z - Constants::SOURCE_TO_ISO_DISTANCE;
	double dsquared = aLocation._x*aLocation._x + zc*zc;

	if (dsquared < _innerRadiusSquared) {
		return regionIndex;
	}
	else if (dsquared >= _innerRadiusSquared && dsquared < _middleRadiusSquared) {
		return _innerLeaves->isWhere(aTimeIndex, aLocation, material);
	}
	else if (dsquared >= _middleRadiusSquared && dsquared < _outerRadiusSquared) {
		int cylinderRegion = _outerLeaves->isWhere(aTimeIndex, aLocation, material);
		return cylinderRegion == -1 ? cylinderRegion : Constants::NLEAF_PAIRS_TOTAL * 4 + cylinderRegion;
	}
	else if (dsquared >= _outerRadiusSquared) {
		return regionIndex;
	}

	return regionIndex;

}
int ViewRayMLC::nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const {
	using namespace Constants;
	if (aCurrentRegion >= 0) {
		int irCyl = aCurrentRegion / (4 * NLEAF_PAIRS_TOTAL);
		int irCylNew = _boundingCylinders.nextRegionAlongPath(aTimeIndex, irCyl + 1, aX, aU, aIntendedDistance, material);
		int leafRegion = aCurrentRegion - 4 * NLEAF_PAIRS_TOTAL * irCyl;
		int leafPair = leafRegion / 4;
		int pairLocalRegion = leafRegion - 4 * leafPair;
		int zplane = 2 * leafPair;
		if (pairLocalRegion < 3) {
			// we are in the area divided by the leafs
			int newPairLocalRegion = _cylinderOfLeaves[irCyl]->_leaf[leafPair]->nextRegionAlongPath(aTimeIndex, pairLocalRegion, aX, aU, aIntendedDistance, material);
			int newZplane = _gaps[irCyl]->nextRegionAlongPath(aTimeIndex, zplane, aX, aU, aIntendedDistance, material);
			if (newZplane < 0) return -1;
			if (newZplane != zplane) {
				int newLeafPair = newZplane > zplane ? leafPair : leafPair - 1;
				material = Air;
				return 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * newLeafPair + 3;
			}
			else if (newPairLocalRegion != pairLocalRegion) {
				if (newPairLocalRegion >= 0) material = newPairLocalRegion == 1 ? Air : Tungsten;
				return newPairLocalRegion >= 0 ? 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * leafPair + newPairLocalRegion : -1;
			}
		}
		else {
			// we are in the air gap between the leaves
			++zplane;
			int newZplane = _gaps[irCyl]->nextRegionAlongPath(aTimeIndex, zplane, aX, aU, aIntendedDistance, material);
			if (newZplane < 0) return -1;
			if (newZplane != zplane) {
				int newLeafPair = newZplane > zplane ? leafPair + 1 : leafPair;
				int newPairLocalRegion = _cylinderOfLeaves[irCyl]->_leaf[newLeafPair]->isWhereNoAirGap(aTimeIndex, aX + aU*aIntendedDistance, material);
				if (newPairLocalRegion >= 0) material = newPairLocalRegion == 1 ? Air : Tungsten;
				return newPairLocalRegion >= 0 ? 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * newLeafPair + newPairLocalRegion : -1;
			}
		}

		// if we are here, we did not intersect any of the planes of the current MLC layer
		if (irCylNew == irCyl + 1) return aCurrentRegion;
		if (irCylNew == 0 || irCylNew == -1) return -1;
		// if we are here, we have entered a new MLC layer
		irCyl = irCylNew - 1;
		Vector xnew(aX + aU*aIntendedDistance);
		zplane = _gaps[irCyl]->isWhere(aTimeIndex, xnew);
		if (zplane < 0) return -1;
		leafPair = zplane / 2;
		if (leafPair * 2 == zplane) {
			// we have entered an area occupied by the leaves
			int newPairLocalRegion = _cylinderOfLeaves[irCyl]->_leaf[leafPair]->isWhere(aTimeIndex, xnew, material);
			if (newPairLocalRegion < 0) return -1;
			material = newPairLocalRegion == 1 ? Air : Tungsten;
			return 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * leafPair + newPairLocalRegion;
		}
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		//if( material ) material = Air;
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		material = Air;
		return 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * leafPair + 3;

	}
	// if here, we are outside.
	int irCyl = _boundingCylinders.isWhere(aTimeIndex, aX, material);
	if (irCyl == 0 || irCyl == -1) {
		if ((irCyl == 0 && aU._z >= 0) || (irCyl == -1 && aU._z <= 0)) return -1;
		double tsave = aIntendedDistance; double tleft = aIntendedDistance; Vector xnew(aX);
		while (1) {
			double tt = tleft;
			int irCylNew = _boundingCylinders.nextRegionAlongPath(aTimeIndex, irCyl, xnew, aU, tt, material);
			if (tt != tt) {
				exitApp("Outside MLC cylinders");
			}
			if (irCylNew < 0) return -1;
			if (irCylNew == irCyl) return aCurrentRegion;
			xnew += aU*tt; tleft -= tt; irCyl = irCylNew;
			if (xnew._z > _maxZ) return -1;
			if (irCyl == 1 || irCyl == 2) {
				int zplane = _gaps[irCyl - 1]->isWhere(aTimeIndex, xnew);
				if (zplane >= 0) {
					aIntendedDistance = tsave - tleft;
					int leafPair = zplane / 2;
					if (leafPair * 2 == zplane) {
						int newPairLocalRegion = _cylinderOfLeaves[irCyl - 1]->_leaf[leafPair]->isWhere(aTimeIndex, xnew, material);
						if (newPairLocalRegion < 0) return -1;
						material = newPairLocalRegion == 1 ? Air : Tungsten;
						return 4 * NLEAF_PAIRS_TOTAL*(irCyl - 1) + 4 * leafPair + newPairLocalRegion;
					}
					else {
						material = Air;
						return 4 * NLEAF_PAIRS_TOTAL*(irCyl - 1) + 4 * leafPair + 3;
					}
				}
				break;
			}
		}
		while (1) {
			double tt = tleft;
			int zplane = _gaps[irCyl - 1]->nextRegionAlongPath(aTimeIndex, -1, xnew, aU, tt, material);
			if (tt != tt) {
				exitApp("Outside MLC second while gaps");
			}
			int irCylNew = _boundingCylinders.nextRegionAlongPath(aTimeIndex, irCyl, xnew, aU, tt, material);
			if (tt != tt) {
				exitApp("Outside MLC second while cylinders");
			}
			if (irCylNew != irCyl) {
				if (irCylNew == 0 || irCylNew < 0) return -1;
				xnew += aU*tt; tleft -= tt; irCyl = irCylNew;
				if (xnew._z > _maxZ) return -1;
				zplane = _gaps[irCyl - 1]->isWhere(aTimeIndex, xnew);
				if (zplane >= 0) {
					aIntendedDistance = tsave - tleft;
					int leafPair = zplane / 2;
					if (leafPair * 2 == zplane) {
						int newPairLocalRegion = _cylinderOfLeaves[irCyl - 1]->_leaf[leafPair]->isWhere(aTimeIndex, xnew, material);
						if (newPairLocalRegion < 0) return -1;
						material = newPairLocalRegion == 1 ? Air : Tungsten;
						return 4 * NLEAF_PAIRS_TOTAL*(irCyl - 1) + 4 * leafPair + newPairLocalRegion;
					}
					else {
						material = Air;
						return 4 * NLEAF_PAIRS_TOTAL*(irCyl - 1) + 4 * leafPair + 3;
					}
				}
			}
			else if (zplane >= 0) {
				aIntendedDistance = tsave - tleft + tt;
				int leafPair = zplane / 2;
				if (leafPair * 2 == zplane) {
					int newPairLocalRegion = _cylinderOfLeaves[irCyl - 1]->_leaf[leafPair]->isWhere(aTimeIndex, xnew, material);
					if (newPairLocalRegion < 0) return -1;
					material = newPairLocalRegion == 1 ? Air : Tungsten;
					return 4 * NLEAF_PAIRS_TOTAL*(irCyl - 1) + 4 * leafPair + newPairLocalRegion;
				}
				else {
					material = Air;
					return 4 * NLEAF_PAIRS_TOTAL*(irCyl - 1) + 4 * leafPair + 3;
				}
			}
			tleft -= tt; irCyl = irCylNew;
			if (tleft < 1e-8) return aCurrentRegion;
		}
	}
	// If here, we are outside but we are actually inside the cylinders with the MLC layers
	// This can be true because we are outside of the leaves, but it can also happen due to numerical roundoff errors.
	// First check for this possibility
	Vector xp(aX);
	double r2 = xp._x*xp._z + xp._y*xp._z;
	if ((r2 < _innerRadiusSquared && aU._z > 0) || (r2 > _outerRadiusSquared && aU._z < 0)) return -1;
	--irCyl;
	int zplane = _gaps[irCyl]->nextRegionAlongPath(aTimeIndex, -1, aX, aU, aIntendedDistance, material);
	if (zplane < 0) return -1;
	if (aIntendedDistance < -1e-3) ZeusLog("ViewRayMLC::nextRegionAlongPath: negative step %g to enter from outside\n", aIntendedDistance);
	// check again this logic !!!!!!!!!!!!!!!!!
	_boundingCylinders.nextRegionAlongPath(aTimeIndex, irCyl + 1, aX, aU, aIntendedDistance, material);
	int leafPair = zplane / 2;
	if (leafPair * 2 == zplane) {
		Vector xnew(aX + aU*aIntendedDistance);
		int newPairLocalRegion = _cylinderOfLeaves[irCyl]->_leaf[leafPair]->isWhere(aTimeIndex, xnew, material);
		if (newPairLocalRegion < 0) return -1;
		material = newPairLocalRegion == 1 ? Air : Tungsten;
		return 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * leafPair + newPairLocalRegion;
	}
	material = Air;
	return 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * leafPair + 3;

}
ViewRayMLC::ViewRayMLC(double aMaxZCoord) : _maxZ(aMaxZCoord) {
	using namespace Constants;
	double sourceToIso = SOURCE_TO_ISO_DISTANCE;

	_innerRadiusSquared = MLC_INNER_RADIUS * MLC_INNER_RADIUS;
	double middleRadius = (MLC_OUTER_RADIUS + MLC_INNER_RADIUS) / 2.0;
	_middleRadiusSquared = middleRadius * middleRadius;
	_outerRadiusSquared = MLC_OUTER_RADIUS * MLC_OUTER_RADIUS;

	int regionIndex = 0;
	Vector innerLeftFocus(0.0, 0.0, sourceToIso);
	Vector innerRightFocus(INNER_RIGHT_FOCUS, 0.0, sourceToIso);
	Vector innerGapFocus(0.0, INNER_GAP_FOCUS, sourceToIso);

	_innerLeaves = new ViewRayLeaves(regionIndex, LEAF_THICKNESS, innerLeftFocus, innerRightFocus, innerGapFocus, INNER_Z_SHIFT);
	regionIndex = NLEAF_PAIRS_TOTAL * 4;

	Vector outerLeftFocus(OUTER_LEFT_FOCUS, 0.0, sourceToIso);
	Vector outerRightFocus(0.0, 0.0, sourceToIso);
	Vector outerGapFocus(0.0, OUTER_GAP_FOCUS, sourceToIso);
	_outerLeaves = new ViewRayLeaves(regionIndex, LEAF_THICKNESS, outerLeftFocus, outerRightFocus, outerGapFocus, OUTER_Z_SHIFT);

	_cylinderOfLeaves.push_back(_innerLeaves);
	_cylinderOfLeaves.push_back(_outerLeaves);

	_gaps.push_back(new ViewRayGaps(innerGapFocus, sourceToIso, INNER_Z_SHIFT));
	_gaps.push_back(new ViewRayGaps(outerGapFocus, sourceToIso, OUTER_Z_SHIFT));

}
ViewRayMLC::~ViewRayMLC() {
	for (size_t i = 0; i < _gaps.size(); i++) delete _gaps[i];
	_gaps.clear();
	for (size_t i = 0; i < _cylinderOfLeaves.size(); i++) delete _cylinderOfLeaves[i];
	_cylinderOfLeaves.clear();
}
void ViewRayMLC::setSegments(const vector<vector<pair<double, double> > > &aListOfSegments) {
	_innerLeaves->setSegments(aListOfSegments);
	_outerLeaves->setSegments(aListOfSegments);
}


class Envelope : public IGeometry {
	/**
	* \brief Class to contain the universe of the simulation.
	*
	* Class contains the base region, which limits the extent of the simulation and marks
	* the boundaries of the inscribed regions.  also then are one or more inscribed regions,
	* which contain the most detailed boundaries.
	**/
public:

	/**
	* \brief Construct the universe.
	*
	* This universe will have one base region, and any number of inscribed regions.
	* \param aBaseRegion - geometry of the base region.
	* \param aInscribedRegions - vector containing the geometry of one or more inscribed regions.
	*/
	Envelope(IGeometry *aBaseRegion, vector<IGeometry *> aInscribedRegions) {

		_baseRegions = aBaseRegion->getNumRegions();
		_baseToInscribedWhere.resize(_baseRegions);
		_baseToInscribedNext.resize(_baseRegions);

		_baseRegion = aBaseRegion;
		_nInscribedRegions = aInscribedRegions.size();
		for (size_t i = 0; i < _nInscribedRegions; ++i) _inscribedRegions.push_back(aInscribedRegions[i]);

		_localStart.push_back(_baseRegions);

		for (size_t i = 0; i < _nInscribedRegions; ++i)
			_localStart.push_back(_localStart[i] + _inscribedRegions[i]->getNumRegions());

	}

	/**
	* \brief Set the relation between a base region and an inscribed region when calling isWhere.
	*
	* If isWhere results in a base region, this relation will yield the appropriate inscribed
	* region to call isWhere.
	* Zero means this base region has no inscribed region.
	* \param aBaseRegion - If isWhere is in this base region, then use index in next argument
	* to determine which inscribed region to call isWhere.
	* \param aInscribedRegion - Index of inscribed region inside this base region.
	*/
	void setBaseToInscribedWhere(int aBaseRegion, int aInscribedRegion) { _baseToInscribedWhere[aBaseRegion] = aInscribedRegion; }

	/**
	* \brief Set the relation between a base region and an inscribed region when calling nextRegionAlongPath.
	*
	* If nextRegionAlongPath results in a base region, this relation will yield the appropriate inscribed
	* region to call nextRegionAlongPath.
	* Zero means this base region has no inscribed region.
	* \param aBaseRegion - If nextRegionAlong is this base region, then use index in next argument
	* to determine which inscribed region to call nextRegionAlongPath.
	* \param aInscribedRegion - Index of inscribed region inside this base region.
	*/
	void setBaseToInscribedNext(int aBaseRegion, int aInscribedRegion) { _baseToInscribedNext[aBaseRegion] = aInscribedRegion; }

	/**
	* \brief Destructor.
	*
	* Clear class memory, which are the vectors maintaining relations between base and inscribed regions.
	*/
	~Envelope() {
		_nInscribedRegions = 0;
		_inscribedRegions.clear(); _baseToInscribedWhere.clear(); _baseToInscribedNext.clear(); _localStart.clear();
	}

	/**
	* \brief Which region of the envelope are we in.
	*
	* Answer can be one of the base regions, or if we are in an inscribed region,
	* then isWhere will be called on the inscribed region.
	* \param aTimeIndex - which MLC segmet to use.
	* \param aX - position of particle.
	* \param material Returned material particle is in.
	* \return Which region particle is in.
	*/
	int isWhere(unsigned int aTimeIndex, const Vector &aX, Medium & material) const {
		int region = _baseRegion->isWhere(aTimeIndex, aX, material);
		if (region > -1) {
			int inscribedRegion = _baseToInscribedWhere[region];
			if (inscribedRegion > 0) {
				int newregion = _inscribedRegions[inscribedRegion - 1]->isWhere(aTimeIndex, aX, material);
				if (newregion != -1) return (_localStart[inscribedRegion - 1] + newregion);
			}
			material = Air;
			return region;
		}
		return -1;
	}

	/**
	* \brief What is the next region the particle will encounter.
	*
	* Answer can be one of the base regions, or if we are in an inscribed region,
	* then nextRegionAlongPath will be called on the inscribed regions.
	* \param aTimeIndex - which MLC segment to use.
	* \param aCurrentRegion - which region are we in.
	* \param aX - position of particle.
	* \param aU - direction of particle.
	* \param aIntendedDistance - intended distance.  This will be reset if we intersect another region in a shorter distance.
	* \param material - Returned material particle is in.
	* \return Index of next region along path of particle.
	*/
	int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, double &aIntendedDistance, Medium& material) const {

		if (aCurrentRegion < _localStart[0]) {

			double distance = aIntendedDistance;
			Medium baseMedium = material;

			int newBaseRegion = _baseRegion->nextRegionAlongPath(aTimeIndex, aCurrentRegion, aX, aU, distance, baseMedium);
			if (newBaseRegion == -1) { aIntendedDistance = distance; material = baseMedium; return newBaseRegion; }

			int inscribedRegion = _baseToInscribedNext[newBaseRegion];
			if (inscribedRegion > 0 && inscribedRegion - 1 < (int)_nInscribedRegions) {
				int newregion = _inscribedRegions[inscribedRegion - 1]->nextRegionAlongPath(aTimeIndex, -1, aX, aU, aIntendedDistance, material);
				if (newregion != -1) return (_localStart[inscribedRegion - 1] + newregion);
			}

			aIntendedDistance = distance;
			material = baseMedium;
			return newBaseRegion;

		}

		// We are in an inscribed region
		if (aCurrentRegion < _localStart[1]) {

			int inscribedRegionIndex = 1; // MLC
			int baseRegion = 2;
			int newInscribedRegion = _inscribedRegions[inscribedRegionIndex - 1]->nextRegionAlongPath(aTimeIndex, aCurrentRegion - _localStart[inscribedRegionIndex - 1], aX, aU, aIntendedDistance, material);

			double distanceBase = 9999;
			Medium baseMedium = material;
			int newBaseRegion = _baseRegion->nextRegionAlongPath(aTimeIndex, baseRegion, aX, aU, distanceBase, baseMedium);

			if (newBaseRegion != baseRegion && distanceBase < aIntendedDistance) {
				if (newBaseRegion < 0) return newBaseRegion;
				int newInscribedRegionIndex = _baseToInscribedNext[newBaseRegion];
				if (newInscribedRegionIndex > 0) {
				}
				aIntendedDistance = distanceBase;
				material = baseMedium;
				return newBaseRegion;
			}

			if (newInscribedRegion < 0) {
				material = baseMedium;
				return newBaseRegion;
			}

			return _localStart[inscribedRegionIndex - 1] + newInscribedRegion;

		}
		else if (aCurrentRegion < _localStart[2]) {

			int inscribedRegionIndex = 2; // Gantry
			int baseRegion = 4;
			int newInscribedRegion = _inscribedRegions[inscribedRegionIndex - 1]->nextRegionAlongPath(aTimeIndex, aCurrentRegion - _localStart[inscribedRegionIndex - 1], aX, aU, aIntendedDistance, material);

			double distanceBase = 9999;
			Medium baseMedium = material;
			int newBaseRegion = _baseRegion->nextRegionAlongPath(aTimeIndex, baseRegion, aX, aU, distanceBase, baseMedium);

			if (newBaseRegion != baseRegion && distanceBase < aIntendedDistance) {
				if (newBaseRegion < 0) return newBaseRegion;
				int newInscribedRegionIndex = _baseToInscribedNext[newBaseRegion];
				if (newInscribedRegionIndex > 0) {
				}
				aIntendedDistance = distanceBase;
				material = baseMedium;
				return newBaseRegion;
			}

			if (newInscribedRegion < 0) {
				material = baseMedium;
				return newBaseRegion;
			}

			return _localStart[inscribedRegionIndex - 1] + newInscribedRegion;

		}

		// 
		double distance = aIntendedDistance;
		Medium baseMedium = material;
		int newBaseRegion = _baseRegion->nextRegionAlongPath(aTimeIndex, -1, aX, aU, distance, baseMedium);
		if (newBaseRegion == -1) return newBaseRegion;
		int inscribedRegion = _baseToInscribedNext[newBaseRegion];
		if (inscribedRegion > 0) {
			int newregion = _inscribedRegions[inscribedRegion - 1]->nextRegionAlongPath(aTimeIndex, -1, aX, aU, aIntendedDistance, material);
			if (newregion != -1) return (_localStart[inscribedRegion - 1] + newregion);
		}
		if (newBaseRegion == aCurrentRegion || newBaseRegion < 0) { material = baseMedium; aIntendedDistance = distance; return newBaseRegion; }

		return -1;
	}

	/**
	* Return maximum number of regions in the envelope.
	*/
	int getNumRegions() { return _localStart[2]; }

private:
	/**
	*  Number of inscribed regions in the envelope
	*/
	size_t _nInscribedRegions;
	/**
	* The geometry of the base region.
	*/
	IGeometry *_baseRegion;
	/**
	* A vector of inscribed region geometries.
	*/
	vector<IGeometry *> _inscribedRegions;
	/**
	*  Given a region code from the base region, which inscribed region
	* to we interrogate to find the region value using isWhere.
	*/
	vector<int> _baseToInscribedWhere;
	/**
	* Given a region code from the base region, which inscribed region
	* to we interrogate to find the region value using nextRegionAlongPath.
	*/
	vector<int> _baseToInscribedNext;

	/**
	* Upper bound on region number.  First entry is the number of base
	* regions, subsequent entries are the number of inscribed regions.
	* Each entry is offset by sum of all previous entries.
	*
	*/
	vector<int> _localStart;

	/**
	* Number of base regions in the envelope.
	*/
	int _baseRegions;

};

template <class DataEntry> class InterpolationData {
	/**
	* \brief A template class for managing interpolation data.
	*/
public:

	/*! \brief Constructor for a not initialized instance. */
	InterpolationData() : _ndata(0), _data(0), _bIndex(0), _aIndex(0), _nmat(0), _isValid(false) { };

	/**
	* \brief Constructor.
	*
	* Read in the cross section data.
	* \param aFileName Which file is the cross section data in.
	* \param aNumMaterials Cross section data for how many materials is in the file.
	* \param useInverseEnergy If set to true, data will be interpolated as a function of 1/E instead of E.
	*/
	InterpolationData(const string &aFileName, const int aNumMaterials, bool useInverseEnergy = false);

	/**
	* \brief Constructor.
	*
	* Construct using the data passed as argument
	* \param[in] aName    The name of the data set
	* \param[in] emin     The minimum energy of the data for each material
	* \param[in] emax     The maximum energy of the data for each material
	* \param[in] nEnergy  The number of energy bins (all data sets have the same number of bins)
	* \param[in] aData    The data. The size of aData vector gives the number of materials, aData[i] is the actual data for material i with nEnergy points
	* Remark: emin.size() must be equal to emax.size() and to aData.size().
	*/
	InterpolationData(const string &aName, const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &aData);
	/**
	* Is the instance valid?
	*/
	bool isValid() const { return _isValid; };
	/**
	* Write all the cross section data to the logger.
	*/
	void dumpData() const;
	/**
	* Destructor.
	* All the data is deleted, so it will have to be loaded again.
	*/
	~InterpolationData();
	/**
	* Query name of cross section file.
	* \return Name of cross section file.
	*/
	const string& fileName() const { return _fileName; };
	/**
	* Interpolate, throwing an exception if x is out of initialized range
	*/
	inline double interpolate(double x, int matid) const {
		int i = (int)(_aIndex[matid] * (x - _bIndex[matid]));
		if (i < 0 || i >= _ndata[matid]) {
			ZeusLog("Spline interpolation error for data %s: index %d out of range. x=%g matid=%d\n", _fileName.c_str(), i, x, matid);
			exitApp("Spline interpolation out of range");
		}
		return _data[matid][i].compute(x);
	};
	/**
	* Interpolate, returning zero if x below range
	*/
	inline double interpolate1(double x, int matid) const {
		double dx = x - _bIndex[matid];
		double result = 1e30;
		if (dx > 0) {
			int i = (int)(_aIndex[matid] * dx);
			if (i >= _ndata[matid]) {
				ZeusLog("Spline interpolation error for data %s: index %d out of range. x=%g matid=%d\n", _fileName.c_str(), i, x, matid);
				exitApp("Spline interpolation out of range");
			}
			result = _data[matid][i].compute(x);
		}
		return result;
	};

	/*! Returns the number of materials */
	inline int numMaterials() const { return _nmat; }

	/**
	* Lower bound of initialized data range
	*/
	inline double lowerBound(int matid) const { return _bIndex[matid]; };

	/**
	* Load data from the data folder \a aDataFolder using the data file \a aFileName and return a InterpolationData instance for the data
	*/
	static InterpolationData* loadData(const string &aDataFolder, const char *aFileName, int aNumMaterials, bool useInverseEnergy = false);

private:
	/**
	* Array of length number of materials, with length of cross section tables for each material.
	*/
	int *_ndata;
	/**
	* Interpolation coefficients.
	* Two dimensional array, first dimension is material.
	*/
	DataEntry** _data;
	/*
	* First entry in table.
	* Array, one entry for each table.
	*/
	double *_bIndex;
	/**
	* Slope of the table values.
	* Array, one entry for each table.
	*/
	double *_aIndex;
	/**
	* Number of materials.
	*/
	int _nmat;
	/**
	* Validity of instance
	*/
	bool _isValid;
	/**
	* Name of file from which data was loaded.
	*/
	string _fileName;
};
template<class DataEntry>
InterpolationData<DataEntry>::InterpolationData(const string &aName, const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &aData) {
	_fileName = aName;
	_nmat = aData.size();
	_ndata = new int[_nmat];
	_bIndex = new double[_nmat];
	_aIndex = new double[_nmat];
	_data = new DataEntry*[_nmat];
	for (int j = 0; j < _nmat; ++j) _data[j] = 0;

	double *xval = new double[nEnergy];

	for (int j = 0; j < _nmat; j++) {
		double deltax = (emax[j] - emin[j]) / (nEnergy - 1);
		for (int i = 0; i < nEnergy; ++i) xval[i] = emin[j] + deltax*i;
		_ndata[j] = nEnergy;
		_data[j] = DataEntry::prepareData(xval, aData[j], _ndata[j]);
		_bIndex[j] = xval[0];
		_aIndex[j] = ((double)_ndata[j] - 1) / (xval[_ndata[j] - 1] - xval[0]);
	}

	delete[] xval;
	_isValid = true;

}

template<class DataEntry>
InterpolationData<DataEntry>::InterpolationData(const string &aFileName, const int aNumMaterials, bool useInverseEnergy) :
_ndata(0), _data(0), _bIndex(0), _aIndex(0), _nmat(aNumMaterials), _isValid(false), _fileName(aFileName) {

	const static char *func = "InterpolationData::InterpolationData";
	FILE *fp = fopen(aFileName.c_str(), "r");
	if (fp == NULL) {
		exitApp("Failed to open material file");
	}

	_ndata = new int[_nmat];
	_bIndex = new double[_nmat];
	_aIndex = new double[_nmat];
	_data = new DataEntry*[_nmat];
	for (int j = 0; j < _nmat; ++j) _data[j] = 0;

	/* skip the header */
	int nitems = fscanf(fp, "%*[^\n] %*1[\n]");
	if (0) ZeusLog("Got %d items\n", nitems); // shut up compiler warning
	nitems = fscanf(fp, "%*[^\n] %*1[\n]");

	/* loop over materials */
	for (int j = 0; j < _nmat; j++) {

		/* skip a line */
		nitems = fscanf(fp, "%*[^\n] %*1[\n]");

		/* read number of data points */
		nitems = fscanf(fp, "%d\n", &_ndata[j]);

		if (_ndata[j] < 2) {
			exitApp(" failed to read number of data points from file");
		}

		double *xval = new double[_ndata[j]];
		double *fval = new double[_ndata[j]];

		/* skip two lines */
		nitems = fscanf(fp, "%*[^\n] %*1[\n]");
		nitems = fscanf(fp, "%*[^\n] %*1[\n]");

		/* loop over data */
		for (int i = 0; i < _ndata[j]; i++) {
			nitems = fscanf(fp, "%lf %lf\n", &xval[i], &fval[i]);
			if (useInverseEnergy) {
				if (xval[i] == 0.0) {
					fclose(fp);
					exitApp("devided by zero!");
				}
				xval[i] = -1.0 / xval[i];
			}
		}

		_data[j] = DataEntry::prepareData(xval, fval, _ndata[j]);

		_bIndex[j] = xval[0];

		_aIndex[j] = ((double)_ndata[j] - 1) / (xval[_ndata[j] - 1] - xval[0]);

		delete[] xval; delete[] fval;

	} // number of materials

	/* close file */
	fclose(fp);
	_isValid = true;

}

template<class DataEntry>
InterpolationData<DataEntry>::~InterpolationData() {
	if (_ndata) delete[] _ndata;
	if (_bIndex) delete[] _bIndex;
	if (_aIndex) delete[] _aIndex;
	if (_data) {
		for (int j = 0; j < _nmat; ++j) if (_data[j]) delete[] _data[j];
		delete[] _data;
	}
}

template<class DataEntry>
InterpolationData<DataEntry>* InterpolationData<DataEntry>::loadData(const string &aDataFolder, const char *aFileName, int aNumMaterials, bool useInverseEnergy) {
	if (!aFileName) return 0;
	string fileName(aDataFolder);
	if (fileName.size() > 0) {
		char lastChar = fileName[fileName.size() - 1];
		if (lastChar != '/' && lastChar != 92) fileName += '/';
	}
	fileName += aFileName;
	return new InterpolationData<DataEntry>(fileName, aNumMaterials, useInverseEnergy);
}

struct LinearDataEntry {
	/*! \brief A structure for the 2 linear interpolation coefficients of an interpolation bin */
	/*! \name Linear interpolation coefficients */
	///@{
	double _a, _b;
	///@}
	/*! \brief Default constructor */
	LinearDataEntry() : _a(0), _b(0) {};
	/*! \brief Construct from input coefficients \a a and \a b. */
	LinearDataEntry(double a, double b) : _a(a), _b(b) {};
	/*! \brief Compute interpolation for \a x. */
	inline double compute(double x) const {
		return _a + _b*x;
	};
	/*! \brief Prepare spline coefficients from the input data \a xval, \a fval with \a ndat grid points. */
	static LinearDataEntry* prepareData(const double *xval, const double *fval, int ndat);
	/*! \brief Write coefficients and \a fileName using  LogProviders::InformationLog(). */
	void logData(const char *fileName) const;
};
LinearDataEntry* LinearDataEntry::prepareData(const double *xval, const double *fval, int ndat) {
	LinearDataEntry *data = new LinearDataEntry[ndat];
	for (int i = 1; i < ndat; i++) {
		double b = (fval[i] - fval[i - 1]) / (xval[i] - xval[i - 1]);
		double a = fval[i] - b*xval[i];
		data[i - 1] = LinearDataEntry(a, b);
	}
	data[ndat - 1] = data[ndat - 2];
	return data;
}

class HeadAttenuation {
	/*! \brief Provides cross section data for the materials in the ViewRay treatment head (MLC and below) */
public:
	/*! \brief Constructor, does nothing. Data has to be loaded using loadData() */
	HeadAttenuation() { };
	/*! \brief Destructor, does nothing. Data has to be unloaded using unloadData() */
	~HeadAttenuation() { };
	/*! \brief Is this instance valid? */
	static bool isValid() { return (_totData.size() > 0 && _compData.size() > 0 && _totData.size() == _compData.size()); };
	/*! \brief Load the data from the files in the folder \a aFolderName. */
	static int loadData(const string & aFolderName);
	/*! \brief Unload the cross section data. */
	static void unloadData();
	/*! \brief Returns the total attenuation coefficient for a photon with energy \a e in material \a matid. */
	inline static double computeAttenuation(double e, int matid) { return _totData[matid]->interpolate(e, 0); };
	/*! \brief Returns the Compton attenuation coefficient for a photon with energy \a e in material \a matid. */
	inline static double computeCompton(double e, int matid) { return _compData[matid]->interpolate(e, 0); };
private:
	static vector<InterpolationData<LinearDataEntry>* >  _totData;    //!< The total attenuation data
	static vector<InterpolationData<LinearDataEntry>* >  _compData;   //!< The Compton attenuation data
	static string                                        _folderName; //!< The folder from which the data was loaded.
};
int HeadAttenuation::loadData(const string & aFolderName) {
	const static char *func = "HeadAttenuation::loadData";
	if (aFolderName == _folderName && isValid()) return 0;

	unloadData();

	_totData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenTungstenTot", 1, false));
	_totData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenAirTot", 1, false));
	_totData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenGradCoilTot", 1, false));
	_totData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenRFShieldTot", 1, false));
	_totData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenTxPCBTot", 1, false));
	_totData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenTxFR4Tot", 1, false));

	_compData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenTungstenInco", 1, false));
	_compData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenAirInco", 1, false));
	_compData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenGradCoilInco", 1, false));
	_compData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenRFShieldInco", 1, false));
	_compData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenTxPCBInco", 1, false));
	_compData.push_back(InterpolationData<LinearDataEntry>::loadData(aFolderName, "vpmc.AttenTxFR4Inco", 1, false));

	_folderName = aFolderName;

	bool isValid = true;
	for (size_t k = 0; k < _totData.size(); ++k) {
		if (!_totData[k]->isValid()) {
			ZeusLog("%s: invalid data set %s\n", func, _totData[k]->fileName().c_str());
			isValid = false;
		}
	}
	for (size_t k = 0; k < _compData.size(); ++k) {
		if (!_compData[k]->isValid()) {
			ZeusLog("%s: invalid data set %s\n", func, _compData[k]->fileName().c_str());
			isValid = false;
		}
	}

	int result = 0;
	if (!isValid) {
		unloadData();
		result = 1;
	}

	return result;

}

void HeadAttenuation::unloadData() {
	if (_totData.size() > 0) {
		for (size_t k = 0; k < _totData.size(); ++k) if (_totData[k]) delete _totData[k];
		_totData.clear();
	}
	if (_compData.size() > 0) {
		for (size_t k = 0; k < _compData.size(); ++k) if (_compData[k]) delete _compData[k];
		_compData.clear();
	}
	_folderName = "";
}

vector<InterpolationData<LinearDataEntry>* >  HeadAttenuation::_totData;
vector<InterpolationData<LinearDataEntry>* >  HeadAttenuation::_compData;
string HeadAttenuation::_folderName;
/***********************************************************************************************************/


/********************************************** HeadTransport **********************************************/
class HeadTransport {

public:
	static HeadTransport* getInstance(const char* dataFolder, const char* geometryFileName);

	/*! \brief Deletes the singleton instance */
	static void deleteSingeltonInstance();

	void transportPhoton(PRNG &rng, int aSegment, double E, Vector &x, Vector &u, int weight, int type, double baseWeight, bool isInFOI,
		int &nParticle, Particle aParticleList[]);

	void setSegments(const vector<vector<pair<double, double> > > &segments);

	/*! \brief Get the Compton splitting parameter _nSplit */
	inline int getNsplit() const { return _nSplit; }

protected:
	HeadTransport(const char* dataFolder, const char* geometryFileName);
	~HeadTransport();

	/*! \brief Initializes the photon cross section run time tables */
	int  initializePhotonCrossSections();

	/*! \brief Transport a photon through the simulation geometry

	This function is essentially the same as the public function with the same name, except that now the geometry region and medium
	index are also passed as input.  It is used by the public function to perform the actual transport.
	*/
	void transportPhoton(PRNG &rng, int aSegment, double E, Vector &x, Vector &u, int weight, int type, double baseWeight, int ireg, Medium imed,
		bool isInFOI, double lastMutr, int &nParticle, Particle aParticleList[]);

	bool doSmartCompton(PRNG &rng, int aSegmentIndex, int ireg, Medium imed, double &E, const Vector &x, Vector &u,
		double baseWeight, double lastMutr, int &nParticle, Particle aParticleList[]);

	static HeadTransport*   _singletonHeadTransport; //!< A pointer to the singleton HeadTransport instance

	IGeometry*              _geometry;       //!< The simulation geometry (i.e., the geometry of the treatment head)

	double*               _weightToLambda; //!< log(weight) tabulated for weight = 2..._nSplit

	string                  _dataFolder;     //!< The cross section data folder used for initializations
	string                  _geometryFile;   //!< The file containing geometry specs used for initialization

	double                _pCut;           //!< Photon transport cutoff

	double                _maxWaterMutr;   //!< Maximum energy transfer coeefficient in water (given at E=1.33 MeV)

	double                _xmin;           //!< The minimum of the rectangle to which to send photons in Compton interactions in x-direction
	double                _xmax;           //!< The maximum of the rectangle to which to send photons in Compton interactions in x-direction
	double                _ymin;           //!< The minimum of the rectangle to which to send photons in Compton interactions in y-direction
	double                _ymax;           //!< The maximum of the rectangle to which to send photons in Compton interactions in y-direction

	int                     _exitRegion;     //!< Particles are collected when they enter this region of the geometry
	int                     _nMedia;         //!< Number of materials in the geometry
	int                     _nSplit;         //!< Number of times to split a Compton interaction

	bool                    _isValid;        //!< Is the source valid?

	bool                    _noScatter;      //!< If true, all scatter is ignored

	/*! If true, scattered photons generated in the simulation are subjected to Russian Roullette game with survival probability
	mutr(E)/mutr(1.33), where mutr(E) is the mass energy absorption coefficient in water
	*/
	bool                    _useMutrRejection;

private:

	/*! \brief A private constructor with no arguments (HeadTransport instances can not be created without the location of the cross section files and
	the geometry specification file).
	*/
	HeadTransport() {};

#ifdef USE_OLD_ATTENUATION_TABLES
	/* \brief Provides attenuation values */
	Attenuation _attenuation;
#endif

	/* \brief Base region of the geometry enevelope. */
	BaseRegion _b;
	/* \brief Base region without gradient coils */
	BaseRegionMlc _bNoGradient;

	/* \brief Model of the gantry coils. */
	ViewRayGantryCoils _gc;

	/*! \brief The MLC geometry */
	ViewRayMLC          _mlcGeometry;

	/*! \brief Envelope containing the mlc and the gantry coils.*/
	Envelope *_e;

	/*! \brief Controls the logging in this class. */
	//LogLevel _logLevel;

	inline double calcAttenuation(Medium aMedium, double aPhotonEnergy);

	inline double calcComptonAttenuation(Medium aMedium, double aPhotonEnergy);
};

HeadTransport* HeadTransport::_singletonHeadTransport = 0;

HeadTransport* HeadTransport::getInstance(const char* dataFolder, const char* geometryFileName) {
	//
	// *** If an instance already exists, check to see of it was initialized using a different data folder or a different spec file. 
	//     If yes, delete it.
	//
	if (_singletonHeadTransport) {
		if (_singletonHeadTransport->_dataFolder != dataFolder || _singletonHeadTransport->_geometryFile != geometryFileName) deleteSingeltonInstance();
	}

	//
	// *** If no singleton instance exist, create one
	//
	if (!_singletonHeadTransport) {
		_singletonHeadTransport = new HeadTransport(dataFolder, geometryFileName);
		// Check if the instance is valid and delete it if not.
		if (!_singletonHeadTransport->_isValid) deleteSingeltonInstance();
	}
	return _singletonHeadTransport;
}

void HeadTransport::deleteSingeltonInstance() {
	if (_singletonHeadTransport) {
		delete _singletonHeadTransport;
		_singletonHeadTransport = 0;
	}
}

void HeadTransport::setSegments(const vector<vector<pair<double, double> > > &segments) {
	// Simply call the setSegments() function of the MLC geometry
	_mlcGeometry.setSegments(segments);
}

HeadTransport::~HeadTransport() {
	if (_weightToLambda) delete[] _weightToLambda;
	HeadAttenuation::unloadData();
	delete _e;
}

HeadTransport::HeadTransport(const char *dataFolder, const char *geometryFile) : _weightToLambda(0) {

	const static char* func = "HeadTransport::HeadTransport";

	//
	// *** Mark as not valid yet. 
	//
	_isValid = false;
	_noScatter = false;

	//
	// *** Remember the data folder and spec file used to initialize
	//
	_dataFolder = dataFolder;
	_geometryFile = geometryFile;

	//
	// *** Parse the spec file
	//
	string path = dataFolder;
	path += geometryFile;
	ConfigFile input(path.c_str());

	//
	// *** Construct the treatment head geometry
	//
	string theHead;
	if (!input.getValue("viewray head model", theHead)) exitApp("head model read failed!");

	vector<IGeometry*> inscribedRegions;
	if (theHead.compare("GradientPlusMlc") == 0) {

		_mlcGeometry.setMaxZCoord(_b.getZmax());
		inscribedRegions.push_back(&_mlcGeometry);
		inscribedRegions.push_back(&_gc);

		// Construct the envelope containing the regions
		_e = new Envelope(&_b, inscribedRegions);

		// Correspondence between base and inscribed
		_e->setBaseToInscribedWhere(0, 0);
		_e->setBaseToInscribedWhere(1, 0);
		_e->setBaseToInscribedWhere(2, 1);
		_e->setBaseToInscribedWhere(3, 0);
		_e->setBaseToInscribedWhere(4, 2);

		_e->setBaseToInscribedNext(0, 0);
		_e->setBaseToInscribedNext(1, 1);
		_e->setBaseToInscribedNext(2, 1);
		_e->setBaseToInscribedNext(3, 2);
		_e->setBaseToInscribedNext(4, 2);

	}
	else if (theHead.compare("Mlc") == 0) {

		_mlcGeometry.setMaxZCoord(_bNoGradient.getZmax());
		inscribedRegions.push_back(&_mlcGeometry);

		// Construct the envelope containing the regions
		_e = new Envelope(&_bNoGradient, inscribedRegions);

		// Correspondence between base and inscribed
		_e->setBaseToInscribedWhere(0, 0);
		_e->setBaseToInscribedWhere(1, 0);
		_e->setBaseToInscribedWhere(2, 1);

		_e->setBaseToInscribedNext(0, 0);
		_e->setBaseToInscribedNext(1, 1);
		_e->setBaseToInscribedNext(2, 1);

	}
	else
	{
		exitApp("missing/bad head model specification");
	}

	_geometry = (IGeometry*)_e;

	//
	// *** From which geometry region do particles exit.
	//
	_exitRegion = 3;
	bool success = input.getValue("exit region", _exitRegion);
	if (!success || _exitRegion < 0) {
		exitApp("missing/bad 'exit region' specification");
	}

	//
	// *** Ignore scatter ?
	//
	int ignore = 0; success = input.getValue("ignore scatter", ignore);
	if (success && ignore) _noScatter = true;

	//
	// *** The photon transport cutoff energy
	//
	_pCut = 0.05; //MeV

	//
	// *** Initialize photon cross sections
	//
	if (initializePhotonCrossSections()) {
		exitApp("%s: failed to initialize the cross sections");
	}

	//
	// *** Set the field of interest in the isocenter plane
	//
	_xmin = -17.5;
	_xmax = 17.5;
	_ymin = -17.5;
	_ymax = 17.5;

	//
	// *** Set the Compton splitting number
	//
	_nSplit = 200;
	input.getValue("head nsplit", _nSplit);
	ZeusLog("%s:Compton splitting number is %d", func, _nSplit);

	//
	// *** Initialize log(weight) in the range 2 ... _nSplit. 
	//     Commented out because it did not improve the speed at all 
	//_weightToLambda = new double [_nSplit-2];
	//for(int i=2; i<_nSplit; ++i) _weightToLambda[i-2] = log((double)i);

	//
	// *** ZeusLog that everything went OK and mark the instance as valid
	//
	_isValid = true;
	ZeusLog("%s: all initializations are OK, source is valid.\n", func);

}

int HeadTransport::initializePhotonCrossSections() {
	int result = 0;
	if (HeadAttenuation::loadData(_dataFolder)) result = 1;
	return result;
}

double HeadTransport::calcAttenuation(Medium aMedium, double aPhotonEnergy) {
	if (aMedium < UnknownMaterial) return HeadAttenuation::computeAttenuation(aPhotonEnergy, (int)aMedium);
	else {
		ZeusLog("Error: calcAttenuation: unknownMaterial"); 
		return 0;
	}
}

double HeadTransport::calcComptonAttenuation(Medium aMedium, double aPhotonEnergy) {
	if (aMedium < UnknownMaterial) return HeadAttenuation::computeCompton(aPhotonEnergy, (int)aMedium);
	else {
		ZeusLog("Error: calcComptonAttenuation: unknownMaterial");
		return 0;
	}
}

void HeadTransport::transportPhoton(PRNG &rng, int aSegmentIndex, double E, Vector &x, Vector &u, int weight, int type,
	double baseWeight, bool isInFOI, int &nParticle, Particle aParticleList[]) {

	if (E <= _pCut) return;
	//
	// *** Just in case, check that the energy makes sense
	//
	//
	// *** Where is the photon?
	//
	Medium medium;
	int ireg = _geometry->isWhere(aSegmentIndex, x, medium);

	if (ireg < 0) {
		//
		// *** The photon is outside the simulation geometry. So see if it will enter
		//
		double t = 1e30; ireg = _geometry->nextRegionAlongPath(aSegmentIndex, -1, x, u, t, medium);
		if (ireg < 0) return;  // no, it doesn't enter
		// yes, it does. Update position to the entrance point. Note that the nextRegionAlongPath() method 
		// of the geometry sets the medium index.
		x += u*t;
	}
	else {
		//
		// *** Already inside the geometry. Get the medium index for this region
		//
	}

	//
	// *** Transport the photon
	//
	transportPhoton(rng, aSegmentIndex, E, x, u, weight, type, baseWeight, ireg, medium, isInFOI, 0, nParticle, aParticleList);
}

void HeadTransport::transportPhoton(PRNG &rng, int aSegmentIndex, double E, Vector &x, Vector &u, int weight, int type,
	double baseWeight, int ireg, Medium imed, bool isInFOI, double lastMutr, int &nParticle, Particle aParticleList[]) {

	Medium newmed = imed;

	if (E <= _pCut || E > 1.34) {
		return;
	}

	//
	// *** Remember the photon position 
	//
	Vector xo(x);

	//
	// *** Set up interpolations for this photon and compute attenuation coefficient
	//

#ifdef USE_OLD_ATTENUATION_TABLES
	double mu = _attenuation.CalcAttenuation(imed, E);
#else 
	double mu = calcAttenuation(imed, E);
#endif

	double mui = 1.0 / mu;
	//
	// *** If the weight is greater than 1 and the photon is going towards the FOI, we follow it until it accumulates lambdaMax = log(weight) mean-free-paths. 
	//     If it happens to interact before, we simulate the interaction using doSmartCompton(). When the photon has traveled lambdaMax MFP, 
	//     we attenuate it by exp(-lambdaMax)=1/weight (so its weight becomes 1) and follow it in the usual way
	//
	if (weight > 1 && isInFOI) {

		//
		// *** Compute lambdaMax. Strangely enough, using the tabulated log(weight) made things slightly slower. I guess, a cache issue
		//
		//double lambdaMax = _weightToLambda[weight-2];           // number of MFPs until weight becomes 1
		double lambdaMax = log((double)weight);

		//
		// *** Sample number of MFP until next interaction and remember it
		//
		double lambda = -log(1 - rng());       // MFPs to next interaction
		double lambda0 = lambda;

		//
		// *** See if the photon will interact before reaching lambdaMax
		//
		bool doInteraction = true;
		if (lambda >= lambdaMax) {
			// No, it will not. So, set the number of MFP to betraveled to lambdaMax.
			doInteraction = false; lambda = lambdaMax;
		}

		double totLambda = 0;  // initialize total number of MFP traveled
		int nstep = 0;           // number of steps taken

		//
		// *** Now trace the photon through the geometry
		//
		while (1) {

			// count steps taken
			++nstep;

			//
			// *** Just in case a particle gets stuck at a boundary, check for the number of steps being taken so far not exceeding 1000
			//
			if (nstep > 1000) {
				printf("HeadTransport::transportPhoton: more than 1000 steps taken through geometry\n");
				return;
			}

			// distance left to travel
			double t = lambda*mui;

			// query geometry
			int inew = _geometry->nextRegionAlongPath(aSegmentIndex, ireg, x, u, t, newmed);

			//
			// *** If inew < 0, the photon has gone out of the geometry => simply return. If t < -0.01, the geometry is confused. This may happen on 
			//     extremely rare occasions when the photon enters the MLC from the side. We solve the issue by discarding the photon in this case. 
			//
			if (t < -0.01 || inew < 0) {
				return;
			}

			// update position and MFP's
			x += u*t; totLambda += mu*t; lambda -= t*mu;

			//
			// *** Check to see if the photon has entered the region where we score particles escaping the treatment head?
			//     If the strategy for sampling high weight photons in the source feeding us particles is correct, this should never happen.
			//
			if (inew == _exitRegion && isInFOI) {
				//
				// *** Yes, it did. We tolerate 'misbehaviour' of up to ln(2) MFP's (so that the weight of the photon recorded in 
				//     aParticleList is not more than twice the weight of other photons.
				//      
				if (lambdaMax - totLambda < 0.693147) {
					//
					// *** Store this photon in aParticleList
					//
					Particle p;

					p.x = x._x;
					p.y = -x._z;
					p.z = x._y;
					p.u = u._x;
					p.v = -u._z;
					p.w = u._y;
					p.E = 1e6*E;
					p.weight = baseWeight*exp(lambdaMax - totLambda);
					p.type = photon;
					if (nParticle >= NSampleStack) exitApp("sampling stack overflow. You have to increase the constant NSampleStack");
					else aParticleList[nParticle] = p;
					++nParticle;
				}
				else {
					//
					// *** Warn about this event. If this happens too many times, whoever is using this class has to redesign the way they 
					//     sample high-weight photons
					//
					ZeusLog("A particle with weight %d is entering scoring region before reaching lambdaMax attenuation. segment = %d\n", weight, aSegmentIndex);
					ZeusLog("E=%g x=(%g,%g,%g) u=(%g,%g,%g) xo=(%g,%g,%g) lambda=%g (%g)\n", E, x._x, x._y, x._z, u._x, u._y, u._z, xo._x, xo._y, xo._z, lambdaMax - totLambda, totLambda);
				}
				return;
			}

			//
			// *** If inew == ireg, we have traveled the required number of MFP's. If doInteraction is true, we have arrived at an interaction site. 
			//     In that case we have to simulate the interaction and continue until we reach lambdaMax MFP. 
			//     If doInteraction is false, we have traveled lambdaMax = log(weight) MFP's, so its time to break out of the loop and 
			//     do normal transport below.
			//
			if (inew == ireg) {
				if (doInteraction) {
					//
					// *** I have given the user the option to turn off scatter, so only simulate the interaction if _noScatter is false.
					//
					if (!_noScatter) {
						//
						// *** Further, take into account Compton interaction probability and Russian Roullette game (if weight < _nSplit). 
						//     The probability to do the Compton interaction becomes Pcomton * weight/_nSplit.
						//
						double comptonSigma = calcComptonAttenuation(imed, E);
						double leftSide = rng()*_nSplit*mu;
						if (leftSide < comptonSigma*weight) {
							//
							// *** OK, we need to do Compton scatter. Perform using doSmartCompton().
							//
							double Enew = E; Vector unew(u);
							if (doSmartCompton(rng, aSegmentIndex, ireg, imed, Enew, x, unew, baseWeight, lastMutr, nParticle, aParticleList)) {
								//
								// *** If doSmartCompton() returned true, the scattered photon is not going towards the FOI => transport it.
								//
								Vector xnew(x);
								transportPhoton(rng, aSegmentIndex, Enew, xnew, unew, _nSplit, 3, baseWeight, ireg, imed, false, lastMutr, nParticle, aParticleList);
							}
						}
					}
					//
					// *** Remember that we no longer need to interact and set lambda to the remaning number of MFP until lambdaMax
					//
					doInteraction = false;
					lambda = lambdaMax - lambda0;
				}
				else break; // i.e., we have reached lambdaMax MFP and its time to do regular transport
			}
			else {
				//
				// *** If here, the photon has entered a new region. Check to see if the medium in that region is different and, if yes, 
				//     update the attenuation coefficient.
				//
				if (newmed != imed) {
					// the photon has entered a region with a different material, so compute attenuation coefficient
					imed = newmed;
#ifdef USE_OLD_ATTENUATION_TABLES
					mu = _attenuation.CalcAttenuation(imed, E);
#else 
					mu = calcAttenuation(imed, E);
#endif
					mui = 1 / mu;
				}
				ireg = inew;
			}
		}
		weight = 1; // if here, the photon has been attenuated to a weight of 1.
	}

	//
	// *** Regular transport section
	//
	while (1) {
		// Pick number of MFP until next interaction
		double lambda = -log(1 - rng());
		int nstep = 0;  // as above, count steps just in case.
		while (1) {
			++nstep;
			// avoid particles stuck on boundaries (or otherwise confused) by discarding them.
			if (nstep > 1000) {
				printf("Warning: HeadTransport::transportPhoton(a): more than 1000 steps taken through geometry\n");
				return;
			}
			// distance to travel for lambda MFP
			double t = lambda*mui;
			// query geometry
			int inew = _geometry->nextRegionAlongPath(aSegmentIndex, ireg, x, u, t, newmed);

			// see comments above about inew < 0 or t < 0
			if (t < -0.01 || inew < 0) {
				return;
			}

			if (newmed == UnknownMaterial) exitApp("encounter unknown material??");
			// update position
			x += u*t;
			// has the photon entered the region where we score particles escaping the treatment head?
			if (inew == _exitRegion) {
				if (weight == 1) {
					// yes, it has. Store the particle into the container provided and terminate transport
					Particle p;
					if (E <= 0 || E > 1.34) {
						ZeusLog("Particle with E=%g arrived to be stored?\n", E);
						return;
					}

					p.x = x._x;
					p.y = -x._z;
					p.z = x._y;
					p.u = u._x;
					p.v = -u._z;
					p.w = u._y;
					p.E = 1e6*E;
					p.weight = baseWeight;
					p.type = photon;
					if (nParticle >= NSampleStack) exitApp("sampling stack overflow! you have to increase NSampleStack");
					else aParticleList[nParticle] = p;
					++nParticle;
				}
				return;
			}
			// inew == ireg means we have reached the interaction site => exit the tracking loop
			if (inew == ireg) break;
			// update lambda
			lambda -= t*mu;
			if (newmed != imed) {
				if (newmed == UnknownMaterial) exitApp("encounter unknown material??");
				// the photon has entered a region with a different material, so compute new attenuation coefficient
				imed = newmed;
#ifdef USE_OLD_ATTENUATION_TABLES
				mu = _attenuation.CalcAttenuation(imed, E);
#else 
				mu = calcAttenuation(imed, E);
#endif
				mui = 1 / mu;
			}
			ireg = inew;
		}

		// 
		//  *** The photon has reached an interaction site. 
		// 
		if (_noScatter) return;// i.e., if the user has turned off scatter, simply terminate.

		//
		// We play RussianRoulette with survival probability of weight/_nSplit. Also, if the interaction is 
		// not Compton scattering, we simply terminate the photon history. Thus, probability to survive RR and to do Compton scattering is 
		// pCompton*weight/_nSplit. 
		//
		double comptonSigma = calcComptonAttenuation(imed, E);
		double leftSide = rng()*_nSplit*mu;
		bool keepGoing = leftSide < comptonSigma*weight ?
			doSmartCompton(rng, aSegmentIndex, ireg, imed, E, x, u, baseWeight, lastMutr, nParticle, aParticleList) : false;
		if (!keepGoing) break;

		// 
		// *** If here, the photon was scattered into a direction not going towards the FOI. 
		//     We need to recompute log(E), setup interpolation and continue transport with the new energy.
		//
		weight = _nSplit;
		type = 3;
#ifdef USE_OLD_ATTENUATION_TABLES
		mu = _attenuation.CalcAttenuation(imed, E);
#else 
		mu = calcAttenuation(imed, E);
#endif
		mui = 1 / mu;
	}

}

bool HeadTransport::doSmartCompton(PRNG &rng, int aSegmentIndex, int ireg, Medium imed, double &E, const Vector &x, Vector &u,
	double baseWeight, double lastMutr, int &nParticle, Particle aParticleList[]) {

	//
	// *** Compute minimum and maximum scattering angle that may make the scattered photon move towards the field of interest.
	//     The interested examiner of this piece of code is encouraged to verify the correctness of this calculation.
	//
	double xmin = _xmin, xmax = _xmax;
	double ymin = _ymin, ymax = _ymax;
	double x1, y1, ctmin, ctmax;
	if (u._z) { double t = -x._z / u._z; x1 = x._x + u._x*t; y1 = x._y + u._y*t; }
	else      { x1 = x._x; y1 = x._y; }
	double xx1, yy1, xx2, yy2;
	if (x1 < xmin) { xx1 = xmin; xx2 = xmax; }
	else if (x1 < xmax) { xx1 = x1; xx2 = xmax - x1 > x1 - xmin ? xmax : xmin; }
	else                 { xx1 = xmax; xx2 = xmin; }
	if (y1 < ymin) { yy1 = ymin; yy2 = ymax; }
	else if (y1 < ymax) { yy1 = y1; yy2 = ymax - y1 > y1 - ymin ? ymax : ymin; }
	else                 { yy1 = ymax; yy2 = ymin; }
	Vector v1(xx1, yy1, 0); v1 -= x; v1.normalizeToUnitLength(); double ct1 = v1*u;
	Vector v2(xx2, yy2, 0); v2 -= x; v2.normalizeToUnitLength(); double ct2 = v2*u;
	if (ct1 < ct2) { ctmin = ct1; ctmax = ct2; }
	else            { ctmin = ct2; ctmax = ct1; }

	//
	// *** Compute an estimate of the probability for scatter into the field of interest
	//     Note: this is an overestimate of the true probability
	//
	// - Minimum and maximum possible energy fraction.
	//    eps1 is the minimum possible energy fraction, eps2 the maximum possible. eps1i=1/eps1, eps2i = 1/eps2
	double ko = E*1.9569341; double broi = 1 + 2 * ko, ko2 = ko*ko;
	double eps1i = 1 + ko*(1 - ctmin), eps2i = 1 + ko*(1 - ctmax);
	if (E <= _pCut*eps2i)  return false;  // i.e., the maximum possible energy towards the field of interest is less than the cuttoff, so no need to simulate.
	double eps1 = 1 / eps1i, eps2 = 1 / eps2i;


	// - The integral of Klein-Nishina over all angles and the Klein-Nishina normalization (fnorm)
	double alpha1_t = log(broi);
	double eps1_t = 1 / broi, eps2_t = 1;
	double w2 = alpha1_t*(ko2 - 2 * ko - 2) + (eps2_t - eps1_t)*(2 * broi + 0.5*ko2*(eps1_t + eps2_t));
	double fnorm = w2 / (ko2*ko);

	// - The maximum of Klein-Nishina in the ctmin...ctmax interval.
	double f2 = (eps2i + eps2 - 1 + ctmax*ctmax)*eps2*eps2;
	double f1 = (eps1i + eps1 - 1 + ctmin*ctmin)*eps1*eps1;
	double fmax = f1 > f2 ? f1 : f2;

	// - The overstimate of the probability to have scattered photons going towards the FOI discussed above
	//   (the probability P mentioned above is called wnew here because it is new compared to an earlier overestimate of P used by 
	//    Kawrakow, Walters and Rogers in their paper on DBS [Med. Phys. 31 (2004) 2883]
	//
	double A = (xmax - xmin)*(ymax - ymin);
	double wnew = A / (2 * PI*x._z*x._z*fnorm)*fmax;

	// - Number of interactions to sample: integer of wnew*_nSplit or, with corresponding probability, 1 + integer of wnew*_nSplit.
	double asample = wnew*_nSplit; int nsample = (int)asample; asample = asample - nsample;
	if (rng() < asample) ++nsample;

	//
	// *** Now simulate nsample Compton interactions that go towards the FOI
	//
	for (int i = 0; i < nsample; ++i) {

		//
		// *** Pick a random point within the FOI in the isocenter plane
		//
		Vector uu; uu._z = 0;
		uu._x = xmin + rng()*(xmax - xmin);
		uu._y = ymin + rng()*(ymax - ymin);

		//
		// *** Compute distance between interaction site and sampled point and set direction uu of scattered photon
		//
		uu -= x; double disti = 1 / uu.getLength(); uu *= disti;

		//
		// *** The cosine of the scattering angle and the corresponding scattered photon energy ( E/epsi )
		//
		double cost = uu*u;
		double epsi = 1 + ko*(1 - cost);

		if (E > _pCut*epsi) {          // i.e., only simulate if scattered photon energy is above the cutoff
			double eps = 1 / epsi;
			double aux = disti*x._z;
			// The condition below is a better way of writing the rejection probability condition avoiding expensive divisions
			if (rng()*fmax < (1 + eps*eps - eps*(1 - cost*cost))*eps*aux*aux*aux) {

				// The scattered photon energy
				double Enew = E*eps;

				// Now trace the scattered photon through the geometry
				Vector xx(x);
				transportPhoton(rng, aSegmentIndex, Enew, xx, uu, 1, 3, baseWeight, ireg, imed, true, lastMutr, nParticle, aParticleList);
			}
		}

	}

	//
	// *** Now simulate one more Compton interaction to any direction and follow the scattered photon if it does not go towards the FOI
	//
	// 1. Sample the energy fraction br of the scattered photon. Experiments show that it is faster to sample the scattering angle uniformely 
	//    between the minimum and maximum allowed kinematically up to about 5.5 MeV. The usual method found in EGS, PENELOPE, and alike, 
	//    is only better (faster) above about 5.5 MeV because of the expensive log, exp and sqrt operatiobns needed. 
	//
	double bro = 1 / broi; double br, temp, sint;
	if (ko < 11) {
		double bro1 = 1 - bro;
		double rejmax = ko2*(broi + bro);
		double br2;
		do {
			br = bro + bro1*rng(); br2 = br*br;
		} while (rng()*br2*rejmax > ko2*br*(br2 + 1) - (1 - br)*(br*broi - 1));
		temp = (1 - br) / (ko*br); sint = temp*(2 - temp);
	}
	else {
		double broi2 = broi*broi;
		double alpha1 = log(broi);
		double alpha2 = ko*(broi + 1)*bro*bro;
		double alphaS = alpha1 + alpha2;
		do {
			br = rng()*alphaS < alpha1 ? exp(alpha1*rng())*bro : sqrt(rng()*(broi2 - 1) + 1)*bro;
			temp = (1 - br) / (ko*br); sint = temp*(2 - temp);
		} while (rng()*(1 + br*br) < br*sint);
	}

	//
	// 2. Now compute direction and see if the scattered photon has energy above _pCut and does not go towards the FOI
	//
	bool keepIt = true; E *= br;
	if (E > _pCut) {
		double cost;
		if (sint > 0) { cost = 1 - temp; sint = sqrt(sint); }
		else { cost = -1; sint = 0; }
		double cphi, sphi;
		rng.getRandomAzimuth(cphi, sphi);
// 		double phi = 2 * PI*rng();
// 		cphi = cos(phi);
// 		sphi = sqrt(1 - cphi*cphi);
		u.rotate(cost, sint, cphi, sphi);
		if (u._z < 0) {
			double t = -x._z / u._z; x1 = x._x + u._x*t; y1 = x._y + u._y*t;
			if (x1 > xmin && x1 < xmax && y1 > ymin && y1 < ymax) keepIt = false;
		}
	}
	else keepIt = false;

	return keepIt;
}

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
SEGMENT::SEGMENT(ConfigFile *cf):primAT(NULL), scatAT(NULL), primMask(NULL), scatMask(NULL)
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
	for (int i = 0; i < 3; ++i)	start = vrs.find("\n", start+1);
	sscanf(vrs.c_str()+start, "     Total Beam-On Time (s): %lf", &onTime);
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
int SEGMENT::sample(PRNG &rng, Particle pars[], int is, CoPhaseSpace *ps, HeadTransport* headTransport, MLCattenuation* mlcAttenuation) //is-> total segments index
{
	//obtain how many particles to sample in one call
	int NP = (int)nSample;
	double RSample = nSample - NP;
	if (rng() < RSample) ++NP;
	if (0 == NP) return 0; // i.e., no particles from this call

	double E;
	Vector X, U;
	int type;
	int nPar = 0; //particle number that reaches the FOI
	for (int i = 0; i < NP; ++i) //sample NP particles from phase space file
	{
		bool isPrimary;
		unsigned short *table;
		char *mask;
		// *** select particle type
		if (rng() < pPrim) // picked a primary photon
		{			
			isPrimary = true; table = primAT; mask = primMask;
		}
		else // picked a scattered photon
		{
			isPrimary = false; table = scatAT; mask = scatMask;
		}
		// *** Sample fluence bin from the position table
		int bin = int(rng()*ps->getFluenceBins());
		if (rng()*65535. > table[2 * bin]) bin = table[2 * bin + 1];
		
		int weight = SourceHead::NReject[(int)mask[bin]]; // this is the corresponding weight
		//generate one particle from phase space file
		ps->sample(rng, isPrimary, bin, E, X, U, type);

		// If we have MLC attenuation data initialized, we use it to further reduce the number of particles transported.
		// weight = _Nrej[2] means that the particle is outside of MLC opening plus margins. The rejection probability in this case (1/_Nrej[2])
		// is based on the minimum attenuation for totally closed MLC (about 1/34). But depending on the position, the actual attenuation 
		// may be much larger. We use the MLC attenuation data to play Russian Roulette with the particle if this is the case.
		if (mlcAttenuation && weight == SourceHead::NReject[2])
		{
			// get 1/min. attenuation for this y-position (rounded to integer)
			int newWeight = mlcAttenuation->getInverseRejection(X._y);
			// get the _nSplit parameter from the HeadTransport instance
			int nmax = headTransport->getNsplit();
			// Restrict 1/min. attenuation to 1/_nSplit.
			if (newWeight > nmax) newWeight = nmax;
			// If newWeight > weight, play RR with the particle with survival probability weight/newWeight
			if (newWeight > weight) 
			{
				if (rng()*newWeight > weight) continue; // i.e., particle was killed in the RR game
				// particle survived, so set weight to newWeight
				weight = newWeight;
			}
		}

		// Now transport this photon through the MLC. The transportPhoton() function of the HeadTransport class will record particles 
		// exiting the last gradient coil into aParticleList and will keep track of the number of such particles in nParticle.
		// As a reminder: is comes in as argument, rng is the random number generator, E,X,U,weight,type were sampled above, 
		//                _phaseSpaceData->getParticleWeight() gives the weight of phase space particles, true means that the particle 
		//                is aimed towards the FOI. nPar and pars would be output
		headTransport->transportPhoton(rng, is, E, X, U, weight, type, ps->getParticleWeight(), true, nPar, pars);
	}
	return nPar;
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

int CoPhaseSpace::sample(PRNG &rng, double x, double y, bool isPrimary, double &E, double &dx, double &dy) {

	//
	// *** Check if position is within the rectangle for which data is available
	//
	if (x < _xmin || x > _xmax || y < _ymin || y > _ymax) return 1;

	//
	// *** Convert position to tile indeces
	//
	double sx = 1; if (x < 0) sx = -1;
	double sy = 1; if (y < 0) sy = -1;
	int ix = (int)(x*sx*_dxti0 + _dxti1); if (ix > _Nx1 - 1) ix = _Nx1 - 1;
	int iy = (int)(y*sy*_dyti0 + _dyti1); if (iy > _Ny1 - 1) iy = _Ny1 - 1;
	int j = ix + iy*_Nx1;

	//
	// *** Set sampling tables depending on whether we are sampling a primary or a scattered particle and also sample the energy
	//
	unsigned short *thetaPhiFlu; double minx, miny, delx, dely; int nbin, Nxb;
	if (isPrimary) {
		// Primary particle. 
		// Pick energy
		E = rng() < _pPrim1 ? 1.1732 : 1.3325;
		// Set the tables (and corresponding bin-widths and ranges)
		thetaPhiFlu = _thetaPhiPflu[j];
		minx = _minxP; miny = _minyP; delx = _delxP; dely = _delyP; nbin = _Nxb*_Nyb; Nxb = _Nxb;
	}
	else {
		// Scattered
		// Sample the energy
		unsigned short *spec = _scatSpectrum[j];
		int iE = (int)(rng()*_nE);
		if (rng()*65535. > spec[2 * iE]) iE = spec[2 * iE + 1];
		E = _Emin + _dE*(rng() + iE);
		// set tables (and corresponding bin-widths and ranges)
		double rnno = rng();
		if (rnno < _pScoarse[j]) {
			// pick photon direction from the fine angular grid
			thetaPhiFlu = _thetaPhiSflu[j];
			minx = _minxS; miny = _minyS; delx = _delxS; dely = _delyS; nbin = _Nxb*_Nyb; Nxb = _Nxb;
		}
		else if (rnno < _pScoarse1[j]) {
			// pick photon direction from the less fine angular grid
			thetaPhiFlu = _thetaPhiSflu1[j];
			minx = _minxS1; miny = _minyS1; delx = _delxS1; dely = _delyS1; nbin = _Nxb1*_Nyb1; Nxb = _Nxb1;
		}
		else {
			// pick photon direction from the coarse angular grid
			thetaPhiFlu = _thetaPhiSfluC[j];
			minx = _minxSc; miny = _minySc; delx = _delxSc; dely = _delySc; nbin = _Nxbc*_Nybc; Nxb = _Nxbc;
		}
	}

	//
	// *** Sample the angular bin
	//
	int aBin = (int)(rng()*nbin);
	if (rng()*65535. > thetaPhiFlu[2 * aBin]) aBin = thetaPhiFlu[2 * aBin + 1];
	int iyb = aBin / Nxb; int ixb = aBin - iyb*Nxb;

	//
	// *** Now sample uniformely within the selected angular bin
	//
	dx = minx + (rng() + ixb)*delx;
	dy = miny + (rng() + iyb)*dely;
	sx *= sx; dy *= sy;
	return 0;
}

int CoPhaseSpace::sample(PRNG &rng, bool isPrimary, int bin, double &E, Vector &x, Vector &u, int &type) {
	//
	// *** Convert fluence bin to a random position within the bin and set position vector
	//
	int iyy = bin / _Nx; int ixx = bin - iyy*_Nx;
	x._x = _xmin + _dx*(rng() + ixx);
	x._y = _ymin + _dy*(rng() + iyy);
	x._z = _phspPos;

	//
	// *** Find the fluence tile in which this position falls
	//
	double sx = 1; if (x._x < 0) sx = -1;
	double sy = 1; if (x._y < 0) sy = -1;
	int ix = (int)(x._x*sx*_dxti0 + _dxti1); if (ix > _Nx1 - 1) ix = _Nx1 - 1;
	int iy = (int)(x._y*sy*_dyti0 + _dyti1); if (iy > _Ny1 - 1) iy = _Ny1 - 1;
	int j = ix + iy*_Nx1;

	//
	// *** Set sampling tables depending on whether we are sampling a primary or a scattered particle and also sample the energy
	//
	unsigned short *thetaPhiFlu; double minx, miny, delx, dely; int nbin, Nxb;
	if (isPrimary) {
		// Primary particle. 
		// Pick energy
		if (rng() < _prob1[bin]) {
			E = 1.1732; type = 0;
		}
		else {
			E = 1.3325; type = 1;
		}
		// Set the tables (and corresponding bin-widths and ranges)
		thetaPhiFlu = _thetaPhiPflu[j];
		minx = _minxP; miny = _minyP; delx = _delxP; dely = _delyP; nbin = _Nxb*_Nyb; Nxb = _Nxb;
	}
	else {
		// Scattered
		// Sample the energy
		unsigned short *spec = _scatSpectrum[j];
		int iE = (int)(rng()*_nE);
		if (rng()*65535. > spec[2 * iE]) iE = spec[2 * iE + 1];
		E = _Emin + _dE*(rng() + iE);
		type = 2;
		// set tables (and corresponding bin-widths and ranges)
		double rnno = rng();
		if (rnno < _pScoarse[j]) {
			// pick photon direction from the fine angular grid
			thetaPhiFlu = _thetaPhiSflu[j];
			minx = _minxS; miny = _minyS; delx = _delxS; dely = _delyS; nbin = _Nxb*_Nyb; Nxb = _Nxb;
		}
		else if (rnno < _pScoarse1[j]) {
			// pick photon direction from the less fine angular grid
			thetaPhiFlu = _thetaPhiSflu1[j];
			minx = _minxS1; miny = _minyS1; delx = _delxS1; dely = _delyS1; nbin = _Nxb1*_Nyb1; Nxb = _Nxb1;
		}
		else {
			// pick photon direction from the coarse angular grid
			thetaPhiFlu = _thetaPhiSfluC[j];
			minx = _minxSc; miny = _minySc; delx = _delxSc; dely = _delySc; nbin = _Nxbc*_Nybc; Nxb = _Nxbc;
		}
	}

	//
	// *** Sample the angular bin
	//
	int aBin = (int)(rng()*nbin);
	if (rng()*65535. > thetaPhiFlu[2 * aBin]) aBin = thetaPhiFlu[2 * aBin + 1];
	int iyb = aBin / Nxb; int ixb = aBin - iyb*Nxb;

	//
	// *** Now sample uniformely within the selected bin
	//
	double dx = minx + (rng() + ixb)*delx;
	double dy = miny + (rng() + iyb)*dely;

	//
	// *** Position projected to the isocenter plane
	//
	u = x; u._z -= _sourcePos; u.normalizeToUnitLength();
	double t = -x._z / u._z;
	Vector xx(x + u*t);

	//
	// *** Add angular spread sampled above
	//
	xx._x -= sx*dx; xx._y -= sy*dy;

	//
	// *** Compute direction
	//
	u = xx - x; u.normalizeToUnitLength();

	return 0;
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
HeadTransport* SourceHead::headTransport = NULL;
MLCattenuation* SourceHead::mlcAttenuation=NULL;
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
	phsp=new CoPhaseSpace(File.c_str());

// 	cf->getValue("Geometry file", File);
// 	headTransport = HeadTransport::getInstance(dir.c_str(), File.c_str());
	headTransport = HeadTransport::getInstance(dir.c_str(), "vrSourceDC.py");
	mlcAttenuation = new MLCattenuation((dir + "vrMLCatt.dat").c_str());
	if (mlcAttenuation&&!mlcAttenuation->isValid())
	{
		delete mlcAttenuation;
		mlcAttenuation = NULL;
	}

	if (!cf->getValue("Beams file", File)||!getViewRayBeams(File)) //try to load from ViewRay's beams file
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
	

	//pass all the segments in vector to headTransport
	vector<vector<pair<double, double> > > theSegments;
	size_t ntotSeg = 0;
	for (size_t b = 0; b < beams.size(); ++b) 
	{
		ntotSeg += beams[b].segments.size();
		for (size_t s = 0; s < beams[b].segments.size(); ++s) theSegments.push_back(beams[b].segments[s].leafPositions);
	}
	if (ntotSeg < 1) exitApp("No segments?");
	headTransport->setSegments(theSegments);

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

int SourceHead::sample(PRNG* rng, Particle pars[])
{
	//sample which segments according to their weight of fluence
	int nSeg = (int)segIndex.size();
	int id = nSeg > 1?segSampler.sampleBin((*rng)()):0;//index in total segments
	if (id >= nSeg)
	{
		id = nSeg - 1;
		//for debug purpose
// 		FILE* fp = fopen("bug.txt", "a");
// 		fprintf(fp, "id >= nSeg happened!\n");
// 		fclose(fp);
		//exit(-2);
	}
	int ib = segIndex[id].first;
	int is = segIndex[id].second;
	int np = beams[ib].segments[is].sample(*rng, pars, id, phsp, headTransport,mlcAttenuation); //sample some particles from this segment
	// this beam may rotate a gantry angle

	//rotate the beam
	if (beams[ib].gantryAngle != 0.0)
	{
		double cosPhi = beams[ib].cosPhi;
		double sinPhi = beams[ib].sinPhi;
		double xp;
		for (int i = 0; i < np; ++i)
		{
			xp = cosPhi*pars[i].x - sinPhi*pars[i].y;
			pars[i].y = sinPhi*pars[i].x + cosPhi*pars[i].y;
			pars[i].x = xp;

			xp = cosPhi*pars[i].u - sinPhi*pars[i].v;
			pars[i].v = sinPhi*pars[i].u + cosPhi*pars[i].v;
			pars[i].u = xp;
		}
	}
	if (kE != 1)
	{
		for (int i = 0; i < np; ++i) pars[i].E *= kE;
	}
	

	return np;
}

bool SourceHead::getViewRayBeams(const string& file) // get beam definition from ViewRay's plan
{
	//let's find the file extension
	string ext = file.substr(file.find_last_of("."),string::npos);
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
	headTransport->deleteSingeltonInstance();
}
/***********************************************************************************************************/

