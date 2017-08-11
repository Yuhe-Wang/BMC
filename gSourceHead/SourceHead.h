#ifndef _SOURCEHEAD_H_
#define _SOURCEHEAD_H_

#include "config.h"

#define ZeusLog printf


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

class MLCattenuation;

class CoPhaseSpace {

public:
	friend class CoPhaseSpace4GPU;
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

	string                  _dataFileName;    //!< The file name of the data file used to initialize this instance.

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
	int                     _nE;              //!< Number of energy bins for the scattered spectrum
	int                     _Nx1;             //!< Number of x-bins for the fluence scored differential in energy and angle
	int                     _Ny1;             //!< Number of y-bins for the fluence scored differential in energy and angle
	int                     _Nxb;             //!< Number of x- angular bins for primary photons and part of the scattered photons
	int                     _Nyb;             //!< Number of y- angular bins for primary photons and part of the scattered photons
	int                     _Nxb1;            //!< Number of x- angular bins for part of the scattered photons
	int                     _Nyb1;            //!< Number of y- angular bins for part of the scattered photons
	int                     _Nxbc;            //!< Number of x- angular bins for part of the scattered photons
	int                     _Nybc;            //!< Number of y- angular bins for part of the scattered photons

	bool                    _isValid;         //!< True, if all initializations went OK

};

class SEGMENT
{
public:
	friend class SEGMENT4GPU;
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
	SourceHead() :isoX(0), isoY(0), isoZ(0), prescriptionDose(0), treatmentFraction(1), kE(1){}
	~SourceHead();
	//SourceHead(ConfigFile *cf){ init(cf); }

	bool init(ConfigFile *cf);
	vector<BEAM> beams;
	BinSampler segSampler;
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
	vector<vector<pair<double, double> > > getListofSegments()
	{
		//pass all the segments in vector to headTransport
		vector<vector<pair<double, double> > > theSegments;
		size_t ntotSeg = 0;
		for (size_t b = 0; b < beams.size(); ++b)
		{
			ntotSeg += beams[b].segments.size();
			for (size_t s = 0; s < beams[b].segments.size(); ++s) theSegments.push_back(beams[b].segments[s].leafPositions);
		}
		if (ntotSeg < 1) exitApp("No segments?");
		return theSegments;
	}
	

	bool getViewRayBeams(const string& file);//get beam config from ViewRay's plan file
	bool exportViewRayBeams(const string& file);
	double maxEnergy() { return kE*1.35e6; } //unit eV
	static CoPhaseSpace* getPhaseSpace(){ return phsp; }

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
	static CoPhaseSpace *phsp;//phase space, only need one instance now

	vector<pair<int, int> > segIndex; //first = beam index; second = segment index
	double isoX, isoY, isoZ;
	double prescriptionDose;
	int treatmentFraction;

	//experiment variable
	double kE; // energy scaling factor
};
#endif