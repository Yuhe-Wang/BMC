#pragma once

#include "zeusConfig.h"

#include <string>
#include <vector>
using namespace std;

#include "zeusExceptions.h"
#include "zeusConstants.h"
#include "zeusLogProviders.h"

namespace Zeus { namespace Random { 
    class RandomGenerator;
}}

namespace Zeus { 

    class VoxelGeometry;

    /**
    * \brief Namespace containing everything related to loading, unloading, and interpolating cross section data.
    *
    * Each cross section is represented by a public class, containing a static attribute class, 
    * with singleton behaviour, which loads and unloads the data.  The data can be loaded once
    * and then many simulations can be done.  You dont ever have to unload the data.
    * 
    * The public class has methods to allow the simulation to query attributes of the cross section 
    * data, and to interpolate between the tables.
    *
    */
    namespace CrossSectionData {

    /*! \brief A structure for the 4 spline coefficients of an interpolation bin */
    struct SplineDataEntry {
        /*! \name Spline coefficients */
        ///@{
        double _a, _b, _c, _d;
        ///@}
        /*! \brief Default constructor */
        SplineDataEntry() : _a(0), _b(0), _c(0), _d(0) {};
        /*! \brief Construct from input coefficients \a a, \a b, \a c and \a d. */
        SplineDataEntry(double a, double b, double c, double d) : _a(a), _b(b), _c(c), _d(d) {};
        /*! \brief Compute interpolation for \a x. */
        inline double compute(double x) const {
            double x2 = x*x;
            return _a + _b*x + _c*x2 + _d*x*x2;
        };
        /*! \brief Prepare spline coefficients from the input data \a xval, \a fval with \a ndat grid points. */
        static SplineDataEntry* prepareData(const double *xval, const double *fval, int ndat);
        /*! \brief Write coefficients and \a fileName using  LogProviders::InformationLog(). */
        void logData(const char *fileName) const;
    };

    /*! \brief A structure for the 2 linear interpolation coefficients of an interpolation bin */
    struct LinearDataEntry {
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

    /**
    * \brief A template class for managing interpolation data.
    */
    template <class DataEntry> class ZEUS_EXPORT InterpolationData {
    public:

        /*! \brief Constructor for a not initialized instance. */
        InterpolationData( ) : _ndata(0), _data(0), _bIndex(0), _aIndex(0), _nmat(0), _isValid(false) { };

        /**
        * \brief Constructor.
        *
        * Read in the cross section data.
        * \param aFileName Which file is the cross section data in.
        * \param aNumMaterials Cross section data for how many materials is in the file.
        * \param useInverseEnergy If set to true, data will be interpolated as a function of 1/E instead of E.
        */
        InterpolationData( const string &aFileName, const int aNumMaterials, bool useInverseEnergy = false );

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
        InterpolationData( const string &aName, const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &aData);
        /**
        * Is the instance valid?
        */
        bool isValid( ) const { return _isValid; };
        /**
        * Write all the cross section data to the logger.
        */
        void dumpData( ) const;
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
            int i = (int) ( _aIndex[matid] * (x - _bIndex[matid]) );
            if ( i < 0 || i >= _ndata[matid] ) {
                Zeus::LogProviders::FatalLog("Spline interpolation error for data %s: index %d out of range. x=%g matid=%d\n",_fileName.c_str(),i,x,matid);
                throw( new Zeus::Exceptions::ComputeException( _fileName, matid, i, 0, _ndata[matid] ) );
            }
            return _data[matid][i].compute(x);
        };
        /**
        * Interpolate, returning zero if x below range 
        */
        inline double interpolate1(double x, int matid) const {
            double dx = x - _bIndex[matid];
            double result = 1e30;
            if( dx > 0 ) {
                int i = (int) ( _aIndex[matid] * dx );
                if ( i >= _ndata[matid] ) {
                    Zeus::LogProviders::FatalLog("Spline interpolation error for data %s: index %d out of range. x=%g matid=%d\n",_fileName.c_str(),i,x,matid);
                    throw( new Zeus::Exceptions::ComputeException( _fileName, matid, i, 0, _ndata[matid] ) );
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
        static InterpolationData* loadData(const string &aDataFolder, const char *aFileName, int aNumMaterials, bool useInverseEnergy = false );

		void outputSpline(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d)
		{
			int n = _ndata[matid];
			a.resize(n);
			b.resize(n);
			c.resize(n);
			d.resize(n);
			aIndex = _aIndex[matid];
			bIndex = _bIndex[matid];
			for (int i = 0; i < n; ++i)
			{
				a[i] = _data[matid][i]._a;
				b[i] = _data[matid][i]._b;
				c[i] = _data[matid][i]._c;
				d[i] = _data[matid][i]._d;
			}
			
		}

		void outputLinear(int matid, double& aIndex, double& bIndex, vector<double>& a, vector<double>& b)
		{
			int n = _ndata[matid];
			a.resize(n);
			b.resize(n);
			aIndex = _aIndex[matid];
			bIndex = _bIndex[matid];
			for (int i = 0; i < n; ++i)
			{
				a[i] = _data[matid][i]._a;
				b[i] = _data[matid][i]._b;
			}

		}
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

    /**
    * \brief Loads and unloads the material data.
    */
    class ZEUS_EXPORT MaterialData {
    public:

        /**
        * Constructor.  This will read in the material data from the disk.
        */
        MaterialData( const string & aFileName );

        /**
        * Write all the cross section data to the logger.
        */
        void dumpData( );

        /**
        * Destructor.  This will delete all the material data that was read in 
        * by the constructor.
        */
        ~MaterialData( );

        /**
        * Is the class valid?
        */
        bool isValid( ) const { return _isValid; };

        /**
        * How many materials are there?
        */
        int numMaterials();

        double _eabs;        //!< Electron absorption threshold
        double _eabsph;      //!< Photon absoprtion threshold
        double _eminph;      //!< Minimum energy in the photon tabulations
        double _emax;        //!< Maximum energy of all tabulations

        double *_zmass;      //!< Effective Z/A for all materials

        double _eloss;       //!< Step length parameter

        double _eRangeCut;   //!< Lowest energy for range calculations
        
        double _wcion;       //!< Production threshold for inelastic collisions.
        double _wcbre;       //!< Production threshold for bremsstrahlung */

        
        int  _csda;          //!< Flag for CSDA calculation 

        /** @name Unused
         *  Inherited from the initial DPM implementation. Was used in subeabs, the function performing sub-threshold transport. 
         */
        ///@{
        double _substp, _subfac, _subden;
        ///@}
   
    private:
        bool _isValid;        //!< Is the data valid
        int _nmat;            //!< Number of materials
        double _emin;         //!< Minimum energy of the data
        double _refden;       //!< Density of the reference material
        double *_matden;      //!< Densities for all materials
        string _fileName;     //!< Name of the file from which the data was loaded
    };

    /**
    * \brief Provides access to material data.
    */
    class ZEUS_EXPORT Material {
    public:
        /*! \brief Constructor. Does nothing. */
        Material( );
        /**
        * Is the class valid?
        */
        bool isValid() const { return _isValid; };
        /**
        * \brief Load from the material cross section file.
        *
        * This load will determine the number of materials stored in the material file.
        * Every other cross section file must have cross sections for this number of materials.
        * \param aFolderName Path to the cross section file.
        * \return true if cross section data loaded successfully.
        */
        static int loadData( const string & aFolderName );
        /**
        * Write all the cross section data to the logger.
        */
        void dumpData( );
        /**
        * Destructor.  
        * All the data is deleted, so it will have to be loaded again.
        */
        static void unloadData( );

        /**
        * Number of materials stored in this class.
        */
        static int numMaterials();
        /*! \brief Returns bremsstrahlung production threshold. */
        double wcbre() const { return _matData->_wcbre; };
        /*! \brief Returns electron absorption threshold */
        double getEabs() const { return _matData->_eabs; };
        /*! \brief Returns photon absorption threshold */
        double getEabsph() const { return _matData->_eabsph; };
        /*! \brief Returns the maximum energy of the tabulations */
        double getEmax() const { return _matData->_emax; };
        /*! \brief Returns the minimum energy of the photon tabulations */
        double getEminph() const { return _matData->_eminph; };
        /*! \brief Returns the step length parameter (which is maximum energy loss per step) */
        double getEloss() const { return _matData->_eloss; };
        /*! \brief Returns the minimum energy of the range tables */
        double getErangeCut() const { return _matData->_eRangeCut; };
        /*! \brief Returns the array of effective Z/A for all materials. */
        const double* getZmass() const { return _matData->_zmass; };
        /*! \brief Returns the electron production threshold */
        double getWcion() const { return _matData->_wcion; };
        /*! \brief Returns true, if the simulation is to be done in the CSDA, false otherwise. */
        int isCSDA() const { return _matData->_csda; };
        /**
        * Destructor.
        */
        ~Material( );

    private:
        static MaterialData *_matData;   //!< The static material data instance
        bool _isValid;                   //!< True, if everythink went OK.
    };


    /**
    * \brief Loads and unloads Bremsstrahling data.
    */
    class ZEUS_EXPORT BremmConstantsData {
    public:
        /**
        * Initialize the class attributes with cross section data.
        * No cross section data is loaded here, as it is all hardcoded inside this class. 
        */
		BremmConstantsData(const string & aFolderName, const int aNumMaterials);
        /**
        * Write all the cross section data to the logger.
        */
        void dumpData( );
        /**
        * Is the class valid?
        */
        bool isValid( ) const { return _isValid; };
        /**
        * Returns and array of coefficients, one per material.
        */
        double *bcb();
        /** 
        * Returns a scalar value for that material.
        */
        double *fo( int materialIndex );
        /**
        * Destructor.
        */
        ~BremmConstantsData();
    private:
        /*! Is the instance valid? */
        bool _isValid; 
        /*! Number of materials */
        int _nmat;
        /*! bremsstrahlung constants for all materials */
        double **_f0, *_bcb;
        /*! The name of the file from which the data was loaded. */
        string _fileName;
    };

    /**
    * \brief Provides access to Bremsstrahling data.
    */
    class ZEUS_EXPORT BremmConstants {
    public:
        /**
        * Constructor.
        */
		BremmConstants();
        /**
        * Is the class valid?
        */
        bool isValid() const { return _isValid; };
        /**
        * Initialize the class holding the cross section parameters.
        * \return true if cross section data loaded successfully.
        */
		static int loadData(const string & aFolderName);
        /**
        * Unload the cross section data.
        */ 
        static void unloadData();
        /**
        * Write all the cross section data to the logger.
        */
        void dumpData( );

        /*! \brief Samples a bremsstrahlung photon energy.
         *
         * \param[in]  e      Energy of the incident electron in eV.
         * \param[in]  matid  Material ID where the interaction occurs.
         * \param[in]  aRng   A pointer to a random number generator.
         * \returns The sampled photon energy.
         */
        double sambre ( double e, int matid, Random::RandomGenerator *aRng );
        /**
        * Destructor.
        */
        ~BremmConstants();
        /**
        * Public method to call this Schiff function, for testing 
        */
        void schiffPublic ( double b, double *f1, double *f2 );

		static BremmConstantsData* getData() { return _bremmData; }
    private:
        /*! \brief Computes the Schiff function */
        void schiff ( double b, double *f1, double *f2 );
        /*! \brief Pointer to the static bremsstrahlung data */
        static BremmConstantsData *_bremmData;
        /*! \brief Pointer to the material data */
        Material *_material;
        /*! \brief Is this instance valid? */
        bool _isValid;
    };

    //class RunTimeTables {
    //public:
    //    RunTimeTables() { };
    //    ~RunTimeTables( ) { };
    //    static bool isValid();
    //    static int loadData( const string & aFolderName, const int aNumMaterials );
    //    static void dumpData( );
    //    static void unloadData( );

    //    inline static double computeKerma(double e, int matid) { return _photonKerma ? _photonKerma->interpolate(e,matid) : 0; };
    //    inline static double computeLambdaPair(double e, int matid) { return _photonPair && e > _photonPair->lowerBound(matid) ? 
    //        _photonPair->interpolate(e,matid) : 0; 
    //    };
    //private:
    //    static InterpolationData<LinearDataEntry>*    _photonPair;
    //    static InterpolationData<LinearDataEntry>*    _photonCompton;
    //    static InterpolationData<LinearDataEntry>*    _photonSigma;
    //    static InterpolationData<LinearDataEntry>*    _photonKerma;
    //    static InterpolationData<LinearDataEntry>*    _electronLambda;
    //    static InterpolationData<LinearDataEntry>*    _electronInvLambda;
    //    static InterpolationData<LinearDataEntry>*    _electronBremsProbability;
    //    // etc.
    //};


    /**
    * \brief Provides access to photon kerma as a function of energy
    */
    class ZEUS_EXPORT Kerma {
    public:
        /*! \brief Constructor. Does nothing. */
        Kerma() {};
        /*! \brief Destructor. Does nothing. */
        ~Kerma() {};
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Returns the Kerma of a photon with energy \a e in material \a matid. */
        inline static double compute(double e, int matid) { return _data ? _data->interpolate(e,matid) : 0; };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to the disk. */
        static void dumpData( );
        /*! \brief Unload the data */
        static void unloadData( );
    private:
        /*! \brief Pointer to the static linear interpolation table. */
        static InterpolationData<LinearDataEntry>* _data;
    };

    /**
    * \brief Provides access to the total number of MFP for electron discrete interactions as a function of energy
    *        i.e., Integral [ (SigmaMoller(E') + SigmaBrems(E')) dE'/ RestrictedStoppingPower(E') , {E' from Eabs to E} ]
    */
    class ZEUS_EXPORT LambdaTot {
    public:
        /*! \brief Constructor. Does nothing. */
        LambdaTot( ) { };
        /*! \brief Destructor. Does nothing. */
        ~LambdaTot( ) { };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Prepare the interpolation data. 
         *
         * \param[in]  emin   Minimum energies of the tabulations to be prepared for all materials.
         * \param[in]  emax   Maximum energies of the tabulations to be prepared for all materials.
         * \param[in]  nEnergy Number of energy bins to be used. 
         * \param[in]  data    Array of pointers to MFP data for all materials.
         */
        static void prepareData(const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &data);
        /*! \brief Returns the total number of MFP's for discrete interactions for an electron of energy \a e in material \a matid */
        inline static double lmbdatot(double e, int matid) { return _data ? _data->interpolate(e,matid) : 0; };
        /*! \brief Unload the data. */
        static void unloadData( );
    private:
        /*! \brief Pointer to the static linear interpolation table. */
        static InterpolationData<LinearDataEntry>* _data;
    };

    /**
    * \brief Provides access to the inverse LambdaTot.
    *
    * That is, the solution of 
    *      lambda = Integral [ (SigmaMoller(E') + SigmaBrems(E')) dE'/ RestrictedStoppingPower(E') , {E' from Eabs to E} ]
    * for E as a function of lambda.
    */
    class ZEUS_EXPORT InvLambdaTot {
    public:
        /*! \brief Constructor, does nothing */
        InvLambdaTot( ) { };
        /*! \brief Destructor, does nothing */
        ~InvLambdaTot( ) { };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Prepare the interpolation data. 
         *
         * \param[in]  emin   Minimum energies of the tabulations to be prepared for all materials.
         * \param[in]  emax   Maximum energies of the tabulations to be prepared for all materials.
         * \param[in]  nEnergy Number of energy bins to be used. 
         * \param[in]  data    Array of pointers to inverse MFP data for all materials.
         */
        static void prepareData(const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &data);
        /*! \brief Returns the energy correpondingto \a lambda MFP's in the material \a matid. */
        inline static double energy(double lambda, int matid) { return _data ? _data->interpolate(lambda,matid) : 0; };
        /*! \brief Unload the data. */
        static void unloadData( );
    private:
        /*! \brief Pointer to the static linear interpolation table. */
        static InterpolationData<LinearDataEntry>* _data;
    };

    /**
    * \brief Provides access to the bremsstrahlung probability 
    */
    class ZEUS_EXPORT BremsProbability {
    public:
        /*! \brief Constructor, does nothing */
        BremsProbability() {};
        /*! \brief Destructor, does nothing */
        ~BremsProbability() {};
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Prepare the interpolation data. 
         *
         * \param[in]  emin   Minimum energies of the tabulations to be prepared for all materials.
         * \param[in]  emax   Maximum energies of the tabulations to be prepared for all materials.
         * \param[in]  nEnergy Number of energy bins to be used. 
         * \param[in]  data    Array of pointers to bremsstrahlung probability data for all materials.
         */
        static void prepareData(const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &data);
        /*! \brief Returns the probability for a bremsstrahlung interaction (as oposed to Moller) for electron with energy \a e in material \a matid */
        inline static double pbrem(double e, int matid) { return _data ? _data->interpolate(e,matid) : 0; };
        /*! \brief Unload the data. */
        static void unloadData( );
    private:
        /*! \brief Pointer to the static linear interpolation table. */
        static InterpolationData<LinearDataEntry>* _data;
    };

    /**
    * \brief Provides access to the Bremsstrahling cross section data.
    */
    class ZEUS_EXPORT LambdaBr {
    public:
        /*! \brief Constructor. Does nothing */
        LambdaBr( ) { };
        /*! \brief Destructor. Does nothing */
        ~LambdaBr( ) { };
        /*! \brief Returns the bremsstrahlung cross section (in cm^2/g) for an incident electron with energy \a e in material \a matid. */
        inline static double lmbdabr(double e, int matid) { return _data ? _data->interpolate(e,matid) : 1e-30; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );
    private:
        /*! \brief Pointer to the static linear interpolation table */
        static InterpolationData<LinearDataEntry>* _data;
    };

    /**
    * \brief Provides access to the Compton cross section.
    */
    class ZEUS_EXPORT LambdaCo {
    public:
        /*! \brief Constructor. Does nothing */
        LambdaCo( ) { };
        /*! \brief Destructor. Does nothing */
        ~LambdaCo( ) { };
        /*! \brief Returns the Compton cross section (in cm^2/g) for an incident photon with energy \a e in material \a matid. */
        inline static double lmbdaco(double e, int matid) { return _data ? _data->interpolate(e,matid) : 0; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );

		static InterpolationData<SplineDataEntry>* getData() { return _data; }
    private:
        /*! \brief Pointer to the static linear interpolation table */
        static InterpolationData<SplineDataEntry>* _data;
    };

    /**
    * \brief Provides access to the Moller cross section. 
    */
    class ZEUS_EXPORT LambdaMoller {
    public:
        /*! \brief Constructor. Does nothing */
        LambdaMoller( ) {}
        /*! \brief Destructor. Does nothing */
        ~LambdaMoller( ) {};
        /*! \brief Returns the Moller cross section (in cm^2/g) for an incident electron with energy \a e in material \a matid. */
        inline static double lmbdamo(double e, int matid) { return _data ? _data->interpolate1(-1/e,matid) : 1e30; };
        /*! \brief Same as above, but now the input argument is -1/e */
        inline static double lmbdamoFast(double ie, int matid) { return _data ? _data->interpolate1(ie,matid) : 1e30; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );
    private:
        /*! \brief Pointer to the static spline interpolation table */
        static InterpolationData<SplineDataEntry>* _data;
    };


    /**
    * \brief Provides access to the pair production cross section.
    */
    class ZEUS_EXPORT LambdaPair {
    public:
        /*! \brief Constructor. Does nothing */
        LambdaPair( ) { };
        /*! \brief Destructor. Does nothing */
        ~LambdaPair( ) { };
        /*! \brief Returns the pair production cross section (in cm^2/g) for an incident photon with energy \a e in material \a matid. */
        inline static double lmbdapp(double e, int matid) { return _data && e > _data->lowerBound(matid) ? _data->interpolate(e,matid) : 0; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );

		static InterpolationData<SplineDataEntry>* getData() { return _data; }
    private:
        /*! \brief Pointer to the static spline interpolation table */
        static InterpolationData<SplineDataEntry>* _data;
    };

    /**
    * \brief Provides access to the total photon cross section, including its maximum across materials (needed for the Woodcock trick).
    */
    class ZEUS_EXPORT LambdaPhoton {
    public:
        /*! \brief Constructor. SetupWoodcock() should be called before the instance can be used in the simulation. */
        LambdaPhoton( ) : _wcka(0), _wckb(0), _wckLen(0) { };
        /*! \brief Destructor. If needed, deletes the Woodcock interpolation table. */
        ~LambdaPhoton( ) { 
            if( _wckLen > 0 ) {
                delete [] _wcka; delete [] _wckb;
            }
        };
        /*! \brief Returns the total cross section (in cm^2/g) for an incident photon with energy \a e in material \a matid. */
        inline static double lmbdaph(double e, int matid) { return _data ? _data->interpolate(e,matid) : 0; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );
        /*! \brief Returns the maximum of the photon cross section (in 1/cm) from all materials and mass densities in the simulation */
        double lmbdawck ( double aEnergy ) const {
            int i = (int) ( _fwck * (aEnergy - _ewck0) );
            if ( i < 0 || i >= _wckLen ) {
                Zeus::LogProviders::FatalLog("lmbdawck error: index out of range for energy %g, _fwck=%g _ewck0=%g",aEnergy,_fwck,_ewck0);
                throw( new Zeus::Exceptions::ComputeException( _data->fileName(), -1, i, 0, _wckLen ) );
            }
            return _wcka[i] + _wckb[i] * aEnergy;
        };
        /*! \brief Prepares the interpolation tables needed for the Woodcock trick and computed by the lmbdawck() function. */
		void SetupWoodcock(int aNumElements, const float *dens, const char *mat, int nvoxels);
		void OutputWoodcock(double& aIndex, double& bIndex, vector<double>& a, vector<double>& b)
		{
			aIndex = _fwck;
			bIndex = _ewck0;
			a.resize(_wckLen);
			b.resize(_wckLen);
			for (int i = 0; i < _wckLen; ++i)
			{
				a[i] = _wcka[i];
				b[i] = _wckb[i];
			}
		}
		static InterpolationData<SplineDataEntry>* getData() { return _data; }
    private:
        /*! \brief Pointer to the static spline interpolation table */
        static InterpolationData<SplineDataEntry>* _data;
        double _ewck0;     //!< Converts energy to interpolation bin
        double _fwck;      //!< Converts energy to interpolation bin
        double *_wcka;     //!< First Woodcock interpolation coefficient.
        double *_wckb;     //!< Second Woodcock interpolation coefficient.
        int _wckLen;       //!< Length of the Woodcock tables.
    };

    /**
    * \brief Provides access to the electron range
    */
    class ZEUS_EXPORT Range {
    public:
        /*! \brief Constructor, does nothing */
        Range ( ) {}
        /*! \brief Destructor, does nothing */
        ~Range ( ) {};
        /*! \brief Returns the range (in g/cm^2) of an electron with energy \a e in material \a matid. */
        inline static double range(double e, int matid) { return _data ? _data->interpolate(e,matid) : 0; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );

		static InterpolationData<LinearDataEntry>* getData() { return _data; }
    private:
        /*! \brief Pointer to the static linear interpolation table */
        static InterpolationData<LinearDataEntry>* _data;
    };

    /**
    * \brief Provides access to the inverse electron range (i.e., energy as a function of range)
    */
    class ZEUS_EXPORT InverseRange {
    public:
        /*! \brief Constructor, does nothing */
        InverseRange ( ) {}
        /*! \brief Destructor, does nothing */
        ~InverseRange ( ) {};
        /*! \brief Returns the energy of the electron that has range \a range (in g/cm^2) in material \a matid */
        inline static double energy(double range, int matid) { return _data ? _data->interpolate(range,matid) : 0; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName, assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );

		static InterpolationData<LinearDataEntry>* getData() { return _data; }
    private:
        /*! \brief Pointer to the static linear interpolation table */
        static InterpolationData<LinearDataEntry>* _data;
    };

    /*! \brief Provides cross section data for the materials in the ViewRay treatment head (MLC and below) */
    class HeadAttenuation {
    public:
        /*! \brief Constructor, does nothing. Data has to be loaded using loadData() */
        HeadAttenuation( ) { };
        /*! \brief Destructor, does nothing. Data has to be unloaded using unloadData() */
        ~HeadAttenuation( ) { };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_totData.size() > 0 && _compData.size() > 0 && _totData.size() == _compData.size()); };
        /*! \brief Load the data from the files in the folder \a aFolderName. */
        static int loadData( const string & aFolderName );
        /*! \brief Unload the cross section data. */
        static void unloadData( );
        /*! \brief Returns the total attenuation coefficient for a photon with energy \a e in material \a matid. */
        inline static double computeAttenuation(double e, int matid) { return _totData[matid]->interpolate(e,0); };
        /*! \brief Returns the Compton attenuation coefficient for a photon with energy \a e in material \a matid. */
        inline static double computeCompton(double e, int matid) { return _compData[matid]->interpolate(e,0); };
    private:
        static vector<InterpolationData<LinearDataEntry>* >  _totData;    //!< The total attenuation data
        static vector<InterpolationData<LinearDataEntry>* >  _compData;   //!< The Compton attenuation data
        static string                                        _folderName; //!< The folder from which the data was loaded.
    };

    /**
    * \brief Loads and unloads the q surface data describing multiple elastic scattering
    *
    *  The data has been calculated using a modified version of the PREDPM code. 
    *  The data is tabulated as a function of -1/E, where E is the electron energy.
    *  In the original version of PREDPM no checks are being made for the step length being longer than the range for this energy. 
    *  Also, it seems that the DPM authors have switched from using the Lewis theory (where elastic scattering moments 
    *  are integrals over the path length, which they said they used in the 2000 PMB paper) to the Goudsmith-Saunderson (GS) theory, 
    *  where moments are computed for the energy and multiplied with 
    *  the step length. In the modified version, the step length is given by range(E) - range(E - dE), where dE is the maximum 
    *  energy loss allowed for electrons with energy E (so dE is either the energy loss parameter, or E itself, if the energy loss 
    *  parameter is the greater than the energy). To compute the elastic scattering momentsi for the GS, the average energy of the step is used. 
    *  The Q-surface is then prepared with the unmodified libelastic.f files from PREDPM, which uses the 
    *  approach described by Kawrakow & Bielajew in NIM B 134 (1998), 325-336.
    */
    class ZEUS_EXPORT QSurfaceData {
    public:
        /**
        * Load the q surface using the old ascii format 
        */
        void asciiLoad( const string & aFileName );
        /**
        * Load the q surface using the binary format
        * Note: this function has been disabled 
        */
        void binaryLoad( const string & aFileName );
        /**
        * Read the qsurface data.  
        * \param aBinaryLoad True to load from binary file, false to load from ascii file.
        * \param aFileName Filename of cross section data.
        * \param aNumMaterials Number of materials in the cross section file.
        */
        QSurfaceData( const bool aBinaryLoad, const string & aFileName, const int aNumMaterials );
        /**
        * Is the class valid?
        */
        bool isValid( ) const { return _isValid; };
        /**
        * Destructor.  
        * All the data is deleted, so it will have to be loaded again.
        */
        ~QSurfaceData();
        /**
        * Write all the cross section data to the logger.
        */
        void dumpData( );
        int *_nuq;     //!< Number of u-values for each material
        int *_ne;      //!< Number of energies for each material
        double *_fq;   //!< First coefficient converging energy to energy bin for each material
        double *_ieq0; //!< Second coefficient converging energy to energy bin for each material
        /** \brief The q-surface data
        *
        * Three dimensional.  Each plane is for different materials, then within
        * each plane, one dimension is energy and one dimension is the angular variable u. 
        */
        double ***_qs;

        /**
        * The maximum of the q surface: _Qmax[matid][ie] is the maximum of q for material matid and energy index ie
        */
        double **_Qmax;

        /**
        * Write the q surface data in a binary chunk.  Method used to create the binary chunk, after loading the ascii
        * file.
        */
        void dumpBinary( const string & aFileName );

		void Output(int matid, double& fq, double& ieq0, vector<vector<double>>& qs)
		{
			fq = _fq[matid];
			ieq0 = _ieq0[matid];
			qs.resize(_ne[matid]);
			for (int i = 0; i < _ne[matid]; ++i)
			{
				qs[i].resize(_nuq[matid]);
				for (int j = 0; j < _nuq[matid]; ++j) qs[i][j] = _qs[matid][i][j];
			}
		}
    private:
        bool _isValid;     //!< Is the Q-surface data valid
        int  _nmat;        //!< Number of materials
        string _fileName;  //!< The file name from which the data was loaded
    };

    /**
    * \brief Provides access to the q surface data describing multiple elastic scattering
    */
    class ZEUS_EXPORT QSurface {
    public:
        /**
        * Interpolates the q surface table. 
        * \param u Angular variable.
        * \param aEnergy Which energy should we interpolate.
        * \param aMaterialId Use cross section tables for which material.
        * \return Interpolated value.
        * This function is depreciated. It was used in the original multiple scattering angle sampling implementation.
        */
        static double qsurf ( double u, double aEnergy, int aMaterialId );

        /*! \brief Set-up energy interpolation of the Q-surface.
         *
         * Given an energy \a aEnergy and a material \a aMaterialId, returns in \a energyIndex the energy index 
         * and in \a probability the probability to use \a energyIndex+1 instead of \a energyIndex.
         */
        static void getEnergyIndex(double aEnergy, int aMaterialId, int &energyIndex, double &probability) {
            double ie = -1/aEnergy;
            double rle = _qSurfaceData->_fq[aMaterialId] * (ie - _qSurfaceData->_ieq0[aMaterialId]);
            if( rle > 0 ) {
                energyIndex = (int) rle;
                probability = rle - energyIndex;
                if( energyIndex >= _qSurfaceData->_ne[aMaterialId] - 1 ) {
                    energyIndex = _qSurfaceData->_ne[aMaterialId] - 1; probability = -1;
                }
            }
            else {
                energyIndex = 0; probability = -1;
            }
        };

        /*! \brief Same as the above function, but now the energy argument is -1/energy instead of energy. */
        static void getEnergyIndex1(double ie, int aMaterialId, int &energyIndex, double &probability) {
            double rle = _qSurfaceData->_fq[aMaterialId] * (ie - _qSurfaceData->_ieq0[aMaterialId]);
            if( rle > 0 ) {
                energyIndex = (int) rle;
                probability = rle - energyIndex;
                if( energyIndex >= _qSurfaceData->_ne[aMaterialId] - 1 ) {
                    energyIndex = _qSurfaceData->_ne[aMaterialId] - 1; probability = -1;
                }
            }
            else {
                energyIndex = 0; probability = -1;
            }
        };

        /*! \brief Interpolates the Q-surface.
         *
         * Unlike the original qsurf(double,double,int) version, here interpolation is only in the u variable. 
         * This is used in the current sampling of multiple scattering angles, where we first set the energy 
         * bin by using the index and probability returned by getEnergyIndex() or getEnergyIndex1() and then 
         * use this much faster version to compute the q-surface value in the sampling loop.
         */
        static double qsurf ( int anEnergyIndex, int aMaterialId, double u ) {
            int nuq = _qSurfaceData->_nuq[aMaterialId] - 1;
            double ru = u * nuq;
            int ju = (int) ru; if( ju > nuq-1 ) ju = nuq-1;
            ru -= ju; const double *q = _qSurfaceData->_qs[aMaterialId][anEnergyIndex];
            return q[ju]*(1-ru) + q[ju+1]*ru;
        };

        /**
        * \brief Load the data from the cross section file.
        * \param aBinaryLoad True to load this file from a binary file.  False will load from the previous ascii file.
        * \param aFolderName Path to the cross section file.
        * \param aNumMaterials How many materials are in the cross section file.
        * \return true if cross section data loaded successfully.
        */
        static int loadData( const bool aBinaryLoad, const string & aFolderName, const int aNumMaterials );
        /**
        * \brief Unload the cross section data.
        */
        static void unloadData( );
        /**
        * \brief Constructor.
        */
        QSurface();
        /**
        * \brief Destructor.
        */
        ~QSurface();
        /**
        * \brief Write all the cross section data to the logger.
        */
        void dumpData( );
        /**
        * Is the class valid?
        */
        bool isValid() const { return _isValid; };

		static QSurfaceData* getData(){ return _qSurfaceData; }
    private:
        bool _isValid;                          //!< True, if the instance is valid
        static QSurfaceData *_qSurfaceData;     //!< Pointer to the static Q-surface data.
    };

    /**
    * \brief Provides access to the restricted stopping power.
    */
    class ZEUS_EXPORT RestrictedStoppingPower {
    public:
        /*! \brief Constructor, does nothing. */
        RestrictedStoppingPower( ) { };
        /*! \brief Destructor, does nothing. */
        ~RestrictedStoppingPower( ) { };
        /*! \brief Returns the restricted stopping power for en electron with energy \a e in material \a matid. 
         *
         * Note: this function is no longer used during the simulation. It is only used when loading the 
         *       cross section data files to prepare the LambdaTot tables. 
         */
        inline static double stpwr(double e, int matid) { return _data ? _data->interpolate(e,matid) : 0; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk using the logger. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );
    private:
        static InterpolationData<LinearDataEntry>* _data;  //!< Pointer to the static data
    };

    /**
    * \brief Provides access to the screening parameter.
    */
    class ZEUS_EXPORT ScreeningParameter {
    public:
        /*! \brief Constructor, does nothing. */
        ScreeningParameter( ) { };
        /*! \brief Destructor, does nothing. */
        ~ScreeningParameter( ) { };
        /*! \brief Returns the effective screening parameter of the multiple scattering distribution for electron with energy \a e in material \a matid.*/
        inline static double bw(double e, int matid) { return _data ? _data->interpolate(-1/e,matid) : 0; };
        /*! \brief Same as bw(), but now the input parameter is -1/e instead of e. */
        inline static double bw1(double ie, int matid) { return _data ? _data->interpolate(ie,matid) : 0; };
        /*! \brief Is this instance valid? */
        static bool isValid() { return (_data != 0); };
        /*! \brief Load the data from the folder \a aFolderName assuming there are \a aNumMaterials materials. */
        static int loadData( const string & aFolderName, const int aNumMaterials );
        /*! \brief Dump the data to disk using the logger. */
        static void dumpData( );
        /*! \brief Unload the data. */
        static void unloadData( );

		static InterpolationData<SplineDataEntry>* getData() { return _data; }
    private:
        static InterpolationData<SplineDataEntry>* _data;   //!< Pointer to the static data.
    };

}} // namespaces

