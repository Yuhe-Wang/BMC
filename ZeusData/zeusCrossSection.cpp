
#include "zeusCrossSection.h"
#include "zeusSpline.h"
#include "zeusRandomGenerator.h"
#include "zeusVoxelGeometry.h"

#include <cassert>
#include <cstdio>

using namespace Zeus::Random;

namespace Zeus { namespace CrossSectionData {

#ifdef MSVC
#pragma region SplineData
#endif

    SplineDataEntry* SplineDataEntry::prepareData(const double *xval, const double *fval, int ndat) {
        double *ya = new double[ndat];
        double *yb = new double[ndat];
        double *yc = new double[ndat];
        double *yd = new double[ndat];
        Spline::prepareSpline(xval,fval,&ya[0],&yb[0],&yc[0],&yd[0],0.0,0.0,ndat);
        ya[ndat-1] = yb[ndat-1] = yc[ndat-1] = yd[ndat-1] = 0.0;
        SplineDataEntry *data = new SplineDataEntry [ndat];
        for ( int i = 0 ; i < ndat ; i++ ) data[i] = SplineDataEntry(ya[i],yb[i],yc[i],yd[i]);
        delete [] ya; delete [] yb; delete [] yc; delete [] yd;
        return data;
    }

    void SplineDataEntry::logData(const char *fileName) const {
        LogProviders::InformationLog( "\t%s a=%g b=%g c=%g d=%g", fileName, _a, _b, _c, _d );
    }

    LinearDataEntry* LinearDataEntry::prepareData(const double *xval, const double *fval, int ndat) {
        LinearDataEntry *data = new LinearDataEntry [ndat];
        for ( int i = 1 ; i < ndat ; i++ ) {
            double b = (fval[i] - fval[i-1])/(xval[i] - xval[i-1]);
            double a = fval[i] - b*xval[i];
            data[i-1] = LinearDataEntry(a,b);
        }
        data[ndat-1] = data[ndat-2];
        return data;
    }

    void LinearDataEntry::logData(const char *fileName) const {
        LogProviders::InformationLog( "\t%s a=%g b=%g", fileName, _a, _b );
    }

    template<class DataEntry>
    InterpolationData<DataEntry>::InterpolationData( const string &aName, const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &aData) {
        _fileName = aName;
        _nmat = (int)aData.size();
        _ndata = new int[_nmat];
        _bIndex = new double[_nmat];
        _aIndex = new double[_nmat];
        _data = new DataEntry* [_nmat];
        for(int j=0; j<_nmat; ++j) _data[j] = 0;

        double *xval = new double [nEnergy];

        for ( int j = 0 ; j < _nmat ; j++ ) {
            double deltax = (emax[j] - emin[j])/(nEnergy-1);
            for(int i=0; i<nEnergy; ++i) xval[i] = emin[j] + deltax*i;
            _ndata[j] = nEnergy;
            _data[j] = DataEntry::prepareData(xval,aData[j],_ndata[j]);
            _bIndex[j] = xval[0];
            _aIndex[j] = ((double)_ndata[j]-1) / (xval[_ndata[j]-1] - xval[0]);
        }

        delete [] xval;
        _isValid = true;

    }

    template<class DataEntry> 
    InterpolationData<DataEntry>::InterpolationData( const string &aFileName, const int aNumMaterials, bool useInverseEnergy ) : 
        _ndata(0), _data(0), _bIndex(0), _aIndex(0), _nmat(aNumMaterials), _isValid(false), _fileName(aFileName) {

        const static char *func = "InterpolationData::InterpolationData";
        FILE *fp = fopen(aFileName.c_str(), "r" );
        if ( fp == NULL ) { 
            LogProviders::FatalLog( "%s: failed to open %s\n", func,aFileName.c_str() ); 
            throw( new Exceptions::FileException( aFileName) );
        }

        _ndata = new int[_nmat];
        _bIndex = new double[_nmat];
        _aIndex = new double[_nmat];
        _data = new DataEntry* [_nmat];
        for(int j=0; j<_nmat; ++j) _data[j] = 0;

        /* skip the header */
        int nitems = fscanf(fp,"%*[^\n] %*1[\n]");
        if( 0 ) LogProviders::InformationLog("Got %d items\n",nitems); // shut up compiler warning
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");

        /* loop over materials */
        for ( int j = 0 ; j < _nmat ; j++ ) {

            /* skip a line */
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");

            /* read number of data points */
            nitems = fscanf(fp,"%d\n",&_ndata[j]);

            if( _ndata[j] < 2 ) {
                LogProviders::FatalLog( "%s: failed to read number of data points from file %s\n",func,aFileName.c_str() );
                throw( new Exceptions::FileException( aFileName) );
            }

            double *xval = new double[_ndata[j]];
            double *fval = new double[_ndata[j]];

            /* skip two lines */
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");

            /* loop over data */
            for ( int i = 0 ; i < _ndata[j] ; i++ ) {
                nitems = fscanf(fp,"%lf %lf\n",&xval[i],&fval[i]);
                if( useInverseEnergy ) {
                    if (xval[i] == 0.0) {
                        fclose( fp );
                        throw( new Exceptions::ZeroDenominatorException( aFileName, j ) );
                    }
                    xval[i] = -1.0 / xval[i];
                }
            }

            _data[j] = DataEntry::prepareData(xval,fval,_ndata[j]);

            _bIndex[j] = xval[0];

            _aIndex[j] = ((double)_ndata[j]-1) / (xval[_ndata[j]-1] - xval[0]);

            delete [] xval; delete [] fval;

        } // number of materials

        /* close file */
        fclose(fp);
        _isValid = true;

    }

    template<class DataEntry> 
    InterpolationData<DataEntry>::~InterpolationData() {
        if( _ndata ) delete [] _ndata;
        if( _bIndex ) delete [] _bIndex;
        if( _aIndex ) delete [] _aIndex;
        if( _data ) {
            for(int j=0; j<_nmat; ++j) if( _data[j] ) delete [] _data[j];
            delete [] _data;
        }
    }

    template<class DataEntry> 
    void InterpolationData<DataEntry>::dumpData( ) const {
        if ( _isValid ) {
            LogProviders::InformationLog( "%s",_fileName.c_str());
            for ( int imat = 0 ; imat < _nmat ; ++imat ) {
                LogProviders::InformationLog( "%s Material %d bIndex %g aIndex %g", _fileName.c_str(), imat, _bIndex[imat], _aIndex[imat] );
                int len = _ndata[ imat ];
                for ( int i = 0 ; i < len ; ++i ) _data[imat][i].logData(_fileName.c_str());
            }
        }
    }

    template<class DataEntry> 
    InterpolationData<DataEntry>* InterpolationData<DataEntry>::loadData(const string &aDataFolder, const char *aFileName, int aNumMaterials, bool useInverseEnergy ) {
        if( !aFileName ) return 0;
        string fileName( aDataFolder );
        if( fileName.size() > 0 ) {
            char lastChar = fileName[ fileName.size()-1 ];
            if( lastChar != '/' && lastChar != 92 ) fileName += '/';
        }
        fileName += aFileName;
        return new InterpolationData<DataEntry>( fileName, aNumMaterials, useInverseEnergy ); 
    }

#ifdef MSVC
#pragma endregion

#pragma region Material
#endif

    MaterialData::MaterialData( const string & aFileName ) {

        FILE *fp = fopen(aFileName.c_str(), "r" );
        if ( fp == NULL ) { LogProviders::FatalLog( "Cant open %s", aFileName.c_str() ); throw( new Exceptions::FileException( aFileName) ); }

        /* skip the .matter header */
        int nitems;
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");	// This file is part of set:
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");	// the data set name
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");	// HEADER section: Material data file for ZeusMC.
        if( 0 ) LogProviders::InformationLog("Got %d items\n",nitems); // shut up compiler warning

        /* Get energy range of data set */
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");   // [Emin_ph,Emin,Emax] (eV): energy interval in which data is to be generated:
        nitems = fscanf(fp,"%lf %lf %lf\n",&_eminph,&_emin,&_emax);

        /* Get particle production threshold energies */
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");   // Wcc & Wcb (eV), cutoffs energies for collision and bremsstrahlung respectively:
        nitems = fscanf(fp,"%lf %lf\n",&_wcion,&_wcbre);

        /* Set transport cutoff energies to be equal to production thresholds */
        //_eabs = _wcion;
        _eabs = _emin + 1e3;
        _eabsph = _wcbre;

        /* Get step length parameters */
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");
        nitems = fscanf(fp,"%lf %lf %d\n",&_eloss,&_eRangeCut,&_csda);

        /* read number of materials */
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");  // No of materials in this file:
        _nmat = 0;
        nitems = fscanf(fp,"%d\n",&_nmat);
        if ( _nmat < 1 ) { LogProviders::FatalLog( "Cant read %s number of materials %d", aFileName.c_str(), _nmat ); return; }

        _matden = new double[_nmat];
        _zmass = new double[_nmat];

        /* skip a line */
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");

        /* loop over material data */
        for ( int i = 0 ; i < _nmat ; i++ ) 
        {

            /* skip a line (material name) */
            nitems = fscanf(fp,"%*[^\n] %*1[\n]]");

            /* material density (g/cm^3) */
            nitems = fscanf(fp,"%lf\n",&_matden[i]);

            /* skip a line (No elements in molecule:)*/
            nitems = fscanf(fp,"%*[^\n] %*1[\n]"); 

            /* number of elements in molecule */
            int k = 0;
            nitems = fscanf(fp,"%d\n",&k);
            for ( int j = 0 ; j < k ; j++ )
                nitems = fscanf(fp,"%*[^\n] %*1[\n]");

            /* skip a line (Total atno(Z), Z^2 & Z(Z+1):)*/
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");

            /* total atno(Z), Z^2 and Z(Z+1) */
            double atno = 0, atno2 = 0, zzp1 = 0;
            nitems = fscanf(fp,"%lf %lf %lf\n",&atno,&atno2,&zzp1);

            /* skip a line (Atomic mass (amu):)*/
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");

            /* atomic mass (amu) */
            double mass = 0;
            nitems = fscanf(fp,"%lf\n",&mass);

            /* skip a line (AtNo AtNo^2 to mass ratios with respect to ref:)*/
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");

            /* atomic no, atomic no^2 ratios with respect to reference */
            double z2mass = 0;
            nitems = fscanf(fp,"%lf %lf %*[^\n] %*1[\n]",&_zmass[i],&z2mass);

        }

        /* close material file */
        fclose(fp);

        /* density of reference material */
        _refden = _matden[0];

        /* density and stopping power for transporting electrons below eabs */
        _subden = 0.1 * _refden;
        _substp = 2.0e6;
        _subfac = _subden / _eabs;

        _fileName = aFileName;
        _isValid = true;

        printf("MaterialData::MaterialData: initialized with %d materials\n",_nmat);

    }

    void MaterialData::dumpData( ) {           
        if ( _isValid ) {
            LogProviders::InformationLog( "%s", _fileName.c_str() );
            LogProviders::InformationLog( "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g",
                    _eabs, _eabsph, _eminph, _emin, _emax, _subfac, _wcion, _wcbre );
            for ( int i = 0 ; i < _nmat ; ++i ) {
                LogProviders::InformationLog( "%g\t%g", _matden[i], _zmass[i] );
            }
        }
    }

    MaterialData::~MaterialData( ) {
        _isValid = false;
        delete [] _matden;
        delete [] _zmass;
    };
    int MaterialData::numMaterials() { return _nmat; }

    Material::Material( ) { };
     int Material::loadData( const string & aFolderName ) { 
        string fileName( aFolderName );
        char lastChar = fileName[ fileName.size()-1 ];
        if( lastChar != '/' && lastChar != 92 ) fileName += '/';
        fileName += "vpmc.matter";
        _matData = new MaterialData( fileName ); 
        return _matData->isValid() ? 0 : 1; 
    }
    void Material::dumpData( ) { if ( _matData ) _matData->dumpData( ); }
     void Material::unloadData( ) { 
        if ( _matData ) delete _matData; _matData = 0; 
    }
    int Material::numMaterials() { return _matData->numMaterials( ); }
    Material::~Material( ) {};

#ifdef MSVC
#pragma endregion 

#pragma region Bremmstrahling
#endif

    //
    // Cross section data, Bremmstrahling
    //

	BremmConstantsData::BremmConstantsData(const string & aFolderName, const int aNumMaterials) {

        _nmat = aNumMaterials;
		
		_fileName = aFolderName;
		if (_fileName.size() > 0) {
			char lastChar = _fileName[_fileName.size() - 1];
			if (lastChar != '/' && lastChar != 92) _fileName += '/';
		}
		_fileName += "vpmc.brecon";
		const static char *func = "BremmConstantsData::BremmConstantsData";
		FILE *fp = fopen(_fileName.c_str(), "r");
		if (fp == NULL) {
			LogProviders::FatalLog("%s: failed to open %s\n", func, _fileName.c_str());
			throw(new Exceptions::FileException(_fileName));
		}

		/* skip the header */
		int nitems = fscanf(fp, "%*[^\n] %*1[\n]");
		if (0) LogProviders::InformationLog("Got %d items\n", nitems); // shut up compiler warning
		nitems = fscanf(fp, "%*[^\n] %*1[\n]");

		/* loop over materials */
		_f0 = new double*[_nmat];
		_bcb = new double[_nmat];

		/* loop over materials */
		for (int j = 0; j < _nmat; j++) {

			/* skip a line */
			nitems = fscanf(fp, "%*[^\n] %*1[\n]");

			_f0[j] = new double[3];

			/* read number of data points */
			fscanf(fp, "%lf %lf %lf\n", &_f0[j][0], &_f0[j][1], &_f0[j][2]);
			fscanf(fp, "%lf\n", &_bcb[j]);

		} // number of materials

		/* close file */
		fclose(fp);
		_isValid = true;
    }
    void BremmConstantsData::dumpData( ) {  
        if ( _isValid ) {
            LogProviders::InformationLog( "%s",_fileName.c_str());
            for ( int i = 0 ; i < _nmat ; ++i ) {
                LogProviders::InformationLog( "%s %g\t%g\t%g\t%g",_fileName.c_str(), _f0[i][0], _f0[i][1], _f0[i][2], _bcb[i] );
            }
        }
    }
    double *BremmConstantsData::bcb() { return _bcb; }
    double *BremmConstantsData::fo( int materialIndex ) { return _f0[materialIndex]; }
    BremmConstantsData::~BremmConstantsData() { delete [] _bcb; for ( int i = 0 ; i < _nmat ; ++i ) delete [] _f0[i]; delete [] _f0; _isValid = false; }


    BremmConstants::BremmConstants( ) : _isValid(false) { _material = new  CrossSectionData::Material(); }
	int BremmConstants::loadData(const string & aFolderName) {
		_bremmData = new BremmConstantsData(aFolderName, Material::numMaterials());
        return _bremmData->isValid() ? 0 : 1; 
    }

    void BremmConstants::unloadData() { if ( _bremmData ) { delete _bremmData; _bremmData = 0; } }

    void BremmConstants::dumpData( ) { if ( _bremmData ) _bremmData->dumpData( ); }

    double BremmConstants::sambre ( double e, int matid, RandomGenerator *aRng )
    {
        double *f0 = _bremmData->fo( matid );

        double etot, gamma, em, ed, ec, ec1, ec2, emecs, xx, f2le;
        double f00, f1, f2, f1c, f2c, f1d, f2d, ed1, bec, bed;
        double eps, eps1, pa, pa0, pa1, pa2, pb, pb2;

        double wcbre = _material->wcbre();

        unsigned int loop = 1;

        etot = e + Constants::ELECTRON_MASS; // >e
        gamma = etot / Constants::ELECTRON_MASS; // >1
        pa = _bremmData->bcb()[matid] * gamma;
        em = e / etot; // in (0,1)
        if ( em > 0.99999 )
            em = 0.99999;

        /* rejection functions */
        ed = (e - 5.0*Constants::ELECTRON_MASS) / etot;// <1 !, >0?
        ec = wcbre / etot;// < 1.0?
        ec2 = ec*ec;
        emecs = em*em - ec2;

        /* low energy correction */
        xx = e / Constants::ELECTRON_MASS;

        f2le = (4.650 - f0[2]*(6.005 - f0[2]*2.946)) / xx -
            (32.42 - f0[2]*(67.08 + f0[2]*3.906)) / pow(1+xx,2) +
            (20.33 + f0[2]*(23.38 - f0[2]*77.42)) / pow(1+xx,3);

        if ( ec < ed ){
            f00 = f0[1] + f2le;
        }
        else{
            f00 = f0[0] + f2le;
        }

        ec1 = 1.0 - ec;
        bec = ec / (pa * ec1);
        schiff(bec,&f1,&f2);

        f1c = f1 + f00;
        f2c = (f2 + f00) * ec1;
        if ( ec < ed ) {
            f00=f0[0] + f2le;
            ed1 = 1.0 - ed;
            bed = ed / (pa * ed1);
            schiff(bed,&f1,&f2);
            f1d = f1 + f00;
            f2d = (f2 + f00) * ec1;
            if ( f1d > f1c ) f1c = f1d;
            if ( f2d > f2c ) f2c = f2d;
        }
        pa1 = emecs * f1c;
        pa2 = 8.0/3.0 * log(em/ec) * f2c;

        while ( loop ) {
            if ( aRng->getUniform()*(pa1+pa2) <= pa1 ) {
                eps = sqrt(ec2 + aRng->getUniform()*emecs);
                pb = eps / (pa * (1.0 - eps));
                f1 = 2.0 - 2.0*log(1.0 + pb*pb);
                if ( pb < 1.0e-8 ){
                    f1 -= pb * 2.0 * Constants::PI;
                }
                else{
                    f1 -= 4.0 * pb * atan2(1.0,pb);
                }
                if ( eps < ed ){
                    f00 = f0[1] + f2le;
                }
                else{
                    f00 = f0[0] + f2le;
                }
                if ( aRng->getUniform()*f1c < f1+f00 ){
                    return eps*etot;
                }
            }
            else {
                eps = ec * pow(em/ec,aRng->getUniform());
                eps1 = 1.0 - eps;
                pb = eps / (pa*eps1);
                pb2 = pb*pb;
                f1 = 2.0 - 2.0*log(1.0+pb2);
                f2 = f1 - 2.0/3.0;
                if ( pb < 1.0e-8 ){
                    f1 -= pb * 2.0 * Constants::PI;
                }
                else {
                    pa0 = 4.0 * pb * atan2(1.0,pb);
                    f1 -= pa0;
                    f2 = f2 + 2.0 * pb2 * (4.0 - pa0 - 3.0*log((1.0+pb2)/pb2));
                }
                f2 = 0.5 * (3.0*f1 - f2);
                if ( eps < ed ){
                    f00 = f0[1] + f2le;
                }
                else{
                    f00 = f0[0] + f2le;
                }
                if ( aRng->getUniform()*f2c <= eps1*(f2+f00) ){
                    return eps*etot;
                }
            }
        }

        return 0.0;

    }
    BremmConstants::~BremmConstants() { delete _material; }

    void BremmConstants::schiffPublic ( double b, double *f1, double *f2 ) {
        schiff( b, f1, f2 );
    }

    void BremmConstants::schiff ( double b, double *f1, double *f2 ) {
        double b2, a0;
        b2 = b*b;
        *f1 = 2.0 - 2.0 * log(1.0+b2);
        *f2 = *f1 - 2.0/3.0;
        if ( b < 1.0e-10 ){
            *f1 = *f1 - b * 2.0*Constants::PI;
        }
        else {
            a0 = 4.0 * b * atan2(1.0,b);
            *f1 -= a0;
            *f2 = *f2 + 2.0 * b2 * (4.0 - a0 - 3.0 * log((1.0+b2)/b2));
        }
        *f2 = 0.5 * (3.0*(*f1) - *f2);
    }

#ifdef MSVC
#pragma endregion 
#endif

    InterpolationData<LinearDataEntry>* LambdaTot::_data = 0;

    void LambdaTot::prepareData(const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &aData) {
        if( _data ) delete _data;
        _data = new InterpolationData<LinearDataEntry>("LambdaTot",emin,emax,nEnergy,aData);
    }

    void LambdaTot::unloadData() {
        if( _data ) { delete _data; _data = 0; }
    }

    InterpolationData<LinearDataEntry>* InvLambdaTot::_data = 0;

    void InvLambdaTot::prepareData(const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &aData) {
        if( _data ) delete _data;
        _data = new InterpolationData<LinearDataEntry>("LambdaTot",emin,emax,nEnergy,aData);
    }

    void InvLambdaTot::unloadData() {
        if( _data ) { delete _data; _data = 0; }
    }

    InterpolationData<LinearDataEntry>* BremsProbability::_data = 0;

    void BremsProbability::prepareData(const vector<double> &emin, const vector<double> &emax, int nEnergy, const vector<double*> &aData) {
        if( _data ) delete _data;
        _data = new InterpolationData<LinearDataEntry>("LambdaTot",emin,emax,nEnergy,aData);
    }

    void BremsProbability::unloadData() {
        if( _data ) { delete _data; _data = 0; }
    }


#define XSECTION_IMPLEMENTATION(ClassName,InerpolationType,dataFile,useInverse) \
    InterpolationData<InerpolationType>* ClassName::_data = 0; \
    int ClassName::loadData( const string & aFolderName, const int aNumMaterials ) { \
        unloadData();\
        _data = InterpolationData<InerpolationType>::loadData(aFolderName,dataFile,aNumMaterials,useInverse);\
        int result = 0;\
        if( !_data->isValid() ) {\
            delete _data; _data = 0; result = 1;\
        }\
        return result;\
    }\
    void ClassName::unloadData() {\
        if( _data ) {\
            delete _data; _data = 0;\
        }\
    }\
    void ClassName::dumpData() {\
        if( _data ) _data->dumpData();\
    }

    XSECTION_IMPLEMENTATION(Kerma,LinearDataEntry,"vpmc.kerma",false);

    XSECTION_IMPLEMENTATION(LambdaBr,LinearDataEntry,"vpmc.lambre",false);

    XSECTION_IMPLEMENTATION(LambdaCo,SplineDataEntry,"vpmc.compt",false);

    XSECTION_IMPLEMENTATION(LambdaMoller,SplineDataEntry,"vpmc.lammo",true);

    XSECTION_IMPLEMENTATION(LambdaPair,SplineDataEntry,"vpmc.pairp",false);

    XSECTION_IMPLEMENTATION(LambdaPhoton,SplineDataEntry,"vpmc.lamph",false);

    XSECTION_IMPLEMENTATION(RestrictedStoppingPower,LinearDataEntry,"vpmc.rstpw",false);

    XSECTION_IMPLEMENTATION(ScreeningParameter,SplineDataEntry,"vpmc.bw",true);

    XSECTION_IMPLEMENTATION(Range,LinearDataEntry,"vpmc.r",false);

    XSECTION_IMPLEMENTATION(InverseRange,LinearDataEntry,"vpmc.ir",false);
   
	void LambdaPhoton::SetupWoodcock(int aNumElements, const float *dens, const char *mat, int nvoxels) {

        if( !_data ) return;
        if( _wckLen > 0 ) {
            delete [] _wcka; delete [] _wckb;
            _wckLen = 0;
        }

        CrossSectionData::Material* _material = new  CrossSectionData::Material(); 

        int nmat = _material->numMaterials();
        double *maxden = new double[nmat];

        _wckLen = aNumElements;
        _wcka = new double[aNumElements];
        _wckb = new double[aNumElements];

        double eminph = _material->getEminph();
        double emax = _material->getEmax( );

        double mine = eminph;
        double maxe = emax;
        for ( int j = 0 ; j < nmat ; j++ ) maxden[j] = 0.0;

//         const float *dens = aGeom->_rho;
//         const char *mat = aGeom->_material;
//         int nvoxels = aGeom->_Nx * aGeom->_Ny * aGeom->_Nz;

        for ( int i = 0 ; i < nvoxels ; i++ ) {
            int matid = mat[i] - 1;
            if ( matid >= 0 && matid < nmat )
            {
                if ( dens[i] > maxden[matid] ) 
                    maxden[matid] = dens[i];
            }
        }

        /* prepare linear interpolation */
        _ewck0 = mine;
        double denom = double(_wckLen);
        double de = (maxe * (1.0 - Constants::EPS) - _ewck0) / denom;
        _fwck = 1.0 / de;

        double ylast = 0;
        for ( int i = 0 ; i <= _wckLen ; i++ ) {
            double e = _ewck0 + de * i;
            double ymax = 0.0;
            for ( int j = 0 ; j < nmat ; j++ ) {
                double ymaybe = lmbdaph(e,j) * maxden[j];
                if ( ymaybe > ymax ) ymax = ymaybe;
            }
            double ymin = 1.0 / ymax;
            if ( i == 0 )
                ylast = ymin;
            else {
                _wckb[i-1] = (ymin - ylast) * _fwck;
                _wcka[i-1] = ylast - (e - de) * _wckb[i-1];
                ylast = ymin;
            }
        }

        delete [] maxden;
        delete _material;

    }

    // Write the q surface data in a binary chunk.
    void QSurfaceData::dumpBinary( const string & aFileName ) {

        string binaryFileName = aFileName;
        binaryFileName += string( ".binary" );
        size_t nsize = 0, ncount = 0, nwritten = 0;

        // Binary write 
        FILE *fp = fopen( binaryFileName.c_str(), "wb" );
        if ( fp == NULL ) { LogProviders::WarningLog( "cant open qsurface binary file for output" ); return; }
        nsize = sizeof(int); ncount = 1; nwritten = fwrite( &_nmat, nsize, ncount, fp ); assert( nwritten == ncount ); 				
        nsize = sizeof(double); ncount = _nmat; nwritten = fwrite( _ieq0, nsize, ncount, fp ); assert( nwritten == ncount );
        nsize = sizeof(double); ncount = _nmat; nwritten = fwrite( _fq, nsize, ncount, fp ); assert( nwritten == ncount );
        nsize = sizeof(int); ncount = _nmat; nwritten = fwrite( _ne, nsize, ncount, fp ); assert( nwritten == ncount );								 		
        nsize = sizeof(int); ncount = _nmat; nwritten = fwrite( _nuq, nsize, ncount, fp ); assert( nwritten == ncount );		
        nsize = sizeof(double);
        for ( int k = 0 ; k < _nmat ; k++ ) 
        {
            for ( int j = 0 ; j < _nuq[k] ; j++ ) 
            {
                ncount = _ne[k]; nwritten = fwrite( _qs[k][j], nsize, ncount, fp ); assert( nwritten == ncount ); 
            }
        }
        fclose(fp);

    }

    // Load the q surface using the old ascii format 
    void QSurfaceData::asciiLoad( const string & aFileName ) {

        _isValid = false;
        FILE *fp = fopen(aFileName.c_str(), "r" );
        if ( fp == NULL ) { LogProviders::FatalLog( "Cant open %s", aFileName.c_str() ); throw( new Exceptions::FileException( aFileName) ); }
        /* skip the header */
        int nitems = fscanf(fp,"%*[^\n] %*1[\n]");
        if( 0 ) LogProviders::InformationLog("Got %d items\n",nitems); // shut up compiler warning
        nitems = fscanf(fp,"%*[^\n] %*1[\n]");
        /* loop over materials */
        _nuq = new int[_nmat];
        _ne = new int[_nmat];
        _ieq0 = new double[_nmat];
        _fq = new double[_nmat];
        _qs = new double**[_nmat];
        _Qmax = new double* [_nmat];

        for ( int k = 0 ; k < _nmat ; k++ ) {
            /* skip a line */
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");
            /* read number of data points */
            int maxj = 0, maxi = 0;
            nitems = fscanf(fp,"%d %d\n",&maxj,&maxi);
            _nuq[k] = maxj;
            /* skip lines */
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");
            /* read qmax */
            double qmax = 0;
            nitems = fscanf(fp,"%lf\n",&qmax);
            if ( qmax == 0.0f) {
                throw( "error: qmax can not be 0.0f" );
            }
            _qs[k] = new double*[maxi];
            for ( int i = 0 ; i < maxi ; i++ ) { 
                _qs[k][i] = new double[maxj];
            }
            /* skip lines */
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");
            nitems = fscanf(fp,"%*[^\n] %*1[\n]");
            /* loop over data */
            _ne[k] = maxi;
            _Qmax[k] = new double [maxi];
            double e0 = 0;
            for ( int i = 0 ; i < maxi ; i++ ) {
                double e1 = 0;
                nitems = fscanf(fp,"%lf\n",&e1);
                if ( e1 == 0.0f) {
                    fclose( fp );
                    throw( "error: e1 can not be 0.0f" );
                }
                if ( i == 0 ) {
                    e0 = e1;
                } else if ( e1 == e0) {				
                    throw("error: e1 can not be same as e0");
                }
                double qmax = 0;
                for ( int j = 0 ; j < maxj ; j++ ) {
                    double u = 0, qval = 0;
                    nitems = fscanf(fp,"%lf %lf\n",&u,&qval);
                    if( qval > qmax ) qmax = qval;
                    _qs[k][i][j] = qval;
                }
                _Qmax[k][i] = qmax;
                double qmaxi = 1/qmax;
                for ( int j = 0 ; j < maxj ; j++ ) _qs[k][i][j] *= qmaxi;
                _ieq0[k]  = -1.0 / e0;
                _fq[k] = (_ne[k] - 1) / (-1.0/e1 - _ieq0[k]);
            }
        }
        /* close file */
        fclose(fp);

        _isValid = true;
        _fileName = aFileName;

        // dumpBinary( _fileName );

    }

    // Load the q surface using the binary format
    void QSurfaceData::binaryLoad( const string & aFileName ) {

        _isValid = false;

        FILE *fp = fopen(aFileName.c_str(), "rb" );
        if ( fp == NULL ) throw( new Exceptions::FileException( aFileName ) );
        size_t nsize = 0, ncount = 0, nread = 0;
        nsize = sizeof(int); ncount = 1; nread = fread( &_nmat, nsize, ncount, fp ); assert( nread == ncount ); 

        _nuq = new int[_nmat];
        _ne = new int[_nmat];
        _ieq0 = new double[_nmat];
        _fq = new double[_nmat];
        _qs = new double**[_nmat];

        nsize = sizeof(double); ncount = _nmat; nread = fread( _ieq0, nsize, ncount, fp ); assert( nread == ncount );
        nsize = sizeof(double); ncount = _nmat; nread = fread( _fq, nsize, ncount, fp ); assert( nread == ncount );
        nsize = sizeof(int); ncount = _nmat; nread = fread( _ne, nsize, ncount, fp ); assert( nread == ncount );								 		
        nsize = sizeof(int); ncount = _nmat; nread = fread( _nuq, nsize, ncount, fp ); assert( nread == ncount );		
        nsize = sizeof(double);
        for ( int k = 0 ; k < _nmat ; k++ ) 
        {
            _qs[k] = new double*[_nuq[k]];
            for ( int j = 0 ; j < _nuq[k] ; j++ ) 
            {
                _qs[k][j] = new double[_ne[k]];
                ncount = _ne[k]; nread = fread( _qs[k][j], nsize, ncount, fp ); assert( nread == ncount ); 
            }
        }
        fclose(fp);

        _isValid = true;
        _fileName = aFileName;

    }

    // Read the qsurface data.  If aBinaryLoad is false, we load the ascii format from vpmc.q.
    // Otherwise we load in binary.
    QSurfaceData::QSurfaceData( const bool aBinaryLoad, const string & aFileName, const int aNumMaterials ) : _isValid(false), _nmat(aNumMaterials) {
        if ( aBinaryLoad == false ) {
            asciiLoad( aFileName ); 
        } else {
            string binaryFileName = aFileName;
            binaryFileName += ".binary";
            binaryLoad( binaryFileName ); 
        }
    }

    QSurfaceData::~QSurfaceData() { 
        _isValid = false; delete [] _ieq0; delete [] _fq; 
        for ( int k = 0 ; k < _nmat ; k++ ) for ( int j = 0 ; j < _nuq[k] ; j++ ) delete [] _qs[k][j];
        for ( int k = 0 ; k < _nmat ; k++ ) { delete [] _qs[k]; delete [] _Qmax[k]; }
        delete [] _qs; delete [] _Qmax;
        delete [] _nuq; delete [] _ne;
    }
    void QSurfaceData::dumpData( ) { 
        if ( _isValid ) {
            LogProviders::InformationLog( "%s",_fileName.c_str());
            for ( int k = 0 ; k < _nmat ; k++ ) {
                LogProviders::InformationLog( "material = %d\n", k );
                LogProviders::InformationLog( "ieq0[%d]=%e fq[%d]=%e nuq[%d]=%d, ne[k]=%d\n", 
                    k, _ieq0[k], k, _fq[k], k, _nuq[k], k, _ne[k] );
                for ( int j = 0 ; j < _nuq[k] ; j++ ) {
                    std::ostringstream streamer;
                    streamer << _fileName << string(" ") << j << string(" ");
                    for ( int i = 0 ; i < _ne[k] ; ++i ) {
                        streamer << _qs[k][j][i];  streamer << string(" ");
                    }
                    LogProviders::InformationLog( "%s", streamer.str().c_str() );
                }
            }
        }
    }

    double QSurface::qsurf ( const double u, const double aEnergy, const int aMaterialId ) {

        double ie = -1.0 / aEnergy;
        int nuq = _qSurfaceData->_nuq[aMaterialId] - 1;
        int ne = _qSurfaceData->_ne[aMaterialId] - 1;

        double ru = u * nuq;
        int ju = (int) ru; if( ju > nuq-1 ) ju = nuq-1;
        ru -= ju; double ru1 = 1 - ru;

        double re = _qSurfaceData->_fq[aMaterialId] * (ie - _qSurfaceData->_ieq0[aMaterialId]);
        int je = (int) re; if( je > ne-1 ) ne = ne-1;
        re -= je; double re1 = 1 - re;

        const double *q1 = _qSurfaceData->_qs[aMaterialId][je];
        const double *q  = _qSurfaceData->_qs[aMaterialId][je+1];
        double qmax1 = _qSurfaceData->_Qmax[aMaterialId][je];
        double qmax  = _qSurfaceData->_Qmax[aMaterialId][je+1];
        double iqmax = qmax > qmax1 ? 1/qmax : 1/qmax1;

        return iqmax*(re1*(q1[ju]*ru1 + q1[ju+1]*ru) + re*(q[ju]*ru1 + q[ju+1]*ru));

    }

     int QSurface::loadData( const bool aBinaryLoad, const string & aFolderName, const int aNumMaterials ) { 
        string fileName( aFolderName );
        char lastChar = fileName[ fileName.size()-1 ];
        if( lastChar != '/' && lastChar != 92 ) fileName += '/';
        fileName += "vpmc.q";
        _qSurfaceData = new QSurfaceData( aBinaryLoad, fileName, aNumMaterials ); 
        return _qSurfaceData->isValid() ? 0 : 1; 
    }
     void QSurface::unloadData( ) { 
        if ( _qSurfaceData ) { delete _qSurfaceData; _qSurfaceData = 0; }
    }

    QSurface::QSurface() : _isValid(false) { }
    QSurface::~QSurface() {}
    void QSurface::dumpData( ) { if ( _qSurfaceData ) _qSurfaceData->dumpData( ); }


    MaterialData* Material::_matData = 0;
    BremmConstantsData* BremmConstants::_bremmData = 0;
    QSurfaceData *QSurface::_qSurfaceData = 0;

    vector<InterpolationData<LinearDataEntry>* >  HeadAttenuation::_totData;
    vector<InterpolationData<LinearDataEntry>* >  HeadAttenuation::_compData;
    string HeadAttenuation::_folderName;

    int HeadAttenuation::loadData( const string & aFolderName ) {
        const static char *func = "HeadAttenuation::loadData";
        if( aFolderName == _folderName && isValid() ) return 0;

        unloadData();

        _totData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenTungstenTot",1,false) );
        _totData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenAirTot",1,false) );
        _totData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenGradCoilTot",1,false) );
        _totData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenRFShieldTot",1,false) );
        _totData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenTxPCBTot",1,false) );
        _totData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenTxFR4Tot",1,false) );
 
        _compData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenTungstenInco",1,false) );
        _compData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenAirInco",1,false) );
        _compData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenGradCoilInco",1,false) );
        _compData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenRFShieldInco",1,false) );
        _compData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenTxPCBInco",1,false) );
        _compData.push_back( InterpolationData<LinearDataEntry>::loadData(aFolderName,"vpmc.AttenTxFR4Inco",1,false) );

        _folderName = aFolderName;

        bool isValid = true;
        for(size_t k=0; k<_totData.size(); ++k) {
            if( !_totData[k]->isValid() ) {
                LogProviders::WarningLog("%s: invalid data set %s\n",func,_totData[k]->fileName().c_str());
                isValid = false;
            }
        }
        for(size_t k=0; k<_compData.size(); ++k) {
            if( !_compData[k]->isValid() ) {
                LogProviders::WarningLog("%s: invalid data set %s\n",func,_compData[k]->fileName().c_str());
                isValid = false;
            }
        }

        int result = 0;
        if( !isValid ) {
            unloadData();
            result = 1;
        }

        return result;

    }

    void HeadAttenuation::unloadData( ) {
        if( _totData.size() > 0 ) {
            for(size_t k = 0; k<_totData.size(); ++k) if( _totData[k] ) delete _totData[k];
            _totData.clear();
        }
        if( _compData.size() > 0 ) {
            for(size_t k = 0; k<_compData.size(); ++k) if( _compData[k] ) delete _compData[k];
            _compData.clear();
        }
        _folderName = "";
    }

}} // namespaces
