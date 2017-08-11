/*!   \file   zeusConfig.h
 *    \brief  Some basic definitions used everywhere
 *    \author Tony Apicella, Iwan Kawrakow
 */
#pragma once

#ifdef WIN32
#define _WINDOWS
#define _CRT_SECURE_NO_WARNINGS
#include "windows.h"
#endif

/*! The floating point type used in the Zeus package */
#ifdef ZEUS_SINGLE
    typedef float  ZFloat;
#else
    typedef double ZFloat;
#endif

/*! To avoid replacing in all KMC-provided classes */
typedef ZFloat KMCFloat;

#ifdef _USE_DLL

	#ifdef _WINDOWS
		#ifdef MAKE_DLL
			#define ZEUS_EXPORT __declspec(dllexport)
		#else
			#define ZEUS_EXPORT __declspec(dllimport)
		#endif
	#else
/*! Macro to be able to build/use a DLL if needed */
		#define ZEUS_EXPORT 
/*! Macro to be able to build/use a DLL if needed */
	#endif

#else

/*! \brief Most classes are declared with this macro, in case we want to build a DLL on Windows. */
	#define ZEUS_EXPORT

#endif

/*! The 64-bit integer type */
#ifdef _WINDOWS
    typedef __int64 ZEUS_I64;
#else
    typedef long long ZEUS_I64;
#endif

/*! \brief Namespace for the Zeus Monte Carlo. */
namespace Zeus {

/*!  \brief Enumeration of the various modalities. This is used in the Image3D teampate.
*
* For this app, the only relevant modalities are Dose and
* Relative Electron Density.
*/
enum Modality { 
    DOSE=2,    //!< The volume is a dose volume
    RED=3      //!< The volume is a relative electron density volume
};

/*! Log level */
enum LogLevel {
    LogLevelNone = 0, 
    LogLevelMedium, 
    LogLevelMaximum,
    LogLevelInfinite,
    LogLevelMoreThanInfinite
};

/*! Windows compilers tend to not have the round function defined in cmath. This is a macro replacing it. */
#define IROUND(x) int( (x > 0.0f) ? floor(x + 0.5f) : ceil(x - 0.5f) );

/**
* \brief A structure to hold a three element vector.  Each component is stored in double precision.
*
* In principle we could remove this structure and use Vector instead, but this is just a simple structure 
* (without all the methods that Vector provides) and is used everywhere during the simulation for 
* particle positions and directions.
*/
struct mcvector {

    /*! \brief Construct from the input \a aX, \a aY and \a aZ. */
    mcvector( double aX, double aY, double aZ ) : x(aX), y(aY), z(aZ) {}

    /*! \brief Dummy constructor. Sets all components to zero. */
    mcvector( ) : x(0), y(0), z(0) {}

    /**
     * X component of vector.
     */
    double x;
    /**
     * Y component of vector.
     */
    double y;
    /**
     * Z component of vector.
     */
    double z;
};

/**
* Status returned from the worker threads.
*/
enum Status {
    OK = 0, 
    ComputeException,
    CalcElectronDirectionException,
    FileException,
    StackException,
    ZeroDenominatorException,
    UnknownMaterialException,
    FloatingPointException,
    HasntRunYet,
    CalculationInterrupted
};

}

