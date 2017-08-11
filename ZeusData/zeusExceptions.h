#pragma once

#include <string>
#include <sstream>
using namespace std;

#include "zeusConfig.h"
#include "zeusVector.h"

namespace Zeus {  
/**
* \brief Namespace containing all the exceptions thrown in ZeusMC.
*/
    namespace Exceptions {

        /** \brief An error has occured in the cross section calculations.
         *
         * Thrown when index into cross section array is out of bounds, or when energy too low.
         */
        class ZEUS_EXPORT ComputeException {
            public:
                /** This exception is thrown when calculations on the cross section data, mostly interpolation between the table entries, will fail because
                 * an array index is out of bounds, or because the energy is too low.
                 * \param aCrossSectionName Which cross section file was this data loaded from.
                 * \param aMaterialId Which material in the cross section file is under consideration.
                 * \param aIndex The index we computed, which is out of bounds.
                 * \param aLowerBound The lower bound of the index
                 * \param aUpperBound The upper bound of the index.
                 */
                ComputeException( const string & aCrossSectionName, const int aMaterialId, const int aIndex, const int aLowerBound, const int aUpperBound );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                /**
                 * Index into the material array.   Most of the arrays are two dimensional, the first index 
                 * being material, and the second index being energy.
                 */
                const int _materialId;
                /**
                 * Name of cross section file that the cross section data was read from.  
                 */
                const string _crossSectionClass;
                /** 
                 * Index into the cross section arrays.
                 */
                const int _index;
                /**
                 * The lower bound on the cross section array.
                 */
                const int _lowerBound;
                /**
                 * The upper bound on the cross section array.
                 */
                const int _upperBound;
        };

        /** \brief Encountered a zero divide.
         *
         * Thrown before the division, when a denominator is either zero or very small.
         */
        class ZEUS_EXPORT ZeroDenominatorException {
            public:
                /**
                 * Constructor
                 * \param aCrossSectionName Which cross section file was this data loaded from.
                 * \param aMaterialId Which material in the cross section file is under consideration.
                 */
                ZeroDenominatorException( const string &aCrossSectionName, const int aMaterialId );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                /**
                 * Index into the cross section arrays for material under consideration.
                 */
                const int _materialId;
                /**
                 * Name of cross section file that the cross section data was read from.  
                 */			
                const string _crossSectionClass;
        };

        /** \brief Cant open a cross section file.
         *
         * Thrown when cant open cross section file.
         */
        class ZEUS_EXPORT FileException {
            public:
                /**
                 * Constructor for file exception.
                 * \param aCrossSectionName The cross section file that could not be opened.
                 */
                FileException( const string & aCrossSectionName );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                /**
                 * Name of cross section file that the cross section data was read from.  
                 */		
                const string _crossSectionClass;
        };

        /** \brief Bad input to CalcElectronDirection.
         *
         * Thrown when the argument to method CalcElectronDirection is out of bounds.
         */
        class ZEUS_EXPORT CalcElectronDirectionException {
            public:
                /** 
                 * Constructor.  
                 * \param aCosTheta Argument to CalcElectronDirection that caused this exception.
                 */
                CalcElectronDirectionException( const double aCosTheta );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                const double _theta;
        };

        /** \brief The stack has overflowed.
         * 
         * Thrown when the stack has overflowed. 
         */
        class ZEUS_EXPORT StackException {
            public:
                /** Constructor.
                 */
                StackException( );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
        };

        /** \brief The material that particle is in is unknown.
         * 
         * Thrown when attenuation must be calculated.  This attenuation is material specific,
         * so if we dont know which material we are in, we cannot proceed.
         */
        class ZEUS_EXPORT UnknownMaterialException {
            public:
                /** Constructor.
                 * Used during transport of particle.  
                 */
                UnknownMaterialException( const string & aLabel, const Vector & aPos, const Vector & aDir, const int aCurRegion, const int aNextRegion, const ZFloat t ) :
                    _label(aLabel), _curRegion(aCurRegion), _nextRegion(aNextRegion), _position(aPos), _direction(aDir), _t(t) {}
                /** Constructor.
                 * Used during calc attenuation.
                 */
                UnknownMaterialException( const string & aLabel ) :
                    _label(aLabel), _curRegion(-1), _nextRegion(-1), _position(Vector(0,0,0)), _direction(Vector(0,0,0)), _t(0) {}
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                /**
                 * Descriptive label, perhaps location in code where this occurred.
                 */
                const string _label;
                /**
                 * Which region are we in.
                 */
                const int _curRegion;
                /**
                 * Next region along path.
                 */
                const int _nextRegion;
                /* 
                 * Particle position.
                 */
                const Vector _position;
                /**
                 * Particle direction.
                 */
                const Vector _direction;
                /**
                 * Intended distance.
                 */
                const ZFloat _t;
        };
        /*! \brief There were too mny geometry errors. 
         *
         * Not really needed anymore. It is here from the early stages where there could be quite a few geometry related errors.
         */
        class ZEUS_EXPORT TooManyGeometryErrorsException {
            public:
                /*! \brief Constructor */
                TooManyGeometryErrorsException(const string &description) : _description(description) {};
                /*! \brief Returns the description of the exception */
                const string whats() { return _description; };
            private:
                /*! The error description string */
                string _description;
        };

        /** \brief Input to spline wasnt monotonic.
         * 
         * Spline function takes two arrays, one for x and one for y.  The values in the x array must be monotonically increasing.
         */
        class ZEUS_EXPORT MonotonicSplineException {
            public:
                /** 
                 * Constructor.
                 */
                MonotonicSplineException( );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
        };

        /** \brief Log linear linear interpolation failed.
         * 
         */
        class ZEUS_EXPORT LogLinLinException {
            public:
                /** 
                 * Constructor.
                 */
                LogLinLinException( const double aInterpLoc, const double aLowerBound, const double aUpperBound );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                const double _interpLoc;
                const double _lowerBound;
                const double _upperBound;
        };

        /** \brief Input to spline did not have enough points.
         * 
         * Spline function takes two arrays, one for x and one for y.  There must be at least 
         * four points in these arrays.
         */
        class ZEUS_EXPORT NotEnoughPointsSplineException {
            public:
                /**
                 * Constructor.
                 * \param aNpts How many points were there in the table.
                 */
                NotEnoughPointsSplineException( int aNpts );
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                /**
                 * How many points were in the spline tables.
                 */
                int _npts;
        };

        /** \brief Variable representation is not floating point.
         * 
         * A variable has a value such as infinite or indeterminant.
         */
        class ZEUS_EXPORT FloatingPointException {
            public:
                /**
                 * Constructor.
                 * \param aLabel Where in the code did this happen.
                 * \param aValue The bad value.
                 */
                FloatingPointException( const string & aLabel, const ZFloat aValue ) : _label(aLabel), _value(aValue) {}
                /**
                 * Generate a string of explanation for this exception.
                 * \return A string describing the exception.
                 */
                const string whats();
            private:
                /**
                 * How many points were in the spline tables.
                 */
                string _label;
                ZFloat _value;
        };

    }} // namespaces

