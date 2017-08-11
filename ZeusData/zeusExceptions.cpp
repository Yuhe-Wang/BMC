
#include "zeusExceptions.h"

namespace Zeus { namespace Exceptions {

    ComputeException::ComputeException( const string & aCrossSectionName, const int aMaterialId, const int aIndex, const int aLowerBound, 
            const int aUpperBound ) : _materialId(aMaterialId), _crossSectionClass( aCrossSectionName ), _index( aIndex ), 
            _lowerBound( aLowerBound ), _upperBound( aUpperBound ) {}

    const string ComputeException::whats() {                    
        std::ostringstream streamer; streamer.str("");
        streamer << "Out of bounds in class " << _crossSectionClass << string(" for material ID ") << _materialId << string( " index=" ) << _index << string( " lower bound=" ) << _lowerBound << string(" upper bound=") << _upperBound << endl;
        return streamer.str();
    }

    ZeroDenominatorException::ZeroDenominatorException( const string & aCrossSectionName, int aMaterialId ) : 
        _materialId(aMaterialId), _crossSectionClass( aCrossSectionName ) {}

    const string ZeroDenominatorException::whats() {                    
        std::ostringstream streamer; streamer.str("");
        streamer << "Zero divide in class " << _crossSectionClass << string(" for material ID ") << _materialId << endl;
        return streamer.str();
    }

    FileException::FileException( const string & aCrossSectionName ) : _crossSectionClass( aCrossSectionName ) {}

    const string FileException::whats() {   
        std::ostringstream streamer; streamer.str("");
        streamer << "Cant open file " << _crossSectionClass << std::ends;
        return streamer.str();
    }

    CalcElectronDirectionException::CalcElectronDirectionException( const double aCosTheta ) : _theta(aCosTheta) {}

    const string CalcElectronDirectionException::whats() {                    
        std::ostringstream streamer; streamer.str("");
        streamer << "CalcElectronDirectionException theta is " << _theta << endl;
        return streamer.str();
    }

    StackException::StackException( ) {}

    const string StackException::whats() {                    
        std::ostringstream streamer; streamer.str("");
        streamer << "Stack overflow" << endl;
        return streamer.str();
    }

    NotEnoughPointsSplineException::NotEnoughPointsSplineException( const int aNpts ) : _npts(aNpts) {}

    const string NotEnoughPointsSplineException::whats() {                    
        std::ostringstream streamer; streamer.str("");
        streamer << "Spline tables only had " << _npts << " points.  Must be 4 or more. " << endl;
        return streamer.str();
    }

    MonotonicSplineException::MonotonicSplineException( ) {}

    const string MonotonicSplineException::whats() {                    
        std::ostringstream streamer; streamer.str("");
        streamer << "Monotonic spline exception" << endl;
        return streamer.str();
    }

    LogLinLinException::LogLinLinException(  const double aInterpLoc, const double aLowerBound, const double aUpperBound ) : _interpLoc(aInterpLoc), _lowerBound( aLowerBound ), _upperBound( aUpperBound ) {}

    const string LogLinLinException::whats() {                    
        std::ostringstream streamer; streamer.str("");
        streamer << "Log lin lin interpolation exception x value " << _interpLoc << " must be within range [" << _lowerBound << _upperBound << "]" << endl;
        return streamer.str();
    }

    const string UnknownMaterialException::whats( ) {
        std::ostringstream streamer; streamer.str("");
        streamer << "Unknown material exception "; 
        streamer << _label << " ";
        streamer << "Current region " << _curRegion;
        streamer << "Next region " << _nextRegion;
        streamer << "Intended distance " << _t;
        streamer << "Particle position X=" << _position._x << " Y=" << _position._y << " Z=" << _position._z;
        streamer << "Particle direction X=" << _direction._x << " Y=" << _direction._y << " Z=" << _direction._z;
        streamer << endl;
        return streamer.str();
    }

    const string FloatingPointException::whats( ) {
        std::ostringstream streamer; streamer.str("");
        streamer << "Floating point exception "; 
        streamer << _label << " value is " << _value;
        streamer << endl;
        return streamer.str();
    }

}} // namespaces
