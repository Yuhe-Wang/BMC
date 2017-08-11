/*!  \file   zeusVector.h
 *   \brief  A simple 3D vector class
 *   \author Iwan Kawrakow, Tony Apicella
 */
#pragma once

#include <cmath>
#include "zeusConfig.h"

namespace Zeus {

/*! \brief A 3D-Vector. */
class Vector {
public:

    double _x;  //!< x-coordinate
    double _y;  //!< y-coordinate
    double _z;  //!< z-coordinate

    /*! \brief Default constructor */
    Vector() : _x(0), _y(0), _z(0) {}

    /*! \brief Construct from input values \a aX, \a aY and \a aZ */
    Vector( double aX, double aY, double aZ ) : _x(aX), _y(aY), _z(aZ) {}

    /*! \brief Assignment operator */
    Vector &operator=(const Vector &v) { _x = v._x; _y = v._y; _z = v._z; return *this; };

    /*! \brief Add \a v to this vector and return the result */
    Vector operator+(const Vector &v) const { return Vector(_x+v._x, _y+v._y, _z+v._z); };

    /*! \brief Add vector \a v and assign the result to this vector */
    Vector &operator+=(const Vector &v) { _x += v._x; _y += v._y; _z += v._z; return *this; };

    /*! \brief Returns the scalar product of vector \a v and this vector */
    double operator*(const Vector &v) const { return _x*v._x + _y*v._y + _z*v._z; }; 

    /*! \brief Multiply all co-ordinates with \a f and assign the result to this vector */
    Vector &operator*=(ZFloat f) { _x *= f; _y *= f; _z *= f; return *this; };

    /*! \brief Multiply \a v with \a f and return the result */
    friend Vector operator*(ZFloat f, Vector &v) { return v*f; };

    /*! \brief Substract the vector \a v from this vector and return the result */
    Vector operator-(const Vector &v) const { return Vector(_x-v._x, _y-v._y, _z-v._z); };

    /*! \brief Substract \a v and assign the result to this vector */
    Vector &operator-=(const Vector &v) { _x -= v._x; _y -= v._y; _z -= v._z; return *this; };

    /*! \brief Return the multiplication of this vector with \a f */
    Vector operator*(ZFloat f) const { return Vector(_x*f,_y*f,_z*f); };

    /*! \brief Get the squared length of the vector */
    double getLengthSquared() const { return _x*_x + _y*_y + _z*_z; };

    /*! \brief Get the length of the vector */
    double getLength() const { return sqrt(getLengthSquared()); };

    /*! \brief Normalize the vector to unit length */
    void normalizeToUnitLength() { 
        double length2 = getLengthSquared(); 
        if( length2 ) { 
            double invLength = 1/sqrt(length2); 
            _x *= invLength; _y *= invLength; _z *= invLength;
        }
    };

    /*! \brief Rotate the vector by the polar angle Theta and the azimutha angle Phi */
    void rotate(ZFloat cosTheta, ZFloat sinTheta, ZFloat cosPhi, ZFloat sinPhi) {
        ZFloat sinz = _x*_x + _y*_y;
        if( sinz > 1e-20 ) {
            sinz = sqrt(sinz);
            ZFloat c1 = sinTheta/sinz;
            ZFloat cphi = _z*cosPhi;
            ZFloat cx = _x*cosTheta, cy = _y*cosTheta;
            ZFloat cx1 = cphi*_x-_y*sinPhi;
            ZFloat cy1 = cphi*_y+_x*sinPhi;
            _x = c1*cx1+cx; _y = c1*cy1+cy;
            _z = _z*cosTheta-sinz*sinTheta*cosPhi;
        }
        else { _x = sinTheta*cosPhi; _y = sinTheta*sinPhi; _z *= cosTheta; }
    };

};

}
