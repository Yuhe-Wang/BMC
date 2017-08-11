#ifndef _COMMON_H_
#define _COMMON_H_

//#define DEPOSIT_ENERGY_IN_PHOTON_INTERACTION 1
//#define USE_INTERPOLATED_LOG 1
#if __CUDA_ARCH__ < 350
#define USE_TEXTURE_CACHE 1
#endif

#define USE_SINGLE_PRECISION

//notice GlueF will operate on float constant or float function like sin, sqrt
#ifdef USE_SINGLE_PRECISION
typedef float ZFloat;
#define GlueF(a) a##f
#else
typedef double ZFloat;
#define GlueF(a) a
#endif

class Vector {
public:

	ZFloat x;  //!< x-coordinate
	ZFloat y;  //!< y-coordinate
	ZFloat z;  //!< z-coordinate

	/*! \brief Default constructor */
	Vector() : x(0), y(0), z(0) {}

	/*! \brief Construct from input values \a aX, \a aY and \a aZ */
	Vector(ZFloat aX, ZFloat aY, ZFloat aZ) : x(aX), y(aY), z(aZ) {}

	/*! \brief Assignment operator */
	Vector &operator=(const Vector &v) { x = v.x; y = v.y; z = v.z; return *this; };

	/*! \brief Add \a v to this vector and return the result */
	Vector operator+(const Vector &v) const { return Vector(x + v.x, y + v.y, z + v.z); };

	/*! \brief Add vector \a v and assign the result to this vector */
	Vector &operator+=(const Vector &v) { x += v.x; y += v.y; z += v.z; return *this; };

	/*! \brief Returns the scalar product of vector \a v and this vector */
	ZFloat operator*(const Vector &v) const { return x*v.x + y*v.y + z*v.z; };

	/*! \brief Returns the cross product of vector \a v and this vector */
	Vector operator&(const Vector &v) const { return Vector(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x); };

	/*! \brief Multiply all co-ordinates with \a f and assign the result to this vector */
	Vector &operator*=(ZFloat f) { x *= f; y *= f; z *= f; return *this; };

	/*! \brief Multiply \a v with \a f and return the result */
	friend Vector operator*(ZFloat f, Vector &v) { return v*f; };

	/*! \brief Subtract the vector \a v from this vector and return the result */
	Vector operator-(const Vector &v) const { return Vector(x - v.x, y - v.y, z - v.z); };

	/*! \brief Subtract \a v and assign the result to this vector */
	Vector &operator-=(const Vector &v) { x -= v.x; y -= v.y; z -= v.z; return *this; };

	/*! \brief Return the multiplication of this vector with \a f */
	Vector operator*(ZFloat f) const { return Vector(x*f, y*f, z*f); };

	/*! \brief Get the squared length of the vector */
	ZFloat getLengthSquared() const { return x*x + y*y + z*z; };

	/*! \brief Get the length of the vector */
	ZFloat getLength() const { return GlueF(sqrt)(getLengthSquared()); };

	/*! \brief Normalize the vector to unit length */
	void normalizeToUnitLength() {
		ZFloat length2 = getLengthSquared();
		if (length2) {
			ZFloat invLength = 1 / GlueF(sqrt)(length2);
			x *= invLength; y *= invLength; z *= invLength;
		}
	};

	/*! \brief Rotate the vector by the polar angle Theta and the azimutha angle Phi */
	void rotate(ZFloat cosTheta, ZFloat sinTheta, ZFloat cosPhi, ZFloat sinPhi) {
		ZFloat sinz = x*x + y*y;
		if (sinz > 1e-20) {
			sinz = GlueF(sqrt)(sinz);
			ZFloat c1 = sinTheta / sinz;
			ZFloat cphi = z*cosPhi;
			ZFloat cx = x*cosTheta, cy = y*cosTheta;
			ZFloat cx1 = cphi*x - y*sinPhi;
			ZFloat cy1 = cphi*y + x*sinPhi;
			x = c1*cx1 + cx; y = c1*cy1 + cy;
			z = z*cosTheta - sinz*sinTheta*cosPhi;
		}
		else { x = sinTheta*cosPhi; y = sinTheta*sinPhi; z *= cosTheta; }
	};
};

struct StoreFace
{
	int iadj; //adjacent tetrahedral. if it's -1, means this face is the boundary
	Vector v[3];
};
struct StoreTet
{
	StoreFace face[4];
	short mid;
	SFloat density;
};
#endif