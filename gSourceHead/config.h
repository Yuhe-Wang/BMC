#ifndef _CONFIG_H_
#define _CONFIG_H_
//This a head file that can be included by both C++ and cu file, just like tools.h. So in this project, we should/can 
//include "tools.h"
//include "config.h"
//in very C++ and cu source file, as well as in gSourceHead.cu.h


#define USE_SINGLE_PRECISION

//notice GlueF will operate on float constant or float function like sin, sqrt
#ifdef USE_SINGLE_PRECISION
typedef float ZFloat;
#define GlueF(a) a##f
#else
typedef double ZFloat;
#define GlueF(a) a
#endif

class ParticleR //store status information of a particle
{
public:
	__device__ void changeByCos(ZFloat cosTheta, ZFloat phi)
	{
		ZFloat cosPhi = GlueF(cos)(phi);
		ZFloat sinPhi = GlueF(sqrt)(1 - cosPhi*cosPhi);
		if (phi > PI) sinPhi = -sinPhi;
		ZFloat dxy = u*u + v*v;
		ZFloat dxyz = dxy + w*w;
		if (GlueF(fabs)(dxyz - GlueF(1.0)) > GlueF(1e-14))
		{
			ZFloat norm = 1 / GlueF(sqrt)(dxyz);
			u *= norm;
			v *= norm;
			w *= norm;
			dxy = u*u + v*v;
		}

		if (dxy > GlueF(1.0e-28))
		{
			ZFloat sinTheta = GlueF(sqrt)((1 - cosTheta*cosTheta) / dxy);
			ZFloat up = u;
			u = u*cosTheta + sinTheta*(up*w*cosPhi - v*sinPhi);
			v = v*cosTheta + sinTheta*(v*w*cosPhi + up*sinPhi);
			w = w*cosTheta - dxy*sinTheta*cosPhi;
		}
		else
		{
			ZFloat sinTheta = GlueF(sqrt)(1 - cosTheta*cosTheta);
			if (w > 0)
			{
				u = sinTheta*cosPhi;
				v = sinTheta*sinPhi;
				w = cosTheta;
			}
			else
			{
				u = -sinTheta*cosPhi;
				v = -sinTheta*sinPhi;
				w = -cosTheta;
			}
		}
	}

	__device__ void copyPosition(ParticleR& r)
	{
		x = r.x; y = r.y; z = r.z;
		ivx = r.ivx; ivy = r.ivy; ivz = r.ivz; iabsv = r.iabsv;
	}
	__device__ void copyDirection(ParticleR& r)
	{
		u = r.u; v = r.v; w = r.w;
	}

	__device__ void normDirection()
	{
		ZFloat length2 = u*u + v*v + w*w;
		if (length2)
		{
			ZFloat invLength = 1 / GlueF(sqrt)(length2);
			u *= invLength; v *= invLength; w *= invLength;
		}
	}

	//data
	ZFloat x, y, z; //position, unit cm
	ZFloat u, v, w; //direction vector
	int ivx, ivy, ivz, iabsv; //voxel index for current particle

	ParticleType type;

	ZFloat E; //energy, unit eV
	ZFloat weight;
};

class ParticleStack
{
	//method
public:
	__device__ void push(ParticleR &par)
	{
		if (cur < 10)
		{
			ss[cur] = par; //copy and store the particle
			++cur;
		}
		else
		{
			printf("exceed the max depth of the GPU stack!!!!!!!!!!!!!!!\n");
		}
	}
	__device__ __forceinline__ bool empty()
	{
		if (cur > 0) return false;
		else return true;
	}
	__device__ __forceinline__ void pop()
	{
		if (cur <= 0) printf("no element to pop in the stack!");
		--cur;
	}
	__device__ __forceinline__ ParticleR& top(){ return ss[cur - 1]; }
	__device__ void init(ParticleR* pp){ cur = 0; ss = pp; }
	//data
private:
	ParticleR* ss;
	int cur; //point to next free space, 0 at initial time
};

class Vector {
public:

	ZFloat _x;  //!< x-coordinate
	ZFloat _y;  //!< y-coordinate
	ZFloat _z;  //!< z-coordinate

	/*! \brief Default constructor */
	__device__ __host__ Vector() : _x(0), _y(0), _z(0) { _x = 0; _y = 0; _z = 0; }

	/*! \brief Construct from input values \a aX, \a aY and \a aZ */
	__device__ __host__ Vector(ZFloat aX, ZFloat aY, ZFloat aZ) : _x(aX), _y(aY), _z(aZ) {}

	/*! \brief Assignment operator */
	__device__ __host__ Vector &operator=(const Vector &v) { _x = v._x; _y = v._y; _z = v._z; return *this; };

	/*! \brief Add \a v to this vector and return the result */
	__device__ __host__ Vector operator+(const Vector &v) const { return Vector(_x + v._x, _y + v._y, _z + v._z); };

	/*! \brief Add vector \a v and assign the result to this vector */
	__device__ __host__ Vector &operator+=(const Vector &v) { _x += v._x; _y += v._y; _z += v._z; return *this; };

	/*! \brief Returns the scalar product of vector \a v and this vector */
	__device__ __host__ ZFloat operator*(const Vector &v) const { return _x*v._x + _y*v._y + _z*v._z; };

	/*! \brief Multiply all co-ordinates with \a f and assign the result to this vector */
	__device__ __host__ Vector &operator*=(ZFloat f) { _x *= f; _y *= f; _z *= f; return *this; };

	/*! \brief Multiply \a v with \a f and return the result */
	__device__ __host__ friend Vector operator*(ZFloat f, Vector &v) { return v*f; };

	/*! \brief Substract the vector \a v from this vector and return the result */
	__device__ __host__ Vector operator-(const Vector &v) const { return Vector(_x - v._x, _y - v._y, _z - v._z); };

	/*! \brief Substract \a v and assign the result to this vector */
	__device__ __host__ Vector &operator-=(const Vector &v) { _x -= v._x; _y -= v._y; _z -= v._z; return *this; };

	/*! \brief Return the multiplication of this vector with \a f */
	__device__ __host__ Vector operator*(ZFloat f) const { return Vector(_x*f, _y*f, _z*f); };

	/*! \brief Get the squared length of the vector */
	__device__ __host__ ZFloat getLengthSquared() const { return _x*_x + _y*_y + _z*_z; };

	/*! \brief Get the length of the vector */
	__device__ __host__ ZFloat getLength() const { return GlueF(sqrt)(getLengthSquared()); };

	/*! \brief Normalize the vector to unit length */
	__device__ __host__ void normalizeToUnitLength() {
		ZFloat length2 = getLengthSquared();
		if (length2) {
			ZFloat invLength = 1 / GlueF(sqrt)(length2);
			_x *= invLength; _y *= invLength; _z *= invLength;
		}
	};

	/*! \brief Rotate the vector by the polar angle Theta and the azimutha angle Phi */
	__device__ __host__ void rotate(ZFloat cosTheta, ZFloat sinTheta, ZFloat cosPhi, ZFloat sinPhi) {
		ZFloat sinz = _x*_x + _y*_y;
		if (sinz > 1e-20) {
			sinz = GlueF(sqrt)(sinz);
			ZFloat c1 = sinTheta / sinz;
			ZFloat cphi = _z*cosPhi;
			ZFloat cx = _x*cosTheta, cy = _y*cosTheta;
			ZFloat cx1 = cphi*_x - _y*sinPhi;
			ZFloat cy1 = cphi*_y + _x*sinPhi;
			_x = c1*cx1 + cx; _y = c1*cy1 + cy;
			_z = _z*cosTheta - sinz*sinTheta*cosPhi;
		}
		else { _x = sinTheta*cosPhi; _y = sinTheta*sinPhi; _z *= cosTheta; }
	};
};


#endif