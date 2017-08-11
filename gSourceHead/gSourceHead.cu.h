#include "gZeus.h"
#include "gSourceHead.h"

/******************************************************Start: GPU memory******************************************************/
__device__ CoPhaseSpace4GPU gps;
__constant__ int nSeg;
__constant__ SEGMENT4GPU* segs;
__constant__ ParticleR* InitPars;
__constant__ int NMaxInitPars;
__constant__ int* NInitPars;

class HeadTransport;
__constant__ HeadTransport* headTransport;
//extern __constant__ GRNG* RNGState;

//define constants that used in head model simulation
/**
* From microsoft, smallest number such that 1.0+FLTEPSILON != 1.0. CAUTION!!!!!!!!
*/
#define FLTEPSILON GlueF(1.192092896e-07)
/**
* Number of leaf pairs in the MLC. Only the middle 30 are real.
*/
#define NLEAF_PAIRS_REAL 30
/**
* Number of leaf pairs, on each side, which are not real.
*/
#define N_LEAF_PAIRS_FAKE 4
/**
* Number of leaf pairs in the MLC. Only the middle 30 are real.
*/
#define NLEAF_PAIRS_TOTAL 38
/**
* Distance in CM from source to isocenter.
*/
#define SOURCE_TO_ISO_DISTANCE GlueF(105.0)
/**
* Leaf thickness in CM.
*/
#define LEAF_THICKNESS GlueF(1.05)
/**
* Leaf length in CM.
*/
#define LEAF_LENGTH GlueF(100.0)
/**
* MLC model, the radius of the cylinder representing the inner leaves.
*/
#define MLC_INNER_RADIUS GlueF(41.0)
/**
* MLC model, the radius of the cylinder representing the inner leaves.
*/
#define MLC_OUTER_RADIUS GlueF(50.0)
/**
* The middle of the MLC, in the Y dimension.
*/
#define MLC_MIDDLE ((MLC_INNER_RADIUS + MLC_OUTER_RADIUS) / GlueF(2.0))
/**
* Minimum X location of MLC.
*/
#define MLC_MINX GlueF(-15.0)
/**
* Maximum X location of MLC.
*/
#define MLC_MAXX GlueF(15.0)
/**
* MLC design constant
*/
#define INNER_RIGHT_FOCUS GlueF(0.03)
/**
* MLC design constant
*/
#define INNER_GAP_FOCUS GlueF(-0.035)
/**
* MLC design constant
*/
#define INNER_Z_SHIFT GlueF(-0.015)
/**
* MLC design constant
*/
#define OUTER_LEFT_FOCUS GlueF(-0.03)
/**
* MLC design constant
*/
#define OUTER_GAP_FOCUS GlueF(0.015)
/**
* MLC design constant
*/
#define OUTER_Z_SHIFT GlueF(0.035)

#define HEAD_MAX_Z GlueF(70.0)

#define HEAD_NO_SCATTER 0 //whether to simulate scattering in head

#define HEAD_N_SPLIT 200

#define HEAD_pCut 0.05

#define HEAD_XMIN -17.5
#define HEAD_XMAX 17.5
#define HEAD_YMIN -17.5
#define HEAD_YMAX 17.5
//end constants definition

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



class Plane {
public:

	/*!
	* \brief Query x component of plane point.
	*
	* Return x component of point used to define this plane.
	*/
	ZFloat X() const { return _point[0]; }
	/*!
	* \brief Query y component of plane point.
	*
	* Return y component of point used to define this plane.
	*/
	ZFloat Y() const { return _point[1]; }
	/*!
	* \brief Query z component of plane point.
	*
	* Return z component of point used to define this plane.
	*/
	ZFloat Z() const { return _point[2]; }

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
	__device__ int nextRegionAlongPath(const Vector &origin, const Vector &direction, ZFloat &distance, int& aNewRegion) const {
		ZFloat d = distance;
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
	__device__ bool hit(const Vector &origin, const Vector &direction, ZFloat &distance) const {

		Vector C(_point[0], _point[1], _point[2]);
		Vector N(_norm[0], _norm[1], _norm[2]);

		Vector X(origin._x, origin._y, origin._z);
		Vector U(direction._x, direction._y, direction._z);

		ZFloat numer = (C._x - X._x) * N._x + (C._y - X._y) * N._y + (C._z - X._z) * N._z;
		ZFloat denom = U._x * N._x + U._y * N._y + U._z * N._z;

		if (GlueF(fabs)(denom)>GlueF(1e-8)) {
			if (GlueF(fabs)(numer)>GlueF(1e-8)) {
				ZFloat t = numer / denom;
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
	ZFloat _point[3];

	/*!
	* \brief Store the normal to this plane.
	*/
	ZFloat _norm[3];

	/*!
	* \brief Return this code if the plane will be hit.
	*/
	int _regionCode;

};

class Cylinder {
public:

	/*! \brief Construct a cylinder with radius \a aRadius and axis passing through \a aCenter */
	__device__ __host__ void init(const ZFloat aRadius, const Vector & aCenter) {
		_radius = aRadius;
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
	__device__ __host__ void setRegionCode(const int aRegionCode) { _regionCode = aRegionCode; }

	/*! \brief Get the region code of this cylinder */
	__device__ __host__ int getRegionCode() const { return _regionCode; }

	/*! \brief Construct a dummy cylinder */
	__device__ __host__ Cylinder() {
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
	__device__ bool nextRegionAlongPath(const Vector & aX, const Vector & aU, ZFloat& aDistance, bool aIsInside) const
	{

		ZFloat xs = aX._x - _center._x, zs = aX._z - _center._z;
		ZFloat ux = aU._x, uz = aU._z;

		//
		// The quadratic equation to solve is A*t^2 + 2*B*t - C, with A, B, C defined below. The discriminant is D = B^2 + C and the solutions, 
		// if they exist, are (-B +/- sqrt(D))/A.
		// There are no solutions when
		//   a) A = 0 (trajectory is parallel to cylinder axis
		//   b) D < 0 Strictly speaking, this can only happen if aIsInside is false (C < 0 in this case)
		//   c) aIsInside is false and B >= 0
		//
		ZFloat A = ux*ux + uz*uz;
		if (A < GlueF(1e-10)) return false; // traveling parallel to cylinder axis
		ZFloat B = ux * xs + uz * zs;
		if (!aIsInside && B >= 0) return false;
		ZFloat C = _radiusSquared - xs*xs - zs*zs;
		ZFloat D = B*B + A*C;
		if (D < 0) return false;

		if (aIsInside) {
			// when we are inside, C >= 0 (except for numerical roundoff), so the discriminant is always > |B|. 
			// Thus, the positive solution is (-B + sqrt(D))/A. To improve precision when B > 0, we multiply and divide the solution by (sqrt(D) + B), which 
			// then simplifies to C/(sqrt(D) + B).
			// To avoid the relatively expensive sqrt and division evaluations, we check if A*aDistance^2 + 2*B*aDistance - C < 0, and if so, simply return false.
			//
			if (aDistance*(A*aDistance + 2 * B) < C) return false;
			ZFloat t = B > 0 ? C / (GlueF(sqrt)(D)+B) : (GlueF(sqrt)(D) - B) / A;
			aDistance = t; return true;
		}
		// If here, we are outside and also B < 0. In that case, the solution (first intersection) is (-B - sqrt(D))/A = -(sqrt(D) + B)/A. To improve precision, 
		// we multiply and divide with (sqrt(D) - B), and the solution simplifies to -C/(sqrt(D) - B)
		ZFloat t = -C / (GlueF(sqrt)(D) - B);
		if (t <= aDistance) { aDistance = t; return true; }
		return false;
	}

	/*! \brief Returns true if the position \a aX is inside the cylinder, false otherwise. */
	bool isInside(const Vector &aX) const;

	/*! \brief Returns the radius of this cylinder. */
	__device__ __host__ ZFloat getRadius() const { return _radius; }

	/*! \brief Returns the radius squared of this cylinder. */
	__device__ __host__ ZFloat getRadiusSquared() const { return _radiusSquared; }

private:
	Vector _center;             //!< A point on the axis of this cylinder
	ZFloat  _radius;             //!< The radius of this cylinder
	ZFloat  _radiusSquared;      //!< The radius squared of this cylinder
	int       _regionCode;         //!< The region code of this cylinder
};

/****************************************** Start: Model of MLC *************************************/
class ViewRayLeaf {
public:
	__device__ int isWhere(const Vector &aLocation, Medium & aMedium) const {
		// Air gap?
		ZFloat normDotLoc = _airGapNormal._x*aLocation._x + _airGapNormal._y*aLocation._y + _airGapNormal._z*aLocation._z;
		if (normDotLoc > _airGapPosition) {
			aMedium = Air;
			return 3;
		}

		ZFloat x1 = _leafNormals[1] * aLocation;
		ZFloat x2 = _leafNormals[2] * aLocation;

		if (x1 > _leafPositions[1] && x2 < _leafPositions[2]) {
			aMedium = Air;
			return (1);
		}
		else if (x1 <= _leafPositions[1]) {
			aMedium = Tungsten;
			return 0;
		}
		else {
			aMedium = Tungsten;
			return (2);
		}
	}
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
	__device__ int isWhereNoAirGap(const Vector &aLocation, Medium & aMedium) const {

		ZFloat x1 = _leafNormals[1] * aLocation;
		ZFloat x2 = _leafNormals[2] * aLocation;

		if (x1 > _leafPositions[1] && x2 < _leafPositions[2]) {
			aMedium = Air;
			return (1);
		}
		else if (x1 <= _leafPositions[1]) {
			aMedium = Tungsten;
			return 0;
		}
		else {
			aMedium = Tungsten;
			return (2);
		}

	}
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
	__device__ int nextRegionAlongPath(int ireg, const Vector &x, const Vector &u, ZFloat &t, Medium& newmed)  const
	{
		if (ireg >= 0) {
			ZFloat up = _leafNormals[ireg] * u;
			int inew = ireg;
			if (up < -FLTEPSILON) {
				ZFloat xp = _leafNormals[ireg] * x;
				if (xp < _leafPositions[ireg] - up*t) {
					ZFloat tt = (_leafPositions[ireg] - xp) / up; inew = ireg - 1;
					if (tt > t) printf("Error: Huh1: %g %g", t, tt);
					t = tt;
					if (inew >= 0 && newmed) newmed = getMedium(inew);
				}
			}
			up = _leafNormals[ireg + 1] * u;
			if (up > FLTEPSILON) {
				ZFloat xp = _leafNormals[ireg + 1] * x;
				if (xp > _leafPositions[ireg + 1] - up*t) {
					ZFloat tt = (_leafPositions[ireg + 1] - xp) / up; inew = ireg + 1;
					if (tt > t) printf("Error: Huh2: %g %g", t, tt);
					t = tt;
					if (inew >= 3) inew = -1;
					if (inew >= 0 && newmed) newmed = getMedium(inew);
				}
			}
			return inew;
		}

		ZFloat xp, up; int inew = ireg;
		xp = _leafNormals[0] * x;
		if (xp < _leafPositions[0]) {
			up = _leafNormals[0] * u;
			if (up > FLTEPSILON) {
				ZFloat tt = (_leafPositions[0] - xp) / up;
				if (tt < t) {
					t = tt; inew = 0;
					if (newmed) newmed = Tungsten;
				}
			}
			return inew;
		}
		xp = _leafNormals[3] * x;
		if (xp > _leafPositions[3]) {
			up = _leafNormals[3] * u;
			if (up < -FLTEPSILON) {
				ZFloat tt = (_leafPositions[3] - xp) / up;
				if (tt < t) {
					t = tt; inew = 3 - 1;
					if (newmed) newmed = Tungsten;
				}
			}
			return inew;
		}
		return inew;
	}
	/**
	* \brief Return the medium this region consists of.
	*/
	__device__ Medium getMedium(int iregion) const {
		if (iregion == 0 || iregion == 2) return Tungsten;
		return Air;
	}

	void setOpening(const pair<double, double>& opening, int ileaf, Vector& leftFocus, Vector& rightFocus, Vector& gapFocus)
	{
		ZFloat z = SOURCE_TO_ISO_DISTANCE;
		ZFloat y = (ileaf + 1 - NLEAF_PAIRS_TOTAL / 2)*LEAF_THICKNESS;
		ZFloat aux = 1/ GlueF(sqrt)(y*y + z*z);
		_airGapNormal = Vector(0, z*aux, y*aux);
		_airGapPosition = gapFocus * _airGapNormal;

		ZFloat leftPos = (ZFloat)opening.first;
		ZFloat rightPos = (ZFloat)opening.second;
		ZFloat minX = leftPos - LEAF_LENGTH;
		ZFloat maxX = rightPos + LEAF_LENGTH;
		Vector a;
		ZFloat zFocus = SOURCE_TO_ISO_DISTANCE;
		a = Vector(zFocus, 0, minX);     a.normalizeToUnitLength(); _leafNormals[0] = a; _leafPositions[0] = a*leftFocus;
		a = Vector(zFocus, 0, leftPos);  a.normalizeToUnitLength(); _leafNormals[1] = a; _leafPositions[1] = a*leftFocus;
		a = Vector(zFocus, 0, rightPos); a.normalizeToUnitLength(); _leafNormals[2] = a; _leafPositions[2] = a*rightFocus;
		a = Vector(zFocus, 0, maxX);     a.normalizeToUnitLength(); _leafNormals[3] = a; _leafPositions[3] = a*rightFocus;
	}
	/**
	* \brief Normal vector to the air gap plane.
	*/
	Vector _airGapNormal;
	/**
	* \brief Position of the air gap plane.
	*/
	ZFloat _airGapPosition;
	/**
	* The normals of the four planes that represent the leaf.
	*/
	Vector _leafNormals[4];
	/**
	* The positions of the four planes that represent the leaf.
	*/
	ZFloat _leafPositions[4];
};

class ViewRayLeaves {
public:	
	void setSegment(const vector<pair<double, double> >  &aSegment, Vector& leftFocus, Vector& rightFocus, Vector& gapFocus)
	{
		// Set the positions of the real leaves.
		for (int leafPair = 0; leafPair < NLEAF_PAIRS_REAL; ++leafPair) {
			_leaf[leafPair + 4].setOpening(aSegment[leafPair], leafPair + 4, leftFocus, rightFocus, gapFocus);
		}

		// Set the fake leaves to closed.
		pair<double, double> outPositions;
		outPositions = pair<double, double>(-20, -20);
		for (int leafPair = 0; leafPair < N_LEAF_PAIRS_FAKE; ++leafPair) {
			_leaf[leafPair].setOpening(outPositions, leafPair, leftFocus, rightFocus, gapFocus);
			_leaf[leafPair + NLEAF_PAIRS_REAL + N_LEAF_PAIRS_FAKE].setOpening(outPositions, leafPair + NLEAF_PAIRS_REAL + N_LEAF_PAIRS_FAKE, leftFocus, rightFocus, gapFocus);
		}
	}
	ViewRayLeaf _leaf[NLEAF_PAIRS_TOTAL];
};

class ViewRayGaps {
public:
	void init(Vector & aGapFocus, ZFloat aZshift)
	{
		_focus = Vector(0, aZshift, 0);
		Vector focus(0, aZshift, SOURCE_TO_ISO_DISTANCE);

		ZFloat yloc = 0, zloc = 0, aux = 0;
		int n = NLEAF_PAIRS_TOTAL;
		ZFloat farLeft = -(n / 2)*LEAF_THICKNESS;

		int ii = 0;
		for (int y = 0; y <= NLEAF_PAIRS_TOTAL; y++) {
			zloc = SOURCE_TO_ISO_DISTANCE;
			yloc = farLeft + ZFloat(y)*LEAF_THICKNESS;
			aux = 1 / GlueF(sqrt)(yloc*yloc + zloc*zloc);
			_airGapNormal[ii] = Vector(0, zloc*aux, yloc*aux);
			_airGapPosition[ii] = (focus * _airGapNormal[ii]);
			ii++;
			if (y != NLEAF_PAIRS_TOTAL) {
				yloc += LEAF_THICKNESS;
				aux = 1 / GlueF(sqrt)(yloc*yloc + zloc*zloc);
				_airGapNormal[ii] = Vector(0, zloc*aux, yloc*aux);
				_airGapPosition[ii] = aGapFocus * _airGapNormal[ii];
				ii++;
			}
		}
	}

	__device__ int isWhere(const Vector &x) const {
		ZFloat dz = SOURCE_TO_ISO_DISTANCE - x._z; ZFloat y = x._y - _focus._y;
		ZFloat minZ = -(NLEAF_PAIRS_TOTAL / 2)*LEAF_THICKNESS;
		int region = (dz <= 0 || 105 * y < dz*minZ || 105 * y >= -dz*minZ) ? -1 :
			(int)((105 * y - dz*minZ) / (dz*LEAF_THICKNESS));
		if (region >= 0) {
			region *= 2;
			if (region >= 76) { return -1; }
			if (_airGapNormal[region + 1] * x > _airGapPosition[region + 1]) ++region;
		}
		return region;
	}

	__device__ Medium getMedium(int aRegion) const { int rmod4 = aRegion % 4; return rmod4 == 1 || rmod4 == 3 ? Air : Tungsten; }

	__device__ int nextRegionAlongPath(int ireg, const Vector &x, const Vector &u, ZFloat &t, Medium& newmed) const {

		if (ireg >= 0) {
			ZFloat up = _airGapNormal[ireg] * u;
			int inew = ireg;
			if (up < -FLTEPSILON) {
				ZFloat xp = _airGapNormal[ireg] * x;
				if (xp < _airGapPosition[ireg] - up*t) {
					ZFloat tt = (_airGapPosition[ireg] - xp) / up; inew = ireg - 1;
					if (tt > t)  printf("Error: Huh1: t=%g tt=%g", t, tt);
					t = tt;
					if (inew >= 0 && newmed) newmed = getMedium(inew);
				}
			}
			up = _airGapNormal[ireg + 1] * u;
			if (up > FLTEPSILON) {
				ZFloat xp = _airGapNormal[ireg + 1] * x;
				if (xp > _airGapPosition[ireg + 1] - up*t) {
					ZFloat tt = (_airGapPosition[ireg + 1] - xp) / up; inew = ireg + 1;
					if (tt > t)  printf("Error: Huh2: t=%g tt=%g", t, tt);
					t = tt;
					if (inew >= 76) inew = -1;
					if (inew >= 0 && newmed) newmed = getMedium(inew);
				}
			}
			return inew;
		}

		ZFloat xp, up; int inew = ireg;
		xp = _airGapNormal[0] * x;
		if (xp < _airGapPosition[0]) {
			up = _airGapNormal[0] * u;
			if (up > FLTEPSILON) {
				ZFloat tt = (_airGapPosition[0] - xp) / up;
				if (tt < t) {
					t = tt; inew = 0;
					if (newmed) newmed = getMedium(0);
				}
			}
			return inew;
		}
		xp = _airGapNormal[76] * x;
		if (xp > _airGapPosition[76]) {
			up = _airGapNormal[76] * u;
			if (up < -FLTEPSILON) {
				ZFloat tt = (_airGapPosition[76] - xp) / up;
				if (tt < t) {
					t = tt; inew = 76 - 1;
					if (newmed) newmed = getMedium(76 - 1);
				}
			}
			return inew;
		}
		return inew;
	}

private:

	/**
	* \brief Normal of the plane representing the air gap.
	*/
	Vector _airGapNormal[2*NLEAF_PAIRS_TOTAL+1];
	/**
	* \brief Position of the plane representing the air gap.
	*/
	ZFloat _airGapPosition[2*NLEAF_PAIRS_TOTAL+1];

	/**
	* \brief The position of the focal point
	*/
	Vector _focus;
};

class ViewRayMlcCylinders {
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
	* Sets the radii of the cylinders modeling the MLC.
	* There are three cylinders, separating the upper and lower
	* leaf halves.
	*
	*/
	ViewRayMlcCylinders()  {
		ZFloat radii[3];
		radii[0] = ZFloat(MLC_INNER_RADIUS);
		radii[1] = ZFloat(MLC_MIDDLE);
		radii[2] = ZFloat(MLC_OUTER_RADIUS);
		for (int i = 0; i < 3; ++i) _radiiSquared[i] = radii[i] * radii[i];
	}

	__device__ int isWhere(const Vector &aLocation, Medium & aMedium) const {

		ZFloat z = aLocation._z - 105;
		ZFloat x = aLocation._x;
		ZFloat d = x*x + z*z;

		aMedium = Air;
		if (d < _radiiSquared[0]) return 0;

		aMedium = Tungsten;
		if (d < _radiiSquared[1]) return 1;
		if (d < _radiiSquared[2]) return 2;

		aMedium = UnknownMaterial;
		return -1;

	}
	// \return Index of next region along path of particle.
	__device__ int nextRegionAlongPath(int aCurrentRegion, const Vector &aX, const Vector &aU, ZFloat &aIntendedDistance, Medium& nextMedium) const {

		Vector xp(aX - Vector(0,0,105));
		ZFloat A, B, C;
		A = aU._x*aU._x + aU._z*aU._z;
		B = xp._x*aU._x + xp._z*aU._z;
		C = xp._x*xp._x + xp._z*xp._z;
		if (A < GlueF(1e-10)) return aCurrentRegion;  // traveling parallel to axis.

		int nextRegion = aCurrentRegion;
		if (aCurrentRegion >= 0) {
			if (B >= 0 || !aCurrentRegion) {  // in the innermost cylinder or traveling outwards
				ZFloat C1 = _radiiSquared[aCurrentRegion] - C;
				if (A*aIntendedDistance*aIntendedDistance + 2 * B*aIntendedDistance >= C1) {
					ZFloat D = B*B + A*C1;
					// if ( D < 0 ) D = 0;
					aIntendedDistance = B > 0 ? C1 / (GlueF(sqrt)(D)+B) : (GlueF(sqrt)(D)-B) / A;
					++nextRegion;
					if (nextRegion >= 3) nextRegion = -1;
				}
			}
			else {
				// traveling inwards
				// -> first check inner cylinder
				ZFloat C1 = _radiiSquared[aCurrentRegion - 1] - C;
				ZFloat D = B*B + A*C1;
				if (D > 0) {
					ZFloat distanceToInner = -C1 / (GlueF(sqrt)(D)-B);
					if (distanceToInner < aIntendedDistance) {
						aIntendedDistance = distanceToInner; --nextRegion;
					}
				}
				else {
					// missing inner cylinder
					C1 = _radiiSquared[aCurrentRegion] - C;
					if (A*aIntendedDistance*aIntendedDistance + 2 * B*aIntendedDistance >= C1) {
						aIntendedDistance = (GlueF(sqrt)(GlueF(fabs)(D)) - B) / A;
						++nextRegion;
						if (nextRegion >= 3) nextRegion = -1;
					}
				}
			}
		}
		else if (B < 0) {
			// if outside, only need to check for intersection if B<0
			ZFloat C1 = _radiiSquared[2] - C;
			ZFloat D = B*B + A*C1;
			if (D > 0) {
				ZFloat distanceToInner = -C1 / (GlueF(sqrt)(D)-B);
				//ZFloat distanceToInner = (sqrt(D) - B)/A;
				if (distanceToInner < aIntendedDistance) {
					aIntendedDistance = distanceToInner; nextRegion = 2;
				}
			}
		}

		if (nextRegion != aCurrentRegion) {
			if (nextMedium && nextRegion >= 0) nextMedium = getMedium(nextRegion);
		}

		return nextRegion;
	}
	
	__device__ Medium getMedium(int aRegion) const {
		if (aRegion == 1 || aRegion == 2) return Tungsten;
		return UnknownMaterial;
	}

private:

	/**
	* \brief The square of the cylinder radii.
	*
	*/
	ZFloat _radiiSquared[3];
};

class ViewRayMLC {
	/**
	* \brief A class to model the ViewRay MLC.
	*
	* A leaf is modeled with an inner and outer half, so there is one model for each leaf for each half.
	*/
public:
	ViewRayMLC()
	{
		// init the MLC gaps
		Vector innerGapFocus(0, INNER_GAP_FOCUS, SOURCE_TO_ISO_DISTANCE);
		_gaps[0].init(innerGapFocus, INNER_Z_SHIFT);
		Vector outerGapFocus(0, OUTER_GAP_FOCUS, SOURCE_TO_ISO_DISTANCE);
		_gaps[1].init(outerGapFocus, OUTER_Z_SHIFT);

		//These must be set to GPU memory pointer later
		_cylinderOfLeaves[0] = NULL;
		_cylinderOfLeaves[1] = NULL;
	}
	void setSegments(const vector<vector<pair<double, double> > > &aListOfSegments, GPUConfig& gc)
	{
		int nSeg = (int)aListOfSegments.size();
		//inner leaves in CPU memory
		Vector innerLeftFocus(0.0, 0.0, SOURCE_TO_ISO_DISTANCE);
		Vector innerRightFocus(INNER_RIGHT_FOCUS, 0.0, SOURCE_TO_ISO_DISTANCE);
		Vector innerGapFocus(0.0, INNER_GAP_FOCUS, SOURCE_TO_ISO_DISTANCE);

		ViewRayLeaves* h_inner = new ViewRayLeaves[nSeg];
		for (int i = 0; i < nSeg; ++i)
			h_inner[i].setSegment(aListOfSegments[i], innerLeftFocus, innerRightFocus, innerGapFocus);

		//outer leaves in CPU memory
		Vector outerLeftFocus(OUTER_LEFT_FOCUS, 0.0, SOURCE_TO_ISO_DISTANCE);
		Vector outerRightFocus(0.0, 0.0, SOURCE_TO_ISO_DISTANCE);
		Vector outerGapFocus(0.0, OUTER_GAP_FOCUS, SOURCE_TO_ISO_DISTANCE);

		ViewRayLeaves* h_outer = new ViewRayLeaves[nSeg];
		for (int i = 0; i < nSeg; ++i)
			h_outer[i].setSegment(aListOfSegments[i], outerLeftFocus, outerRightFocus, outerGapFocus);

		//allocate memory in GPU
		cudaSetDevice(gc.id);
		//set inner MLC
		cudaMalloc(&_cylinderOfLeaves[0], nSeg*sizeof(ViewRayLeaves));
		gc.addGPUPointer(_cylinderOfLeaves[0]); //remember to free when exiting
		cudaMemcpy(_cylinderOfLeaves[0], h_inner, nSeg*sizeof(ViewRayLeaves), cudaMemcpyHostToDevice);

		//set outer MLC
		cudaMalloc(&_cylinderOfLeaves[1], nSeg*sizeof(ViewRayLeaves));
		gc.addGPUPointer(_cylinderOfLeaves[1]); //remember to free when exiting
		cudaMemcpy(_cylinderOfLeaves[1], h_outer, nSeg*sizeof(ViewRayLeaves), cudaMemcpyHostToDevice);
		//release the CPU memory
		delete[] h_inner;
		delete[] h_outer;
	}

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
	__device__ int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, ZFloat &aIntendedDistance, Medium& material) const {
		if (aCurrentRegion >= 0) {
			int irCyl = aCurrentRegion / (4 * NLEAF_PAIRS_TOTAL);
			int irCylNew = _boundingCylinders.nextRegionAlongPath(irCyl + 1, aX, aU, aIntendedDistance, material);
			int leafRegion = aCurrentRegion - 4 * NLEAF_PAIRS_TOTAL * irCyl;
			int leafPair = leafRegion / 4;
			int pairLocalRegion = leafRegion - 4 * leafPair;
			int zplane = 2 * leafPair;
			if (pairLocalRegion < 3) {
				// we are in the area divided by the leafs
				int newPairLocalRegion = _cylinderOfLeaves[irCyl][aTimeIndex]._leaf[leafPair].nextRegionAlongPath(pairLocalRegion, aX, aU, aIntendedDistance, material);
				int newZplane = _gaps[irCyl].nextRegionAlongPath(zplane, aX, aU, aIntendedDistance, material);
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
				int newZplane = _gaps[irCyl].nextRegionAlongPath(zplane, aX, aU, aIntendedDistance, material);
				if (newZplane < 0) return -1;
				if (newZplane != zplane) {
					int newLeafPair = newZplane > zplane ? leafPair + 1 : leafPair;
					int newPairLocalRegion = _cylinderOfLeaves[irCyl][aTimeIndex]._leaf[newLeafPair].isWhereNoAirGap(aX + aU*aIntendedDistance, material);
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
			zplane = _gaps[irCyl].isWhere(xnew);
			if (zplane < 0) return -1;
			leafPair = zplane / 2;
			if (leafPair * 2 == zplane) {
				// we have entered an area occupied by the leaves
				int newPairLocalRegion = _cylinderOfLeaves[irCyl][aTimeIndex]._leaf[leafPair].isWhere(xnew, material);
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
		int irCyl = _boundingCylinders.isWhere(aX, material);
		if (irCyl == 0 || irCyl == -1) {
			if ((irCyl == 0 && aU._z >= 0) || (irCyl == -1 && aU._z <= 0)) return -1;
			ZFloat tsave = aIntendedDistance; ZFloat tleft = aIntendedDistance; Vector xnew(aX);
			while (1) {
				ZFloat tt = tleft;
				int irCylNew = _boundingCylinders.nextRegionAlongPath(irCyl, xnew, aU, tt, material);
				if (tt != tt) {
					printf("Fatal error: Outside MLC cylinders");
				}
				if (irCylNew < 0) return -1;
				if (irCylNew == irCyl) return aCurrentRegion;
				xnew += aU*tt; tleft -= tt; irCyl = irCylNew;
				if (xnew._z > HEAD_MAX_Z) return -1;
				if (irCyl == 1 || irCyl == 2) {
					int zplane = _gaps[irCyl - 1].isWhere(xnew);
					if (zplane >= 0) {
						aIntendedDistance = tsave - tleft;
						int leafPair = zplane / 2;
						if (leafPair * 2 == zplane) {
							int newPairLocalRegion = _cylinderOfLeaves[irCyl - 1][aTimeIndex]._leaf[leafPair].isWhere(xnew, material);
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
				ZFloat tt = tleft;
				int zplane = _gaps[irCyl - 1].nextRegionAlongPath(-1, xnew, aU, tt, material);
				if (tt != tt) {
					printf("Fatal error:Outside MLC second while gaps");
				}
				int irCylNew = _boundingCylinders.nextRegionAlongPath(irCyl, xnew, aU, tt, material);
				if (tt != tt) {
					printf("Fatal error:Outside MLC second while cylinders");
				}
				if (irCylNew != irCyl) {
					if (irCylNew == 0 || irCylNew < 0) return -1;
					xnew += aU*tt; tleft -= tt; irCyl = irCylNew;
					if (xnew._z > HEAD_MAX_Z) return -1;
					zplane = _gaps[irCyl - 1].isWhere(xnew);
					if (zplane >= 0) {
						aIntendedDistance = tsave - tleft;
						int leafPair = zplane / 2;
						if (leafPair * 2 == zplane) {
							int newPairLocalRegion = _cylinderOfLeaves[irCyl - 1][aTimeIndex]._leaf[leafPair].isWhere(xnew, material);
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
						int newPairLocalRegion = _cylinderOfLeaves[irCyl - 1][aTimeIndex]._leaf[leafPair].isWhere(xnew, material);
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
		ZFloat r2 = xp._x*xp._z + xp._y*xp._z;
		if ((r2 < MLC_INNER_RADIUS * MLC_INNER_RADIUS && aU._z > 0) || (r2 > MLC_OUTER_RADIUS * MLC_OUTER_RADIUS && aU._z < 0)) return -1;
		--irCyl;
		int zplane = _gaps[irCyl].nextRegionAlongPath(-1, aX, aU, aIntendedDistance, material);
		if (zplane < 0) return -1;
		if (aIntendedDistance < -1e-3) printf("ViewRayMLC::nextRegionAlongPath: negative step %g to enter from outside\n", aIntendedDistance);
		// check again this logic !!!!!!!!!!!!!!!!!
		_boundingCylinders.nextRegionAlongPath(irCyl + 1, aX, aU, aIntendedDistance, material);
		int leafPair = zplane / 2;
		if (leafPair * 2 == zplane) {
			Vector xnew(aX + aU*aIntendedDistance);
			int newPairLocalRegion = _cylinderOfLeaves[irCyl][aTimeIndex]._leaf[leafPair].isWhere(xnew, material);
			if (newPairLocalRegion < 0) return -1;
			material = newPairLocalRegion == 1 ? Air : Tungsten;
			return 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * leafPair + newPairLocalRegion;
		}
		material = Air;
		return 4 * NLEAF_PAIRS_TOTAL*irCyl + 4 * leafPair + 3;

	}
	/**
	* \brief Return the number of regions in this class.
	*
	*/
	int getNumRegions() { return NLEAF_PAIRS_TOTAL * 4 * 2; }

private:

	/**
	* \brief
	*
	* Leaves are modeled with an upper and lower half.
	* Each half fits inside a cylinder.
	*/
	ViewRayLeaves* _cylinderOfLeaves[2];
	/**
	* \brief A vector of air gaps between the leaves.
	*
	* Each air gap is one class inside this vector.
	*
	*/
	ViewRayGaps _gaps[2];

	/**
	* \brief Two cylinders lining the inside and the outside of the MLC.
	*
	*/
	ViewRayMlcCylinders _boundingCylinders;
};
/******************************************  End: Model of MLC  *************************************/

/****************************************** Start: Model of GantryCoil *************************************/
class ViewRayGantryCoils {
	/*! \brief Implements the IGeometry interface for the gradient coils of the ViewRay delivery system.
	*
	* This class contains the coils around the patient.
	* Gradient, RF, TxPCB and TxFr.
	*/
public:
	ViewRayGantryCoils(){

		Vector center(0, 0, 0);

		ZFloat radii[8];
		radii[0] = GlueF(35.000000);
		radii[1] = GlueF(35.200000);
		radii[2] = GlueF(35.501854);
		radii[3] = GlueF(35.550000);
		radii[4] = GlueF(39.600000);
		radii[5] = GlueF(39.634540);
		radii[6] = GlueF(40.000000);
		radii[7] = GlueF(40.500000);

		for (int i = 0; i < 8; ++i)
		{
			_cyls[i].init(radii[i], center);
			_radiiSquared[i] = radii[i] * radii[i];
		}

		_media[0] = (Air);
		_media[1] = (TxFR4);
		_media[2] = (Air);
		_media[3] = (TxPCB);
		_media[4] = (Air);
		_media[5] = (RFShield);
		_media[6] = (Air);
		_media[7] = (GradientCoil);
	}

	// \return Which region particle is in.
	__device__ int isWhere(const Vector &aX, Medium& medium) const {
		ZFloat distance = (aX._x ) * (aX._x);
		distance += (aX._z) * (aX._z);

		if (distance >= _radiiSquared[8 - 1]) { medium = Air; return 9; }

		for (int i = 0; i < 8; i++) {
			if (distance < _radiiSquared[i]) { medium = _media[i]; return i; }
		}
		return 9;
	}

	__device__ int nextRegionAlongPath(int aCurrentRegion, const Vector &aX, const Vector &aU, ZFloat &aIntendedDistance, Medium& medium) const {

		int nextRegion = -1;

		const Vector& shifted = aX;

		if (aCurrentRegion == -1) aCurrentRegion = isWhere(aX, medium);

		if (aCurrentRegion == 0) {

			if (_cyls[0].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
				medium = _media[1];
				return 1;
			}
			else
			{
				return 0;
			}

		}
		else if (aCurrentRegion == 9) {

			ZFloat B = shifted._x*aU._x + shifted._z*aU._z;
			if (B < 0) { // traveling inward
				if (_cyls[8 - 1].nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) {
					medium = _media[8 - 1];
					return (8 - 1);
				}
			}

		}
		else if (aCurrentRegion <(8 - 1)) {

			ZFloat B = shifted._x*aU._x + shifted._z*aU._z;

			if (B >= 0) { // traveling outward

				if (_cyls[aCurrentRegion].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					medium = _media[aCurrentRegion + 1];
					return aCurrentRegion + 1;
				}
				else
				{
					return aCurrentRegion;
				}

			}
			else { // traveling inward

				if (_cyls[aCurrentRegion - 1].nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) {
					medium = _media[aCurrentRegion - 1];
					return aCurrentRegion - 1;
				}
				if (_cyls[aCurrentRegion].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					medium = _media[aCurrentRegion];
					return aCurrentRegion;
				}

				return aCurrentRegion;

			}

		}
		else if (aCurrentRegion == 8 - 1) {

			ZFloat B = shifted._x*aU._x + shifted._z*aU._z;

			if (B >= 0) { // traveling outward

				if (_cyls[aCurrentRegion].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					return -1;
				}
				else
				{
					return (8 - 1);
				}

			}
			else { // traveling inward

				if (_cyls[aCurrentRegion - 1].nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) {
					medium = _media[aCurrentRegion - 1];
					return aCurrentRegion - 1;
				}
				if (_cyls[aCurrentRegion].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) {
					medium = _media[aCurrentRegion];
					return aCurrentRegion;
				}

				return aCurrentRegion;

			}

		}

		return nextRegion;
	}

private:

	ZFloat _radiiSquared[8];
	Medium _media[8];
	Cylinder _cyls[8];
};
/******************************************  End:  Model of GantryCoil *************************************/
class CylinderPair {

public:

	/*! \brief Construct the pair of cylinders.
	*
	* @param[in]  aRadii  The two radii of the cylinder pair
	* @param[in]  aCenter A point on the axis of the cylinder pair
	* @param[in]  aInnerRegion The inner region index (i.e., the region inside the smaller cylinder)
	* @param[in]  aMiddleRegion The middle region index (i.e., the region between the two cylinders)
	* @param[in]  aOuterRegion The outer region index (i.e., the region outside the larger cylinder)
	*/
	void init(ZFloat r0, ZFloat r1, const Vector aCenter, const int aInnerRegion, const int aMiddleRegion, const int aOuterRegion) {

		_innerRegion = aInnerRegion;
		_middleRegion = aMiddleRegion;
		_outerRegion = aOuterRegion;

		_radiiSquared[0] = r0*r0;
		_radiiSquared[1] = r1*r1;

		_center._x = aCenter._x;
		_center._y = aCenter._y;
		_center._z = aCenter._z;

		_cyls[0].init(r0, _center);
		_cyls[1].init(r1, _center);

	}

	/*! \brief Implements the IGeometry::getNumRegions() interface */
	int getNumRegions() { return 3; }

	/*! \brief Implements the IGeometry::isWhere() interface */
	__device__ int isWhere(unsigned int, const Vector &aX, Medium&) const {
		ZFloat distance = (aX._x - _center._x) * (aX._x - _center._x);
		distance += (aX._z - _center._z) * (aX._z - _center._z);
		if (distance <= _radiiSquared[0]) return _innerRegion;
		if (distance <= _radiiSquared[1]) return _middleRegion;
		return _outerRegion;
	}

	/*! \brief Implements the IGeometry::nextRegionAlongPath() interface */
	__device__ int nextRegionAlongPath(int aCurrentRegion, const Vector &aX, const Vector &aU, ZFloat &aIntendedDistance, Medium& aMedium) const {

		int nextRegion = -1;

		Vector shifted(aX);
		shifted -= _center;

		if (aCurrentRegion == _innerRegion) {

			if (_cyls[0].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) { aMedium = Air; return _middleRegion; }

		}
		else if (aCurrentRegion == _middleRegion) {

			ZFloat B = shifted._x*aU._x + shifted._z*aU._z;
			if (B >= 0) { // traveling outward
				if (_cyls[1].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) { aMedium = Air; return _outerRegion; }
			}
			else { // traveling inward
				if (_cyls[0].nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) { aMedium = Air; return _innerRegion; }
				if (_cyls[1].nextRegionAlongPath(aX, aU, aIntendedDistance, true) == true) { aMedium = Air; return _outerRegion; }
			}


		}
		else if (aCurrentRegion == _outerRegion) {

			ZFloat B = shifted._x*aU._x + shifted._z*aU._z;
			if (B < 0) { // traveling inward
				if (_cyls[1].nextRegionAlongPath(aX, aU, aIntendedDistance, false) == true) { aMedium = Air; return _middleRegion; }
			}

		}

		return nextRegion;

	}

private:

	ZFloat _radiiSquared[2];      //!< The two radii squared
	Vector _center;                          //!< A point on the cylinder pair axis
	int _innerRegion;                           //!< The inner region index
	int _middleRegion;                          //!< The middle region index
	int _outerRegion;                           //!< The outer region index
	Cylinder _cyls[2];      //!< The actual implementation of the geometry methods.
};

class BaseRegion {
	/**
	* \brief  Base region for MLC and gradient
	*
	* Class describes the geometry of the base region, a large parallelpiped marking the boundaries
	* of the universe, with cylinders for the MLC and gradient.
	*/
public:

	BaseRegion()
	{
		Plane _top(Vector(0, 0, 70), Vector(0, 0, -1));
		Plane _bot(Vector(0, 0, -5), Vector(0, 0, +1));
		Plane _lef(Vector(-42, 0, 0), Vector(+1, 0, 0));
		Plane _rig(Vector(+42, 0, 0), Vector(-1, 0, 0));
		Plane _bac(Vector(0, +42, 0), Vector(0, -1, 0));
		Plane _fro(Vector(0, -42, 0), Vector(0, +1, 0));
		_bot.setRegionCode(-1);
		_top.setRegionCode(-1);
		_lef.setRegionCode(-1);
		_rig.setRegionCode(-1);
		_fro.setRegionCode(-1);
		_bac.setRegionCode(-1);

		_xplanes[0] = _lef;
		_xplanes[1] = _rig;
		_yplanes[0] = _fro;
		_yplanes[1] = _bac;
		_zplanes[0] = _bot;
		_zplanes[1] = _top;

		_mlc.init(41.0, 50.0, Vector(0, 0, SOURCE_TO_ISO_DISTANCE), 1, 2, 0);
		_gradient.init(35.0, 40.5, Vector(0, 0, 0), 3, 4, 0);
	}


	/**
	* Return the base region that the particle is in.
	* \param aX position of the particle.
	* \return Which region particle is in.
	*/
	__device__ int isWhere(const Vector &aX, Medium&) const {

		Medium medium;

		if (aX._x > _xplanes[0].X() && aX._x < _xplanes[1].X()) {
			if (aX._y > _yplanes[0].Y() && aX._y < _yplanes[1].Y()) {
				if (aX._z > _zplanes[0].Z() && aX._z < _zplanes[1].Z()) {

					// Check if inside MLC
					int mlcReg = _mlc.isWhere(0, aX, medium);
					if (mlcReg != 0) return mlcReg;

					// Check if inside gantry.
					int gradReg = _gradient.isWhere(0, aX, medium);
					if (gradReg != 0) return gradReg;

					// Inside universe.
					return 0;

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
	__device__ int nextRegionAlongPath(int aCurrentRegion, const Vector &aX, const Vector &aU, ZFloat &aIntendedDistance, Medium& aMedium) const {

		int newRegion = -1;

		if (aCurrentRegion == 0) {

			for (int i = 0; i < 2; i++) _xplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			ZFloat distance = aIntendedDistance;
			int tempMlcRegion = _mlc.nextRegionAlongPath(0, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempMlcRegion;
				aIntendedDistance = distance;
				if (newRegion == 0) newRegion = 0;
			}

			distance = aIntendedDistance;
			int tempGradientRegion = _gradient.nextRegionAlongPath(0, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempGradientRegion;
				aIntendedDistance = distance;
				if (newRegion == 0) newRegion = 0;
			}

		}
		else if (aCurrentRegion == 1 || aCurrentRegion == 2) {

			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			ZFloat distance = aIntendedDistance;
			int tempMlcRegion = _mlc.nextRegionAlongPath(aCurrentRegion, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempMlcRegion;
				aIntendedDistance = distance;
			}

		}
		else if (aCurrentRegion == 3 || aCurrentRegion == 4) {

			for (int i = 0; i < 2; i++) _yplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);
			for (int i = 0; i < 2; i++) _zplanes[i].nextRegionAlongPath(aX, aU, aIntendedDistance, newRegion);

			ZFloat distance = aIntendedDistance;
			int tempGradientRegion = _gradient.nextRegionAlongPath(aCurrentRegion, aX, aU, distance, aMedium);
			if (distance < aIntendedDistance) {
				newRegion = tempGradientRegion;
				aIntendedDistance = distance;
			}

		}

		return newRegion;

	}

private:

	/**
	* The left and right planes are stored in this vector.
	*/
	Plane _xplanes[2];
	/**
	* The top and bottom planes are stored in this vector.
	*/
	Plane _yplanes[2];
	/**
	* The front and back planes are stored in this vector.
	*/
	Plane _zplanes[2];

	/**
	* Two cylinders marking inner and outer MLC boundaries.
	*/
	CylinderPair _mlc;

	/**
	* Two cylinders marking inner and outer gradient boundaries.
	*/
	CylinderPair _gradient;
};

class Envelope {
	/**
	* \brief Class to contain the universe of the simulation.
	*
	* Class contains the base region, which limits the extent of the simulation and marks
	* the boundaries of the inscribed regions.  also then are one or more inscribed regions,
	* which contain the most detailed boundaries.
	**/
public:
	Envelope()
	{
		_localStart[0] = 5;
		_localStart[1] = _localStart[0] + NLEAF_PAIRS_TOTAL * 8;
		_localStart[2] = _localStart[1] + 8;
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
	__device__ int nextRegionAlongPath(unsigned int aTimeIndex, int aCurrentRegion, const Vector &aX, const Vector &aU, ZFloat &aIntendedDistance, Medium& material) const {

		if (aCurrentRegion < 5) {

			ZFloat distance = aIntendedDistance;
			Medium baseMedium = material;

			int newBaseRegion = _baseRegion.nextRegionAlongPath(aCurrentRegion, aX, aU, distance, baseMedium);
			if (newBaseRegion == -1) { aIntendedDistance = distance; material = baseMedium; return newBaseRegion; }

			int inscribedRegion = newBaseRegion % 2 + newBaseRegion / 2;
			if (inscribedRegion > 0 && inscribedRegion <3) {
				int newregion = 0;
				if (inscribedRegion == 1) newregion = _inscribedRegion1.nextRegionAlongPath(aTimeIndex, -1, aX, aU, aIntendedDistance, material);
				else newregion = _inscribedRegion2.nextRegionAlongPath(-1, aX, aU, aIntendedDistance, material);
				if (newregion != -1) return (_localStart[inscribedRegion - 1] + newregion);
			}

			aIntendedDistance = distance;
			material = baseMedium;
			return newBaseRegion;

		}

		// We are in an inscribed region
		if (aCurrentRegion < 5 + NLEAF_PAIRS_TOTAL * 8) {

			int inscribedRegionIndex = 1; // MLC
			int baseRegion = 2;
			int newInscribedRegion = _inscribedRegion1.nextRegionAlongPath(aTimeIndex, aCurrentRegion - _localStart[inscribedRegionIndex - 1], aX, aU, aIntendedDistance, material);

			ZFloat distanceBase = 9999;
			Medium baseMedium = material;
			int newBaseRegion = _baseRegion.nextRegionAlongPath(baseRegion, aX, aU, distanceBase, baseMedium);

			if (newBaseRegion != baseRegion && distanceBase < aIntendedDistance) {
				if (newBaseRegion < 0) return newBaseRegion;
				int newInscribedRegionIndex = newBaseRegion % 2 + newBaseRegion / 2;
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
		else if (aCurrentRegion < 5 + NLEAF_PAIRS_TOTAL * 8 + 8) {

			int inscribedRegionIndex = 2; // Gantry
			int baseRegion = 4;
			int newInscribedRegion = _inscribedRegion2.nextRegionAlongPath(aCurrentRegion - _localStart[inscribedRegionIndex - 1], aX, aU, aIntendedDistance, material);

			ZFloat distanceBase = 9999;
			Medium baseMedium = material;
			int newBaseRegion = _baseRegion.nextRegionAlongPath(baseRegion, aX, aU, distanceBase, baseMedium);

			if (newBaseRegion != baseRegion && distanceBase < aIntendedDistance) {
				if (newBaseRegion < 0) return newBaseRegion;
				int newInscribedRegionIndex = newBaseRegion % 2 + newBaseRegion / 2;
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
		ZFloat distance = aIntendedDistance;
		Medium baseMedium = material;
		int newBaseRegion = _baseRegion.nextRegionAlongPath(-1, aX, aU, distance, baseMedium);
		if (newBaseRegion == -1) return newBaseRegion;
		int inscribedRegion = newBaseRegion % 2 + newBaseRegion / 2;
		if (inscribedRegion > 0) {
			int newregion = 0;
			if (inscribedRegion == 1) newregion = _inscribedRegion1.nextRegionAlongPath(aTimeIndex, -1, aX, aU, aIntendedDistance, material);
			else newregion = _inscribedRegion2.nextRegionAlongPath(-1, aX, aU, aIntendedDistance, material);
			if (newregion != -1) return (_localStart[inscribedRegion - 1] + newregion);
		}
		if (newBaseRegion == aCurrentRegion || newBaseRegion < 0) { material = baseMedium; aIntendedDistance = distance; return newBaseRegion; }

		return -1;
	}
	ViewRayMLC& getViewRayMLC(){ return _inscribedRegion1; }
private:

	/**
	* The geometry of the base region.
	*/
	BaseRegion _baseRegion;
	/**
	* A vector of inscribed region geometries.
	*/
	ViewRayMLC _inscribedRegion1;
	ViewRayGantryCoils _inscribedRegion2;

	/**
	* Upper bound on region number.  First entry is the number of base
	* regions, subsequent entries are the number of inscribed regions.
	* Each entry is offset by sum of all previous entries.
	*
	*/
	int _localStart[3];
};

class HeadAttenuation
{
public:
	void loadData(const char *aFileName)
	{
		FILE* fp = fopen(aFileName, "rb");
		if (!fp) exitApp("Cannot find the head attenuation data file");
		double elem = 0;
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 2000; ++j)
			{
				fread(&elem, sizeof(double), 1, fp);
				total[i][j] = ZFloat(elem);
			}
		}
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 2000; ++j)
			{
				fread(&elem, sizeof(double), 1, fp);
				compton[i][j] = ZFloat(elem);
			}
		}
		fclose(fp);
	}
	__device__ ZFloat calcAttenuation(Medium aMedium, ZFloat aPhotonEnergy)
	{
		const ZFloat invDE = ZFloat(1 / 0.06495e-2);
		ZFloat r = (aPhotonEnergy - GlueF(0.05)) *invDE;
		int i = (int)r;
		r -= i;
		return total[aMedium][i] * (1 - r) + total[aMedium][i + 1] * r;
	}

	__device__ ZFloat calcComptonAttenuation(Medium aMedium, ZFloat aPhotonEnergy)
	{
		const ZFloat invDE = ZFloat(1 / 0.06495e-2);
		ZFloat r = (aPhotonEnergy - GlueF(0.05))*invDE;
		int i = (int)r;
		r -= i;
		return compton[aMedium][i] * (1 - r) + compton[aMedium][i + 1] * r;
	}
	//energy range E0 = 5e-2, dE = 0.06495
	ZFloat total[6][2000];
	ZFloat compton[6][2000];
};

class HeadTransport {

public:

	/*! \brief Transport a photon through the simulation geometry

	This function is essentially the same as the public function with the same name, except that now the geometry region and medium
	index are also passed as input.  It is used by the public function to perform the actual transport.
	*/
	__device__ void transportPhoton(GRNG &rng, int aSegment, ZFloat E, Vector &x, Vector &u, int weight, ZFloat baseWeight, int ireg, Medium imed, bool isInFOI,
		ZFloat lastMutr, int &nParticle, ParticleR aParticleList[]);

	void setSegments(const vector<vector<pair<double, double> > > &aListOfSegments, GPUConfig& gc);

protected:

	__device__ bool doSmartCompton(GRNG &rng, int aSegmentIndex, int ireg, Medium imed, ZFloat &E, const Vector &x, Vector &u,
		ZFloat baseWeight, ZFloat lastMutr, int &nParticle, ParticleR aParticleList[]);

private:

	/*! \brief Envelope containing the mlc and the gantry coils.*/
	Envelope _e;
	HeadAttenuation _headAttenuation;

	__device__ __forceinline__ ZFloat calcAttenuation(Medium aMedium, ZFloat aPhotonEnergy){ return _headAttenuation.calcAttenuation(aMedium, aPhotonEnergy); }

	__device__ __forceinline__ ZFloat calcComptonAttenuation(Medium aMedium, ZFloat aPhotonEnergy){ return _headAttenuation.calcComptonAttenuation(aMedium, aPhotonEnergy); }
};

//__device__ HeadTransport headTransport;

void HeadTransport::setSegments(const vector<vector<pair<double, double> > > &aListOfSegments, GPUConfig& gc)
{
	_headAttenuation.loadData("Head.binary");

	_e.getViewRayMLC().setSegments(aListOfSegments, gc);
}

__device__ void HeadTransport::transportPhoton(GRNG &rng, int aSegmentIndex, ZFloat E, Vector &x, Vector &u, int weight, 
	ZFloat baseWeight, int ireg, Medium imed, bool isInFOI, ZFloat lastMutr, int &nParticle, ParticleR aParticleList[]) {

	Medium newmed = imed;

	if (E <= HEAD_pCut || E > 1.34) {
		return;
	}

	//
	// *** Remember the photon position 
	//
	Vector xo(x);

	//
	// *** Set up interpolations for this photon and compute attenuation coefficient
	//
	ZFloat mu = calcAttenuation(imed, E);

	ZFloat mui = 1.0 / mu;
	//
	// *** If the weight is greater than 1 and the photon is going towards the FOI, we follow it until it accumulates lambdaMax = log(weight) mean-free-paths. 
	//     If it happens to interact before, we simulate the interaction using doSmartCompton(). When the photon has traveled lambdaMax MFP, 
	//     we attenuate it by exp(-lambdaMax)=1/weight (so its weight becomes 1) and follow it in the usual way
	//
	if (weight > 1 && isInFOI) {

		//
		// *** Compute lambdaMax. Strangely enough, using the tabulated log(weight) made things slightly slower. I guess, a cache issue
		//
		//ZFloat lambdaMax = _weightToLambda[weight-2];           // number of MFPs until weight becomes 1
		ZFloat lambdaMax = GlueF(log)((ZFloat)weight);

		//
		// *** Sample number of MFP until next interaction and remember it
		//
		ZFloat lambda = -GlueF(log)(1 - rng());       // MFPs to next interaction
		ZFloat lambda0 = lambda;

		//
		// *** See if the photon will interact before reaching lambdaMax
		//
		bool doInteraction = true;
		if (lambda >= lambdaMax) {
			// No, it will not. So, set the number of MFP to betraveled to lambdaMax.
			doInteraction = false; lambda = lambdaMax;
		}

		ZFloat totLambda = 0;  // initialize total number of MFP traveled
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
				printf("Fatal error: HeadTransport::transportPhoton: more than 1000 steps taken through geometry\n");
			}

			// distance left to travel
			ZFloat t = lambda*mui;

			// query geometry
			int inew = _e.nextRegionAlongPath(aSegmentIndex, ireg, x, u, t, newmed);

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
			if (inew == 3 && isInFOI) {
				//
				// *** Yes, it did. We tolerate 'misbehaviour' of up to ln(2) MFP's (so that the weight of the photon recorded in 
				//     aParticleList is not more than twice the weight of other photons.
				//      
				if (lambdaMax - totLambda < 0.693147) {
					//
					// *** Store this photon in aParticleList
					//
					if (nParticle >= NMaxInitPars) printf("Fatal error: sampling stack overflow. You have to increase the constant NSampleStack");
					ParticleR& p = aParticleList[nParticle];
					p.x = x._x;
					p.y = -x._z;
					p.z = x._y;
					p.u = u._x;
					p.v = -u._z;
					p.w = u._y;
					p.E = 1e6*E;
					p.weight = baseWeight*GlueF(exp)(lambdaMax - totLambda);
					p.type = photon;

					++nParticle;
				}
				else {
					//
					// *** Warn about this event. If this happens too many times, whoever is using this class has to redesign the way they 
					//     sample high-weight photons
					//
					printf("A particle with weight %d is entering scoring region before reaching lambdaMax attenuation. segment = %d\n", weight, aSegmentIndex);
					printf("E=%g x=(%g,%g,%g) u=(%g,%g,%g) xo=(%g,%g,%g) lambda=%g (%g)\n", E, x._x, x._y, x._z, u._x, u._y, u._z, xo._x, xo._y, xo._z, lambdaMax - totLambda, totLambda);
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
					if (!HEAD_NO_SCATTER) {
						//
						// *** Further, take into account Compton interaction probability and Russian Roullette game (if weight < HEAD_N_SPLIT). 
						//     The probability to do the Compton interaction becomes Pcomton * weight/HEAD_N_SPLIT.
						//
						ZFloat comptonSigma = calcComptonAttenuation(imed, E);
						ZFloat leftSide = rng()*HEAD_N_SPLIT*mu;
						if (leftSide < comptonSigma*weight) {
							//
							// *** OK, we need to do Compton scatter. Perform using doSmartCompton().
							//
							ZFloat Enew = E; Vector unew(u);
							if (doSmartCompton(rng, aSegmentIndex, ireg, imed, Enew, x, unew, baseWeight, lastMutr, nParticle, aParticleList)) {
								//
								// *** If doSmartCompton() returned true, the scattered photon is not going towards the FOI => transport it.
								//
								Vector xnew(x);
								transportPhoton(rng, aSegmentIndex, Enew, xnew, unew, HEAD_N_SPLIT, baseWeight, ireg, imed, false, lastMutr, nParticle, aParticleList);
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
					mu = calcAttenuation(imed, E);
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
		ZFloat lambda = -GlueF(log)(1 - rng());
		int nstep = 0;  // as above, count steps just in case.
		while (1) {
			++nstep;
			// avoid particles stuck on boundaries (or otherwise confused) by discarding them.
			if (nstep > 1000) {
				printf("Fatal error: HeadTransport::transportPhoton(a): more than 1000 steps taken through geometry\n");
			}
			// distance to travel for lambda MFP
			ZFloat t = lambda*mui;
			// query geometry
			int inew = _e.nextRegionAlongPath(aSegmentIndex, ireg, x, u, t, newmed);

			// see comments above about inew < 0 or t < 0
			if (t < -0.01 || inew < 0) {
				return;
			}

			if (newmed == UnknownMaterial) printf("Fatal error: encounter unknown material??");
			// update position
			x += u*t;
			// has the photon entered the region where we score particles escaping the treatment head?
			if (inew == 3) {
				if (weight == 1) {
					// yes, it has. Store the particle into the container provided and terminate transport
					
					if (E <= 0 || E > 1.34) {
						printf("Particle with E=%g arrived to be stored?\n", E);
						return;
					}
					if (nParticle >= NMaxInitPars) printf("Fatal error: sampling stack overflow! you have to increase NSampleStack");
					ParticleR& p = aParticleList[nParticle];
					p.x = x._x;
					p.y = -x._z;
					p.z = x._y;
					p.u = u._x;
					p.v = -u._z;
					p.w = u._y;
					p.E = 1e6*E;
					p.weight = baseWeight;
					p.type = photon;

					++nParticle;
				}
				return;
			}
			// inew == ireg means we have reached the interaction site => exit the tracking loop
			if (inew == ireg) break;
			// update lambda
			lambda -= t*mu;
			if (newmed != imed) {
				if (newmed == UnknownMaterial) printf("Fatal error: encounter unknown material??");
				// the photon has entered a region with a different material, so compute new attenuation coefficient
				imed = newmed;
				mu = calcAttenuation(imed, E);
				mui = 1 / mu;
			}
			ireg = inew;
		}

		// 
		//  *** The photon has reached an interaction site. 
		// 
		if (HEAD_NO_SCATTER) return;// i.e., if the user has turned off scatter, simply terminate.

		//
		// We play RussianRoulette with survival probability of weight/HEAD_N_SPLIT. Also, if the interaction is 
		// not Compton scattering, we simply terminate the photon history. Thus, probability to survive RR and to do Compton scattering is 
		// pCompton*weight/HEAD_N_SPLIT. 
		//
		ZFloat comptonSigma = calcComptonAttenuation(imed, E);
		ZFloat leftSide = rng()*HEAD_N_SPLIT*mu;
		bool keepGoing = leftSide < comptonSigma*weight ?
			doSmartCompton(rng, aSegmentIndex, ireg, imed, E, x, u, baseWeight, lastMutr, nParticle, aParticleList) : false;
		if (!keepGoing) break;

		// 
		// *** If here, the photon was scattered into a direction not going towards the FOI. 
		//     We need to recompute log(E), setup interpolation and continue transport with the new energy.
		//
		weight = HEAD_N_SPLIT;
		mu = calcAttenuation(imed, E);
		mui = 1 / mu;
	}

}

__device__ bool HeadTransport::doSmartCompton(GRNG &rng, int aSegmentIndex, int ireg, Medium imed, ZFloat &E, const Vector &x, Vector &u,
	ZFloat baseWeight, ZFloat lastMutr, int &nParticle, ParticleR aParticleList[]) {

	//
	// *** Compute minimum and maximum scattering angle that may make the scattered photon move towards the field of interest.
	//     The interested examiner of this piece of code is encouraged to verify the correctness of this calculation.
	//
	ZFloat xmin = HEAD_XMIN, xmax = HEAD_XMAX;
	ZFloat ymin = HEAD_YMIN, ymax = HEAD_YMAX;
	ZFloat x1, y1, ctmin, ctmax;
	if (u._z) { ZFloat t = -x._z / u._z; x1 = x._x + u._x*t; y1 = x._y + u._y*t; }
	else      { x1 = x._x; y1 = x._y; }
	ZFloat xx1, yy1, xx2, yy2;
	if (x1 < xmin) { xx1 = xmin; xx2 = xmax; }
	else if (x1 < xmax) { xx1 = x1; xx2 = xmax - x1 > x1 - xmin ? xmax : xmin; }
	else                 { xx1 = xmax; xx2 = xmin; }
	if (y1 < ymin) { yy1 = ymin; yy2 = ymax; }
	else if (y1 < ymax) { yy1 = y1; yy2 = ymax - y1 > y1 - ymin ? ymax : ymin; }
	else                 { yy1 = ymax; yy2 = ymin; }
	Vector v1(xx1, yy1, 0); v1 -= x; v1.normalizeToUnitLength(); ZFloat ct1 = v1*u;
	Vector v2(xx2, yy2, 0); v2 -= x; v2.normalizeToUnitLength(); ZFloat ct2 = v2*u;
	if (ct1 < ct2) { ctmin = ct1; ctmax = ct2; }
	else            { ctmin = ct2; ctmax = ct1; }

	//
	// *** Compute an estimate of the probability for scatter into the field of interest
	//     Note: this is an overestimate of the true probability
	//
	// - Minimum and maximum possible energy fraction.
	//    eps1 is the minimum possible energy fraction, eps2 the maximum possible. eps1i=1/eps1, eps2i = 1/eps2
	ZFloat ko = E*1.9569341; ZFloat broi = 1 + 2 * ko, ko2 = ko*ko;
	ZFloat eps1i = 1 + ko*(1 - ctmin), eps2i = 1 + ko*(1 - ctmax);
	if (E <= HEAD_pCut*eps2i)  return false;  // i.e., the maximum possible energy towards the field of interest is less than the cuttoff, so no need to simulate.
	ZFloat eps1 = 1 / eps1i, eps2 = 1 / eps2i;


	// - The integral of Klein-Nishina over all angles and the Klein-Nishina normalization (fnorm)
	ZFloat alpha1_t = GlueF(log)(broi);
	ZFloat eps1_t = 1 / broi, eps2_t = 1;
	ZFloat w2 = alpha1_t*(ko2 - 2 * ko - 2) + (eps2_t - eps1_t)*(2 * broi + 0.5*ko2*(eps1_t + eps2_t));
	ZFloat fnorm = w2 / (ko2*ko);

	// - The maximum of Klein-Nishina in the ctmin...ctmax interval.
	ZFloat f2 = (eps2i + eps2 - 1 + ctmax*ctmax)*eps2*eps2;
	ZFloat f1 = (eps1i + eps1 - 1 + ctmin*ctmin)*eps1*eps1;
	ZFloat fmax = f1 > f2 ? f1 : f2;

	// - The overestimate of the probability to have scattered photons going towards the FOI discussed above
	//   (the probability P mentioned above is called wnew here because it is new compared to an earlier overestimate of P used by 
	//    Kawrakow, Walters and Rogers in their paper on DBS [Med. Phys. 31 (2004) 2883]
	//
	ZFloat A = (xmax - xmin)*(ymax - ymin);
	ZFloat wnew = A / (2 * PI*x._z*x._z*fnorm)*fmax;

	// - Number of interactions to sample: integer of wnew*HEAD_N_SPLIT or, with corresponding probability, 1 + integer of wnew*HEAD_N_SPLIT.
	ZFloat asample = wnew*HEAD_N_SPLIT; int nsample = (int)asample; asample = asample - nsample;
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
		uu -= x; ZFloat disti = 1 / uu.getLength(); uu *= disti;

		//
		// *** The cosine of the scattering angle and the corresponding scattered photon energy ( E/epsi )
		//
		ZFloat cost = uu*u;
		ZFloat epsi = 1 + ko*(1 - cost);

		if (E > HEAD_pCut*epsi) {          // i.e., only simulate if scattered photon energy is above the cutoff
			ZFloat eps = 1 / epsi;
			ZFloat aux = disti*x._z;
			// The condition below is a better way of writing the rejection probability condition avoiding expensive divisions
			if (rng()*fmax < (1 + eps*eps - eps*(1 - cost*cost))*eps*aux*aux*aux) {

				// The scattered photon energy
				ZFloat Enew = E*eps;

				// Now trace the scattered photon through the geometry
				Vector xx(x);
				transportPhoton(rng, aSegmentIndex, Enew, xx, uu, 1, baseWeight, ireg, imed, true, lastMutr, nParticle, aParticleList);
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
	ZFloat bro = 1 / broi; ZFloat br, temp, sint;
	if (ko < 11) {
		ZFloat bro1 = 1 - bro;
		ZFloat rejmax = ko2*(broi + bro);
		ZFloat br2;
		do {
			br = bro + bro1*rng(); br2 = br*br;
		} while (rng()*br2*rejmax > ko2*br*(br2 + 1) - (1 - br)*(br*broi - 1));
		temp = (1 - br) / (ko*br); sint = temp*(2 - temp);
	}
	else {
		ZFloat broi2 = broi*broi;
		ZFloat alpha1 = GlueF(log)(broi);
		ZFloat alpha2 = ko*(broi + 1)*bro*bro;
		ZFloat alphaS = alpha1 + alpha2;
		do {
			br = rng()*alphaS < alpha1 ? exp(alpha1*rng())*bro : sqrt(rng()*(broi2 - 1) + 1)*bro;
			temp = (1 - br) / (ko*br); sint = temp*(2 - temp);
		} while (rng()*(1 + br*br) < br*sint);
	}

	//
	// 2. Now compute direction and see if the scattered photon has energy above HEAD_pCut and does not go towards the FOI
	//
	bool keepIt = true; E *= br;
	if (E > HEAD_pCut) {
		ZFloat cost;
		if (sint > 0) { cost = 1 - temp; sint = GlueF(sqrt)(sint); }
		else { cost = -1; sint = 0; }
		ZFloat cphi, sphi;
		randomAzimuth(rng(),cphi, sphi);
		// 		ZFloat phi = 2 * PI*rng();
		// 		cphi = cos(phi);
		// 		sphi = sqrt(1 - cphi*cphi);
		u.rotate(cost, sint, cphi, sphi);
		if (u._z < 0) {
			ZFloat t = -x._z / u._z; x1 = x._x + u._x*t; y1 = x._y + u._y*t;
			if (x1 > xmin && x1 < xmax && y1 > ymin && y1 < ymax) keepIt = false;
		}
	}
	else keepIt = false;

	return keepIt;
}



/******************************************************Start: GPU memory******************************************************/
__device__ int CoPhaseSpace4GPU::sample(GRNG &rng, bool isPrimary, int bin, ZFloat& E, Vector& x, Vector& u)
{
	//
	// *** Convert fluence bin to a random position within the bin and set position vector
	//
	int iyy = bin / 240; int ixx = bin - iyy * 240;
	x._x = -6 + GlueF(0.05)*(rng() + ixx);
	x._y = -6 + GlueF(0.05)*(rng() + iyy);
	x._z = GlueF(65.39561);

	//
	// *** Find the fluence tile in which this position falls
	//
	ZFloat sx = 1; if (x._x < 0) sx = -1;
	ZFloat sy = 1; if (x._y < 0) sy = -1;
	int ix = (int)(x._x*sx*GlueF(2.5)); if (ix > 15 - 1) ix = 15 - 1;
	int iy = (int)(x._y*sy*GlueF(2.5)); if (iy > 15 - 1) iy = 15 - 1;
	int j = ix + iy * 15;

	//
	// *** Set sampling tables depending on whether we are sampling a primary or a scattered particle and also sample the energy
	//
	unsigned short *thetaPhiFlu; ZFloat minx, miny, delx, dely; int nbin, Nxb;
	if (isPrimary) {
		// Primary particle. 
		// Pick energy
		if (rng() < _prob1[bin]) {
			E = GlueF(1.1732);
		}
		else {
			E = GlueF(1.3325);
		}
		// Set the tables (and corresponding bin-widths and ranges)
		thetaPhiFlu = _thetaPhiPflu[j];
		minx = -GlueF(2.1); miny = -GlueF(2.1); delx = GlueF(0.01640625); dely = GlueF(0.01640625); nbin = 256 * 256; Nxb = 256;
	}
	else {
		// Scattered
		// Sample the energy
		unsigned short *spec = _scatSpectrum[j];
		int iE = (int)(rng() * 254);
		if (rng()*65535 > spec[2 * iE]) iE = spec[2 * iE + 1];
		E = GlueF(0.05) + GlueF(0.0051181102362204724)*(rng() + iE);
		// set tables (and corresponding bin-widths and ranges)
		ZFloat rnno = rng();
		if (rnno < _pScoarse[j]) {
			// pick photon direction from the fine angular grid
			thetaPhiFlu = _thetaPhiSflu[j];
			minx = -GlueF(2.4); miny = -GlueF(2.4); delx = GlueF(0.01875); dely = GlueF(0.01875); nbin = 256 * 256; Nxb = 256;
		}
		else if (rnno < _pScoarse1[j]) {
			// pick photon direction from the less fine angular grid
			thetaPhiFlu = _thetaPhiSflu1[j];
			minx = -GlueF(7.2); miny = -GlueF(7.2); delx = GlueF(0.1125); dely = GlueF(0.1125); nbin = 128 * 128; Nxb = 128;
		}
		else {
			// pick photon direction from the coarse angular grid
			thetaPhiFlu = _thetaPhiSfluC[j];
			minx = -16; miny = -16; delx = GlueF(0.5); dely = GlueF(0.5); nbin = 64 * 64; Nxb = 64;
		}
	}

	//
	// *** Sample the angular bin
	//
	int aBin = (int)(rng()*nbin);
	if (rng()*65535 > thetaPhiFlu[2 * aBin]) aBin = thetaPhiFlu[2 * aBin + 1];
	int iyb = aBin / Nxb; int ixb = aBin - iyb*Nxb;

	//
	// *** Now sample uniformly within the selected bin
	//
	ZFloat dx = minx + (rng() + ixb)*delx;
	ZFloat dy = miny + (rng() + iyb)*dely;

	//
	// *** Position projected to the iso-center plane
	//
	u._x = x._x;
	u._y = x._y;
	u._z = x._z - 105;
	u.normalizeToUnitLength();
	ZFloat t = -x._z / u._z;

	u._x = u._x*t - sx*dx;
	u._y = u._y*t - sy*dy;
	u._z = u._z*t;
	u.normalizeToUnitLength();

	return 0;
}

__device__ int SEGMENT4GPU::sample(GRNG &rng, int iseg, ParticleR pars[])
{
	//obtain how many particles to sample in one call
	int NP = (int)nSample;
	ZFloat RSample = nSample - NP;
	if (rng() < RSample) ++NP;
	if (0 == NP) return 0; // i.e., no particles from this call

	int nPar = 0; //particle number that reaches the FOI
	for (int i = 0; i < NP; ++i) //sample NP particles from phase space file
	{
		bool isPrimary;
		unsigned short *table;
		char* mask;
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
		int bin = int(rng() * 240 * 240);
		if (rng()*65535 > table[2 * bin]) bin = table[2 * bin + 1];
		int weight = 1;
		if (mask[bin] == 1) weight = 10;
		if (mask[bin] == 2) weight = 34;

		ZFloat E = 0;
		Vector X, U;
		//generate one particle from phase space file
		gps.sample(rng, isPrimary, bin, E, X, U);
		headTransport->transportPhoton(rng, iseg, E, X, U, weight, GlueF(3.9843677959043429e-005), 1, Air, true, 0, nPar, pars);
	}


	//rotate the beam
	if (gantryAngle != 0.0)
	{
		for (int i = 0; i < nPar; ++i)
		{
			ZFloat xp = cosPhi*pars[i].x - sinPhi*pars[i].y;
			pars[i].y = sinPhi*pars[i].x + cosPhi*pars[i].y;
			pars[i].x = xp;

			xp = cosPhi*pars[i].u - sinPhi*pars[i].v;
			pars[i].v = sinPhi*pars[i].u + cosPhi*pars[i].v;
			pars[i].u = xp;
		}
	}
	return nPar;
}

__device__ int sample(GRNG &rng, ParticleR pars[])
{
	__shared__ int iss[16]; //assume block dim <= 512, which is paractical for current GPUs
	int iwarp = threadIdx.x / 32;
	if (threadIdx.x % 32 == 0)
	{
		int is = 0;
		if (nSeg > 1)
		{
			ZFloat as = rng()*nSeg;
			is = (int)as;
			as -= is;
			if (as >= segs[is].wfirst) is = segs[is].wsecond;
			if (is >= nSeg)	is = nSeg - 1;
		}
		iss[iwarp] = is;
	}
	iwarp = iss[iwarp];//iwarp now is the selected segment index
	return segs[iwarp].sample(rng, iwarp, pars);
}

__global__ void gSample()
{
	unsigned int it = blockIdx.x * blockDim.x + threadIdx.x; //thread index
	GRNG rng = RNGState[it];//read the rng status
	__shared__ int iss[16]; //assume block dim <= 512, which is paractical for current GPUs
	int iwarp = threadIdx.x / 32;
	if (threadIdx.x % 32 == 0)
	{
		int is = 0;
		if (nSeg > 1)
		{
			ZFloat as = rng()*nSeg;
			is = (int)as;
			as -= is;
			if (as >= segs[is].wfirst) is = segs[is].wsecond;
			if (is >= nSeg)	is = nSeg - 1;
		}
		iss[iwarp] = is;
	}
	__syncthreads();
	iwarp = iss[iwarp];//iwarp now is the selected segment index
	NInitPars[it] = segs[iwarp].sample(rng, iwarp, InitPars + it*NMaxInitPars);

// 	int np = NInitPars[it] - NMaxInitPars / 2;
// 	for (int i = 0; i < np; ++i) //copy from the end to the front
// 	{
// 		InitPars[it*NMaxInitPars + i] = InitPars[it*NMaxInitPars + i + NMaxInitPars / 2];
// 	}
// 
// 	//BufferSize == NMaxInitPars/2, the rest half is for safe access
// 	while (np < NMaxInitPars / 2)
// 	{
// 		np += segs[iwarp].sample(rng, iwarp, InitPars + it*NMaxInitPars + np);
// 	}
// 	NInitPars[it] = np;

// 	int is = 0;
// 	if (nSeg > 1)
// 	{
// 		ZFloat as = rng()*nSeg;
// 		is = (int)as;
// 		as -= is;
// 		if (as >= segs[is].wfirst) is = segs[is].wsecond;
// 		if (is >= nSeg)	is = nSeg - 1;
// 	}
//	NInitPars[it] = segs[is].sample(rng, is, InitPars + it*NMaxInitPars);
	RNGState[it] = rng; //write back the rng state
}

void SourceHead4GPU_Init(ConfigFile* cf, void* vgc)
{
	CoPhaseSpace4GPU* hps = NULL;
	int hnSeg = 0;
	SEGMENT4GPU* hseg = NULL;
	int nMaxPhoton = 0;
	vector<vector<pair<double, double> > > aListofSegments;
	getHostData(cf, hps, hseg, hnSeg, nMaxPhoton, aListofSegments);
	nMaxPhoton *= 10; // the actual sampled number can be much bigger
	HeadTransport* hhead = new HeadTransport;
	
	//copy them to GPU
	vector<GPUConfig>& gc = *((vector<GPUConfig>*)vgc);
	
	int NGPU = (int)gc.size();
	for (int i = 0; i < NGPU; ++i)
	{
		cudaSetDevice(gc[i].id);
		cudaMemcpyToSymbol(gps, hps, sizeof(CoPhaseSpace4GPU)); //copy the phase space data
		void* gp = NULL;
		cudaMalloc(&gp, hnSeg*sizeof(SEGMENT4GPU));
		gc[i].addGPUPointer(gp);//remember to free when exiting
		cudaMemcpyToSymbol(segs, &gp, sizeof(SEGMENT4GPU*));//copy the pointer to GPU
		cudaMemcpy(gp, hseg, hnSeg*sizeof(SEGMENT4GPU), cudaMemcpyHostToDevice); //copy the segment data to GPU
		cudaMemcpyToSymbol(nSeg, &hnSeg, sizeof(int)); //copy segment number

		int NT = gc[i].BlockSize*gc[i].NBlock;//number of thread in GPU

		cudaMalloc(&gp, nMaxPhoton*NT*sizeof(ParticleR));
		gc[i].addGPUPointer(gp); //remember to free when exiting
		cudaMemcpyToSymbol(InitPars, &gp, sizeof(ParticleR*));//copy the pointer to GPU

		cudaMemcpyToSymbol(NMaxInitPars, &nMaxPhoton, sizeof(int)); //copy the max initial photon number

		cudaMalloc(&gp, NT*sizeof(int)); //allocate memory to record how many photons were provided for each thread
		gc[i].addGPUPointer(gp); //remember to free when exiting
		cudaMemcpyToSymbol(NInitPars, &gp, sizeof(int*));//copy the pointer to GPU

		hhead->setSegments(aListofSegments, gc[i]);
		//copy the instance of the class HeadTranport to GPU
		cudaMalloc(&gp, sizeof(HeadTransport));
		gc[i].addGPUPointer(gp); //remember to free when exiting
		cudaMemcpy(gp, hhead, sizeof(HeadTransport), cudaMemcpyHostToDevice);
		cudaMemcpyToSymbol(headTransport, &gp, sizeof(HeadTransport*));//copy the pointer to GPU
	}

	//clean up
	delete hps;
	delete[] hseg;
	delete hhead;
}