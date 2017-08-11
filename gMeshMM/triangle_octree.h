#ifndef  _TRIANGLE_OCTREE_H
#define  _TRIANGLE_OCTREE_H
#include "octree_common.h"

///////////////////////////////////////////////////////////////////////////////////////////


// ** Macros used in function "bool triangle_octree::intersect":

// -Macro to calculate the cross product of two vectors ('##' pastes together two tokens):
#define CROSS(dest,v1,v2)                   \
						        { dest##_X = v1##_Y*v2##_Z - v1##_Z*v2##_Y; \
          dest##_Y = v1##_Z*v2##_X - v1##_X*v2##_Z; \
          dest##_Z = v1##_X*v2##_Y - v1##_Y*v2##_X; }

// -Macro to calculate the dot product of two vectors:
#define DOT(v1,v2) (v1##_X*v2##_X + v1##_Y*v2##_Y + v1##_Z*v2##_Z)

// -Calculate the vector that goes from v1 to v2:
#define SUB(dest,v1,v2)         \
						        { dest##_X = v1##_X - v2##_X; \
          dest##_Y = v1##_Y - v2##_Y; \
          dest##_Z = v1##_Z - v2##_Z; }

/** Declaration of a triangle object that will be ray-traced in the Monte Carlo simulation. */
class triangle_octree
{
    public:

      // --Public variables:
       double vert0_X;  /**< X coordinate of the first triangle vertex.  */
       double vert0_Y;  /**< Y coordinate of the first triangle vertex.  */
       double vert0_Z;  /**< Z coordinate of the first triangle vertex.  */

       double vert1_X;  /**< X coordinate of the second triangle vertex. */
       double vert1_Y;  /**< Y coordinate of the second triangle vertex. */
       double vert1_Z;  /**< Z coordinate of the second triangle vertex. */

       double vert2_X;  /**< X coordinate of the third triangle vertex.  */
       double vert2_Y;  /**< Y coordinate of the third triangle vertex.  */
       double vert2_Z;  /**< Z coordinate of the third triangle vertex.  */

	   double nv[3]; // normal vector
	   MatInfo mat;


      // --Public functions:

       /** Implicit constructor for the triangle objects.
        *  The triangle vertices are set to the origin and the material and organ numbers to -1.
        */
        triangle_octree();

        /** Explicit class constructor.        */
        triangle_octree(double new_vert0_X,double new_vert0_Y,double new_vert0_Z,
                        double new_vert1_X,double new_vert1_Y,double new_vert1_Z,
                        double new_vert2_X,double new_vert2_Y,double new_vert2_Z,
                        short int new_organ, short int new_material);


       /** Function that sets the value of all the class members. */
       void set_triangle_members (double new_vert0_X,double new_vert0_Y,double new_vert0_Z,
               double new_vert1_X,double new_vert1_Y,double new_vert1_Z,
               double new_vert2_X,double new_vert2_Y,double new_vert2_Z,
               short int new_organ, short int new_material);


       /** Calculate whether the input ray intersects this triangle or not using the Moller-Trumbore
        *  algorithm. The original algorithm and code was published in  J. Graphics Tools 2, p 21-28
        *  (1997) and is available at http://www.acm.org/jgt/papers/MollerTrumbore97/ .
        *
        *  We have adapted the algorithm to our problem:
        *    - A negative or zero intersection distance is considered no intersection.
        *    - This "bool" function returns 'true' when there is intersection and 'false' when not.
        *    - The ray position and directions are not arrays.
        *    - The vertices are read as member variables and not entered as input parameters.
        *    - The barycentric coordinates of the intersection point are not calculated nor returned.
        *    - All the triangles can be intersected through both sides (non culling version).
        *
        *    @return The returned integer value is 1 if there is an intersection or 0 otherwise.
        *            A negative or zero intersection distance is considered no intersection.
        *    @param[in]  orig_X    X coordinate of the incoming particle. (/track/ x).
        *    @param[in]  orig_Y    Y coordinate of the incoming particle. (/track/ y).
        *    @param[in]  orig_Z    Z coordinate of the incoming particle. (/track/ z).
        *    @param[in]  dir_X     X coordinate of the particle direction (/track/ u).
        *    @param[in]  dir_Y     Y coordinate of the particle direction (/track/ v).
        *    @param[in]  dir_Z     Z coordinate of the particle direction (/track/ w).
        *    @param[out] distance  Distance to the intersection point (if it exists).
        *
        */
	   /**
	   *
	   *  Calculate whether the input ray intersects this triangle or not using the Moller-Trumbore
	   *  algorithm. The original algorithm and code was published in  J. Graphics Tools 2, p 21-28
	   *  (1997) and is available at http://www.acm.org/jgt/papers/MollerTrumbore97/ .
	   *
	   *  The algorithm used here has been adapted to our problem:
	   *
	   *    - A negative or zero intersection distance is considered no intersection!!!
	   *
	   *    - This "bool" function returns 'true' when there is intersection and 'false' when not.
	   *    - The ray position and directions are not arrays.
	   *    - The vertices are read as member variables and not entered as input parameters.
	   *    - The barycentric coordinates of the intersection point are not calculated nor returned.
	   *    - All the triangles can be intersected through both sides (non culling version).
	   *
	   *    @return The returned integer value is 1 if there is an intersection or 0 otherwise.
	   *    @param[in]  orig_X    X coordinate of the incoming particle. (/track/ x).
	   *    @param[in]  orig_Y    Y coordinate of the incoming particle. (/track/ y).
	   *    @param[in]  orig_Z    Z coordinate of the incoming particle. (/track/ z).
	   *    @param[in]  dir_X     X coordinate of the particle direction (/track/ u).
	   *    @param[in]  dir_Y     Y coordinate of the particle direction (/track/ v).
	   *    @param[in]  dir_Z     Z coordinate of the particle direction (/track/ w).
	   *    @param[out] distance  Distance to the intersection point (if it exists).
	   *
	   */
	   __device__ __host__ bool intersect(double &orig_X, double &orig_Y, double &orig_Z,
		   double &dir_X, double &dir_Y, double &dir_Z,
		   double &distance)
	   {
		   // !!DeBuG!! const double EPSILON = 1.0e-6;
		   const double EPSILON = 1.0e-8;
		   double edge1_X, edge1_Y, edge1_Z;  // Triangle edge 1.
		   double edge2_X, edge2_Y, edge2_Z;  // Triangle edge 2.
		   double tvec_X, tvec_Y, tvec_Z;   // Distance from triangle vertex 0 to ray origin.
		   double pvec_X, pvec_Y, pvec_Z;   // [dir * edge2]
		   double qvec_X, qvec_Y, qvec_Z;   // [tvec X edge2]
		   double det, inv_det;               // [edge1*pvec]
		   double uu, vv;                     // Intersection point in the triangle barycentric coordinates.

		   // find vectors for two edges sharing vert0
		   SUB(edge1, vert1, vert0);
		   SUB(edge2, vert2, vert0);

		   // begin calculating determinant - also used to calculate U parameter
		   CROSS(pvec, dir, edge2);

		   // if determinant is near zero, ray lies in plane of triangle
		   det = DOT(edge1, pvec);

		   if (det > -EPSILON && det < EPSILON)
			   return false;
		   inv_det = 1.0 / det;

		   // calculate distance from vert0 to ray origin
		   SUB(tvec, orig, vert0);

		   // calculate U parameter and test bounds
		   uu = DOT(tvec, pvec) * inv_det;  // I don't return the barycentric coordinates, but I need the correct sign of uu!

		   if (uu < 0.0 || uu > 1.0)
			   return false;

		   // prepare to test V parameter
		   CROSS(qvec, tvec, edge1);

		   // calculate V parameter and test bounds
		   vv = DOT(dir, qvec) * inv_det;
		   if (vv < 0.0 || (uu + vv) > 1.0)
			   return false;

		   // calculate the distance ray intersects triangle
		   distance = DOT(edge2, qvec) * inv_det;

		   if (distance > 0.0) return true;

		   return false;  // Negative or zero intersection distance: no real intersection occurs.
	   }


       /** Calculate whether the present triangle is inside or overlaps the input box (octree node).
        *  This function returns 'true' when the triangle is in the box (overlap) and 'false' when
        *  the triangle is outside.
        *
        *  The triangle--Box intersection is calculated using a version of the algorithm and code
        *  from:  Tomas Akenine-Möller, 2001, "Fast 3D Triangle-Box Overlap Testing",
        *         Journal of Graphics Tools, 6, pp. 29-33.
        *         (Code available at: http://jgt.akpeters.com/papers/AkenineMoller01/)
        *
        *    @return Return 'true' when the triangle is in the box.
        *    @param[in]  x_min   X coordinate of the box lower corner.
        *    @param[in]  y_min   Y coordinate of the box lower corner.
        *    @param[in]  z_min   Z coordinate of the box lower corner.
        *    @param[in]  x_max   X coordinate of the box upper corner.
        *    @param[in]  y_max   Y coordinate of the box upper corner.
        *    @param[in]  z_max   Z coordinate of the box upper corner.
        */
       bool inside_box(double x_min, double y_min, double z_min,
                       double x_max, double y_max, double z_max);
};





#endif

