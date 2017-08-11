//VC++ cannot recognize the following cuda systex, so we need to redefine it in cpp files
#define __device__ 
#define __host__

#include "triangle_octree.h"
#include <math.h>

#define CROSS_PRODUCT(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT_PRODUCT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
// Implicit class constructor: the triangle vertices are set to the origin and the material and organ numbers to -1.
triangle_octree::triangle_octree()
{
  set_triangle_members(0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0,
                      -1);
}

// Explicit class constructor:
triangle_octree::triangle_octree(double new_vert0_X,double new_vert0_Y,double new_vert0_Z,
                                 double new_vert1_X,double new_vert1_Y,double new_vert1_Z,
                                 double new_vert2_X,double new_vert2_Y,double new_vert2_Z,
                                 int new_idt)
{
  set_triangle_members(new_vert0_X, new_vert0_Y, new_vert0_Z,
                       new_vert1_X, new_vert1_Y, new_vert1_Z,
                       new_vert2_X, new_vert2_Y, new_vert2_Z,
                       new_idt);
}

// Function that sets the value of all the class members:
void triangle_octree::set_triangle_members
                          (double new_vert0_X,double new_vert0_Y,double new_vert0_Z,
                           double new_vert1_X,double new_vert1_Y,double new_vert1_Z,
                           double new_vert2_X,double new_vert2_Y,double new_vert2_Z,
						   int new_idt)
{
  vert0_X = new_vert0_X;
  vert0_Y = new_vert0_Y;
  vert0_Z = new_vert0_Z;

  vert1_X = new_vert1_X;
  vert1_Y = new_vert1_Y;
  vert1_Z = new_vert1_Z;

  vert2_X = new_vert2_X;
  vert2_Y = new_vert2_Y;
  vert2_Z = new_vert2_Z;

  idt = new_idt;

  //calculate the normal vector
  double v1[3], v2[3];
  v1[0] = vert1_X - vert0_X;
  v1[1] = vert1_Y - vert0_Y;
  v1[2] = vert1_Z - vert0_Z;

  v2[0] = vert2_X - vert0_X;
  v2[1] = vert2_Y - vert0_Y;
  v2[2] = vert2_Z - vert0_Z;

  CROSS_PRODUCT(nv, v1, v2)

  
  double isr = 1.0/sqrt(DOT_PRODUCT(nv,nv));
  nv[0] *= isr;
  nv[1] *= isr;
  nv[2] *= isr;
}



/*
/  Calculate whether the present triangle is inside or overlaps the input box (octree node).
/  This function returns 'true' when the triangle is in the box (overlap) and 'false' when 
/  the triangle is outside.
/  The triangle-Box intersection is calculated using a version of the algorithm and code
/  from:  Tomas Akenine-Möller, 2001, "Fast 3D Triangle-Box Overlap Testing",
/         Journal of Graphics Tools, 6, pp. 29-33.
/         (Code available at: http://jgt.akpeters.com/papers/AkenineMoller01/)
/
/        ********************************************************
/        * AABB-triangle overlap test code                      *
/        * by Tomas Akenine-Möller                              *
/        * Function: int triBoxOverlap(float boxcenter[3],      *
/        *          float boxhalfsize[3],float triverts[3][3]); *
/        * History:                                             *
/        *   2001-03-05: released the code in its first version *
/        *   2001-06-18: changed the order of the tests, faster *
/        *                                                      *
/        * Acknowledgement: Many thanks to Pierre Terdiman for  *
/        * suggestions and discussions on how to optimize code. *
/        * Thanks to David Hunt for finding a ">="-bug!         *
/        ********************************************************
*/
/*======== Defines for the "inside_box" function (triangle-Box intersection): =====*/
#define X 0
#define Y 1
#define Z 2


#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;
/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)              \
     p0 = a*v0[Y] - b*v0[Z];                         \
     p2 = a*v2[Y] - b*v2[Z];                         \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
     rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
     if(min>rad || max<-rad) return 0;
#define AXISTEST_X2(a, b, fa, fb)               \
     p0 = a*v0[Y] - b*v0[Z];                       \
     p1 = a*v1[Y] - b*v1[Z];                         \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
     rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
     if(min>rad || max<-rad) return 0;
/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)              \
     p0 = -a*v0[X] + b*v0[Z];                   \
     p2 = -a*v2[X] + b*v2[Z];                        \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
     rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
     if(min>rad || max<-rad) return 0;
#define AXISTEST_Y1(a, b, fa, fb)               \
     p0 = -a*v0[X] + b*v0[Z];                   \
     p1 = -a*v1[X] + b*v1[Z];                        \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
     rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
     if(min>rad || max<-rad) return 0;
/*======================== Z-tests ========================*/
#define AXISTEST_Z12(a, b, fa, fb)              \
     p1 = a*v1[X] - b*v1[Y];                       \
     p2 = a*v2[X] - b*v2[Y];                         \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
     rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
     if(min>rad || max<-rad) return 0;
#define AXISTEST_Z0(a, b, fa, fb)               \
     p0 = a*v0[X] - b*v0[Y];                    \
     p1 = a*v1[X] - b*v1[Y];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
     rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
     if(min>rad || max<-rad) return 0;
/*=========================================================*/
bool triangle_octree::inside_box(double x_min, double y_min, double z_min,
                                 double x_max, double y_max, double z_max)
{
  /*    Use separating axis theorem to test overlap between triangle and box       */
  /*    need to test for overlap in these directions:                              */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
  /*       we do not even need to test these)                                      */
  /*    2) normal of the triangle                                                  */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin)                           */
  /*       this gives 3x3=9 more tests                                             */

  // (Variables defines as "double" instead of the original "float", Andreu B)
   double v0[3],v1[3],v2[3];
   double min,max,d,p0,p1,p2,rad,fex,fey,fez;
   double normal[3],e0[3],e1[3],e2[3];
   double boxhalfsize[3] = { (0.50*(x_max - x_min)),
                             (0.50*(y_max - y_min)),
                             (0.50*(z_max - z_min)) };
   /* This is the fastest branch on Sun */
   /* move everything so that the boxcenter is in (0,0,0) */
   v0[0]= vert0_X - 0.50*(x_max + x_min);
   v0[1]= vert0_Y - 0.50*(y_max + y_min);
   v0[2]= vert0_Z - 0.50*(z_max + z_min);
   v1[0]= vert1_X - 0.50*(x_max + x_min);
   v1[1]= vert1_Y - 0.50*(y_max + y_min);
   v1[2]= vert1_Z - 0.50*(z_max + z_min);
   v2[0]= vert2_X - 0.50*(x_max + x_min);
   v2[1]= vert2_Y - 0.50*(y_max + y_min);
   v2[2]= vert2_Z - 0.50*(z_max + z_min);

   /* compute triangle edges */
   e0[0]= v1[0]-v0[0];      /* tri edge 0 */
   e0[1]= v1[1]-v0[1];
   e0[2]= v1[2]-v0[2];
   e1[0]= v2[0]-v1[0];      /* tri edge 1 */
   e1[1]= v2[1]-v1[1];
   e1[2]= v2[2]-v1[2];
   e2[0]= v0[0]-v2[0];      /* tri edge 2 */
   e2[1]= v0[1]-v2[1];
   e2[2]= v0[2]-v2[2];

   /* Bullet 3:  */
   /*  test the 9 tests first (this was faster) */
   fex = fabs(e0[X]);
   fey = fabs(e0[Y]);
   fez = fabs(e0[Z]);
   AXISTEST_X01(e0[Z], e0[Y], fez, fey);
   AXISTEST_Y02(e0[Z], e0[X], fez, fex);
   AXISTEST_Z12(e0[Y], e0[X], fey, fex);

   fex = fabs(e1[X]);
   fey = fabs(e1[Y]);
   fez = fabs(e1[Z]);
   AXISTEST_X01(e1[Z], e1[Y], fez, fey);
   AXISTEST_Y02(e1[Z], e1[X], fez, fex);
   AXISTEST_Z0(e1[Y], e1[X], fey, fex);

   fex = fabs(e2[X]);
   fey = fabs(e2[Y]);
   fez = fabs(e2[Z]);
   AXISTEST_X2(e2[Z], e2[Y], fez, fey);
   AXISTEST_Y1(e2[Z], e2[X], fez, fex);
   AXISTEST_Z12(e2[Y], e2[X], fey, fex);

   /* Bullet 1: */
   /*  first test overlap in the {x,y,z}-directions */
   /*  find min, max of the triangle each direction, and test for overlap in */
   /*  that direction -- this is equivalent to testing a minimal AABB around */
   /*  the triangle against the AABB */

   /* test in X-direction */
   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
   if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return false;

   /* test in Y-direction */
   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
   if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return false;

   /* test in Z-direction */
   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
   if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return false;

   /* Bullet 2: */
   /*  test if the box intersects the plane of the triangle */
   /*  compute plane equation of triangle: normal*x+d=0 */
   CROSS_PRODUCT(normal,e0,e1);
   d=-DOT_PRODUCT(normal,v0);  /* plane eq: normal.x+d=0 */
   
   
   // Original code: int planeBoxOverlap(float normal[3],float d, float maxbox[3])
      int q, planeBoxOverlap = 0;
      double vmin[3],vmax[3];
      for(q=X;q<=Z;q++)
      {
         if(normal[q]>0.0)
         {
            vmin[q]=-boxhalfsize[q];
            vmax[q]= boxhalfsize[q];
         }
         else
         {
            vmin[q]= boxhalfsize[q];
            vmax[q]=-boxhalfsize[q];
         }
      }
      if(DOT_PRODUCT(normal,vmin)+d> 0.0) planeBoxOverlap = 0;
      if(DOT_PRODUCT(normal,vmax)+d>=0.0) planeBoxOverlap = 1;
    
   if(!planeBoxOverlap) return false;

   // - If the execution arrives here there is not a plane separating the triangle and the box and,
   //   for the "Separating Axis Theorem" we know that the box and the triangle must overlap:
   
   return true;   /* box and triangle overlaps!! */
}

