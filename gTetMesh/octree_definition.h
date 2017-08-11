
/* -- Program documented for DOXYGEN: www.stack.nl/~dimitri/doxygen/
   JAVADOC_AUTOBRIEF has to be set to YES in the configuration file: comment blocks
   with ' / * * ' start a brief description which ends at the first dot followed by a
   space or new line. */


#ifndef  _OCTREE_DEFINITION_H
#define  _OCTREE_DEFINITION_H


/** Constant that limits the octree maximum recursion level.
*  Despite the private function termination_condition() controls the octree recursion
*  level (depending on the number of triangles), the octree should not be subdivided
*  more than MAX_LEVEL times (or we might easily run out of memory).
*/
#define MAX_LEVEL 11

/** Parameter that controls the termination condition of the octree nodes. If the parameter is
*  negative the octree is subdivided until the number of triangles is smaller than the recursion
*  level (this will produce many small nodes); otherwise the node is subdivided while there are
*  more triangles in the node than value of this parameter.
*/
#define MAX_TRIANGLES  9

/** A large distance treated as infinity**/
#define INFINIT_LENGTH (1e9)


#include "triangle_octree.h"
#include "../Tools/Tools.h"

class octree_node;

/** if maxTriangles<0,  the octree is subdivided until the number of triangle
*  is smaller than the recursion level
*/


/** Class describing a node of the octree structure.
 *  Each node of the octree, sorting the triangles in space, is an instance of this class.
 *  The root node (the triangle geometry bounding box) has level 0 and contains all the other
 *  sub-nodes. The leaves are the ending sub-nodes (highest recursion level) and contain
 *  pointers to the triangles that are located inside the node volume.
 *
 *  In this code the word 'node' will be considered a masculine word, as it is in
 *  Spanish, and the sub-nodes will be refered to as 'father' or 'son'.
 *
 */
class octree_node
{
    private:
         /** Check the condition that decides whether the node is a leaf (final node)
         *  or it has to be further sub-divided.
         *
         *    @return Return 'true' when the node is a leaf or 'false' if it has to be sub-divided.
         */
		 inline bool termination_condition(int &max_level, int &maxTriangles);


         /** Set the 'triangle_list' class member variable assigning pointers to the triangles
         *  inside each node. This function uses the static class variables 'triangle_mesh' and
         *  'triangle_mesh_size', and also the input node --the root node-- to get the triangles
         *  in the father node. The father's 'triangle_list' variable can be deleted after its
         *  triangles have been distributed between the 8 sons calling 'clean_octree'.
         *  The number of triangles assigned to the node is stored in the variable 'num_triangles'.
         *
         *     @param[in]  root_node Pointer to the root node (i.e., the triangle geometry bounding box).
         */
        void sort_triangles(octree_node *root_node);

    public:

        double x_min,   /**< Node lower corner: minimum X. */
               y_min,   /**< Node lower corner: minimum Y. */
               z_min;   /**< Node lower corner: minimum Z. */

        double x_max,   /**< Node upper corner: maximum X. */
               y_max,   /**< Node upper corner: maximum Y. */
               z_max;   /**< Node upper corner: maximum Z. */

        /** Recursion level of this node in the octree hierarchy (0 for the root node). */
        char  level;

        /** Pointers to the 8 sons of the node (NULL for leaves). */
        octree_node* son[8];


        /** Array with pointers to the 6 neighbor nodes (NULL for root).
         *  When all the octree structure has been constructed we can calculate
         *  the nodes that are found behind the 6 node walls.
         *  Some of the neighbors will simply be the "brothers" of the current,
         *  node but other neighbors will have to be found locating points from
         *  the root node. The neighbors can not have a higher level than the
         *  current node (that is, they can not have a smaller size).
         *  The information of the neighbors is essential to easily move the
         *  particles across the octree.
         *  This variable is NULL only for the root node (level=0).
         */
        octree_node* neighbor[6];

        /** List of pointers that point to the triangles that are found inside this sub-node.
         *  The intersection with the triangles only has to be checked at the leaves of the
         *  octree; therefore the triangle_list array will be NULL for all the nodes
         *  that are not leaves.
         */
        triangle_octree** triangle_list;

        /** Number of triangles that are inside --or intersect-- this node.
         *  This variable gives the size of the array 'triangle_list' which is dynamically
         *  allocated after the triangles are read.
         */
        int num_triangles;

         /** Static variable that stores the collection of triangles that define the simulation
          *  geometry. This variable should be initialized by the function that reads the
          *  files with the triangles. */
        static triangle_octree* triangle_mesh;

        /** Size of the triangle mesh (ie, number of triangles). */
        static int triangle_mesh_size;


        /** Amount of nodes for each octree level.
         *  This static variable is updated in the class constructor and stores the
         *  the amount of nodes that have been created for each octree level.
         *  The information in this variable can be used to optimize the octree structure.
         */
        static int amount_nodes[MAX_LEVEL+1];

        /** Amount of leaves in each level of the octree.
         *  This static variable stores the total number of leaves that are found in
         *  each level of the octree.
         *  This information is also useful to optimize the octree structure.
         */
        static int amount_leaves[MAX_LEVEL+1];

        // Public functions:
		static octree_node* create_octree(vector<triangle_octree>& triangles, int input_max_level, int maxTriangles, vector<double> bound = vector<double>());

		static void delete_octree(octree_node* root);

        /** Default class constructor.
         *  The node lower and upper corners and its level are initialized to 0. The
         *  function 'set_node_parameters' has to be called to correctly initialize the node data.
         */
        octree_node();


        /** Explicit class constructor.
         *      @param[in] x0  X coordinate of the node lower corner.
         *      @param[in] y0  Y coordinate of the node lower corner.
         *      @param[in] z0  Z coordinate of the node lower corner.
         *      @param[in] x0  X coordinate of the node upper corner.
         *      @param[in] y0  Y coordinate of the node upper corner.
         *      @param[in] z0  Z coordinate of the node upper corner.
         *      @param[in] level0   Level of the created sub-node.
         */
        octree_node(double x0,double y0,double z0,
                    double x1,double y1,double z1, char level0);

        /** Class destructor.
         *  It recursively deletes all the sons of the node. */
        ~octree_node();

        /** Sets the value of the node members.
         *      @param[in] x0  X coordinate of the node lower corner.
         *      @param[in] y0  Y coordinate of the node lower corner.
         *      @param[in] z0  Z coordinate of the node lower corner.
         *      @param[in] x0  X coordinate of the node upper corner.
         *      @param[in] y0  Y coordinate of the node upper corner.
         *      @param[in] z0  Z coordinate of the node upper corner.
         *      @param[in] level0   Level of the created sub-node.
         */
        void set_node_parameters(double x0,double y0,double z0,
                                 double x1,double y1,double z1, char level0);

		bool projectIn(double& px, double& py, double& pz, double& pu, double& pv, double& pw);

		// Project the particle into the octree in a straight line, return the corresponding tetrahedral index.
		// If it doesn't collide any tetrahedral, return -1. Call it when this node is root.
		int lineIn(double& x, double& y, double& z, double& u, double& v, double& w);

        /** Function to search from the input node (the octree root usually) for the sub-node that
         *  contains the input point and that has the maximum level specified (this is needed in
         *  function 'set_neighbors' because the stored neighbors can not have a higher level, ie
         *  smaller size). This functions checks if the particle is outside the node.
         *
         *     @param[in]  px  X coordinate of the point where the node has to be searched.
         *     @param[in]  py  Y coordinate of the point where the node has to be searched.
         *     @param[in]  pz  Z coordinate of the point where the node has to be searched.
         *     @param[in]  level_limit  Maximum level of the returned node (the sons of
         *                    the node at this level are not checked). The default value is -1,
         *                    which means that the limit will have no effect.
         *      @return  The return value is a pointer to the sub-node that contains the input
         *               point, or NULL if the point is not inside the current node.
         */
        octree_node* get_node(double &px,double &py,double &pz, char level_limit=char(-1));


        /** Function to search from the input node for the sub-node that contains the input point.
         *  This version of 'get_node' does not set a limit for the maximum level of the returned
         *  sub-node and it does not check if the point is actually inside the node. This will speed up
         *  a little the part of the simulation where particles pass from one octree leave its neighbor.
         *
         *     @param[in]  px  X coordinate of the point where the node has to be searched.
         *     @param[in]  py  Y coordinate of the point where the node has to be searched.
         *     @param[in]  pz  Z coordinate of the point where the node has to be searched.
         *      @return  The return value is a pointer to the sub-node that contains the input point.
         *               If the point is not inside the current node the returned pointer will be wrong.
         */
        octree_node* get_node_fast(double &px,double &py,double &pz);


        /** Generate the octree structure and distribute the triangles.
         *  This subroutine finds the triangles that intersect the node (using the
         *  private function 'sort_triangles'), assigns the triangles to the member
         *  variable 'triangle_list', checks the 'termination_condition', and --if the node
         *  has to be subdivided-- creates the instances of the 8 sub-nodes and repeat
         *  the previous operations for each sub-node.
         *
         *  Typically this function is only called once for the "root" octree_node instance
         *  and then the complete octree structure is generated below this root node.
         *
         *  The triangles must be stored in the static class variable *triangle_mesh
         *  before calling this function.
         *
         *      @param[in]  input_max_level  Maximum number of subdivisions allowed for the
         *                             octree structure. This value is read from an input file
         *                             and passed to the function 'termination_condition'.
         *      @param[in]  root_node  Pointer to the root node (i.e., the triangle geometry bounding box).
         */
		void create_subnodes_with_triangles(int &input_max_level, int &maxTriangles, octree_node *root);



        /** Delete the unnecessary instances of 'triangle_list' in the nodes that are not leaves. */
        /// (NOTE: I could simplify the octree structure merging equivalent neighbor nodes with few triangles.)
        void clean_octree();


        /** Calculate the distance and neighbor number of the nearest node,
         *  for the input position and inverse value of the direction cosines.
         *
         *      @param[in]  x   X coordinate of the particle position.
         *      @param[in]  y   Y coordinate of the particle position.
         *      @param[in]  z   Z coordinate of the particle position.
         *      @param[in]  invers_vx  inverse of the X coordinate of the particle direction vector.
         *      @param[in]  invers_vy  inverse of the Y coordinate of the particle direction vector.
         *      @param[in]  invers_vz  inverse of the Z coordinate of the particle direction vector.
         *      @param[out] num_neighbor  Neighbor number (1->6) of the nearest node.
         *                      This value will be used to find the next node using the
         *                      member variable octree_node* neighbor[6].
         *      @return     Distance to the nearest node.
         */
        inline double get_wall_distance(double &x,double &y,double &z,
                        double &vx,double &vy,double &vz,
                        int &num_neighbor);

        /** Set the 'neighbor' class variable: pointers to the 6 neighbor nodes of each node.
         *  The neighbor nodes are ordered as: X++, Y+, Z ("front,back,right,left,top,bottom").
         *
         *     @param[in]  root_node  Pointer to the root node (i.e., the triangle geometry bounding box).
         */
        void set_neighbors(octree_node* root_node);

		int get_triangle_intersection(double &x, double &y, double &z, double& u, double& v, double& w);

		int step_octree(double& x, double& y, double& z, double u, double v, double w);


        // Friend functions:

        /** This friend function of the octree_node class draws the octree structure at the z plane given
         *  in the input file, in encapsulated postscript format (one square drawn for
         *  each octree leaf on the plane). The output image size is approximately a vertical letter page
         *  and it is automatically scaled to the input node size (usually the root node).
         *  If the input z plane is exactly z==1000, the full 3D octree structure is drawn
         *  in EPS format, one image for each plane. This option can be time-consuming and require
         *  large storage space for all the EPS files. The generated EPS images can be easily
         *  converted to an animated GIF using the free software ImageMagick (in Linux):
         *     - convert -delay 25 octree_3D_*.eps octree_3D.gif".
         *  After converting the EPS to PNG format with ImageMagick, it is also simple to create
         *  an AVI video of the octree structure with the mencoder program:
         *     - mencoder "mf://*.png" -mf fps=18 -o octree_3D_Level9_mencoder.avi -ovc lavc -lavcopts vcodec=mpeg4
         *
         *     @param[in] octree   Pointer to the octree root node.
         *     @param[in] octree_max_level Maximum octree subdivision level: useful information
         *                                 written in the EPS header.
         *     @param[in] z_plane  Z plane where the octree is to be plotted
         *     @param[in] file_name  Name of the output postscript file.
         *     @param[in] first_call  Optional variable used internally to detect the first
         *                            call to the function (the function is called recursively
         *                            and only the first call must write the PS header).
         */
         static void draw_octree_ps(octree_node *octree, int &octree_max_level, double &z_plane,
                                   const char *file_ps, bool first_call=true);


        /** Calculate the mean quantity (and standard deviation) of triangles inside the leaves with
         *  maximum level, and the amount of empty leaves. This statistical information gives an idea
         *  of how the triangles are sorted in the octree. The best situation is a few triangles in
         *  the final leaves and a lot of empty leaves (no triangle intersection there).
         *
         *     @param[in]  octree       Pointer to the current octree node.
         *     @param[in]  octree_max_level Input maximum octree subdivision level.
         *     @param[out] mean_tri     Mean quantity of triangles inside the final leaf.
         *     @param[out] sigma_tri    Standard deviation of the number of triangles.
         *     @param[out] leaves_full  Amount of leaves with triangles inside in the full octree.
         *     @param[out] leaves_empty Amount of empty leaves in the full octree.
         */
         static void triangles_per_leaf(octree_node *octree, int &octree_max_level, double &mean_tri,
                       double &sigma_tri, int &leaves_full, int &leaves_empty, bool first_call=true);
};



// /////////////////////////////////////////////////////////////////////////////

// --------- Inline function definition: ----------------------

//  Function to calculate the distance (return value) and neighbor number (output argument) of the
//   nearest node, for the input position and inverse value of the direction cosines.
//   (the neighbor nodes are ordered as: X++, Y+, Z ("front,back,right,left,top,bottom"))
inline double octree_node::get_wall_distance(double &x,double &y,double &z,
                             double &vx,double &vy,double &vz,
                             int &num_neighbor)
{
    double dist, dist_tmp;

    if (vx > 0.0)
    {
        dist=((x_max-x)/vx);
        num_neighbor=0;
    }
    else
    {
        dist=((x_min-x)/vx);
        num_neighbor=1;
    }

    if (vy > 0.0)
    {
        dist_tmp=((y_max-y)/vy);
        if (dist_tmp<dist)
        {
            dist=dist_tmp;
            num_neighbor=2;
        }
    }
    else
    {
        dist_tmp=((y_min-y)/vy);
        if (dist_tmp<dist)
        {
            dist=dist_tmp;
            num_neighbor=3;
        }
    }

    if (vz > 0.0)
    {
        dist_tmp=((z_max-z)/vz);
        if (dist_tmp<dist)
        {
            dist=dist_tmp;
            num_neighbor=4;
        }
    }
    else
    {
        dist_tmp=((z_min-z)/vz);
        if (dist_tmp<dist)
        {
            dist=dist_tmp;
            num_neighbor=5;
        }
    }

	if (dist < 0) //for debug
	{
		cout << "dis < 0 !!!" << endl;
		getchar();
		exit(1);
	}
    return dist;
}


/* Check the condition that decides whether the node is a leaf (final node)
*  or it has to be further sub-divided.
*
*    @return Return 'true' when the node is a leaf or 'false' if it has to be sub-divided.
*/
inline bool octree_node::termination_condition(int &input_max_level, int& maxTriangles)
{
	if (maxTriangles >= MAX_TRIANGLES) maxTriangles = MAX_TRIANGLES;
	if (input_max_level >= MAX_LEVEL) input_max_level = MAX_LEVEL;

    // ** Do not sub-divide nodes which have no triangles or when the input maximum level is found:
     if (level >= input_max_level || num_triangles == 0 )  return true;

     if(maxTriangles<0)
     {                  
      // ** "T=L condition": subdivide if there are less triangles than the node level:
       if (num_triangles <= level) return true;
     }
     else
     {
      // ** Don't sub-divide the node if it has = or < triangles than the limit set above:
       if (num_triangles <= maxTriangles) return true;
     }

     return false;  // If all the conditions are FALSE, return FALSE ==> sub-divide the node!
}


#endif

