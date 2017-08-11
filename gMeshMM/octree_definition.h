
/* -- Program documented for DOXYGEN: www.stack.nl/~dimitri/doxygen/
   JAVADOC_AUTOBRIEF has to be set to YES in the configuration file: comment blocks
   with ' / * * ' start a brief description which ends at the first dot followed by a
   space or new line. */


#ifndef  _OCTREE_DEFINITION_H
#define  _OCTREE_DEFINITION_H

#include "triangle_octree.h"
#include "octree_common.h"

class octree_node;

struct OBJ3D //class to read 3d mesh object. Currently supported format: *.OBJ
{
	struct VERTEX
	{
		double x, y, z;
		int obj_id; //object index to retrieve OBJINFO
	};
	struct FACEINDEX
	{
		int v[3];
	};
	struct OBJINFO
	{
		//the name format in object files must be name_matid(_priority). If you leave _priority empty, it will be uniquely assigned to a number >1000
		string name;
		int mat_id; //the material index defined in PENELOPE data base
		//In case the objects overlap, the mat_id will be the one has the highest priority.
		int priority;
	};
	vector<VERTEX> vs; //vertex array
	vector<FACEINDEX> fs; //triangle array
	vector<OBJINFO> info; //object information array
	VERTEX maxv, minv; //the boundary coordinates
	bool load(string fname)
	{
		clear();
		ifstream ifs(fname);//cube bunny Eight
		string s;
		string head;
		int default_priority = 1001;
		bool hasVT = false, hasVN = false;
		int vertexStart = 0;
		while (getline(ifs, s))
		{
			if (s.length() < 2) continue;
			if (s[0] == 'v')
			{
				if ('t' == s[1]) hasVT = true;
				else if ('n' == s[1]) hasVN = true;
				else //the vertex
				{
					std::istringstream in(s);
					VERTEX v;
					in >> head >> v.x >> v.y >> v.z;
					vs.push_back(v);
				}
			}
			else if ('f' == s[0] && ' ' == s[1])
			{
				for (int k = (int)s.size() - 1; k >= 0; k--) //replace all '/' to ' '
				{
					if (s[k] == '/') s[k] = ' ';
				}
				FACEINDEX fid;
				int useless;

				std::istringstream in(s);
				in >> head;
				for (int i = 0; i < 3; ++i)
				{
					in >> fid.v[i];
					--fid.v[i];
					if (hasVT) in >> useless;
					if (hasVN) in >> useless;
				}
				fs.push_back(fid);
			}
			else if ('g' == s[0] && ' ' == s[1])
			{
				std::istringstream in(s);
				in >> head; //pass the head "s "
				in >> head;
				//parse the name
				size_t pos = head.find('_');
				int NUnderline = 0;
				while (pos != string::npos)
				{
					head[pos] = ' ';
					pos = head.find('_');
					++NUnderline;
				}
				if (NUnderline != 1 && NUnderline != 2)
				{
					cout << "Error: The object name is not correctly assigned" << endl;
					getchar();
					return false;
				}
				std::istringstream inf(head);
				OBJINFO objinf;
				inf >> objinf.name;
				inf >> s;

				try
				{
					objinf.mat_id = std::stoi(s);

					if (2 == NUnderline)
					{
						inf >> s;
						objinf.priority = std::stoi(s);
					}
					else
					{
						objinf.priority = default_priority;
						++default_priority;
					}
				}
				catch (const std::invalid_argument& ia)
				{
					std::cerr << "Invalid argument: " << ia.what() << '\n';
					getchar();
					return false;
				}

				info.push_back(objinf);
				//assign the points to the corresponding object
				int vertexEnd = (int)vs.size();
				for (int i = vertexStart; i < vertexEnd; ++i) vs[i].obj_id = (int)info.size() - 1;
				vertexStart = vertexEnd;
			}
		}

		//calculate the boundary
		maxv.x = minv.x = vs[0].x;
		maxv.y = minv.y = vs[0].y;
		maxv.z = minv.z = vs[0].z;
		int nv = (int)vs.size();
		for (int i = 1; i < nv; ++i)
		{
			maxv.x = max(maxv.x, vs[i].x);
			minv.x = min(minv.x, vs[i].x);
			maxv.y = max(maxv.y, vs[i].y);
			minv.y = min(minv.y, vs[i].y);
			maxv.z = max(maxv.z, vs[i].z);
			minv.z = min(minv.z, vs[i].z);
		}
		//enlarge the boundary by a little bit
		const double EPS = 1e-7;
		maxv.x += EPS;
		maxv.y += EPS;
		maxv.z += EPS;
		minv.x -= EPS;
		minv.y -= EPS;
		minv.z -= EPS;

		//check if there are overlapped priority numbers
		int nobj = (int)info.size();
		for (int i = 0; i < nobj; ++i)
			for (int j = i + 1; j < nobj; ++j)
			{
				if (info[i].priority == info[j].priority)
				{
					cout << "Error: The objects " << info[i].name << " and " << info[j].name << " have a duplicated priority!" << endl;
					getchar();
					return false;
				}
			}
		return true;
	}
	void clear()
	{
		vs.resize(0);
		fs.resize(0);
		info.resize(0);
	}
};

/** if maxTriangles<0,  the octree is subdivided until the number of triangle
*  is smaller than the recursion level
*/
octree_node* creatOctree(const string& fname, int input_max_level, int maxTriangles, vector<int> matList = vector<int>());

void deleteOctree(octree_node* root);

struct OctreeData
{
	//distance regarding the crossed triangles
	CrossedTriangle ds_triangles[MAX_TRIANGLES];
	int ntri;
	double ds_min_triangle;
	int ind_min_triangle;

	//wall distance
	double ds_wall;
	int num_neighbor;

	MatInfo matIDs[MAX_OBJ_OVERLAP];

	void resetMatIDs()
	{
		for (int i = 0; i < MAX_OBJ_OVERLAP; ++i) matIDs[i].id = -1;
	}
	void update(octree_node* node, double& x, double& y, double& z, double& u, double& v, double& w);

	void update();
};

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

inline	bool projectIn(double& px, double& py, double& pz, double& pu, double& pv, double& pw)
		{
			const double Delta = 1e-5;
			double sx, sy, sz;
			double ds = 0;
			if (px < x_min || px >= x_max || py < y_min || py >= y_max || pz < z_min || pz >= z_max) //initial position lays outside the phantom
			{
				if (px < x_min && pu > 0)
				{
					sx = x_min;
					ds = (x_min - px) / pu;
					sy = py + ds*pv;
					sz = pz + ds*pw;
					if (y_min <= sy&&sy < y_max && z_min <= sz&&sz < z_max)
					{
						px = sx;
						py = sy;
						pz = sz;
						return true;
					}
				}
				else if (px >= x_max && pu < 0)
				{
					sx = x_max - Delta;
					ds = (sx - px) / pu;
					sy = py + ds*pv;
					sz = pz + ds*pw;
					if (y_min <= sy &&sy < y_max && z_min <= sz&&sz < z_max)
					{
						px = sx;
						py = sy;
						pz = sz;
						return true;
					}
				}
				else if (py < y_min&&pv>0)
				{
					sy = y_min;
					ds = (sy - py) / pv;
					sx = px + ds*pu;
					sz = pz + ds*pw;
					if (x_min <= sx &&sx < x_max && z_min <= sz&&sz < z_max)
					{
						px = sx;
						py = sy;
						pz = sz;
						return true;
					}
				}
				else if (py >= y_max && pv < 0)
				{
					sy = y_max - Delta;
					ds = (sy - py) / pv;
					sx = px + ds*pu;
					sz = pz + ds*pw;
					if (x_min <= sx &&sx < x_max && z_min <= sz&&sz < z_max)
					{
						px = sx;
						py = sy;
						pz = sz;
						return true;
					}
				}
				else if (pz < z_min && pw > 0)
				{
					sz = z_min;
					ds = (sz - pz) / pw;
					sy = py + ds*pv;
					sx = px + ds*pu;
					if (y_min <= sy &&sy < y_max && x_min <= sx&&sx < x_max)
					{
						px = sx;
						py = sy;
						pz = sz;
						return true;
					}
				}
				else if (pz >= z_max && pw < 0)
				{
					sz = z_max - Delta;
					ds = (sz - pz) / pw;
					sy = py + ds*pv;
					sx = px + ds*pu;
					if (y_min <= sy &&sy < y_max && x_min <= sx&&sx < x_max)
					{
						px = sx;
						py = sy;
						pz = sz;
						return true;
					}
				}
				return false;
			}
			else return true;
		}

		// project the particle into the octree in a straight line, return the corresponding node
		// if it doesn't collide with the octree, return NULL. Call it when this node is root.
		octree_node* lineIn(double& x, double& y, double& z, double& u, double& v, double& w);

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



        /** Return the distance and a pointer to the nearest triangle intersected.
         *  All the triangles from the node 'triangle_list[]' are checked.
         *  A pointer to the crossed triangle is also stored in the private class static variable
         *  'last_triangles'. This triangle will not be checked in the next step to avoid intersecting
         *  it forever. Before transporting a new particle the variable has to be cleared (it is
         *  automatically done when calling 'init_flight_log' or 'get_flight_log').
         *
         *      @param[in]  x   X coordinate of the particle position.
         *      @param[in]  y   Y coordinate of the particle position.
         *      @param[in]  z   Z coordinate of the particle position.
         *      @param[in]  u   X coordinate of the particle direction vector.
         *      @param[in]  v   Y coordinate of the particle direction vector.
         *      @param[in]  w   Z coordinate of the particle direction vector.
         *      @param[out] triangle_array  The array storing the crossed triangles.
         *      @param[out] triangle_num  The number of crossed triangles.
         *      @return  minimum distance
         */
		int get_triangle_intersection(double &x, double &y, double &z, double& u, double& v, double& w,
			CrossedTriangle* triangle_array, int& triangle_num);


        /** Moves the current particle (in common /track/) across the octree and the triangular mesh.
         *  The particle is stopped at the interaction point or at a triangle interface or whenever the
         *  particle goes out of the octree. The particle position, material, and organ (ibody) are
         *  updated before returning.
         *
         *     @param[in]  ds     Sampled distance to the next interaction.
         *     @param[out] dsef   Effective distance jumped.
         *     @param[out] ncross Switch that is 0 when no interface (triangle) has been crossed, or 1 otherwise.
         *     @return     A pointer to the ending node is returned (NULL when the particle escapes the octree).
         *                 Multiple nodes may be traversed in a single leap.
         *
         */
		 octree_node* step_octree(double ds, double& x, double& y, double& z, double u, double v, double w, OctreeData& od, int& cMatID, bool& dirChanged);


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
    // ** Never create nodes below the programmed MAX_LEVEL (too many levels will need need too much memory):
     if (level == MAX_LEVEL) return true;   // TRUE ==> do not sub-divide.

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

