#ifndef  _OCTREE_GPU_H
#define  _OCTREE_GPU_H
/** This head file will be included by gPENELOPE **/
#include "triangle_octree.h"

class goctree_node
{
public:
	// the box position
	double x_min, y_min, z_min;
	double x_max, y_max, z_max;

	short  level;

	/** Pointers to the 8 sons of the node (NULL for leaves). */
	goctree_node* son[8];

	goctree_node* neighbor[6];

	triangle_octree* triangle_list; //the triangle array belong to this box

	int num_triangles; // size of above array

	/*********************************simulation functions************************************************/
	__device__ bool projectIn(double& px, double& py, double& pz, double& pu, double& pv, double& pw);

	__device__ goctree_node* get_node_fast(double &px, double &py, double &pz);

	// project the particle into the octree in a straight line, return the corresponding node
	// if it doesn't collide with the octree, return NULL. Call it when this node is root.
	__device__ goctree_node* lineIn(double& x, double& y, double& z, double& u, double& v, double& w);

	__device__ double get_wall_distance(double &x, double &y, double &z, double &vx, double &vy, double &vz, int &num_neighbor);

	__device__ double get_triangle_intersection(double &x, double &y, double &z, double& u, double& v, double& w, int& ind);

	__device__ goctree_node* step_octree(double ds, double& x, double& y, double& z, double u, double v, double w, MatInfo matIDs[MAX_OBJ_OVERLAP], int& cMatID);
	
	__device__ goctree_node* get_node(double &px, double &py, double &pz, short level_limit);
	
	__device__ void set_neighbors(goctree_node* root_node);

	__device__ void deleteOctree(goctree_node* node);
};

// instructions: call createCPUOctree first, then call createGPUOctree in each GPU context to copy the data to GPU,
// and deleteCPUOctree to delete the temporary data. In the end, you need to call deleteGPUOctree to release the resource.

void createCPUOctree(const string& fname, int input_max_level, int maxTriangles, vector<int> matList);

void deleteCPUOctree();

goctree_node* creatGPUOctree();

void deleteGPUOctree(goctree_node* root);
#endif