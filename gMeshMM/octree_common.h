#ifndef  _OCTREE_COMMON_H
#define  _OCTREE_COMMON_H

using namespace std;
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

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

/** Maximum number of organs that may overlap (one inside the other) inside the triangular phantom. */
#define MAX_OBJ_OVERLAP 4

/** A large distance treated as infinity**/
#define INFINIT_LENGTH (1e9)

struct CrossedTriangle
{
	double dist;
	int ind;
};

struct MatInfo
{
	int id;
	int priority;
};

#endif