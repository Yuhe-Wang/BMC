#include "octree_gpu.h"
#include "octree_definition.h"

octree_node* octree_root = NULL;

goctree_node* copyOctree(octree_node* node) //note this function doesn't establish the neighor relation
{
	//resize memory for current node
	goctree_node* cnode = NULL;
	cudaMalloc(&cnode, sizeof(goctree_node));
	if (cnode == NULL) printf("fail to resize GPU memory for goctree_node");
	//copy the data
	goctree_node node_h;
	node_h.x_min = node->x_min;
	node_h.x_max = node->x_max;
	node_h.y_min = node->y_min;
	node_h.y_max = node->y_max;
	node_h.z_min = node->z_min;
	node_h.z_max = node->z_max;

	node_h.level = node->level;

	if (node->triangle_list != NULL && node->num_triangles > 0)
	{
		node_h.num_triangles = node->num_triangles;
		//resize memory in CPU for the triangles
		triangle_octree* ptri = new triangle_octree[node->num_triangles];
		for (int i = 0; i < node->num_triangles; ++i)
		{
			ptri[i] = *(node->triangle_list[i]);
		}
		//resize memory in GPU for the triangles
		cudaMalloc(&node_h.triangle_list, node->num_triangles * sizeof(triangle_octree));
		if (node_h.triangle_list == NULL)
		{
			printf("fail to resize GPU memory for triangle_list\n");
			getchar();
			exit(1);
		}
		//copy the content
		cudaMemcpy(node_h.triangle_list, ptri, node->num_triangles * sizeof(triangle_octree), cudaMemcpyHostToDevice);
		delete[] ptri; //release the CPU memory
	}

	//need to create the son
	if (node->son[0] != NULL)
	{
		for (int i = 0; i < 8; ++i) node_h.son[i] = copyOctree(node->son[i]);
	}
	else
	{
		for (int i = 0; i < 8; ++i) node_h.son[i] = NULL;
	}

	//copy the content of node_h to GPU
	cudaMemcpy(&node_h, cnode, sizeof(goctree_node), cudaMemcpyHostToDevice);
	return cnode;
}

/************************************  Stsrt: kernel called by one thread **********************************/
__global__ void kernel_setNeighbors(goctree_node* groot) //GPU thread to set the octree neighbor pointers
{
	groot->set_neighbors(groot);
}

__global__ void kernel_deleteOctree(goctree_node* groot) //GPU thread to set the octree neighbor pointers
{
	groot->deleteOctree(groot);
}
/************************************* End: kernel called by one thread ************************************/

/********************************* Start: interface functions ************************************************/
void createCPUOctree(const string& fname, int input_max_level, int maxTriangles, vector<int> matList)
{
	if (octree_root == NULL) octree_root = creatOctree(fname, input_max_level, maxTriangles, matList);
}

void deleteCPUOctree()
{
	if (octree_root != NULL) deleteOctree(octree_root);
	octree_root = NULL;
}

goctree_node* creatGPUOctree()
{
	goctree_node* groot = copyOctree(octree_root);
	kernel_setNeighbors << <1, 1 >> >(groot); //set the neighbor pointers
	cudaDeviceSynchronize();
	return groot;
}

void deleteGPUOctree(goctree_node* groot)
{
	kernel_deleteOctree<< <1, 1 >> >(groot);
	cudaDeviceSynchronize();
}

/********************************* End: interface functions ************************************************/


/************************************* SIMULATION¡¡FUNCTIONS ******************************************/
__device__ bool goctree_node::projectIn(double& px, double& py, double& pz, double& pu, double& pv, double& pw)
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

__device__ goctree_node* goctree_node::get_node_fast(double &px, double &py, double &pz)
{
	goctree_node* nd = this;
	while (true)
	{
		if (nd->son[0] == NULL)
		{
			return nd;     //  Return a pointer to the current leaf
		}
		else
		{
			int num_son = 0;
			if (pz < (0.5*(nd->z_min + nd->z_max)))
			{
				if (py < (0.5*(nd->y_min + nd->y_max)))
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 0;
					else
						num_son = 1;
				}
				else
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 2;
					else
						num_son = 3;
				}
			}
			else
			{
				if (py < (0.5*(nd->y_min + nd->y_max)))
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 4;
					else
						num_son = 5;
				}
				else
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 6;
					else
						num_son = 7;
				}
			}
			nd = nd->son[num_son];
		}
	}
}

__device__ goctree_node* goctree_node::lineIn(double& px, double& py, double& pz, double& pu, double& pv, double& pw)
{
	bool in = projectIn(px, py, pz, pu, pv, pw);
	if (in) return get_node_fast(px, py, pz);
	else return NULL;
}

__device__ double goctree_node::get_wall_distance(double &x, double &y, double &z, double &vx, double &vy, double &vz, int &num_neighbor)
{
	double dist, dist_tmp;

	if (vx > 0.0)
	{
		dist = ((x_max - x) / vx);
		num_neighbor = 0;
	}
	else
	{
		dist = ((x_min - x) / vx);
		num_neighbor = 1;
	}

	if (vy > 0.0)
	{
		dist_tmp = ((y_max - y) / vy);
		if (dist_tmp < dist)
		{
			dist = dist_tmp;
			num_neighbor = 2;
		}
	}
	else
	{
		dist_tmp = ((y_min - y) / vy);
		if (dist_tmp < dist)
		{
			dist = dist_tmp;
			num_neighbor = 3;
		}
	}

	if (vz > 0.0)
	{
		dist_tmp = ((z_max - z) / vz);
		if (dist_tmp < dist)
		{
			dist = dist_tmp;
			num_neighbor = 4;
		}
	}
	else
	{
		dist_tmp = ((z_min - z) / vz);
		if (dist_tmp < dist)
		{
			dist = dist_tmp;
			num_neighbor = 5;
		}
	}

	if (dist < 0) //for debug
	{
		printf("dis < 0 !!!\n");
		//asm("trap;");
	}
	return dist;
}

__device__ double goctree_node::get_triangle_intersection(double &x, double &y, double &z, double& u, double& v, double& w, int& ind)
{
	double dist, min_dist = INFINIT_LENGTH;
	ind = -1;
	//get the interaction distance array and find the nearest one
	for (int i = 0; i < num_triangles; i++)
	{
		if (triangle_list[i].intersect(x, y, z, u, v, w, dist))
		{
			if (dist < min_dist)
			{
				min_dist = dist;
				ind = i;
			}
		}
	}
	return min_dist;
}

__device__ goctree_node* goctree_node::step_octree(double ds, double& x, double& y, double& z, double u, double v, double w, MatInfo matIDs[MAX_OBJ_OVERLAP], int& cMatID)
{
	goctree_node* cNode = this;
	//calculate the ds_min_triangle and ds_wall
	int num_neighbor, ind_triangle;
	double ds_min_triangle;
	double ds_wall;

	while (true)
	{
		ds_min_triangle = cNode->get_triangle_intersection(x, y, z, u, v, w, ind_triangle);
		ds_wall = cNode->get_wall_distance(x, y, z, u, v, w, num_neighbor);

		if (ds < ds_min_triangle && ds < ds_wall)
		{
			// An interaction occurs inside this node:
			x += u*ds;
			y += v*ds;
			z += w*ds;
			return cNode; // call knock function next
		}
		else // A wall or a triangle is found before the next interaction.
		{
			if (ds_wall < ds_min_triangle) // wall crossed: move to the neighbor node and continue the jump:
			{
				if (cNode->neighbor[num_neighbor] != NULL)
				{
					// Update the particle position on the new node wall:
					x += ds_wall*u;
					y += ds_wall*v;
					z += ds_wall*w;

					// Discount distance jumped to this wall and continue jump in neighbor node:
					ds -= ds_wall;

					// Update the current node pointer (I can use the fast version because I already
					// know that the particle will be in the neighbor node or inside one of his sons).
					cNode = cNode->neighbor[num_neighbor]->get_node_fast(x, y, z); //it's actually questionable...
				}
				else  // Particle escaped the octree:
				{
					ds_wall += ds_wall*1e-6; //slightly increase ds_wall so that the particle will get off the octree
					x += ds_wall*u;
					y += ds_wall*v;
					z += ds_wall*w;
					return NULL;
				}
			}
			else // Triangle intersection found inside this node
			{
				ds_min_triangle += ds_min_triangle*1e-6; //artificially make the distance a little bit further
				x += ds_min_triangle*u;
				y += ds_min_triangle*v;
				z += ds_min_triangle*w;
				ds -= ds_min_triangle;

				//judge going in or out, and get the current material index
				triangle_octree& tri = cNode->triangle_list[ind_triangle];
				double cp = u*tri.nv[0] + v*tri.nv[1] + w*tri.nv[2];
				int i = 0;
				if (cp < 0) //entering a new object
				{
					for (; i < MAX_OBJ_OVERLAP; ++i)
					{
						if (matIDs[i].id < 0)
						{
							matIDs[i] = tri.mat;
							break;
						}
					}
					if (i == MAX_OBJ_OVERLAP)
					{
						printf("\nError: Too many objects overlap!\n");
						//asm("trap;");
					}
				}
				else //leaving an object
				{
					for (; i < MAX_OBJ_OVERLAP; ++i)
					{
						if (matIDs[i].id == tri.mat.id) //should find a match
						{
							matIDs[i].id = -1;
							break;
						}
					}
					if (i == MAX_OBJ_OVERLAP)
					{
						printf("\nError: Leaving an object that hasn't been entered!\n");
						//asm("trap;");
					}
				}

				//get the current material id
				int ic = 0;
				for (int i = 1; i < MAX_OBJ_OVERLAP; ++i)
				{
					if (matIDs[i].id >= 0 && matIDs[i].priority > matIDs[ic].priority) ic = i;
				}
				int newMatID = matIDs[ic].id;

				if (newMatID == cMatID || newMatID < 0)
				{
					cMatID = newMatID;
					if (newMatID < 0) ds = INFINIT_LENGTH; //make a long way to arrive at a new material or get off the octree
				}
				else //newMatID != cMatID && newMatID >= 0
				{
					cMatID = newMatID;
					return cNode; //arriving at a new material and needs to re-sample ds
				}
			}
		}
	}
}


__device__ goctree_node* goctree_node::get_node(double &px, double &py, double &pz, short level_limit)
{
	goctree_node* nd = this;
	while (true)
	{
		if (px < nd->x_min || px >= nd->x_max ||
			py < nd->y_min || py >= nd->y_max ||
			pz < nd->z_min || pz >= nd->z_max)
		{
			return NULL;    //  Point outside the node
		}
		else if (nd->son[0] == NULL || nd->level == level_limit)
		{
			return nd;     //  Return a pointer to the current leaf
		}
		else
		{
			int num_son = 0;
			if (pz < (0.5*(nd->z_min + nd->z_max)))
			{
				if (py < (0.5*(nd->y_min + nd->y_max)))
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 0;
					else
						num_son = 1;
				}
				else
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 2;
					else
						num_son = 3;
				}
			}
			else
			{
				if (py < (0.5*(nd->y_min + nd->y_max)))
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 4;
					else
						num_son = 5;
				}
				else
				{
					if (px < (0.5*(nd->x_min + nd->x_max)))
						num_son = 6;
					else
						num_son = 7;
				}
			}
			nd = nd->son[num_son];
		}
	}
}

// non-recursive version
__device__ void goctree_node::set_neighbors(goctree_node* root_node)
{
	const double EPSILON = 1.0e-6;
	double px, py, pz;
	int ind[MAX_LEVEL]; // variable stack
	goctree_node* ps[MAX_LEVEL];
	int cur = 0; // current stack index
	int i; // variable can be altered widely
	goctree_node* cn = root_node; 
	
ENTRY:
	if (cn->son[0] != NULL)
	{
		i = 0;
		for (; i < 8; )
		{
			//push the i and cnode on stack
			ind[cur] = i;
			ps[cur] = cn;
			++cur;
			cn = cn->son[i];
			goto ENTRY;

REENTER:    //wait for the function to goto this point
			++i;
		}
	}
	else
	{
		px = cn->x_max + EPSILON;        // 1- Set the neighbor node behind the plane x=x_max
		py = 0.5*(cn->y_min + cn->y_max);
		pz = 0.5*(cn->z_min + cn->z_max);
		cn->neighbor[0] = root_node->get_node(px, py, pz, cn->level);
			
		px = cn->x_min - EPSILON;        // 2- Set the neighbor node before the plane x=x_min
		cn->neighbor[1] = root_node->get_node(px, py, pz, cn->level);

		px = 0.5*(cn->x_min + cn->x_max);    // 3- Set the neighbor node behind the plane y=y_max
		py = cn->y_max + EPSILON;
		cn->neighbor[2] = root_node->get_node(px, py, pz, cn->level);

		py = cn->y_min - EPSILON;
		cn->neighbor[3] = root_node->get_node(px, py, pz, cn->level);

		py = 0.5*(cn->y_min + cn->y_max);
		pz = cn->z_max + EPSILON;
		cn->neighbor[4] = root_node->get_node(px, py, pz, cn->level);

		pz = cn->z_min - EPSILON;
		cn->neighbor[5] = root_node->get_node(px, py, pz, cn->level);
	}

	//got the neighbor pointers for cn
	--cur;
	if (cur > 0)
	{
		i = ind[cur];
		cn = ps[cur];
		goto REENTER;
	} //else: end the call
}

// non-recursive version
__device__ void goctree_node::deleteOctree(goctree_node* root_node) //delete the octree recursively
{
	int ind[MAX_LEVEL]; // variable stack
	goctree_node* ps[MAX_LEVEL];
	int cur = 0; // current stack index
	int i; // variable can be altered widely
	goctree_node* cn = root_node;

ENTRY:
	if (cn->son[0] != NULL)
	{
		i = 0;
		for (; i < 8; ++i)
		{
			//push the stack, and save the current variable
			ind[cur] = i;
			ps[cur] = cn;
			++cur;
			cn = cn->son[i];
			goto ENTRY;

REENTER:    //wait for the function to goto this point
			++i;
		}
	}

	if (cn->num_triangles>0 && cn->triangle_list != NULL) free(cn->triangle_list);
	free(cn);

	//pop the stack, and restore the current variable
	--cur;
	if (cur > 0)
	{
		i = ind[cur];
		cn = ps[cur];
		goto REENTER;
	} //else: end the call
}

/*
// recursive version
__device__ void goctree_node::set_neighbors(goctree_node* root_node)
{
	const double EPSILON = 1.0e-6;
	double px, py, pz;

	if (son[0] != NULL)
	{
		for (int i = 0; i < 8; i++)
		{
			son[i]->set_neighbors(root_node);
		}
	}
	else
	{  
		// It is a leaf: set the neighbors scanning a point near the node wall
		// from the lowest level input node (usually the octree root):

		px = x_max + EPSILON;        // 1- Set the neighbor node behind the plane x=x_max
		py = 0.5*(y_min + y_max);
		pz = 0.5*(z_min + z_max);
		neighbor[0] = root_node->get_node(px, py, pz, level);
		// (The nodes on the octree borders will have NULL neighbors)

		px = x_min - EPSILON;        // 2- Set the neighbor node before the plane x=x_min
		neighbor[1] = root_node->get_node(px, py, pz, level);

		px = 0.5*(x_min + x_max);    // 3- Set the neighbor node behind the plane y=y_max
		py = y_max + EPSILON;
		neighbor[2] = root_node->get_node(px, py, pz, level);

		py = y_min - EPSILON;
		neighbor[3] = root_node->get_node(px, py, pz, level);

		py = 0.5*(y_min + y_max);
		pz = z_max + EPSILON;
		neighbor[4] = root_node->get_node(px, py, pz, level);

		pz = z_min - EPSILON;
		neighbor[5] = root_node->get_node(px, py, pz, level);
	}
}

// recursive version
__device__ void goctree_node::deleteOctree(goctree_node* root_node) //delete the octree recursively
{
	if (root_node->son[0] != NULL)
	{
		for (int i = 0; i < 8; ++i) deleteOctree(root_node->son[i]);
	}
	if (root_node->num_triangles > 0 && root_node->triangle_list != NULL) free(root_node->triangle_list);
	free(root_node);
}
*/