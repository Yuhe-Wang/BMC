//VC++ cannot recognize the following cuda systex, so we need to redefine it in cpp files
#define __device__ 
#define __host__

#include "octree_definition.h"   
#include "triangle_octree.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iomanip> // (Required to use "precision (11)")


// -- Define the static class members: they will be initialized to 0 and NULL by default.
int octree_node::amount_nodes[MAX_LEVEL+1];
int octree_node::amount_leaves[MAX_LEVEL+1];
int octree_node::triangle_mesh_size=0;
triangle_octree* octree_node::triangle_mesh = NULL;

// -- Definition of the class constructor:

// Default class constructor:
// The node lower and upper corners and its level are initialized to 0. The
// function 'set_node_parameters' has to be called to correctly initialize
// the node data.
octree_node::octree_node()
{
    set_node_parameters(0.0,0.0,0.0, 0.0,0.0,0.0, 0);    // Set the node corners and level

    // Initialize all the other node pointers to NULL:
    triangle_list = NULL;
    for (register int i=0 ; i<6 ; i++)
      neighbor[i] = NULL;
    for (register int j=0 ; j<8 ; j++)
      son[j] = NULL;
}


// Explicit class constructor:
octree_node::octree_node(double x0,double y0,double z0,
                         double x1,double y1,double z1,
                         char level0)
{
    set_node_parameters(x0,y0,z0, x1,y1,z1, level0);     // Set the node corners and level

    // Initialize all the other node pointers to NULL:
    triangle_list = NULL;
    for (register int i=0 ; i<6 ; i++)
      neighbor[i] = NULL;
    for (register int j=0 ; j<8 ; j++)
      son[j] = NULL;    
}

// -- Definition of the class destructor.
//    The memory of all the node sons and grand sons is also freed
octree_node::~octree_node()
{
    for (int i=0 ; i<8 ; i++)
       if (son[i]!=NULL)                 // Delete the son object, if it exists
          delete son[i];
    if (triangle_list!=NULL)
    {
       delete [] triangle_list;  // Delete the triangle list, if it exists.
    }
}


// -- Class function definition:

// Sets the value of the node lower and upper corners and the node level.
void octree_node::set_node_parameters(double x0,double y0,double z0,
                                      double x1,double y1,double z1,
                                      char level0)
{
    x_min = x0;        // - Lower corner
    y_min = y0;
    z_min = z0;
    x_max = x1;        // - Upper corner
    y_max = y1;
    z_max = z1;
    level = level0;    // - Node level
}

/* Generate the octree structure and distribute the triangles.
 *  This subroutine finds the triangles that intersect the node (using the
 *  private function sort_triangles), assigns the triangles to the member
 *  variable triangle_list, checks the termination_condition, and --if the node
 *  has to be subdivided-- creates the instances of the 8 sub-nodes and repeat
 *  the previous operations for each sub-node.
 *
 *  Typically this function is only called once for the "root" octree_node instance
 *  (global variable), and then the complete octree structure is generated below
 *  this root node.
 *
 *  The triangles must be stored in the global variable triangle_octree *triangle_mesh
 *  before calling this function.
 */
void octree_node::create_subnodes_with_triangles(int &input_max_level, int &maxTriangles, octree_node *root)
{
   // Assign the triangles to this node, using the triangle_mesh static class variable or the node father's triangle_list:
    this->sort_triangles(root);

   // Check if the present node is a leaf and doesn't have to be subdivided:
    if (termination_condition(input_max_level, maxTriangles)) 
    {   
       // --This is a leaf!
        (amount_leaves[int(level)])++;
    }
    else
    {   
       // -- Regular node: recursively sub-divide it:        
        (amount_nodes[int(level)])++ ;

        double x_mid = 0.5*(x_max+x_min);    // Node middle point
        double y_mid = 0.5*(y_max+y_min);
        double z_mid = 0.5*(z_max+z_min);
        char level_son = level + (char)(1);  // Sub-node level

		try
		{
			//  The sub-nodes are ordered as: X++, Y+, Z (the X coordinate is changed first and the Z last)
			//    -Sub-node 0: (X,Y,Z)=(0,0,0)
			son[0] = new octree_node(x_min, y_min, z_min, x_mid, y_mid, z_mid, level_son);
			//    -Sub-node 1: (1,0,0)
			son[1] = new octree_node(x_mid, y_min, z_min, x_max, y_mid, z_mid, level_son);
			//    -Sub-node 2: (0,1,0)
			son[2] = new octree_node(x_min, y_mid, z_min, x_mid, y_max, z_mid, level_son);
			//    -Sub-node 3: (1,1,0)
			son[3] = new octree_node(x_mid, y_mid, z_min, x_max, y_max, z_mid, level_son);
			//    -Sub-node 4: (0,0,1)
			son[4] = new octree_node(x_min, y_min, z_mid, x_mid, y_mid, z_max, level_son);
			//    -Sub-node 5: (1,0,1)
			son[5] = new octree_node(x_mid, y_min, z_mid, x_max, y_mid, z_max, level_son);
			//    -Sub-node 6: (0,1,1)
			son[6] = new octree_node(x_min, y_mid, z_mid, x_mid, y_max, z_max, level_son);
			//    -Sub-node 7: (1,1,1)
			son[7] = new octree_node(x_mid, y_mid, z_mid, x_max, y_max, z_max, level_son);
		}
		catch (std::bad_alloc&)
		{
			cout << "\n ERROR in {create_subnodes_with_triangles}!! The memory for sub nodes could not be allocated!\n" << endl;
			exit(1);
		}

        for (int i=0 ; i<8 ; i++)
        {
            // Create the 8 sons, and the grandsons, and the grand-grandsons... :
			son[i]->create_subnodes_with_triangles(input_max_level, maxTriangles, root);
        }
    }
}

/* Delete the unnecessary instances of 'triangle_list' in the nodes that are not leaves. */
// I could also simplify the octree structure merging neighbor nodes with same size and few triangles.  !!! TO DO: !!DeBuG!!
void octree_node::clean_octree()
{
   if (this->son[0] != NULL)
   {
     // The node is not a leaf! Delete the triangle list and step into the son sub-nodes:
      delete [] this->triangle_list;
      this->triangle_list = NULL;
      
      for (int i=0 ; i<8 ; i++)
        this->son[i]->clean_octree();
   }
}


/*  Draw the octree structure at the input z plane in postscript format (one square for
*  each leaf in the plane). The output fills approximately a vertical letter page and it
*  is automatically scaled to the input node size (usually the root node).
*
*     @param[in] octree   Pointer to the octree root node.
*     @param[in] octree_max_level Maximum octree subdivision level: useful information written in the EPS header.
*     @param[in] z_plane  Z plane where the octree is to be plotted
*     @param[in] file_name  Name of the output postscript file.
*     @param[in] first_call  Optional variable used internally to detect the first
*                           call to the function (the function is called recursively
*                           and only the first call must write the PS header).
*/
void octree_node::draw_octree_ps(octree_node *octree, int &octree_max_level, double &z_plane, const char *file_ps, bool first_call)
{
   ofstream file;

   if (first_call)  // Only the first call will write the header and the end of the postscript file:
   {
      file.open(file_ps);

      double page_side_short  = 21.0,                 // Maximum image size in cm
                                                      // ** US letter: 8.5 x 11 inch == 21.59 x 27.94 cm
                                                      // ** A4 page:    21 x 29.7 cm
             page_side_long   = 27.94,                // (==  11 in)
             page_margin      =  0.5,                 // Page margin (for the 4 sides) in cm.
             cmTOinchTOpoints = 28.3465,              // Convert cm->inches (1/2.54 cm)->PS points (1/72 inch) [72/2.54=28.3465]
             width_x = octree->x_max - octree->x_min, // Octree size
             width_y = octree->y_max - octree->y_min,
             scale;                                   // Scale factor to fit the image in the paper.
    
      if( width_x > width_y)
      {
         scale = (page_side_short-2.0*page_margin)/width_x;  // Scale to the paper long axis
         if ((width_y*scale)>page_side_long) scale = (page_side_long-2.0*page_margin)/width_x; // Scale to the short axis
      }
      else
      {
         scale = (page_side_long-2.0*page_margin)/width_y;  // Scale to the paper short axis
         if ((width_x*scale)>page_side_short) scale = (page_side_short-2.0*page_margin)/width_y; // Scale to the long axis
      }
      file<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
      file<<"%%BoundingBox: 0 0 "<<ceil((width_x*scale+2.0*page_margin)*cmTOinchTOpoints)
                            <<' '<<ceil((width_y*scale+2.0*page_margin)*cmTOinchTOpoints)<<endl;
      file<<"%%Title: "<<file_ps<<endl;
      file<<"%%Creator: penMesh (Andreu Badal)"<<endl;
      file<<"%%Description: Drawing of the octree plane "<<z_plane<<" using a square box for each leaf."<<endl;
      file<<"%%             Maximum octree subdivision level is "<<octree_max_level<<endl;
      file<<"%%             For information on postscript drawing see: "<<endl;
      file<<"%%               www.cs.indiana.edu/docproject/programming/postscript/postscript.html"<<endl;
      file<<"%%Page: 1 1 \n"<<endl;

      file<<"/cm {"<<cmTOinchTOpoints<<" mul} def   % Convert cm->inches (1/2.54 cm)-> PS points (1/72 inch)"<<endl;
      file<<page_margin<<" cm "<<page_margin<<" cm translate   % Translate the page origin to have a nice margins"<<endl;
      file<<scale<<' '<<scale<<" scale   % Scale the phantom coordinates to fit the page"<<endl;
      file<<(-octree->x_min)<<" cm "<<(-octree->y_min)
                            <<" cm translate   % Translate the phantom origin to the page origin\n"<<endl;

      file<<"% --Aliased used to reduce the file size:"<<endl;
      file<<"/N {newpath} bind def"<<endl;
      file<<"/M {moveto} bind def"<<endl;
      file<<"/L {lineto} bind def"<<endl;
      file<<"/LS {lineto stroke} bind def \n"<<endl;
      // file<<"/LCS {lineto closepath stroke} bind def \n"<<endl;

      file<<"% --Octree bounding box: "<<endl;
      file<<"newpath"<<' '
        <<octree->x_min<<" cm "<<octree->y_min<<" cm moveto"<<' '
        <<octree->x_max<<" cm "<<octree->y_min<<" cm lineto"<<' '
        <<octree->x_max<<" cm "<<octree->y_max<<" cm lineto"<<' '
        <<octree->x_min<<" cm "<<octree->y_max<<" cm lineto"<<' '
        <<"closepath"<<' '<<"stroke"<<endl;
      file<<"% --Octree leaves (using alias): "<<endl;
      
      file.close();
   }  // End code for first_call
   
   if (octree->son[0] != NULL)
   {
      for (int i=0; i<8; i++)
        draw_octree_ps(octree->son[i], octree_max_level, z_plane, file_ps, false);   // Calling the function recursively!
   }
   else
   {
      // This is a leaf: If zplane is in the node, draw the corresponding square in PostScript format:
      if (octree->z_min <= z_plane && octree->z_max > z_plane)
      {
        file.open(file_ps, ofstream::app);
        
        // - Draw the image with 3 decimal positions (requires: #include <iomanip>):
        // file << fixed << setprecision(3);
        
        // - Draw the image using alias to reduce the size of the output file:
        file <<"N "<<octree->x_min<<" cm "<<octree->y_min<<" cm M "
                   <<octree->x_max<<" cm "<<octree->y_min<<" cm L "
                   <<octree->x_max<<" cm "<<octree->y_max<<" cm LS"<<endl;
          // - Old: drawing the 4 box sides and without aliases:
          //   file<<"newpath"<<' '<<octree->x_min<<" cm "<<octree->y_min<<" cm moveto"<<' '<<octree->x_max<<" cm "<<octree->y_min<<" cm lineto"<<' '<<octree->x_min<<" cm "<<octree->y_max<<" cm lineto"<<' '<<"closepath"<<' '<<"stroke"<<endl;

        file.close();
      }
   }

   if (first_call)
   {
      // Write the postscript trailer at the end of the EPS, for the first call:
      file.open(file_ps, ofstream::app);
      file<<'\n';
      file<<"showpage"<<endl;
      file<<'\n';
      file<<"%%Trailer"<<endl;
      file<<"%%EOF"<<endl;
      file<<'\n';
      file.close();
   }
}


/* Calculate the mean quantity (and standard deviation) of triangles inside the leaves with
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
void octree_node::triangles_per_leaf(octree_node *octree, int &octree_max_level, double &mean_tri,
                       double &sigma_tri, int &leaves_full, int &leaves_empty, bool first_call)
{
   if (first_call)
   {
      mean_tri  = 0.0;
      sigma_tri = 0.0;
      leaves_full  = 0;
      leaves_empty = 0;
   }


   if (octree->son[0] != NULL)
   {
      for (int i=0; i<8; i++)
      // -Call the function recursively through all the octree structure:
      triangles_per_leaf(octree->son[i], octree_max_level, mean_tri, sigma_tri, leaves_full, leaves_empty, false);
   }
   else
   {
    // This is a leaf:
     if (octree->num_triangles>0)
     {
       leaves_full++;  // -Leaf with triangles inside
       if (octree->level==octree_max_level)
       {
         // - This leaf has the maximum recursion level: score the amount of triangles:
          mean_tri += octree->num_triangles;
          sigma_tri += pow(double(octree->num_triangles),2);
       }
    }
    else
      leaves_empty++;
   }

   if (first_call)
   {
      if(leaves_full!=0)
      {
         mean_tri  = mean_tri/double(leaves_full);
         sigma_tri = sqrt(sigma_tri/double(leaves_full) - pow(mean_tri,2));
      }
      else
      {
         mean_tri  = 0.0;
         sigma_tri = 0.0;
      }
   }
}

/* Set the neighbor[6] class variable: pointers to the 6 neighbor nodes.
*  The neighbor nodes are ordered as: X++, Y+, Z ("front,back,right,left,top,bottom")
*
*     @param[in]  root_node Pointer to the root node (i.e., the triangle geometry bounding box).
*/
void octree_node::set_neighbors(octree_node* root_node)
{
	const double EPSILON = 1.0e-6;
	double px, py, pz;

	if (son[0] != NULL)
	{ // ** The node is not a leaf. Go recursively into the sons to set the leaves neighbors:
		for (int i = 0; i < 8; i++)
		{
			son[i]->set_neighbors(root_node);
		}
	}
	else
	{  // ** It is a leaf: set the neighbors scanning a point near the node wall
		//    from the lowest level input node (usually the octree root):

		px = x_max + EPSILON;        // 1- Set the neighbor node behind the plane x=x_max
		py = 0.5*(y_min + y_max);
		pz = 0.5*(z_min + z_max);
		neighbor[0] = root_node->get_node(px, py, pz, level);
		// (The nodes on the octree borders will have NULL neighbors)

		px = x_min - EPSILON;        // 2- Set the neighbor node before the plane x=x_min
		//py = 0.5*(y_min+y_max);
		//pz = 0.5*(z_min+z_max);
		neighbor[1] = root_node->get_node(px, py, pz, level);

		px = 0.5*(x_min + x_max);    // 3- Set the neighbor node behind the plane y=y_max
		py = y_max + EPSILON;
		//pz = 0.5*(z_min+z_max);
		neighbor[2] = root_node->get_node(px, py, pz, level);

		//px = 0.5*(x_min+x_max);  // 4- Set the neighbor node before the plane y=y_min
		py = y_min - EPSILON;
		//pz = 0.5*(z_min+z_max);
		neighbor[3] = root_node->get_node(px, py, pz, level);

		//px = 0.5*(x_min+x_max);  // 5- Set the neighbor node behind the plane z=z_max
		py = 0.5*(y_min + y_max);
		pz = z_max + EPSILON;
		neighbor[4] = root_node->get_node(px, py, pz, level);

		//px = 0.5*(x_min+x_max);  // 6- Set the neighbor node before the plane z=z_min
		//py = 0.5*(y_min+y_max);
		pz = z_min - EPSILON;
		neighbor[5] = root_node->get_node(px, py, pz, level);
	}
}

/* Set the 'triangle_list' class member variable assigning pointers to the triangles
*  inside each node. This function uses the static class variables 'triangle_mesh' and
*  'triangle_mesh_size', and also the input node --the root node-- to get the triangles
*  in the father node. The father's 'triangle_list' variable can be deleted after its
*  triangles have been distributed between the 8 sons calling 'clean_octree'.
*  The number of triangles assigned to the node is stored in the variable 'num_triangles'.
*
*     @param[in]  root_node  Pointer to the root node (i.e., the triangle geometry bounding box).
*/
void octree_node::sort_triangles(octree_node *root_node)
{
	if (triangle_list != NULL)
	{
		cout << "\n ERROR in sort_triangles!! Trying to set the triangles for a node that has already been initiated (triangle_list!=NULL) \n" << endl;
		exit(1);
	}

	int i, triangle_counter = 0;

	triangle_octree *triangle_ptr, **tmp_list;

	if (level == 0)  // --> processing root node
	{
		if (triangle_mesh_size == 0)
		{
			cout << "\n\n\n WARNING!! {sort_triangles} No triangles found inside the input bounding box: triangle geometry not used!!" << endl;
			// exit(1);
		}
		// -Root node; check all the triangles from triangle_mesh:
		try
		{
			tmp_list = new triangle_octree*[triangle_mesh_size];  // Generate the temporal list of pointers.
		}
		catch (const bad_alloc&)
		{
			cout << "\n ERROR in sort_triangles!! The temporal memory for the temporal triangle_list could not be allocated for the root node!\n" << endl;
			exit(1);
		}
		triangle_ptr = triangle_mesh;
		for (i = 0; i < triangle_mesh_size; i++)
		{
			if (triangle_ptr->inside_box(x_min, y_min, z_min, x_max, y_max, z_max))
			{ // Return true if the triangle is inside the node:
				tmp_list[triangle_counter] = triangle_ptr;
				triangle_counter++;
			}
			triangle_ptr++;   // Point to the next triangle in the mesh.
		}
	}
	else
	{

		// Find the father node searching the root node for the node at the current position and "level-1":
		double x_center = 0.5*(x_min + x_max), y_center = 0.5*(y_min + y_max), z_center = 0.5*(z_min + z_max);
		octree_node *node_father = root_node->get_node(x_center, y_center, z_center, (level - 1));

		if (node_father->num_triangles == 0)
		{
			// - This sub-node can not have any triangle inside!
			num_triangles = 0;
			return;
		}
		// -New sub-node; check only the father triangles:
		try
		{
			tmp_list = new triangle_octree*[node_father->num_triangles];
		}
		catch (const bad_alloc&)
		{
			cout << "\n ERROR in sort_triangles!! The temporal memory for the temporal triangle_list could not be allocated!\n" << endl;
			exit(1);
		}
		for (i = 0; i < node_father->num_triangles; i++)
		{
			triangle_ptr = node_father->triangle_list[i];  // Point to the next father triangle.

			if (triangle_ptr->inside_box(x_min, y_min, z_min, x_max, y_max, z_max))
			{ // Return true if the triangle is inside the node:
				tmp_list[triangle_counter] = triangle_ptr;
				triangle_counter++;
			}
		}
	}

	// -Set the node member variable triangle_list to point to every inner triangle:
	try
	{
		triangle_list = new triangle_octree*[triangle_counter];
	}
	catch (const bad_alloc&)
	{
		cout << "\n ERROR in sort_triangles!! The temporal memory for the triangle_list could not be allocated!\n" << endl;
		exit(1);
	}


	for (register int ii = 0; ii < triangle_counter; ii++)
	{
		triangle_list[ii] = tmp_list[ii];
	}

	// -Set the node member variable num_triangles:
	num_triangles = triangle_counter;

	// -Free the space of the temporal list of pointers:
	delete[] tmp_list;
}

/* Function to search from the input node (the octree root usually) for the sub-node that
*  contains the input point and that has the maximum level specified (this is needed in
*  function 'set_neighbors' because the stored neighbors can not have a higher level, ie
*  smaller size).
*
*      @param[in]  px   X coordinate of the point.
*      @param[in]  py   Y coordinate of the point.
*      @param[in]  pz   Z coordinate of the point.
*      @param[in]  level_limit  Maximum level of the returned node (the sons of
*                     the node at this level are not checked). The default value is -1,
*                     which means that the limit will have no effect.
*      @return  The return value is a pointer to the sub-node that contains the input
*               point, or NULL if the point is not inside the current node.
*/
octree_node* octree_node::get_node(double &px,double &py,double &pz, char level_limit)
{
	octree_node* nd = this;
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

/* Function to search from the input node for the sub-node that contains the input point.
*  This version of 'get_node' does not set a limit for the maximum level of the returned
*  sub-node and it does not check if the point is actually inside the node. This will speed up
*  a little the part of the simulation where particles pass from one octree leave its neighbor.
*
*      @param[in]  px   X coordinate of the point.
*      @param[in]  py   Y coordinate of the point.
*      @param[in]  pz   Z coordinate of the point.
*      @return  The return value is a pointer to the sub-node that contains the input point.
*               If the point is not inside the current node the returned pointer will be wrong.
*/
octree_node* octree_node::get_node_fast(double &px,double &py,double &pz)
{
	octree_node* nd = this;
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


//advance the particle in the octree space. If it encounters a new object, it will stop on this surface and give newMatID
//It returns NULL if the particle exit the octree, otherwise it returns the node where the particle locates.
octree_node* octree_node::step_octree(double ds, double& x, double& y, double& z, double u, double v, double w, OctreeData& od, int& cMatID, bool& dirChanged)
{
	octree_node* cNode = this;

	//get the interaction distance array and find the nearest one
	
	double& ds_min_triangle = od.ds_min_triangle;
	double& ds_wall = od.ds_wall;
	int& ind_min_triangle = od.ind_min_triangle;
	MatInfo* matIDs = od.matIDs;

	if(dirChanged) od.update(cNode, x, y, z, u, v, w);

	while (true)
	{
		if (ds < ds_min_triangle && ds < ds_wall)
		{
			// ** An interaction occurs inside this node:
			x += u*ds;
			y += v*ds;
			z += w*ds;
			dirChanged = true;
			return cNode; // call knock function next
		}
		else // ** A wall or a triangle is found before the next interaction.
		{
			if (ds_wall < ds_min_triangle) // ** Node wall crossed: move to the neighbor node and continue the jump:
			{
				if (cNode->neighbor[od.num_neighbor] != NULL)
				{
					// ** Update the particle position on the new node wall:
					x += ds_wall*u;
					y += ds_wall*v;
					z += ds_wall*w;

					// ** Discount distance jumped to this wall and continue jump in neighbor node:
					ds -= ds_wall;

					// Update the current node pointer (I can use the fast version because I already
					// know that the particle will be in the neighbor node or inside one of his sons).
					cNode = cNode->neighbor[od.num_neighbor]->get_node_fast(x, y, z); //it's actually questionable...

					od.update(cNode, x, y, z, u, v, w);
				}
				else  // ** Particle escaped the octree:
				{
					ds_wall += ds_wall*1e-6; //slightly increase ds_wall so that the particle will get off the octree
					x += ds_wall*u;
					y += ds_wall*v;
					z += ds_wall*w;
					return NULL;
				}
			}
			else // ** Triangle intersection found before the end of the node:
			{
				x += ds_min_triangle*u;
				y += ds_min_triangle*v;
				z += ds_min_triangle*w;
				ds -= ds_min_triangle;

				//judge going in or out.
				triangle_octree* ptri = cNode->triangle_list[od.ds_triangles[od.ind_min_triangle].ind];
				double cp = u*ptri->nv[0] + v*ptri->nv[1] + w*ptri->nv[2];
				if (cp < 0) //entering a new object
				{
					int i = 0;
					for (; i < MAX_OBJ_OVERLAP; ++i)
					{
						if (matIDs[i].id < 0)
						{
							matIDs[i] = ptri->mat;
							break;
						}
					}
					if (i == MAX_OBJ_OVERLAP || i == 1)
					{
						cout << "\nToo many objects overlap!" << endl;
						getchar();
						exit(1);
					}
				}
				else //leaving an object
				{
					for (int i = 0; i < MAX_OBJ_OVERLAP; ++i)
					{
						if (matIDs[i].id == ptri->mat.id)
						{
							matIDs[i].id = -1;
							break;
						}
					}
				}

				//get the current material id
				int ic = 0;
				for (int i = 1; i < MAX_OBJ_OVERLAP; ++i)
				{
					if (matIDs[i].id >= 0 && matIDs[i].priority > matIDs[ic].priority) ic = i;
				}

				int newMatID = matIDs[ic].id;
				od.update();
				if (newMatID == cMatID || newMatID < 0)
				{
					cMatID = newMatID;
					if (newMatID < 0) ds = INFINIT_LENGTH; //make a long way to arrive at a new material or get off the octree
				}
				else //newMatID != cMatID && newMatID >= 0, which means arriving at a new material and needs to resample ds
				{
					cMatID = newMatID;
					dirChanged = false;
					return cNode;
				}
			}
		}
	}
	
	return NULL;
}

int octree_node::get_triangle_intersection(double &x, double &y, double &z, double& u, double& v, double& w, 
	CrossedTriangle* triangle_array, int& triangle_num)
{
	triangle_num = 0;
	double dist = 0, min_dist = INFINIT_LENGTH;
	int ind_min_dist = -1;
	
	//get the interaction distance array and find the nearest one
	for (int i = 0; i < num_triangles; i++)
	{
		if (triangle_list[i]->intersect(x, y, z, u, v, w, dist))
		{
			if (dist < min_dist)
			{
				min_dist = dist;
				ind_min_dist = triangle_num;
			}
			else if (dist == min_dist) dist += dist*1e-6; //slightly increase the distance to ensure that they're not identical
			triangle_array[triangle_num].dist = dist;
			triangle_array[triangle_num].ind = i;
			++triangle_num;
		}
	}

	return ind_min_dist;
}

octree_node* octree_node::lineIn(double& px, double& py, double& pz, double& pu, double& pv, double& pw)
{
	bool in = projectIn(px, py, pz, pu, pv, pw);
	if (in) return get_node_fast(px, py, pz);
	else return NULL;
}

void OctreeData::update(octree_node* node, double& x, double& y, double& z, double& u, double& v, double& w)
{
	ds_min_triangle = 1e9;
	ind_min_triangle = node->get_triangle_intersection(x, y, z, u, v, w, ds_triangles, ntri);
	if (ind_min_triangle >= 0)  ds_min_triangle = ds_triangles[ind_min_triangle].dist;
	ds_wall = node->get_wall_distance(x, y, z, u, v, w, num_neighbor);
}
void OctreeData::update()
{
	//recalculate ds_wall
	ds_wall -= ds_min_triangle;
	//find the new ds_min_triangle
	double ds_min = INFINIT_LENGTH;
	ind_min_triangle = 0;
	for (int i = 0; i < ntri; ++i)
	{
		ds_triangles[i].dist -= ds_min_triangle;
		if (ds_triangles[i].dist == 0) ds_triangles[i].dist = INFINIT_LENGTH;
		if (ds_triangles[i].dist < ds_min)
		{
			ds_min = ds_triangles[i].dist;
			ind_min_triangle = i;
		}
		else if (ds_triangles[i].dist == ds_min) ds_triangles[i].dist += 1e-6*ds_triangles[i].dist; //to make them unique
	}
	ds_min_triangle = ds_min;
}

octree_node* creatOctree(const string& fname, int input_max_level, int maxTriangles, vector<int> matList)
{
	OBJ3D obj;
	obj.load(fname);
	octree_node* root = new octree_node;
	root->set_node_parameters(obj.minv.x, obj.minv.y, obj.minv.z, obj.maxv.x, obj.maxv.y, obj.maxv.z, char(0));
	//copy triangles to the static triangle array of octree_node
	octree_node::triangle_mesh_size = (int)obj.fs.size();
	octree_node::triangle_mesh = new triangle_octree[octree_node::triangle_mesh_size];
	int numMat = (int)matList.size();
	for (int i = 0; i < octree_node::triangle_mesh_size; ++i)
	{
		int i0 = obj.fs[i].v[0];
		int i1 = obj.fs[i].v[1];
		int i2 = obj.fs[i].v[2];
		int oid = obj.vs[i0].obj_id;
		//change mat_id to array index j. 
		int j = 0;
		if (numMat > 0) //if matList is not provided, use j = 0 as default
		{
			for (; j < numMat; ++j)
			{
				if (matList[j] == obj.info[oid].mat_id) break;
			}
			if (j == numMat)
			{
				cout << "Error: the geometry contains more materials than provided!" << endl;
				getchar();
				exit(1);
			}
		}
		octree_node::triangle_mesh[i].set_triangle_members(obj.vs[i0].x, obj.vs[i0].y, obj.vs[i0].z, obj.vs[i1].x, obj.vs[i1].y, obj.vs[i1].z,
			obj.vs[i2].x, obj.vs[i2].y, obj.vs[i2].z, j, obj.info[oid].priority);
	}

	root->create_subnodes_with_triangles(input_max_level, maxTriangles, root);   // The root node will be used to find the father node of each leaf.
	root->set_neighbors(root);
	root->clean_octree();

	cout << "\nOctree successfully created! Amount of nodes and leaves for each level:" << endl;
	int i, count_nodes = 0, count_leaves = 0;
	for (i = 0; i < (MAX_LEVEL + 1); i++)
	{
		cout << "\t\t\t" << i << ":  \t\t" << root->amount_nodes[i]
			<< "  \t\t" << root->amount_leaves[i] << endl;
		count_nodes += root->amount_nodes[i];
		count_leaves += root->amount_leaves[i];
	}

	cout << "\n        -Input octree subdivision limit: " << input_max_level << endl;
	cout << "        -Input octree max triangles per node: ";
	if (maxTriangles < 0) cout << " = Node Level" << endl;
	else cout << maxTriangles << endl;
	cout << "        -Built-in octree subdivision limit:  MAX_LEVEL = " << MAX_LEVEL << endl;
	cout << "        -Total number of nodes:  " << count_nodes << endl;
	cout << "        -Total number of leaves: " << count_leaves << endl;
	cout << "        -Number of nodes+leaves: " << (count_nodes + count_leaves) << endl;

	// -Calculate the mean number of triangles in the leaf with maximum level:
	double mean_num_triangles, sigma_num_triangles;
	int leaves_full, leaves_empty;
	octree_node::triangles_per_leaf(root, input_max_level, mean_num_triangles, sigma_num_triangles, leaves_full, leaves_empty);
	cout << "        -Number of leaves with triangles inside = " << leaves_full << " ; Empty leaves = " << leaves_empty << endl;
	cout << "        -Triangles per occupied leaf with maximum level = "
		<< mean_num_triangles << " +- " << 2 * sigma_num_triangles << endl;
	cout << "         (average +- 2sigma)" << endl;

	return root;
}

void deleteOctree(octree_node* root)
{
	delete root;
	delete[] octree_node::triangle_mesh;
	octree_node::triangle_mesh = NULL;
	octree_node::triangle_mesh_size = 0;
	for (int i = 0; i < MAX_LEVEL + 1; ++i)
	{
		octree_node::amount_nodes[i] = 0;
		octree_node::amount_leaves[i] = 0;
	}
}