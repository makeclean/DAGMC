// MCNP5/dagmc/SurfaceMeshTally.cpp

#include <iostream>
#include <sstream>
#include <cmath> 
#include <set>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/OrientedBoxTreeTool.hpp"
#include "moab/Skinner.hpp"
#include "moab/CN.hpp"

#include <cassert>
#include "SurfaceMeshTally.hpp"

namespace moab {
//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
SurfaceMeshTally::SurfaceMeshTally(const TallyInput& input)
    : MeshTally(input),
      mb (new moab::Core())
{
   std::cout << "Creating surface mesh tally" << input.tally_id 
            << ", input: " << input_filename 
            << ", output: " << output_filename << std::endl;

   parse_tally_options(); // parse the tally options
 
   set_tally_meshset(); // set the meshset

   // reduce the loaded MOAB mesh set to include only 2D elements
   Range facets;
   ErrorCode rval = reduce_meshset_to_2D(mb, tally_mesh_set, facets);  
   if(rval != MB_SUCCESS)
     {
       std::cout << "Failed to reduce meshset to 2D" << std::endl;
       exit(1);
     }
   assert (rval == MB_SUCCESS);

   // initialize MeshTally::tally_points to include all mesh cells
   set_tally_points(all_tets);

   // precompute the surface normals
   rval = compute_surface_normals(facets);
   // precompute the adjacency data

   assert (rval == MB_SUCCESS);
  
   // Perform tasks
   rval = setup_tags( mb );
   assert (rval == MB_SUCCESS);

   // build the kd-trees
   build_trees(all_tets);

}
//---------------------------------------------------------------------------//
// DESTRUCTOR
//---------------------------------------------------------------------------//
SurfaceMeshTally::~SurfaceMeshTally()
{
  delete mb;
}

// need a write_data()

// clone of steves parse_tally_options
void SurfaceMeshTally::parse_tally_options()
{
  // does nothing currently
  // will be used to set if flux or current tally
  return;
}

//
// Computes the score for the given event, 
// 
void SurfaceMeshTally::compute_score(const TallyEvent &event)
{
  std::vector<EntityHandle> hit_facets;
  unsigned int num_hits;
  // would like to record surface crossing event
  if (event.type != TallyEvent::TRACK) return;

  hit_facets = get_intersections( kdtree, kdtree_root,event.direction,event.position,
				 event.track_length);
  num_hits = hit_facets.size();
  if ( num_hits == 0 )
    return;

  for ( i = 0 ; i < num_hits ; i++ )
    {
      int index = get_surface_index(hit_facets[i]);
      double angle = get_angle(surface_normals[index],event.direction);
      // score is weight*cos(theta)
      data->add_score_to_tally(index,event.weight*angle,energy_bin);
    }

  return;
}

/*
* Wrapper function, return the surface normal of the triangle 
*/
CartVect SurfaceMeshTally::surface_normal(MBEntityHandle triangle)
{
  const MBEntityHandle *connectivity; 
  int number_nodes = 0;
  MBErrorCode result = MBI()->get_connectivity(triangle,connectivity,number_nodes);
 
  //get the coordinates of each node
  MBCartVect nodes[3],normal;
  for ( int i = 0 ; i <= 2 ; i++ )
    {
      result = MBI()->get_coords(&(connectivity[i]),1,nodes[i].array());
    }
  
  MBCartVect V1 = nodes[1]-nodes[0];
  MBCartVect V2 = nodes[2]-nodes[1];
  MBCartVect CrossProd = V1*V2;
  //    return CrossProd.normalize();
  MBCartVect normalize(CrossProd); 
  return CrossProd;
}
  
double SurfaceMeshTally::get_angle( const MBCartVect u, const MBCartVect v )
{
  double angle=(u % v) / std::sqrt((u % u) * (v % v));
  if (angle>1.) angle=1.;
  if (angle<-1.) angle = -1.;
  return angle;
}
  

// compute the surface normals
ErrorCode SurfaceMeshTally::compute_surface_normals(const Range &all_facets)
{
  ErrorCode rval;
  // iterate over the facets
  int num_facets = all_facets.size();
  std::cerr << "  There are << " num_facets << " triangles in this mesh tally." << std::endl;

  if (num_facets != 0 )
    {
      // resize the surface normal array
      surface_normals.resize(num_facets);
    }

  for( Range::const_iterator i=all_facets.begin(); i!=all_facets.end(); ++i)
    {
      CartVect normal = surface_normal(*i); // get the surface normal
      surface_normals.push_back(normal);
    }
  return MB_SUCCESS;
}

// return the intersection data, could return more than one of there are multiple surfaces
// in each meshset
std::vector<MBEntityHandle> SurfaceMeshTally::get_intersections(MBAdaptiveKDTree* kdtree,
				MBEntityHandle kdtree_root, double dir[3], double box_center[3],
				double ray_length)
{
  // we do not keep the intersection distances
  std::vector<MBEntityHandle> ret_triangles;
  std::vector<double> ret_intersections;
  
  // get the intersections along the ray
  MBErrorCode result = kdtree->ray_intersect_triangles( kdtree_root, 1.0e-6, 
							dir, box_center,
							ret_triangles, ret_intersections, 0, 
							ray_length );      
  
  if ( result != MB_SUCCESS)
    gen::error("failed to intersect ray with triangles");
  
  if ( ret_intersections.size() == 0 )
    {
      if(DEBUG)    
	std::cout << "ray misses mesh" << std::endl;
      return;
    }
  
  return ret_triangles;
}
  

}
