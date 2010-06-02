#include "DagMC.hpp"
#include "MBTagConventions.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/Core.hpp"
#include "moab/GeomUtil.hpp"
#include "FileOptions.hpp"

#ifdef CGM
#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "CubitDefines.h"
#include "GeometryQueryTool.hpp"
#include "CubitVector.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#ifndef M_PI  /* windows */
# define M_PI 3.14159265358979323846
#endif

#define MB_OBB_TREE_TAG_NAME "OBB_TREE"
#define CATEGORY_TAG_LENGTH 32

#define MAT_GROUP 0
#define BC_SPEC 1
#define BC_WHITE 2
#define IMP_ZERO 3
#define TALLY_GROUP 4

#define MATERIAL_TAG_NAME "DAGMC_MATERIAL_ID"
#define DENSITY_TAG_NAME  "DAGMC_MATERIAL_DENSITY"
#define BC_TAG_NAME       "DAGMC_BOUNDARY_CONDITION"
#define IMP_TAG_NAME      "DAGMC_IMPORTANCE"
#define TALLY_TAG_NAME    "DAGMC_TALLY"


#define MBI moab_instance()

namespace moab {

DagMC *DagMC::instance_ = NULL;

const bool debug = false;

unsigned int DagMC::interface_revision()
{
  unsigned int result = 0;
  sscanf( DAGMC_INTERFACE_REVISION, "$Rev: %u $\n", &result );
  return result;
}
  

void DagMC::create_instance(Interface *mb_impl) 
{
  if (NULL == mb_impl) mb_impl = new Core();
  instance_ = new DagMC(mb_impl);
}

void DagMC::set_useCAD( bool use_cad ){
  useCAD = use_cad;
  if( useCAD ){
    if( !have_cgm_geom ){
      std::cerr << "Warning: CAD-based ray tracing not avaiable, because CGM has no data." << std::endl;
      std::cerr << "         your input file was probably not a CAD format." << std::endl;
      useCAD = false;
    }

#ifndef HAVE_CGM_FIRE_RAY
    {
      std::cerr << "Warning: use_cad = 1 not supported with this build of CGM/DagMC." << std:: endl;
      std::cerr << "         Required ray-fire query not available. (Cubit-based CGM?)" <<  std::endl;
      useCAD = false;
    }
#endif
  }
}

DagMC::DagMC(Interface *mb_impl) 
  : mbImpl(mb_impl), obbTree(mb_impl), have_cgm_geom(false)
{
    // This is the correct place to uniquely define default values for the dagmc settings
  options[0] = Option( "source_cell",        "source cell ID, or zero if unknown", "0" );
  options[1] = Option( "discard_distance_tolerance", "positive real value", "1e-8" );
  options[2] = Option( "use_distance_limit", "one to enable distance limit optimization, zero otherwise", "0" );
  options[3] = Option( "use_cad", "one to ray-trace to cad, zero for just facets", "0" );
  options[4] = Option( "faceting_tolerance", "graphics faceting tolerance", "0.001" );
  options[5] = Option( "add_distance_tolerance", "positive real value", "1e-6" );

    // call parse settings to initialize default values for settings from options
  parse_settings();

  memset( specReflectName, 0, NAME_TAG_SIZE );
  strcpy( specReflectName, "spec_reflect" );
  memset( whiteReflectName, 0, NAME_TAG_SIZE );
  strcpy( whiteReflectName, "white_reflect" );
  memset( implComplName, 0, NAME_TAG_SIZE );
  strcpy( implComplName , "impl_complement" );

  distanceLimit = std::numeric_limits<double>::max();

  nameTag = get_tag(NAME_TAG_NAME, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, NULL, false);
  
  idTag = get_tag( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER );
  
  geomTag = get_tag( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER );

  obbTag = get_tag( MB_OBB_TREE_TAG_NAME, sizeof(EntityHandle), MB_TAG_DENSE, MB_TYPE_HANDLE );
  
    // get sense of surfaces wrt volumes
  senseTag = get_tag( "GEOM_SENSE_2", 2*sizeof(EntityHandle), MB_TAG_DENSE, MB_TYPE_HANDLE );

  int matid = 0;
  const void *def_matid = &matid;
  matTag   = get_tag(MATERIAL_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, def_matid );
  densTag  = get_tag(DENSITY_TAG_NAME, sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE );
  bcTag    = get_tag(BC_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER);

  double imp_one = 1;
  const void *def_imp = &imp_one;
  impTag   = get_tag(IMP_TAG_NAME, sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, def_imp );
  tallyTag = get_tag(TALLY_TAG_NAME, sizeof(int), MB_TAG_SPARSE, MB_TYPE_INTEGER);

}

ErrorCode DagMC::ray_fire(const EntityHandle vol, const EntityHandle last_surf_hit, 
                             const int num_pts,
                             const double uuu, const double vvv, const double www,
                             const double xxx, const double yyy, const double zzz,
                             const double huge_val,
                             double &dist_traveled, EntityHandle &next_surf_hit, 
                             OrientedBoxTreeTool::TrvStats* stats ) 
{
  if (debug) {
    std::cout << "Vol " << id_by_index(3, index_by_handle(vol)) << ", xyz = " 
              << xxx << " " << yyy << " " << zzz 
              << ", uvw = " 
              << uuu << " " << vvv << " " << www << std::endl;
  }
  
  assert(vol - setOffset < rootSets.size());
  
  EntityHandle root = rootSets[vol - setOffset];
  ErrorCode rval;
    // delcare some stuff static so we don't need to re-created
    // it for every call
  std::vector<double> &distances = distList;
  std::vector<EntityHandle> &surfaces = surfList;
  std::vector<EntityHandle> &facets = facetList;
  
  
    // do ray fire
  const double point[] = {xxx, yyy, zzz};
  const double dir[] = {uuu, vvv, www};
  distances.clear();
  surfaces.clear();
  facets.clear();
  double len = use_dist_limit() ? distance_limit() : huge_val;
  unsigned min_tolerance_intersections = 1000;

  rval = obbTree.ray_intersect_sets( distances,
                                     surfaces,
                                     facets, 
                                     root, 
                                     add_dist_tol(),
                                     min_tolerance_intersections,
                                     point, dir,
                                     &len, stats);
  assert( MB_SUCCESS == rval );


  // if useCAD is true at this point, then we know we can call CGM's ray casting function.
  if (useCAD) {
    rval = CAD_ray_intersect(point, dir, huge_val,
                             distances, surfaces, len);
    if (MB_SUCCESS != rval) return rval;
  }

  
    // Find smallest intersection
  if (distances.empty()) {
    next_surf_hit = 0;
    dist_traveled = (use_dist_limit() ? len*10.0 : huge_val);
    if (debug) {
      std::cout << "next surf hit = " << 0 << ", dist = (huge_val)" << std::endl;
    }
    return MB_SUCCESS;
  }
  int smallest = std::min_element( distances.begin(), distances.end() ) - distances.begin();

    // Sometimes a ray hits the boundary of two triangles.  In the next ray_fire the 
    // two hits are at zero distance.  An infinite loop occurs between the two neighboring
    // triangles.  This fix rejects surfaces closer than discard_dist_tol to avoid 
    // infinite looping at zero distance.  It required changing min_tolerance_intersections
    // to a very large number to make sure that we get at least one hit beyond add_dist_tol.
    // Last_surf_hit!=0 ensures that source particles close to a surface are
    // not rejected.
  while (   ( last_surf_hit!=0 )    &&    ( distances[smallest] < discard_dist_tol() )   ) {
    distances.erase( distances.begin() + smallest );
    surfaces.erase(  surfaces.begin()  + smallest );
  
      // Find smallest intersection
    if (distances.empty()) {
      next_surf_hit = 0;
      dist_traveled = (use_dist_limit() ? distance_limit()*10.0 : huge_val);
      if (debug) {
        std::cout << "next surf hit = " << 0 << ", dist = (huge_val)" << std::endl;
      }
      return MB_SUCCESS;
    }
    smallest = std::min_element( distances.begin(), distances.end() ) - distances.begin();
  }

    // return results
  dist_traveled = distances[smallest];
  next_surf_hit = surfaces[smallest];
  
  if (debug) {
    std::cout << "next surf hit = " <<  id_by_index(2, index_by_handle(next_surf_hit)) 
              << ", dist = " << dist_traveled << std::endl;
  }

  return MB_SUCCESS;
}

// **** HANDLING FILES TO INPUT DAGMC SETTINGS ****
// there is a movement away from using a settings file and, instead, incorporating this
// information into the native input file of the calling application 

void DagMC::read_settings( const char* filename )
{
  int num_opt = sizeof(options) / sizeof(options[0]);
  FILE* file;
  
  if (filename && (file = fopen( filename, "r" ))) {
    int line = 0;
    char buffer[256];
    bool havenl = true;
    while ( fgets( buffer, sizeof(buffer), file ) ) {
      if (havenl) {
        havenl = false;
        ++line;
      }
      int len = strlen(buffer);
      if (len && buffer[len-1] == '\n')
        havenl = true;
        
        // chop off any thing after the '#' comman indicator
      char* c = strchr( buffer, '#' );
      if (c) *c = '\0';
        // skip leading white space
      char* p = buffer;
      while (isspace(*p)) ++p;
        // if empty line, done
      if (!*p)
        continue;
        // find '='
      c = strchr( p, '=' );
      if (!c) {
        fprintf(stderr, "Invalid option at line %d of config file '%s'\n", line, filename );
        exit(2);
      }
        // get rid of white space around '='
      char* v = c+1;
      *c = ' ';
      while (c > p && isspace(*c)) {
        *c = '\0';
        --c;
      }
      while (isspace(*v))
        ++v;
      
        // search for option
      bool found = false;
      for (int i = 0; i < num_opt; ++i) {
        if (options[i].name == p) {
          found = true;
          options[i].value = v;
          options[i].user_set = true;
        }
      }
      if (!found) 
        fprintf(stderr,"Warning: unknown option at line %d of '%s': %s\n", line, filename, p );
    } // while( fgets() )
    
    fclose(file);
  } // file(file)
  
  // always do this, even if no file was found.
  // it will either read the default values, or re-parse
  // the values previously read from some other file.
  parse_settings();
}

  // this method sets the default values from options[] if none were read from file
void DagMC::parse_settings() {
  sourceCell = atoi( options[0].value.c_str() );
  if (sourceCell < 0) {
    std::cerr << "Invalid source_cell = " << sourceCell << std::endl;
    exit(2);
  }

  addDistTol = strtod( options[5].value.c_str(), 0 );
  if (addDistTol <= 0 || addDistTol > 1) {
    std::cerr << "Invalid add_distance_tolerance = " << addDistTol << std::endl;
    exit(2);
  }

  discardDistTol = strtod( options[1].value.c_str(), 0 );
  if (discardDistTol <= 0 || discardDistTol > 1) {
    std::cerr << "Invalid discard_distance_tolerance = " << discardDistTol << std::endl;
    exit(2);
  }

  useDistLimit = !!atoi( options[2].value.c_str() );

  set_useCAD( !!atoi( options[3].value.c_str() ) );

  facetingTolerance = strtod( options[4].value.c_str(), 0 );
  if (facetingTolerance <= 0) {
    std::cerr << "Invalid faceting_tolerance = " << facetingTolerance << std::endl;
    exit(2);
  }

}

void DagMC::write_settings( FILE* fptr, bool desc ) {
  int num_opt = sizeof(options) / sizeof(options[0]);
  if (desc) for (int i = 0; i < num_opt; ++i) 
    fprintf( fptr, "%s = %s  # %s\n", options[i].name.c_str(), options[i].value.c_str(), options[i].desc.c_str() );
  else for (int i = 0; i < num_opt; ++i) 
    fprintf( fptr, "%s = %s\n", options[i].name.c_str(), options[i].value.c_str() );
}

// **** PASS SETTINGS FROM CALLING APPLICATION ****
// This is the prefered way to set DAGMC settings

void DagMC::set_settings(int source_cell, int use_cad, int use_dist_limit,
			 double add_distance_tolerance, 
			 double discard_distance_tolerance) {
  sourceCell = source_cell;
  if (sourceCell < 0) {
    std::cerr << "Invalid source_cell = " << sourceCell << std::endl;
    exit(2);
  }

  std::cout << "Set Source Cell = " << sourceCell << std::endl;

  addDistTol = add_distance_tolerance;
  if (addDistTol <= 0 || addDistTol > 1) {
    std::cerr << "Invalid add_distance_tolerance = " << addDistTol << std::endl;
    exit(2);
  }

  std::cout << "Set add distance tolerance = " << addDistTol << std::endl;

  discardDistTol = discard_distance_tolerance;
  if (discardDistTol <= 0 || discardDistTol > 1) {
    std::cerr << "Invalid discard_distance_tolerance = " << discardDistTol << std::endl;
    exit(2);
  }

  std::cout << "Set discard distance tolerance = " << discardDistTol << std::endl;

  useDistLimit = !!(use_dist_limit);

  std::cout << "Turned " << (useDistLimit?"ON":"OFF") << " distance limit." << std::endl;

  set_useCAD( use_cad );

  std::cout << "Turned " << (useCAD?"ON":"OFF") << " ray firing on full CAD model." << std::endl;


}

void DagMC::get_settings(int *source_cell, int *use_cad, int *use_dist_limit,
			 double *discard_distance_tolerance, double *facet_tol) {

  *source_cell = sourceCell;

  *discard_distance_tolerance = discardDistTol;

  *use_dist_limit = !!(useDistLimit);

  *use_cad = !!(useCAD);

  *facet_tol = facetingTolerance;

}



// point_in_volume_slow, including poly_solid_angle helper subroutine
// are adapted from "Point in Polyhedron Testing Using Spherical Polygons", Paulo Cezar 
// Pinto Carvalho and Paulo Roma Cavalcanti, _Graphics Gems V_, pg. 42.  Original algorithm
// was described in "An Efficient Point In Polyhedron Algortihm", Jeff Lane, Bob Magedson, 
// and Mike Rarick, _Computer Visoin, Graphics, and Image Processing 26_, pg. 118-225, 1984.

// helper function for point_in_volume_slow.  calculate area of a polygon 
// projected into a unit-sphere space
ErrorCode DagMC::poly_solid_angle( EntityHandle face, const CartVect& point, double& area )
{
  ErrorCode rval;
  
    // Get connectivity
  const EntityHandle* conn;
  int len;
  rval = MBI->get_connectivity( face, conn, len, true );
  if (MB_SUCCESS != rval)
    return rval;
  
    // Allocate space to store vertices
  CartVect coords_static[4];
  std::vector<CartVect> coords_dynamic;
  CartVect* coords = coords_static;
  if ((unsigned)len > (sizeof(coords_static)/sizeof(coords_static[0]))) {
    coords_dynamic.resize(len);
    coords = &coords_dynamic[0];
  }
  
    // get coordinates
  rval = MBI->get_coords( conn, len, coords->array() );
  if (MB_SUCCESS != rval)
    return rval;
  
    // calculate normal
  CartVect norm(0.0), v1, v0 = coords[1] - coords[0];
  for (int i = 2; i < len; ++i) {
    v1 = coords[i] - coords[0];
    norm += v0 * v1;
    v0 = v1;
  }
  
    // calculate area
  double s, ang;
  area = 0.0;
  CartVect r, n1, n2, b, a = coords[len-1] - coords[0];
  for (int i = 0; i < len; ++i) {
    r = coords[i] - point;
    b = coords[(i+1)%len] - coords[i];
    n1 = a * r;
    n2 = r * b;
    s = (n1 % n2) / (n1.length() * n2.length());
    ang = s <= -1.0 ? M_PI : s >= 1.0 ? 0.0 : acos(s);
    s = (b * a) % norm;
    area += s > 0.0 ? M_PI - ang : M_PI + ang;
    a = -b;
  }
  
  area -= M_PI * (len - 2);
  if ((norm % r) > 0)
    area = -area;
  return MB_SUCCESS;
}
  
  
// use spherical area test to determine inside/outside of a polyhedron.
ErrorCode DagMC::point_in_volume_slow( EntityHandle volume,
                                  double x, double y, double z,
                                  int& result )
{
  ErrorCode rval;
  Range faces;
  std::vector<EntityHandle> surfs;
  std::vector<int> senses;
  double sum = 0.0;
  const CartVect point(x,y,z);
  
  rval = MBI->get_child_meshsets( volume, surfs );
  if (MB_SUCCESS != rval)
    return rval;
  
  senses.resize( surfs.size() );
  rval = surface_sense( volume, surfs.size(), &surfs[0], &senses[0] );
  if (MB_SUCCESS != rval)
    return rval;
  
  for (unsigned i = 0; i < surfs.size(); ++i) {
    if (!senses[i])  // skip non-manifold surfaces
      continue;
     
    double surf_area = 0.0, face_area;
    faces.clear();
    rval = MBI->get_entities_by_dimension( surfs[i], 2, faces );
    if (MB_SUCCESS != rval)
      return rval;
    
    for (Range::iterator j = faces.begin(); j != faces.end(); ++j) {
      rval = poly_solid_angle( *j, point, face_area );
      if (MB_SUCCESS != rval)
        return rval;
      
      surf_area += face_area;
    }
    
    sum += senses[i] * surf_area;
  }
  
  result = fabs(sum) > 2.0*M_PI;
  return MB_SUCCESS;
}

// If point is on boundary, then this function is called to 
// discriminate cases in which the ray is entering or leaving.
// result= 1 -> inside volume or entering volume
// result= 0 -> outside volume or leaving volume
// result=-1 -> on boundary with null or tangent uvw
ErrorCode DagMC::boundary_case( EntityHandle volume, int& result, 
				  double u, double v, double w,
				  EntityHandle facet,
				  EntityHandle surface)
{
  ErrorCode rval;

  // test to see if uvx is provided
  if ( u <= 1.0 && v <= 1.0 && w <= 1.0 ) {

    const CartVect ray_vector(u, v, w);
    CartVect coords[3], normal(0.0);
    const EntityHandle *conn;
    int len, sense_out;
   
    rval = mbImpl->get_connectivity( facet, conn, len );
      assert( MB_SUCCESS == rval );
      assert( 3 == len );
  
    rval = mbImpl->get_coords( conn, 3, coords[0].array() );
      assert(MB_SUCCESS == rval);
   
    rval = surface_sense( volume, surface, sense_out );
      assert( MB_SUCCESS == rval);

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal = sense_out * (coords[1] * coords[2]);

    double sense = ray_vector % normal;

    if ( sense < 0.0 ) {
      result = 1;     // inside or entering
    } else  if ( sense > 0.0 ) {
      result = 0;     // outside or leaving
    } else  if ( sense == 0.0 ) {
      result = -1;    // tangent, therefore on boundary
    } else {
      result = -1;    // failure
      return MB_FAILURE;
    }
  
  // if uvw not provided, return on_boundary.
  } else {
    result = -1;      // on boundary
    return MB_SUCCESS;
  
  }

  return MB_SUCCESS;
}

// Fast point_in_volume code.  This function assumes that there is an 
// OBB tree containing only manifold surfaces in the volume, such that
// the closest facet to the input position found using the OBB tree is
// guaranteed to be on the boundary of the volume.  
//
// This function will find the closest point on the volume boundary.
// If that point is closest to the interior of a facet, the direction
// of the surface normal of that facet can be used to determine if the
// point is inside or outside of the volume.  If the point is closest to
// an edge joining two facets, then the point is inside the volume if 
// that edge is at a concavity wrt the volume.  If the point is closest
// to a vertex in the facetting, it relies on the fact that the closest
// vertex cannot be a saddle: it must be a strict concavity or convexity.
ErrorCode DagMC::point_in_volume( EntityHandle volume, 
                             double x, double y, double z,
                             int& result,
			     double u, double v, double w)
{
  ErrorCode rval;
  const double epsilon = discardDistTol;
  const double boundary_epsilon = epsilon;

    // Get OBB Tree for volume
  assert(volume - setOffset < rootSets.size());
  EntityHandle root = rootSets[volume - setOffset];

    // Get closest point in triangulation
  const CartVect point(x,y,z);
  EntityHandle facet, surface;
  CartVect closest, diff;
  rval = obbTree.closest_to_location( point.array(), root, closest.array(), facet, &surface );
  if (MB_SUCCESS != rval)  return rval;
  
    // Check for on-boundary case
  diff = closest - point;
  if (diff%diff <= boundary_epsilon*boundary_epsilon) {
    return boundary_case( volume, result, u, v, w, facet, surface );
  }
  
    // Get triangles at closest point
  std::vector<EntityHandle> &tris = triList, &surfs = surfList;
  tris.clear();
  surfs.clear();
  rval = obbTree.sphere_intersect_triangles( closest.array(), epsilon, root, tris, &surfs );
  if (MB_SUCCESS != rval) return rval;
  
    // One-triangle case : determine if point is above or below triangle
  const EntityHandle* conn;
  int len, sense;
  CartVect coords[3], normal;
  if (tris.size() == 1) {
    rval = MBI->get_connectivity( tris.front(), conn, len );
    if (MB_SUCCESS != rval || len != 3) 
      return MB_FAILURE;
    rval = MBI->get_coords( conn, len, coords[0].array() );
    if (MB_SUCCESS != rval) return rval;
    
    rval = surface_sense( volume, surfs.front(), sense );
    if (MB_SUCCESS != rval) return rval;

      // get triangle normal
    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal = sense * coords[1] * coords[2];
      // compare relative sense
    result = (normal % (point - closest) < 0.0);
    return MB_SUCCESS;
  }
  
    // many triangle case : determine if closest point is convexity or concavity
  
  
    // use algorithm from:
    // (referance)
    // to determine inside vs. outside.
  
    // first find single triangle for which all other triangles are to
    // one side.
  bool some_above, some_below;
  for (unsigned i = 0; i < tris.size(); ++i) {
    some_above = some_below = false;

    rval = MBI->get_connectivity( tris[i], conn, len );
    if (MB_SUCCESS != rval || len != 3) 
      return MB_FAILURE;
    rval = MBI->get_coords( conn, len, coords[0].array() );
    if (MB_SUCCESS != rval) return rval;
    
    rval = surface_sense( volume, surfs[i], sense );
    if (MB_SUCCESS != rval) return rval;

      // get triangle normal
    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal = sense * coords[1] * coords[2];
    closest = coords[0];
    const double norm_len_sqr = normal % normal;

    for (unsigned j = 0; j < tris.size(); ++j) {
      if (j == i)
        continue;
      
      rval = MBI->get_connectivity( tris[j], conn, len );
      if (MB_SUCCESS != rval || len != 3)
        return MB_FAILURE;
      rval = MBI->get_coords( conn, len, coords[0].array() );
      if (MB_SUCCESS != rval)
        return rval;
      
      for (int k = 0; k < len; ++k) {
        const double dot = normal % (coords[k] - closest);
        if (dot * dot / norm_len_sqr > epsilon * epsilon) {
          if (dot > 0.0)
            some_above = true;
          else
            some_below = true;
        }
      }
    }
    
      // All triangles are roughly co-planar: if we're inside of
      // one then we're inside of all.
    if (!some_above && !some_below) {
      result = (normal % (closest - point)) > 0.0;
      return MB_SUCCESS;
    }
      // All other triangles to one side of this triangle:
    else if (some_above != some_below) {
        // If all other triangles are above this one (some_above == true),
        // then the vertex is a concavity so the input point is inside.
        // If all other triangles are below this one (some_above == false),
        // then the vertex is a convexity and the input point must be outside.
      result = some_above;
      return MB_SUCCESS;
    }
  }
  
    // If we got here, then something's wrong somewhere.
    // We appear to be closest to a "saddle" vertex in the
    // triangulation.  That shoudn't be possible (must be 
    // closer to at least one of the adjacent triagles than
    // to the saddle vertex.)
  std::cout << "point_in_volume fast test failure: xyz= " << x << " " << y << " " << z 
            << " uvw= " << u << " " << v << " " << w 
            << " vol=" << id_by_index(3, index_by_handle(volume)) << std::endl; 
  // For now, proceed with the slow test.
  return point_in_volume_slow( volume, x, y, z, result );
}

// detemine distance to nearest surface
ErrorCode DagMC::closest_to_location( EntityHandle volume, double* coords, double& result)
{
    // Get OBB Tree for volume
  assert(volume - setOffset < rootSets.size());
  EntityHandle root = rootSets[volume - setOffset];

    // Get closest triangles in volume
  const CartVect point(coords);
  CartVect nearest;
  EntityHandle facet_out;
  ErrorCode rval = obbTree.closest_to_location( point.array(), root, nearest.array(), facet_out );
  if (MB_SUCCESS != rval) return rval;

  // calculate distance between point and nearest facet
  result = (point-nearest).length();
  
  return MB_SUCCESS;

}




// calculate volume of polyhedron
ErrorCode DagMC::measure_volume( EntityHandle volume, double& result )
{
  ErrorCode rval;
  std::vector<EntityHandle> surfaces, surf_volumes;
  result = 0.0;
  
   // don't try to calculate volume of implicit complement
  if (volume == impl_compl_handle) {
    result = 1.0;
    return MB_SUCCESS;
  }

    // get surfaces from volume
  rval = MBI->get_child_meshsets( volume, surfaces );
  if (MB_SUCCESS != rval) return rval;
  
    // get surface senses
  std::vector<int> senses( surfaces.size() );
  rval = surface_sense( volume, surfaces.size(), &surfaces[0], &senses[0] );
  if (MB_SUCCESS != rval) {
    std::cerr << "ERROR: Surface-Volume relative sense not available. "
              << "Cannot calculate volume." << std::endl;
    return rval;
  }
  
  for (unsigned i = 0; i < surfaces.size(); ++i) {
      // skip non-manifold surfaces
    if (!senses[i])
      continue;
    
      // get triangles in surface
    Range triangles;
    rval = MBI->get_entities_by_dimension( surfaces[i], 2, triangles );
    if (MB_SUCCESS != rval) 
      return rval;
    if (!triangles.all_of_type(MBTRI)) {
      std::cout << "WARNING: Surface " << id_by_index(2, index_by_handle(surfaces[i]))
                << " contains non-triangle elements. Volume calculation may be incorrect." 
                << std::endl;
      triangles.clear();
      rval = MBI->get_entities_by_type( surfaces[i], MBTRI, triangles );
      if (MB_SUCCESS != rval) return rval;
    }
    
      // calculate signed volume beneath surface (x 6.0)
    double surf_sum = 0.0;
    const EntityHandle *conn;
    int len;
    CartVect coords[3];
    for (Range::iterator j = triangles.begin(); j != triangles.end(); ++j) {
      rval = MBI->get_connectivity( *j, conn, len, true );
      if (MB_SUCCESS != rval) return rval;
      assert(3 == len);
      rval = MBI->get_coords( conn, 3, coords[0].array() );
      if (MB_SUCCESS != rval) return rval;
    
      coords[1] -= coords[0];
      coords[2] -= coords[0];
      surf_sum += (coords[0] % (coords[1] * coords[2]));
    }
    result += senses[i] * surf_sum;
  }
  
  result /= 6.0;
  return MB_SUCCESS;
}

// sum area of elements in surface
ErrorCode DagMC::measure_area( EntityHandle surface, double& result )
{
    // get triangles in surface
  Range triangles;
  ErrorCode rval = MBI->get_entities_by_dimension( surface, 2, triangles );
  if (MB_SUCCESS != rval) 
    return rval;
  if (!triangles.all_of_type(MBTRI)) {
    std::cout << "WARNING: Surface " << id_by_index(2, index_by_handle(surface))
              << " contains non-triangle elements. Area calculation may be incorrect." 
              << std::endl;
    triangles.clear();
    rval = MBI->get_entities_by_type( surface, MBTRI, triangles );
    if (MB_SUCCESS != rval) return rval;
  }

    // calculate sum of area of triangles
  result = 0.0;
  const EntityHandle *conn;
  int len;
  CartVect coords[3];
  for (Range::iterator j = triangles.begin(); j != triangles.end(); ++j) {
    rval = MBI->get_connectivity( *j, conn, len, true );
    if (MB_SUCCESS != rval) return rval;
    assert(3 == len);
    rval = MBI->get_coords( conn, 3, coords[0].array() );
    if (MB_SUCCESS != rval) return rval;

    coords[1] -= coords[0];
    coords[2] -= coords[0];
    coords[0] = coords[1] * coords[2];
    result += coords[0].length();
  }
  result *= 0.5;
  return MB_SUCCESS;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense( EntityHandle volume, 
                           int num_surfaces,
                           const EntityHandle* surfaces,
                           int* senses_out )
{

  /* The sense tags do not reference the implicit complement handle.
     All surfaces that interact with the implicit complement should have
     a null handle in the direction of the implicit complement. */
  if (volume == impl_compl_handle)
    volume = (EntityHandle) 0;

  std::vector<EntityHandle> surf_volumes( 2*num_surfaces );
  ErrorCode rval = MBI->tag_get_data( sense_tag(), surfaces, num_surfaces, &surf_volumes[0] );
  if (MB_SUCCESS != rval)  return rval;
  
  const EntityHandle* end = surfaces + num_surfaces;
  std::vector<EntityHandle>::const_iterator surf_vols = surf_volumes.begin();
  while (surfaces != end) {
    EntityHandle forward = *surf_vols; ++surf_vols;
    EntityHandle reverse = *surf_vols; ++surf_vols;
    if (volume == forward) 
      *senses_out = (volume != reverse); // zero if both, otherwise 1
    else if (volume == reverse)
      *senses_out = -1;
    else 
      return MB_ENTITY_NOT_FOUND;
    
    ++surfaces;
    ++senses_out;
  }
  
  return MB_SUCCESS;
}

// get sense of surface(s) wrt volume
ErrorCode DagMC::surface_sense( EntityHandle volume, 
                                  EntityHandle surface,
                                  int& sense_out )
{
  /* The sense tags do not reference the implicit complement handle.
     All surfaces that interact with the implicit complement should have
     a null handle in the direction of the implicit complement. */
  if (volume == impl_compl_handle)
    volume = (EntityHandle) 0;

    // get sense of surfaces wrt volumes
  EntityHandle surf_volumes[2];
  ErrorCode rval = MBI->tag_get_data( sense_tag(), &surface, 1, surf_volumes );
  if (MB_SUCCESS != rval)  return rval;
  
  if (surf_volumes[0] == volume)
    sense_out = (surf_volumes[1] != volume); // zero if both, otherwise 1
  else if (surf_volumes[1] == volume) 
    sense_out = -1;
  else
    return MB_ENTITY_NOT_FOUND;
  
  return MB_SUCCESS;
}

ErrorCode DagMC::load_file(const char* cfile,
			   const double facet_tolerance)
{
  ErrorCode rval;
  
#ifdef CGM
  // cgm must be initialized so we can check it for CAD data after the load
  InitCGMA::initialize_cgma();
#endif

    // override default value of facetingTolerance with passed value
  if (facet_tolerance > 0 )
    facetingTolerance = facet_tolerance;

  char facetTolStr[16];

  sprintf(facetTolStr,"%g",facetingTolerance);

  // what if we are using default faceting tolerance???
  char options[120] = "CGM_ATTRIBS=yes;FACET_DISTANCE_TOLERANCE=";
  strcat(options,facetTolStr);
  
  EntityHandle file_set;
  rval = MBI->create_meshset( MESHSET_SET, file_set );
  if (MB_SUCCESS != rval)
    return rval;

  rval = MBI->load_file(cfile, &file_set, options, NULL, 0, 0);
  
  if( MB_UNHANDLED_OPTION == rval ){
    // Some options were unhandled; this is common for loading h5m files.
    // Print a warning if an option was unhandled for a file that does not end in '.h5m'
    std::string filename(cfile);
    if( filename.length() < 4 || filename.substr(filename.length()-4) != ".h5m"){
      std::cerr << "DagMC warning: unhandled file loading options." << std::endl;
    }
  }
  else if (MB_SUCCESS != rval) {
    std::cerr << "DagMC Couldn't read file " << cfile << std::endl;
    return rval;
  }

#ifdef CGM  
  // check to see if CGM has data; if so, assume it corresponds to the data we loaded in.
  if( GeometryQueryTool::instance()->num_ref_volumes() > 0 ){
    have_cgm_geom = true;
  }
#endif


  // search for a tag that has the faceting tolerance
  Tag faceting_tol_tag;
  Range tagged_sets;
  double facet_tol_tagvalue = 0;
  bool other_set_tagged = false, root_tagged = false;
  rval = MBI->tag_create( "FACETING_TOL", sizeof(double), MB_TAG_SPARSE,
			  MB_TYPE_DOUBLE, faceting_tol_tag, 0, true );
  // get list of entity sets that are tagged with faceting tolerance 
  // (possibly empty set)
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &faceting_tol_tag,
					    NULL, 1, tagged_sets );
  // if NOT empty set
  if (MB_SUCCESS == rval && !tagged_sets.empty()) {
    rval = MBI->tag_get_data( faceting_tol_tag, &(*tagged_sets.begin()), 1, &facet_tol_tagvalue );
    if (MB_SUCCESS != rval) return rval;
    other_set_tagged = true;
  }
  else if (MB_SUCCESS == rval) {
    // check to see if interface is tagged
    rval = MBI->tag_get_data( faceting_tol_tag, 0, 0, &facet_tol_tagvalue );
    if (MB_SUCCESS == rval) root_tagged = true;
    else rval = MB_SUCCESS;
  }

  if ( (root_tagged || other_set_tagged) && facet_tol_tagvalue > 0) {
    facetingTolerance = facet_tol_tagvalue;
  }

  std::cout << "Using faceting tolerance: " << facetingTolerance << std::endl;

  return MB_SUCCESS;

}

ErrorCode DagMC::init_OBBTree() 
{

  ErrorCode rval;

  Range surfs, vols;
  const int three = 3;
  const void* const three_val[] = {&three};
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, 
				       three_val, 1, vols );
  if (MB_SUCCESS != rval)
    return rval;

  const int two = 2;
  const void* const two_val[] = {&two};
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, 
					    two_val, 1, surfs );
  if (MB_SUCCESS != rval)
    return rval;

    // If it doesn't already exist, create implicit complement
    // Create data structures for implicit complement
  rval = get_impl_compl();
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to find or create implicit complement handle." << std::endl;
    return rval;
  }

  // Build OBB trees for everything, but only if we only read geometry
  // Changed to build obb tree if tree does not already exist. -- JK
  if (!have_obb_tree()) {
    rval = build_obbs(surfs, vols);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to build obb." << std::endl;
      return rval;
    }
  }

    // build the various index vectors used for efficiency
  rval = build_indices(surfs, vols);
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to build surface/volume indices." << std::endl;
    return rval;
  }
  
  return MB_SUCCESS;
}


ErrorCode DagMC::write_mesh(const char* ffile,
			      const int flen)
{
  ErrorCode rval;
  
    // write out a mesh file if requested
  if (ffile && 0 < flen) {
    rval = MBI->write_mesh(ffile);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to write mesh to " << ffile << "." << std::endl;
      return rval;
    }
  }
 
  return MB_SUCCESS;
 
}

bool DagMC::have_obb_tree()
{
  Range entities;
  ErrorCode rval = mbImpl->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                           &obbTag, 0, 1,
                                                           entities );
  return MB_SUCCESS == rval && !entities.empty();
}                                                    

EntityHandle DagMC::entity_by_id( int dimension, int id )
{
  assert(0 <= dimension && 3 >= dimension);
  const Tag tags[] = { idTag, geomTag };
  const void* const vals[] = { &id, &dimension };
  ErrorCode rval;
  
  Range results;
  rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, results );
  if (MB_SUCCESS != rval || results.empty())
    return 0;
  
  return results.front();
}

int DagMC::id_by_index( int dimension, int index )
{
  EntityHandle h = entity_by_index( dimension, index );
  if (!h)
    return 0;
  
  int result = 0;
  MBI->tag_get_data( idTag, &h, 1, &result );
  return result;
}

// functionality in write_log
//
// 1. for each group
//    a) check for material assignment -> get list of member vols and set mat/dens array
//    b) check for b.c. -> get list of surfaces and set bc array (incl. periodic???)
//    c) check for graveyard -> set importance
//    d) check for tallies
//       * record tally ID
// 2. write cell cards
// 3. write surf cards
// 4. if tallies
//    a) loop through all groups to find tallies
//       * parse tally type, etc
//       * get list of cells or surfaces
//       * write tally cards



ErrorCode DagMC::write_mcnp(std::string ifile, const bool overwrite) 
{
  std::map<std::string,int> tallyKeywords;
  
  tallyKeywords["surf.current"] = 1;
  tallyKeywords["surf.flux"]    = 2;
  tallyKeywords["cell.flux"]    = 4;
  tallyKeywords["cell.heating"] = 6;
  tallyKeywords["cell.fission"] = 7;
  tallyKeywords["pulse.height"] = 8;

  std::vector<EntityHandle>::iterator iter;

  if (!overwrite) {
      // if not overwriting, test for file, and if it exists, just return;
      // this replaces the old inquir function
    std::ifstream testfile;
    testfile.open(ifile.c_str());

    if (!testfile.fail()) return MB_SUCCESS;
  }
  
  const char *cifile = ifile.c_str();
  std::ofstream cgmfile;
  cgmfile.open(cifile);

  std::cout << "  # cells: " << vol_handles().size()-1 << std::endl;
  std::cout << "  # surfs: " << surf_handles().size()-1 << std::endl;
  std::cout << " # groups: " << group_handles().size()-1 << std::endl;
  std::cout << "# tallies: " << tallyList.size() << std::endl;

  Range grp_sets, grp_ents;
  ErrorCode rval;

  Range surfs, vols;


  // write cell information (skip first entry)             
  for (iter = vol_handles().begin()+1; iter != vol_handles().end(); iter++) {
    
    int cellid,matid;
    double imp, density;

    rval = MBI->tag_get_data( idTag , &(*iter), 1, &cellid   );
    if (MB_SUCCESS != rval) return rval;
    rval = MBI->tag_get_data(impTag , &(*iter), 1, &imp      );
    if (MB_SUCCESS != rval) return rval;
    rval = MBI->tag_get_data(matTag , &(*iter), 1, &matid    );
    if (MB_SUCCESS != rval) return rval;
    if (0 == matid) {
      cgmfile << cellid << " " << matid << " imp:n=" << imp ;
    } else {
      rval = MBI->tag_get_data(densTag, &(*iter), 1, &density  );
      if (MB_SUCCESS != rval) return rval;
      cgmfile << cellid << " " << matid << " " << density << " imp:n=" << imp ;
    }      
    if (*iter == impl_compl_handle)
      cgmfile << "   $ implicit complement ";
    cgmfile << std::endl;
  }

  // skip a line                                                               
  cgmfile << std::endl;

  // write surface info (skip first entry)
  for (iter = surf_handles().begin()+1; iter != surf_handles().end(); iter++ ) {
    
    int surfid, bc_id;

    rval = MBI->tag_get_data( idTag, &(*iter), 1, &surfid );
    if (MB_SUCCESS != rval) return rval;
    rval = MBI->tag_get_data( bcTag, &(*iter), 1, &bc_id  );
    if (MB_SUCCESS != rval && MB_TAG_NOT_FOUND != rval) return rval;

    if (MB_TAG_NOT_FOUND != rval) {
      if (BC_SPEC == bc_id)
	cgmfile << "*";
      if (BC_WHITE == bc_id)
	cgmfile << "+";
    }
    cgmfile << surfid << std::endl;

  }

  // add a final blank line                                                    
  cgmfile << std::endl;
  
  std::vector<std::string> grp_names;
  std::vector<std::string> tokens;

  // write tally info
  for (std::vector<int>::iterator grp = tallyList.begin(); grp != tallyList.end(); grp++) {
    
    std::stringstream tallyCard;
    
    EntityHandle group = group_handles()[*grp];
    grp_names.clear();
    bool success = get_group_names(group, grp_names);
    if (!success)
      return MB_FAILURE;
    if (grp_names.empty()) continue;
    
    // get sets associated with this group
    grp_sets.clear();
    rval = MBI->get_entities_by_type(group, MBENTITYSET, grp_sets);
    if (MB_SUCCESS != rval) continue;

    tokens.clear();
    tokenize(grp_names[0],tokens,"_");
    
    // Get user number for tally
    int talNum = atoi(tokens[1].c_str());
    char talMod = ' ';
    
    // Get tally modifier 
    if ( 'e' == tokens[2][0] )
      talMod = '*';
    if ( 'q' == tokens[2][0] )
      talMod = '+';
    if (talMod != ' ')
      tokens[2].erase(0,1);  // strip special character for modified tallies
    
    // Get tally type
    if ( tallyKeywords.count(tokens[2]) == 0) {
      std::cout << "Invalid tally defined in geometry: " << grp_names[0] << std::endl;
      continue;
    }
    
    // Create MCNP number for tally by combining user number and MCNP type code
    talNum = talNum*10 + tallyKeywords[tokens[2]];
    
    // get tally particles
    std::string partList = "n";
    if (tokens.size()>3)
      {
	partList = tokens[3].substr(0,1);
	for (std::string::iterator pos = tokens[3].begin()+1; pos != tokens[3].end(); pos++)
	  {
	    partList += ",";
	    partList += *pos;
	  } 
      }
    
    // write tally type info
    tallyCard << talMod << "f" << talNum << ":" << partList ;
    
    int dim;
    if (tallyKeywords[tokens[2]] < 3) 
      dim = 2;
    else
      dim = 3;

    // get entities of correct dimension
    const void* dim_val[] = {&dim};
    Range ents;
    rval = MBI->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, dim_val, 1, ents );
    if (MB_SUCCESS != rval)
      return rval;
    Range  grp_ents = intersect(grp_sets,ents);

    for (Range::iterator ent = grp_ents.begin(); ent != grp_ents.end(); ent++)
      tallyCard << " " << get_entity_id(*ent);

    tallyCard << " T";

    std::string tallyCardStr = tallyCard.str();
    tallyCard.str("");

    // write tallyCardStr with 80 char limit per line
    
    while ( tallyCardStr.length() > 72 )
      {
	std::string::size_type pos = tallyCardStr.rfind(' ',72);
	cgmfile << tallyCardStr.substr(0,pos) << " &" << std::endl;
	tallyCardStr.erase(0,pos);

	cgmfile << "     ";
      }
    cgmfile << tallyCardStr << std::endl;

  }

  cgmfile.close();
  return MB_SUCCESS;
}



// desired functionality
//
// 1. for each group
//    a) check for material assignment -> get list of member vols and tag mat & dens
//    b) check for b.c. -> get list of surfaces and tag b.c.
//    c) check for graveyard -> tag imp

ErrorCode DagMC::parse_metadata()
{
  std::vector<std::string> grp_names;
  std::vector<std::string> tokens;
  std::map<std::string,int> keywords;

  EntityHandle group;
  Range grp_sets, grp_ents;
  ErrorCode rval;

  Range vols, surfs;

  std::copy(vol_handles().begin()+1,  vol_handles().end(),
	    range_inserter(vols));
  std::copy(surf_handles().begin()+1, surf_handles().end(),
	    range_inserter(surfs));

  keywords["mat"]           = MAT_GROUP;
  keywords["spec.reflect"]  = BC_SPEC;
  keywords["white.reflect"] = BC_WHITE;
  keywords["graveyard"]     = IMP_ZERO;
  keywords["outside.world"] = IMP_ZERO;
  keywords["rest.of.world"] = IMP_ZERO;
  keywords["tally"]         = TALLY_GROUP;

  for (unsigned int grp=0; grp < group_handles().size(); grp++) {

    // get group handle & names
    group = group_handles()[grp];
    grp_names.clear();
    bool success = get_group_names(group, grp_names);
    if (!success)
      return MB_FAILURE;
    if (grp_names.empty()) continue;

    // get sets associated with this group
    grp_sets.clear();
    rval = MBI->get_entities_by_type(group, MBENTITYSET, grp_sets);
    if (MB_SUCCESS != rval) continue;

    tokens.clear();
    tokenize(grp_names[0],tokens,"_");

    // is this a keyword?
    if (keywords.count(tokens[0]) == 0 ) {
      std::cout << "Ignoring group with name: " << grp_names[0] << std:: endl;
      continue;
    }

    int matid, bc_id;
    std::vector<int> matid_list, bc_id_list;
    double density, imp_zero = 0;
    std::vector<double> density_list, imp_zero_list;

    // actions for each keyword
    switch (keywords[tokens[0]]) {
    case MAT_GROUP:
      // missing density
      if ( tokens[2].compare("rho") != 0 ) {
	std::cout << "Specified material with no density " << grp_names[0] << std::endl;
	continue;
      }

        matid = atoi(tokens[1].c_str());
      density = atof(tokens[3].c_str());
      
      // check to see if this is the implicit complement
      if ( tokens.size() > 4 && tokens[4].compare("comp") == 0 ) {
        std::cout<<"parse_metadata: complement material and density specified" << std::endl;
	// TAG IMPL COMPL WITH MAT & DENS
 	MBI->tag_set_data(matTag ,&impl_compl_handle,1,&matid);
 	MBI->tag_set_data(densTag,&impl_compl_handle,1,&density);
      } else {
	// apply tag to range of volumes
	grp_ents.clear();
	grp_ents = intersect(grp_sets,vols);
	  matid_list.assign(grp_ents.size(),matid);
	density_list.assign(grp_ents.size(),density);
	MBI->tag_set_data(matTag ,grp_ents,&(*matid_list.begin()));
	MBI->tag_set_data(densTag,grp_ents,&(*density_list.begin()));
      }
      break;
    case BC_SPEC:
    case BC_WHITE:
      bc_id = keywords[tokens[0]];
      grp_ents.clear();
      grp_ents = intersect(grp_sets,surfs);
      bc_id_list.assign(grp_ents.size(),bc_id);
      MBI->tag_set_data(bcTag ,grp_ents, &(*bc_id_list.begin()));
      break;
    case IMP_ZERO:
      grp_ents.clear();
      grp_ents = intersect(grp_sets,vols);
      imp_zero_list.assign(grp_ents.size(),imp_zero);
      MBI->tag_set_data(impTag ,grp_ents, &(*imp_zero_list.begin()));
      break;
    case TALLY_GROUP:
      // make list of groups that are tallies
      grp_ents.clear();
      grp_ents = intersect(grp_sets,vols);
      tallyList.push_back(grp);
      break;
    }
    
  }

  return MB_SUCCESS;
}

void DagMC::tokenize( const std::string& str,
		      std::vector<std::string>& tokens,
		      const char* delimiters )
{
  std::string::size_type last = str.find_first_not_of( delimiters, 0 );
  std::string::size_type pos  = str.find_first_of( delimiters, last );
  if ( std::string::npos == pos )
    tokens.push_back(str);
  else
    while (std::string::npos != pos && std::string::npos != last) {
      tokens.push_back( str.substr( last, pos - last ) );
      last = str.find_first_not_of( delimiters, pos );
      pos  = str.find_first_of( delimiters, last );
      if(std::string::npos == pos)
	pos = str.size();
    }
}

std::string DagMC::itos(int ival) {
  std::stringstream s;
  s << ival;
  return s.str();
}

bool DagMC::get_group_names(EntityHandle group_set, 
                             std::vector<std::string> &grp_names) 
{
    // get names
  char name0[NAME_TAG_SIZE];
  std::fill(name0, name0+NAME_TAG_SIZE, '\0');
  ErrorCode result = MBI->tag_get_data(name_tag(), &group_set, 1,
                                                     &name0);
  if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) return false;

  if (MB_TAG_NOT_FOUND != result) grp_names.push_back(std::string(name0));

  int extra_num = 0;
  while (true) {
    sprintf(name0, "%s%s%d", "EXTRA_", NAME_TAG_NAME, extra_num);
    extra_num++;
    Tag this_tag = get_tag(name0, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, 
                             false);
    if (0 == this_tag) break;
    std::fill(name0, name0+NAME_TAG_SIZE, '\0');
    result = MBI->tag_get_data(this_tag, &group_set, 1, &name0);
    if (MB_SUCCESS != result && MB_TAG_NOT_FOUND != result) return false;
    if (MB_TAG_NOT_FOUND == result) break;
    else grp_names.push_back(std::string(name0));
  }
  
  return true;
}

Tag DagMC::get_tag( const char* name, int size, TagType store, 
		      DataType type, const void* def_value,
		      bool create_if_missing) 
{
  Tag retval = 0;
  ErrorCode result = MBI->tag_create(name, size, store, type,
                                                   retval, def_value, create_if_missing);
  if (create_if_missing && MB_SUCCESS != result) 
    std::cerr << "Couldn't find nor create tag named " << name << std::endl;
  
  return retval;
}

int DagMC::get_entity_id(EntityHandle this_ent) 
{
  int id = 0;
  ErrorCode result = MBI->tag_get_data(idTag, &this_ent, 1, &id);
  if (MB_TAG_NOT_FOUND == result)
    id = MBI->id_from_handle(this_ent);
    
  return id;
}

ErrorCode DagMC::get_angle(EntityHandle surf, 
                              double xxx, double yyy, double zzz, double *ang)
{
  EntityHandle root = rootSets[surf - setOffset];
  
  const double in_pt[] = { xxx, yyy, zzz };
  std::vector<EntityHandle> &facets = triList;
  facets.clear();
  ErrorCode rval = obbTree.closest_to_location( in_pt, root, add_dist_tol(), facets );
  assert(MB_SUCCESS == rval);
  
  CartVect coords[3], normal(0.0);
  const EntityHandle *conn;
  int len;
  for (unsigned i = 0; i < facets.size(); ++i) {
    rval = mbImpl->get_connectivity( facets[i], conn, len );
    assert( MB_SUCCESS == rval );
    assert( 3 == len );
  
    rval = mbImpl->get_coords( conn, 3, coords[0].array() );
    assert(MB_SUCCESS == rval);
    
    coords[1] -= coords[0];
    coords[2] -= coords[0];
    normal += coords[1] * coords[2];
  }
  
  normal.normalize();
  normal.get( ang );

  return MB_SUCCESS;
}

ErrorCode DagMC::CAD_ray_intersect(const double *point, 
                                      const double *dir, 
                                      const double huge_val,
                                      std::vector<double> &distances,
                                      std::vector<EntityHandle> &surfaces, 
                                      double &len) 
{
#ifdef CGM
#ifndef HAVE_CGM_FIRE_RAY
  return MB_NOT_IMPLEMENTED;
#else
  std::vector<double>::iterator dit = distances.begin();
  std::vector<EntityHandle>::iterator sit = surfaces.begin();
  static DLIList<double> ray_params;
  
  for (; dit != distances.end(); dit++, sit++) {
      // get the RefFace
    RefEntity *this_face = geomEntities[*sit - setOffset];
      // get the ray distance to this face
    ray_params.clean_out();
    int result = GeometryQueryTool::instance()->
      fire_ray(dynamic_cast<RefFace*>(this_face), CubitVector(point),
               CubitVector(dir), ray_params);
    assert(CUBIT_SUCCESS == result);
    if (ray_params.size() != 0) {
      ray_params.reset();
      *dit = ray_params.get();
    }
    else *dit = huge_val;
  }
  
    // now bubble sort list
  bool done = false;
  while (!done) {
    dit = distances.begin();
    sit = surfaces.begin();
    done = true;
    for (; dit != distances.end(); dit++, sit++) {
      if (dit+1 != distances.end() && *dit > *(dit+1)) {
        double tmp_dist = *dit;
        *dit = *(dit+1);
        *(dit+1) = tmp_dist;
        EntityHandle tmp_hand = *sit;
        *sit = *(sit+1);
        *(sit+1) = tmp_hand;
        done = false;
      }
    }
  }

  if (!distances.empty()) len = distances[0];
  
  return MB_SUCCESS;
#endif
#else
  return MB_FAILURE;
#endif
}

ErrorCode DagMC::build_indices(Range &surfs, Range &vols)
{
  ErrorCode rval = MB_SUCCESS;

    // surf/vol offsets are just first handles
  setOffset = (*surfs.begin() < *vols.begin() ? *surfs.begin() : *vols.begin());
  setOffset = (impl_compl_handle < setOffset ? impl_compl_handle : setOffset);
    // max
  EntityHandle tmp_offset = (surfs.back() > vols.back() ? 
                               surfs.back() : vols.back());
  tmp_offset = (impl_compl_handle > tmp_offset ? 
                impl_compl_handle : tmp_offset);
    // set size
  rootSets.resize(tmp_offset - setOffset + 1);
  entIndices.resize(rootSets.size());

    // store surf/vol handles lists (surf/vol by index) and
    // index by handle lists
  surf_handles().resize( surfs.size() + 1 );
  std::vector<EntityHandle>::iterator iter = surf_handles().begin();
  *(iter++) = 0;
  std::copy( surfs.begin(), surfs.end(), iter );
  int idx = 1;
  for (Range::iterator rit = surfs.begin(); rit != surfs.end(); rit++)
    entIndices[*rit-setOffset] = idx++;
  
  vol_handles().resize( vols.size() + 1 );
  iter = vol_handles().begin();
  *(iter++) = 0;
  std::copy( vols.begin(), vols.end(), iter );
  vol_handles().push_back(impl_compl_handle);
  idx = 1;
  int max_id = -1;
  for (Range::iterator rit = vols.begin(); rit != vols.end(); rit++)    {
    entIndices[*rit-setOffset] = idx++;
    int result=0;
    MBI->tag_get_data( idTag, &*rit, 1, &result );
    max_id = (max_id > result ? max_id : result);
  }
    // add implicit complement to entity index
  entIndices[impl_compl_handle-setOffset] = idx++ ;

    // assign ID to implicit complement
  max_id++;
  MBI->tag_set_data(idTag, &impl_compl_handle, 1, &max_id);
  

 
#ifdef CGM
  if ( have_cgm_geom ) {
    // TODO: this block should only execute if the user has explicitly requested useCAD for ray firing.
    // however, this function curently executes before we know if useCAD will be specified, so do it every time.

    geomEntities.resize(rootSets.size());
      // get geometry entities by id and cache in this vector
    std::vector<int> ids;
    ids.resize(surfs.size());
    rval = MBI->tag_get_data(id_tag(), surfs, &ids[0]);
    if (MB_SUCCESS != rval) return MB_FAILURE;
    int i = 0;
    Range::iterator rit = surfs.begin();
    for (; rit != surfs.end(); rit++, i++) {
      RefEntity *this_surf = GeometryQueryTool::instance()->
        get_ref_face(ids[i]);
      assert(NULL != this_surf);
      geomEntities[*rit - setOffset] = this_surf;
    }
    ids.resize(vols.size());
    rval = MBI->tag_get_data(id_tag(), vols, &ids[0]);
    if (MB_SUCCESS != rval) return MB_FAILURE;
    i = 0;
    rit = vols.begin();
    for (; rit != vols.end(); rit++, i++) {
      RefEntity *this_vol = GeometryQueryTool::instance()->
        get_ref_volume(ids[i]);
      assert(NULL != this_vol);
      geomEntities[*rit - setOffset] = this_vol;
    }
  }
#endif  

    // get group handles
  Tag category_tag = get_tag(CATEGORY_TAG_NAME, CATEGORY_TAG_LENGTH, 
                               MB_TAG_SPARSE, MB_TYPE_OPAQUE);
  char group_category[CATEGORY_TAG_SIZE];
  std::fill(group_category, group_category+CATEGORY_TAG_SIZE, '\0');
  sprintf(group_category, "%s", "Group");
  const void* const group_val[] = {&group_category};
  Range groups;
  rval = MBI->get_entities_by_type_and_tag(0, MBENTITYSET, &category_tag, 
					   group_val, 1, groups);
  if (MB_SUCCESS != rval)
    return rval;
  group_handles().resize(groups.size()+1);
  group_handles()[0] = 0;
  std::copy(groups.begin(), groups.end(), &group_handles()[1]);

    // populate root sets vector
  std::vector<EntityHandle> rsets;
  rsets.resize(surfs.size());
  rval = MBI->tag_get_data(obb_tag(), surfs, &rsets[0]);
  if (MB_SUCCESS != rval) return MB_FAILURE;
  Range::iterator rit;
  int i;
  for (i = 0, rit = surfs.begin(); rit != surfs.end(); rit++, i++)
    rootSets[*rit-setOffset] = rsets[i];

  rsets.resize(vols.size());
  rval = MBI->tag_get_data(obb_tag(), vols, &rsets[0]);
  if (MB_SUCCESS != rval) return MB_FAILURE;
  for (i = 0, rit = vols.begin(); rit != vols.end(); rit++, i++)
    rootSets[*rit-setOffset] = rsets[i];

  // add implicit complement root set
  rsets.resize(1);
  rval = MBI->tag_get_data(obb_tag(), &impl_compl_handle, 1, 
                                       &rsets[0]);
  if (MB_SUCCESS != rval) return MB_FAILURE;
  rootSets[impl_compl_handle-setOffset] = rsets[0];

  return MB_SUCCESS;
}

ErrorCode DagMC::build_obbs(Range &surfs, Range &vols) 
{
  ErrorCode rval = MB_SUCCESS;
  
  for (Range::iterator i = surfs.begin(); i != surfs.end(); ++i) {
    EntityHandle root;
    Range tris;
    rval = MBI->get_entities_by_dimension( *i, 2, tris );
    if (MB_SUCCESS != rval) 
      return rval;
    if (tris.empty()) 
      std::cerr << "WARNING: Surface " << get_entity_id(*i) << " has no facets." << std::endl;
    rval = obbTree.build( tris, root );
    if (MB_SUCCESS != rval) 
      return rval;
    rval = MBI->add_entities( root, &*i, 1 );
    if (MB_SUCCESS != rval)
      return rval;
    rval = MBI->tag_set_data( obbTag, &*i, 1, &root );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  for (Range::iterator i = vols.begin(); i != vols.end(); ++i) {
      // get all surfaces in volume
    Range tmp_surfs;
    rval = MBI->get_child_meshsets( *i, tmp_surfs );
    if (MB_SUCCESS != rval)
      return rval;
    
      // get OBB trees for each surface
    EntityHandle root;
    Range trees;
    for (Range::iterator j = tmp_surfs.begin();  j != tmp_surfs.end(); ++j) {
        // skip any surfaces that are non-manifold in the volume
        // because point containment code will get confused by them
      int sense;
      rval = surface_sense( *i, *j, sense );
      if (MB_SUCCESS != rval) {
        std::cerr << "Surface/Volume sense data missing." << std::endl;
        return rval;
      }
      if (!sense)
        continue;
      
      rval = MBI->tag_get_data( obbTag, &*j, 1, &root );
      if (MB_SUCCESS != rval || !root)  
        return MB_FAILURE;
      trees.insert( root );
    }
    
      // build OBB tree for volume
    rval = obbTree.join_trees( trees, root );
    if (MB_SUCCESS != rval)
      return rval;
    
    rval = MBI->tag_set_data( obbTag, &*i, 1, &root );
    if (MB_SUCCESS != rval)
      return rval;

  }

  rval = build_obb_impl_compl(surfs);
  if (MB_SUCCESS != rval) {
    std::cerr << "Unable to build OBB tree for implicit complement." << std::endl;
    return rval;
  }

  return MB_SUCCESS;
}

ErrorCode DagMC::build_obb_impl_compl(Range &surfs) 
{
  EntityHandle comp_root, surf_obb_root;
  Range comp_tree;
  ErrorCode rval;
  std::vector<EntityHandle> parent_vols;
  
    // search through all surfaces  
  for (Range::iterator surf_i = surfs.begin(); surf_i != surfs.end(); ++surf_i) {
    
    parent_vols.clear();
      // get parents of each surface
    rval = MBI->get_parent_meshsets( *surf_i, parent_vols );
    if (MB_SUCCESS != rval)
      return rval;

      // if only one parent, get the OBB root for this surface
    if (parent_vols.size() == 1 ) {
      rval = MBI->tag_get_data( obbTag, &*surf_i, 1, &surf_obb_root );
      if (MB_SUCCESS != rval)
        return rval;
      if (!surf_obb_root)
        return MB_FAILURE;
      
        // add obb root to list of obb roots
      comp_tree.insert( surf_obb_root );

      // add this surf to the topology of the implicit complement volume
      rval = MBI->add_parent_child(impl_compl_handle,*surf_i);
      if (MB_SUCCESS != rval)
	return rval;
    }  
  }
  
    // join surface trees to make OBB tree for implicit complement
  rval = obbTree.join_trees( comp_tree, comp_root );
  if (MB_SUCCESS != rval)
    return rval;
  
    // tag the implicit complement handle with the handle for its own OBB tree
  rval = MBI->tag_set_data( obbTag, &impl_compl_handle, 1, &comp_root );
  if (MB_SUCCESS != rval)
    return rval;
  
  return MB_SUCCESS;

}

ErrorCode DagMC::getobb(EntityHandle volume, double minPt[3],
                          double maxPt[3])
{
  double center[3], axis1[3], axis2[3], axis3[3];
 
    // get center point and vectors to OBB faces
  ErrorCode rval = getobb(volume, center, axis1, axis2, axis3);
  if (MB_SUCCESS != rval)
    return rval;
     
    // compute min and max verticies
  for (int i=0; i<3; i++) 
  {
    double sum = fabs(axis1[i]) + fabs(axis2[i]) + fabs(axis3[i]);
    minPt[i] = center[i] - sum;
    maxPt[i] = center[i] + sum;
  }
  return MB_SUCCESS;
}

ErrorCode DagMC::getobb(EntityHandle volume, double center[3], double axis1[3],
                          double axis2[3], double axis3[3])
{
    //find EntityHandle node_set for use in box
  EntityHandle root = rootSets[volume - setOffset];
  
    // call box to get center and vectors to faces
  return obbTree.box(root, center, axis1, axis2, axis3);
  
}

ErrorCode DagMC::get_impl_compl()
{
  Range entities;
  const void* const tagdata[] = {implComplName};
  ErrorCode rval = mbImpl->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                           &nameTag, tagdata, 1,
                                                           entities );
  // query error
  if (MB_SUCCESS != rval) {
    std::cerr << "Unable to query for implicit complement." << std::endl;
    return rval;
  }

  // found too many
  if (entities.size() > 1) {
    std::cerr << "Too many implicit complement sets." << std::endl;
    return MB_MULTIPLE_ENTITIES_FOUND;
  }

  // found none
  if (entities.empty()) {
    rval= MBI->create_meshset(MESHSET_SET,impl_compl_handle);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to create mesh set for implicit complement." << std::endl;
      return rval;
    }
      // tag this entity with name for implicit complement
    rval = MBI->tag_set_data(nameTag,&impl_compl_handle,1,&implComplName);
    if (MB_SUCCESS != rval) {
      std::cerr << "Failed to tag new entity as implicit complement." << std::endl;
    }

    return rval;

  } else {
    // found a single implicit complement
    impl_compl_handle = entities.front();
    return MB_SUCCESS;
  }
  
  
}

} // namespace moab

