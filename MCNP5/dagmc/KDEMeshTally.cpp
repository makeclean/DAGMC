// KDEMeshTally.cpp

#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "moab/AdaptiveKDTree.hpp"
#include "moab/CartVect.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/Types.hpp"

#include "KDEKernel.hpp"
#include "KDEMeshTally.hpp"
#include "KDENeighborhood.hpp"
#include "TallyEvent.hpp"

// initialize static variables
moab::CartVect KDEMeshTally::default_bandwidth( 1, 1, 1 );
bool KDEMeshTally::seed_is_set = false;

// quadrature points and weights for the integrate_path_kernel function
const double quad_points[4] = {0.339981043585, -0.339981043585,
                               0.861136311594, -0.861136311594};

const double quad_weights[4] = {0.652145154863, 0.652145154863,
                                0.347854845137, 0.347854845137};

//-----------------------------------------------------------------------------
static double parse_bandwidth_param( const std::string& key,
                                     const std::string& val, 
                                     double default_value = 1 )
{
  char* end;
  double ret = strtod( val.c_str(), &end );
  if( val.c_str() == end || ret <= 0 ){
    // parsing failed or value is invalid
    std::cerr << "Error parsing bandwidth param " << key << ": '" << val << "' is invalid" << std::endl;
    std::cerr << "      Using default value " << key << " = " << default_value << std::endl;
    ret = default_value;
  }
  return ret;
}
//-----------------------------------------------------------------------------
static KDEKernel::KernelType set_kernel_type( const std::string& key,
                                              const std::string& val )
{

  KDEKernel::KernelType k;

  if ( val == "epanechnikov" || val == "e" )
    k = KDEKernel::EPANECHNIKOV;
  else if ( val == "uniform" || val == "u" )
    k = KDEKernel::UNIFORM;
  else {

    k = KDEKernel::EPANECHNIKOV;
    std::cerr << "\nWarning: " << val << " is not a valid kernel function" << std::endl;
        std::cerr << "      Using default " << key << " = epanechnikov\n" << std::endl;

  }

  return k;

}
//-----------------------------------------------------------------------------
KDEMeshTally* KDEMeshTally::setup( const MeshTallyInput& fmesh, 
                                   moab::Interface* mbi, TallyType type )
{

  bool use_dagmc_mesh = false; // true if mesh data should be pulled from DagMC object
  moab::CartVect bandwidth = default_bandwidth;
  KDEKernel::KernelType kernel = KDEKernel::EPANECHNIKOV;
  unsigned int subtracks = 0;

  const MeshTallyInput::TallyOptions& params = fmesh.options;
  int id = fmesh.tally_id;

  for( MeshTallyInput::TallyOptions::const_iterator i = params.begin(); i != params.end(); ++i )
  {
    std::string key = (*i).first, val = (*i).second;
    if( key == "out" ) continue;  // processed in MeshTally
    else if( key == "hx" ) bandwidth[0] = parse_bandwidth_param( key, val );
    else if( key == "hy" ) bandwidth[1] = parse_bandwidth_param( key, val );
    else if( key == "hz" ) bandwidth[2] = parse_bandwidth_param( key, val );
    else if( key == "kernel" ) kernel = set_kernel_type( key, val );
    else if ( key == "seed" && type == SUB_TRACK)
    {
        // override random number seed if requested by user
        int seed = atoi(val.c_str());
        srand(seed);
        seed_is_set = true;
    }
    else if( key == "subtracks" && type != COLLISION ) { 
      subtracks = atoi( val.c_str() );
      if( subtracks == 0 ) {
        std::cerr << "\nWarning: number of subtracks requested is invalid" << std::endl;
        std::cerr << "      Using default value " << key << " = 3\n" << std::endl;
        subtracks = 3;
      }
    }
    else{
      std::cerr << "Warning: KDE tally's FC" << id << " card has unknown key '" << key << "'" << std::endl;
    }
  }

  // set random number seed if it has not already been set by another instance
  if (type == SUB_TRACK && !seed_is_set)
  {
    srand(time(NULL));
    seed_is_set = true;
  }

  moab::EntityHandle moab_set = NULL ; // TODO: this should be queried from DagMC
  if( !use_dagmc_mesh ){
   
    moab::ErrorCode rval;
    rval = mbi->create_meshset( moab::MESHSET_SET, moab_set );
    assert( rval == moab::MB_SUCCESS );

    rval = mbi->load_file( fmesh.input_filename.c_str(), &moab_set );

    if( rval != moab::MB_SUCCESS ){
      std::cerr << "Error: could not load KDE tally mesh file " << fmesh.input_filename << std::endl;
      exit( EXIT_FAILURE );
    }

  }

  KDEMeshTally *kde = new KDEMeshTally( fmesh, mbi, moab_set, bandwidth,
                                        type, kernel, subtracks );

  std::cout << "Creating KDE ";

  if ( type == COLLISION ) std::cout << "collision ";
  else if ( type == INTEGRAL_TRACK ) std::cout << "track ";
  else std::cout << "subtrack ";
 
  std::cout << "fmesh" << id 
            << ", input: " << (use_dagmc_mesh ? "(pre-loaded DagMC data)" : fmesh.input_filename.c_str())  
            << ", output: " << kde->output_filename << std::endl;
  
  std::cout << "   using the " << KDEKernel::kernel_names[kernel] << " kernel";
  std::cout << " with bandwidth = " << bandwidth;

  if ( subtracks != 0 )
    std::cout << "\n   and splitting tracks into " << subtracks << " subtracks" << std::endl;
  else
    std::cout << std::endl;

  // create a tally set that contains only the 3D mesh cells (i.e. hexes/tets)
  moab::ErrorCode rval;
  rval = mbi->create_meshset( moab::MESHSET_SET, kde->tally_set );
  assert( rval == moab::MB_SUCCESS );

  moab::Range mesh_cells;
  rval = mbi->get_entities_by_dimension( moab_set, 3, mesh_cells );
  assert( rval == moab::MB_SUCCESS );
   
  rval = mbi->add_entities( kde->tally_set, mesh_cells );
  assert( rval == moab::MB_SUCCESS );

  return kde;

}
//-----------------------------------------------------------------------------
KDEMeshTally::KDEMeshTally( const MeshTallyInput& settings,
                            moab::Interface* moabMesh,
                            moab::EntityHandle moabSet,
                            moab::CartVect bandwidth,
                            KDEMeshTally::TallyType type,
                            KDEKernel::KernelType k,
                            unsigned int subtracks )
: MeshTally(settings), mb( moabMesh ), bandwidth( bandwidth ),
  kde_tally( type ), kernel( new KDEKernel(k) )
{

  build_tree(moabSet); 

  // initialize running variance variables
  max_collisions = false;
  numCollisions = 0;
  Mn = moab::CartVect( 0, 0, 0 );
  Sn = moab::CartVect( 0, 0, 0 );

  // check numSubtracks is a valid parameter for SUBTRACK tallies
  if (subtracks == 0 && kde_tally == SUB_TRACK)
  {

    std::cerr << "\nWarning: number of subtracks must be non-zero for KDE subtrack tallies" << std::endl;
    std::cerr << "      Using default value subtracks = 3\n" << std::endl;
    num_subtracks = 3;

  }
  else
    num_subtracks = subtracks;
}
//-----------------------------------------------------------------------------
KDEMeshTally::~KDEMeshTally()
{

  delete tree;
  delete kernel;
  
}
//-----------------------------------------------------------------------------
void KDEMeshTally::compute_score(const TallyEvent& event, int ebin)
{
    // initialize common weighting factor for this tally event
    double weight = event.get_weighting_factor();

    // set up tally event based on KDE mesh tally type
    TrackData track;
    CollisionData collision;
    std::vector<moab::CartVect> subtrack_points;
    bool event_is_set = false;

    if (kde_tally == INTEGRAL_TRACK || kde_tally == SUB_TRACK)
    {
        event_is_set = event.get_track_data(track);

        if (kde_tally == SUB_TRACK && event_is_set == true)
        {
            // multiply weight by track length and set up sub-track points
            weight *= track.track_length;
            subtrack_points = choose_points(track, num_subtracks);
        }
    }
    else // kde_tally == COLLISION
    {
        event_is_set = event.get_collision_data(collision);

        if (event_is_set)
        {
            // divide weight by cross section and update optimal bandwidth
            weight /= collision.total_cross_section;
            update_variance(collision.collision_point);
        }
    }

    // check that a valid tally event has been set for this KDE mesh tally
    if (!event_is_set)
    {
        std::cerr << "\nError: Tally event is not valid for KDE mesh tally ";
        std::cerr << input_data.tally_id << std::endl;
        exit(EXIT_FAILURE);
    }

    // create the neighborhood region and find all of the calculations points
    KDENeighborhood region(event, bandwidth, *tree, tree_root);

    std::vector<moab::EntityHandle> calculation_points;
    moab::ErrorCode rval = region.get_points(calculation_points);

    assert(moab::MB_SUCCESS == rval);

    // iterate through the calculation points
    std::vector<moab::EntityHandle>::iterator i;
    double coords[3];

    for (i = calculation_points.begin(); i != calculation_points.end(); ++i)
    {
        // get coordinates of this point
        moab::EntityHandle point = *i;
        rval = mb->get_coords(&point, 1, coords);

        assert(moab::MB_SUCCESS == rval);

        // compute the final contribution to the tally for this point
        double score = weight;

        if (kde_tally == INTEGRAL_TRACK)
        {
            score *= integral_track_score(track, moab::CartVect(coords));
        }
        else if (kde_tally == SUB_TRACK)
        {
            score *= subtrack_score(subtrack_points, moab::CartVect(coords));
        }
        else // kde_tally == COLLISION
        {
            score *= collision_score(collision, moab::CartVect(coords));
        }

        // add score to KDE mesh tally for the current history
        add_score_to_tally(point, score, ebin);
    }
}
//-----------------------------------------------------------------------------
void KDEMeshTally::end_history()
{
  
  std::set<moab::EntityHandle>::iterator i;
 
  // add results from current history to the tally for each calculation point
  for ( i = visited_this_history.begin() ; i != visited_this_history.end() ; ++i ) {
    
    for( unsigned int j = 0; j < num_energy_bins; ++j ){
      double& history = get_data( temp_tally_data, *i, j );
      double& tally =   get_data( tally_data, *i, j );
      double& error =   get_data( error_data, *i, j );
      
      tally += history;
      error += ( history * history );
      
      // reset temp_tally for the next particle history
      history = 0;
      
    }
  }
  visited_this_history.clear();

}
//-----------------------------------------------------------------------------
void KDEMeshTally::print( double sp_norm, double fmesh_fact )
{

  // tags tally/error results to the nodes and writes mesh to output file
  write_results( sp_norm, fmesh_fact );

}
//-----------------------------------------------------------------------------
void KDEMeshTally::write_results( double sp_norm, double fmesh_fact )
{

  double tally = 0;
  double error = 0, rel_err = 0;
  
  moab::ErrorCode rval = moab::MB_SUCCESS;

  // print the optimal bandwidth if it was computed
  if ( kde_tally == COLLISION ) {
  
    std::cout << std::endl << "optimal bandwidth for " << numCollisions;
    std::cout  << " collisions is: " << get_optimal_bandwidth() << std::endl;

  }

  // tag tally and relative error results to the mesh for each entity
  moab::Range::iterator i;
  
  for ( i = tally_points.begin() ; i != tally_points.end() ; ++i ) {

    moab::EntityHandle point = *i;

    for ( unsigned int j = 0; j < num_energy_bins; ++ j){

      tally = get_data( tally_data, point, j);
      error = get_data( error_data, point, j );
      
      // compute relative error for the tally
      // Use 0 as the rel_err value if nothing has been computed for this tally point;
      // this reflects MCNP's approach to avoiding a divide-by-zero situation.
      rel_err = 0; 
      if( error != 0 ){
        rel_err = sqrt( error / ( tally * tally ) - 1.0 / sp_norm );
      }
      
      // normalizing mesh tally results by the number of source particles
      tally /= sp_norm;
      
      // applying the fmesh multiplication FACTOR to the mesh tally results
      tally *= fmesh_fact;
      
      // set tally and error tag values for this entity
      rval = mb->tag_set_data( tally_tags[j], &point, 1, &tally );
      assert( moab::MB_SUCCESS == rval );
      
      rval = mb->tag_set_data( error_tags[j], &point, 1, &rel_err );
      assert( moab::MB_SUCCESS == rval ); 
      
    } 
    
  }

  // create a global tag to store the bandwidth value
  moab::Tag bandwidth_tag;
  rval = mb->tag_get_handle( "BANDWIDTH_TAG", 3, moab::MB_TYPE_DOUBLE, bandwidth_tag,
                             moab::MB_TAG_MESH|moab::MB_TAG_CREAT );
  assert( moab::MB_SUCCESS == rval );

  // add bandwidth tag to the root set
  moab::EntityHandle bandwidth_set = mb->get_root_set();
  rval = mb->tag_set_data( bandwidth_tag, &bandwidth_set, 1, &bandwidth );
  assert( moab::MB_SUCCESS == rval );

  // define list of tags to include and write mesh to output file
  std::vector<moab::Tag> output_tags = tally_tags;
  output_tags.insert( output_tags.end(), error_tags.begin(), error_tags.end() );
  output_tags.push_back( bandwidth_tag );  

  rval = mb->write_file( output_filename.c_str(),
                         NULL, NULL, &tally_set, 1, &(output_tags[0]), output_tags.size() );
  assert( moab::MB_SUCCESS == rval );
  
}
//-----------------------------------------------------------------------------
void KDEMeshTally::build_tree( moab::EntityHandle meshset )
{

  // Obtain all of the calculation points in the mesh and store into tally_points
  moab::EntityType type = moab::MBVERTEX;
  
  moab::ErrorCode rval = mb->get_entities_by_type( meshset, type, tally_points );
  assert( moab::MB_SUCCESS == rval );  

  resize_data_arrays( tally_points.size() );
    
  // Measure the number of divisions in the moab::Range used to represent the tally points
  // If there are many divisions (for some rather arbitrary definition of "many"), print
  // a warning about performance compromise
  int psize = tally_points.psize();
  std::cout << "   Tally range has psize: " << psize << std::endl;
  if( psize > 4 ){
    std::cerr << "Warning: large tally range psize " << psize 
              << ", may reduce performance." << std::endl;
  }
  // Build a KD-tree using all of the calculation points in the mesh
  moab::AdaptiveKDTree::Settings settings;
  settings.maxEntPerLeaf = 10;
  settings.maxTreeDepth = 30;

  tree = new moab::AdaptiveKDTree( mb );
  rval = tree->build_tree( tally_points, tree_root, &settings );
  assert( moab::MB_SUCCESS == rval );  

  rval = setup_tags( mb, "KDE_" );

}
//-----------------------------------------------------------------------------
void KDEMeshTally::update_variance(const moab::CartVect& collision_point)
{
 
  if ( numCollisions != LLONG_MAX ) {

    ++numCollisions;
    
    // obtain previous value for the mean
    moab::CartVect Mn_prev = Mn;
  
    // compute new values for the mean and variance
    if ( numCollisions == 1 )
      Mn = collision_point;
    else {

      Mn += (collision_point - Mn_prev) / numCollisions;
    
      for ( int i = 0 ; i < 3 ; ++i )
        Sn[i] += (collision_point[i] - Mn_prev[i]) * (collision_point[i] - Mn[i]);
    
    }

  }
  else if ( !max_collisions ) {
  
    std::cerr << "/nWarning: number of collisions exceeds maximum/n"
              << "  optimal bandwidth will be based on " << numCollisions
              << " collisions./n/n";

    max_collisions = true;

  }

}
//-----------------------------------------------------------------------------
moab::CartVect KDEMeshTally::get_optimal_bandwidth()
{
  
  double stdev = 0;
  moab::CartVect optimal_bandwidth;
  
  for ( int i = 0 ; i < 3 ; ++i ) {

    stdev = sqrt( Sn[i] / ( numCollisions - 1 ) );
    optimal_bandwidth[i] = 0.968625 * stdev * pow( numCollisions, -1.0/7.0 );

  }

  return optimal_bandwidth;

}
//-----------------------------------------------------------------------------
void KDEMeshTally::add_score_to_tally( moab::EntityHandle mesh_point,
                                       double score,
                                       int ebin )
{

  get_data( temp_tally_data, mesh_point, ebin ) += score;

  // tally the total energy bin if requested
  if ( input_data.total_energy_bin )
    get_data( temp_tally_data, mesh_point, (num_energy_bins-1) ) += score;

  visited_this_history.insert( mesh_point );

}
//-----------------------------------------------------------------------------
// NOTE: integral_track_estimator uses the 4-point gaussian quadrature method
double KDEMeshTally::integral_track_score(const TrackData& data,
                                          const moab::CartVect& X)
{
    // determine the limits of integration
    std::pair<double, double> limits;
    
    bool valid_limits = set_integral_limits(data, X, limits);

    // compute value of the integral only if valid limits exist
    if (valid_limits)
    {
        // define scaling constants
        double c1 = 0.5 * (limits.second - limits.first);
        double c2 = 0.5 * (limits.second + limits.first);

        // sum contributions for all quadrature points
        double sum = 0;

        for (int i = 0; i < 4; ++i)
        {
            // define scaled quadrature point
            double s = c1 * quad_points[i] + c2;

            // compute the value of the kernel function K(X, s)
            double kernel_value = 1;

            for (int j = 0; j < 3; ++j)
            {
                double u = X[j] - data.start_point[j] - s * data.direction[j];
                u /= bandwidth[j];
                kernel_value *= kernel->evaluate(u) / bandwidth[j];
            }
        
            // multiply by quadrature weight and add to sum
            sum += quad_weights[i] * kernel_value;
        }

        // return value of the integral
        return c1 * sum;
    }
    else
    {
        // integration limits are not valid so no score is computed
        return 0;
    }
}
//-----------------------------------------------------------------------------
bool KDEMeshTally::set_integral_limits(const TrackData& data,
                                       const moab::CartVect& X,
                                       std::pair<double, double>& limits)
{
    bool valid_limits = false;

    // set initial integral limits to the full track length (default values)
    limits = std::make_pair(0, data.track_length);

    // check limits against the valid path length interval for each dimension
    for (int i = 0; i < 3; ++i)
    {
        double path_min = limits.first;
        double path_max = limits.second;

        // compute valid path length interval Si = [path_min, path_max]
        if (data.direction[i] > 0)
        {
            path_min = X[i] - data.start_point[i] - bandwidth[i];
            path_min /= data.direction[i];

            path_max = X[i] - data.start_point[i] + bandwidth[i];
            path_max /= data.direction[i];
        }
        else if (data.direction[i] < 0)
        {
            path_min = X[i] - data.start_point[i] + bandwidth[i];
            path_min /= data.direction[i];

            path_max = X[i] - data.start_point[i] - bandwidth[i];
            path_max /= data.direction[i];
        }

        // set lower limit to highest minimum
        if (path_min > limits.first)
        {
            limits.first = path_min;
        }

        // set upper limit to lowest maximum
        if (path_max < limits.second)
        {
            limits.second = path_max;
        }
    }
  
    // limits are only valid if upper limit is greater than lower limit
    if (limits.first < limits.second)
    {
        valid_limits = true;
    }

    return valid_limits;
}
//-----------------------------------------------------------------------------
double KDEMeshTally::subtrack_score(const std::vector<moab::CartVect>& points,
                                    const moab::CartVect& X)
{
    // iterate through the sub-track points
    std::vector<moab::CartVect>::const_iterator i;
    double score = 0;

    for (i = points.begin(); i != points.end(); ++i)
    {
        // Compute the value of the kernel function K(X)
        double kernel_value = 1;

        for (int j = 0; j < 3; ++j)
        {
            double u = (X[j] - (*i)[j]) / bandwidth[j];
            kernel_value *= kernel->evaluate(u) / bandwidth[j];
        }

        // add kernel contribution for sub-track point to sum
        score += kernel_value;
    }

    // normalize by the total number of sub-track points
    score /= points.size();

    return score;
}
//-----------------------------------------------------------------------------
std::vector<moab::CartVect> KDEMeshTally::choose_points(const TrackData& data,
                                                        int p)
{
    // make sure the number of sub-tracks is valid
    assert(p > 0);

    // compute sub-track length, assumed to be equal for all sub-tracks
    double sub_track_length = data.track_length / p;

    // set the starting point to the beginning of the track segment
    moab::CartVect start_point = data.start_point;

    // choose a random position along each sub-track
    std::vector<moab::CartVect> random_points;

    for (int i = 0; i < p; ++i)
    {
        double path_length = rand() * sub_track_length / RAND_MAX;
        
        // add the coordinates of the corresponding point
        random_points.push_back(start_point + path_length * data.direction);

        // shift starting point to the next sub-track
        start_point += sub_track_length * data.direction;
    }
 
    return random_points;
}
//-----------------------------------------------------------------------------
double KDEMeshTally::collision_score(const CollisionData& data,
                                     const moab::CartVect& X)
{
    // compute the value of the kernel function K(X)
    double score = 1;

    for (int i = 0; i < 3; ++i)
    {
        double u = (X[i] - data.collision_point[i]) / bandwidth[i];
        score *= kernel->evaluate(u) / bandwidth[i];
    }

    return score;
}
//-----------------------------------------------------------------------------
