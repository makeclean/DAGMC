#include "ray_urchin.hpp"
#include "moab/ProgOptions.hpp"
#include <fstream>

// constructor
Ray_Urchin::Ray_Urchin(const std::string cad_file,
                       const std::string ray_file,
                       const moab::CartVect start_pt_in,
                       const int vol_id) {
  
  dagmc = new moab::DagMC();
  dmd = new dagmcMetaData(dagmc);
  
  dagmc_file = cad_file;
  urchin_file = ray_file;
  start_position = start_pt_in;
  start_vol_id = vol_id;
}

// destructor
Ray_Urchin::~Ray_Urchin() {
  delete dagmc;
  delete dmd;
}

// setup the DAGMC part of the input
moab::ErrorCode Ray_Urchin::Setup() {
  moab::ErrorCode ret = moab::MB_FAILURE;

  // load the file
  ret = dagmc->load_file( dagmc_file.c_str() );
  MB_CHK_SET_ERR_CONT(ret,"Failed to load file " + dagmc_file);
  
  // initalise dagmc
  ret = dagmc->init_OBBTree();
  MB_CHK_SET_ERR_CONT(ret, "Failed to build OBB tree and/or indices.");
  
  // parse all the properties
  dmd->load_property_data();

  // read rayfile
  ret = rays_load_and_init(ray_nums, ray_list);
  MB_CHK_ERR_CONT(ret);

  // get the handle of the start volume
  if (start_vol_id != -1 ) {
    start_vol_handle = dagmc->entity_by_id(3,start_vol_id);
  } else {
    // given the start point 
    ret = find_start_vol();
    MB_CHK_SET_ERR_CONT(ret,"Can't find starting volume for provided starting point.");
  }
  // find graveyard handle
  find_graveyard();
  return moab::MB_SUCCESS;
}

// function to open a text file and read a set of unit vectors
moab::ErrorCode Ray_Urchin::rays_load_and_init( int (&ray_nums)[3], std::vector< moab::CartVect > &ray_list)
{
  std::ifstream ray_file_stream;
  ray_file_stream.open(urchin_file);

  if (!ray_file_stream)
    MB_SET_ERR(moab::MB_FAILURE,"Unable to open ray file " + urchin_file);

  // read the total number of rays as well as the longitudinal and lattitudinal numbers
  ray_file_stream >> ray_nums[0] >> ray_nums[1] >> ray_nums[2];

  // populate with the up/down unit vectors
  ray_list.push_back(moab::CartVect(0,0,1));
  ray_list.push_back(moab::CartVect(0,0,-1));
  
  // read unit vectors until the file ends
  while (!ray_file_stream.eof())
    {
      double u, v, w;
      ray_file_stream >> u >> v >> w;
      ray_list.push_back(moab::CartVect(u,v,w));
    }

  ray_file_stream.close();

  return moab::MB_SUCCESS;
}

// function to find the volume that contains the starting point,
//    first checking a provided volume ID
// input:
//     DagMC& dagmc - reference to DagMC instance
//     CartVect start_p - location of starting point
//     int start_vol_id - volume ID of possible starting volume
// output:
//     EntityHandle& start_vol_handle - entity handle of volume with starting point
moab::ErrorCode Ray_Urchin::find_start_vol() {

  moab::ErrorCode ret;
  int inside = 0;

  // if vol id already set
  if (start_vol_id > 0 ) {
    ret = dagmc->point_in_volume(start_vol_handle,start_position.array(),inside);
    MB_CHK_ERR(ret);
    return moab::MB_SUCCESS;
  }
  // otherwise search through them all
  if(!inside) MB_SET_ERR_CONT("Starting point is not in provided starting volume.");
  moab::EntityHandle vol_handle;
  int i = 0;
  int num_vols = dagmc->num_entities(3);
  for ( i = 1 ; i < num_vols ; i++) {
    vol_handle = dagmc->entity_by_index(3,i);
    ret = dagmc->point_in_volume(vol_handle,start_position.array(), inside);
    MB_CHK_ERR(ret);
    if(inside) break;
  }
  // if point was inside 
  if(inside) {
    start_vol_id = dagmc->id_by_index(3,i);
    start_vol_handle = vol_handle;  
    return moab::MB_SUCCESS;
  } else { // point outside all volumes
    std::cout << "Start coordinate not found to be inside any volume" << std::endl;
    return moab::MB_FAILURE;
  }
}

// function to find the graveyard volume
// modifies class memeber variable to set graveyard
// entity dmd will fail if no graveyard found
void Ray_Urchin::find_graveyard() {
  int i = 0;
  int num_vols = dagmc->num_entities(3);
  for ( i = 1 ; i < num_vols ; i++) {
    moab::EntityHandle vol_handle = dagmc->entity_by_index(3,i);
    std::string mat_name = dmd->volume_material_property_data_eh[vol_handle];
    if(mat_name.find("Graveyard") != std::string::npos) {
      graveyard_handle = vol_handle;
      return;
    }
  }
  return;
}

// function to repeatedly fire a ray from a point to the boundary and
//    accumulate a list of pairs for each volume traversed,
//    with each pair containing the volume handle and the distance traversed
// input:
//    CartVect start_pt - the starting point of the rays
//    EntityHandle start_vol_handle - the handle for the starting volume
//    CartVect dir - the unit vector for the current ray
//    EntityHandle graveyard_handle - handle for the graveyard volume
// output:
//    std::set< EntityHandle >& include_vols - list of volume entities touched so far
//    (return) Ray_History& - ray history from start_pt to boundary
Ray_History Ray_Urchin::get_ray_hist(const moab::CartVect dir,
                         std::set< moab::EntityHandle >& include_vols){

  Ray_History ray_hist;
  moab::ErrorCode ret = moab::MB_FAILURE;
  moab::EntityHandle vol = start_vol_handle;
  double dist = 0.0;
  moab::EntityHandle new_vol = vol;
  moab::EntityHandle surf = 0;
  
  //   fire ray to get all intersections to graveyard : {dir => {cell => track length}}
  while (vol != graveyard_handle){
    ret = dagmc->ray_fire( vol, start_position.array(), dir.array(), surf, dist);
    MB_CHK_ERR_CONT(ret);
    ray_hist.push_back(std::make_pair(vol,dist));
    include_vols.insert(vol);

    ret = dagmc->next_vol(surf,vol,new_vol);
    MB_CHK_ERR_CONT(ret);
    vol = new_vol;
  }
  
  // pop all the void regions off the end of the list
  while (dmd->volume_material_property_data_eh[ray_hist.end()->first].find("Vacuum") != std::string::npos) {
    include_vols.erase(ray_hist.end()->first);
    ray_hist.pop_back();
  }
  
  return ray_hist;
}

// loop over the rays read in and process each ray
moab::ErrorCode Ray_Urchin::process_rays() {
  // for each ray
  for (std::vector< moab::CartVect >::iterator dir = ray_list.begin();
       dir != ray_list.end();
       dir++){
    int count = dir - ray_list.begin();
    // it makes sense to use cartvec as key in map
    // but result in compile error 
    // cartvect would need to define a < operator
    ray_hist[count] = get_ray_hist(*dir, include_vols);  
  }
  return moab::MB_SUCCESS;
}

// function to write out all ray data in the HZETRN2015 format
// input:
//    int[3]& ray_nums - total number of rays and lat/long rays
//    EntityHandle start_vol_handle - handle for starting volume
//    std::set< EntityHandle >& include_vols - superset of all volumes touched by alll rays
//    std::vector< CartVect >& ray_list - list of ray unit directions
//    std::map< CartVect, Ray_History >& ray_hist - map of directions to ray histories
void Ray_Urchin::write_hzetrn2015() {

  // - write header
  std::cout << "***GEOMETRY_DEFN***" << std::endl;
  std::cout << ray_nums[0] << "\t" << ray_nums[1] << "\t" << ray_nums[2] << std::endl << std::endl;
  
  // - write 'raytrace' flag
  std::cout << "raytrace" << std::endl << std::endl;
  
  // - write start point & start cell
  std::cout << start_position[0] << "\t" << start_position[1] << "\t" << start_position[2] << "\t" << start_vol_id << std::endl << std::endl;
  
  // - write cell list
  std::cout << include_vols.size() << std::endl;
  /* originally this said to iterate through the sorted list but std set 
  is already sorted */
  for (std::set< moab::EntityHandle>::iterator vol_id = include_vols.begin();
       vol_id != include_vols.end();
       vol_id++) {
    //   get cell ID, cell name, material
    std::cout << dagmc->index_by_handle(*vol_id) << "\t";
    std::cout << dagmc->get_entity_id(*vol_id) << "\t";
    std::cout << dmd->get_volume_property("material",*vol_id) << std::endl;
  }

  std::cout << std::endl << ray_list.size() << std::endl;
  
  // - write number of rays
  for (std::vector< moab::CartVect >::iterator dir = ray_list.begin();
       dir != ray_list.end();
       dir++) {
    
    
    // get the index of the current ray
    int count = dir - ray_list.begin();

    //    - write ray direction
    moab::CartVect direction = *dir;
    std::cout << direction << "\t" << ray_hist[count].size() << std::endl;
    //    - write track-map for that ray
    for (Ray_History_iter slab = ray_hist[count].begin();
         slab != ray_hist[count].end();
         slab++) {
      std::cout << dagmc->get_entity_id(slab->first) << "\t" << slab->second << std::endl;
    }
  }  
}

// write sector thicknesses
void Ray_Urchin::write_sectors() {

  // - write header
  std::cout << "# " << ray_nums[0] << "\t" << ray_nums[1] << "\t" << ray_nums[2] << std::endl << std::endl;
  
  // - write start point & start cell
  std::cout "# " << start_position[0] << "\t" << start_position[1] << "\t" << start_position[2] << "\t" << start_vol_id << std::endl << std::endl;
  
  
  // - write number of rays
  for (std::vector< moab::CartVect >::iterator dir = ray_list.begin();
       dir != ray_list.end();
       dir++) {
    
    
    // get the index of the current ray
    int count = dir - ray_list.begin();

    //  - write ray direction
    moab::CartVect direction = *dir;
    std::cout << direction << "\t";
    //    - write track-map for that ray
    double total_thickness = 0.0;
    for (Ray_History_iter slab = ray_hist[count].begin();
         slab != ray_hist[count].end();
         slab++) {
      total_thickness += slab->second;
    }
    std::cout << total_thickness << std::endl;
  }  
}

// turn string into vector of double values assuming space delimiting
std::vector<double> string_to_double_vec(std::string input) {
  // turn string into vector of strings
  std::stringstream ss(input);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings(begin, end);
  // 
  std::vector<double> double_vector;
  // for each string push back to array of doubles
  for ( int i = 0 ; i < vstrings.size() ; i++ ) {
    double_vector.push_back(std::stod(vstrings[i]));
  }
  return double_vector;
}

int main(int argc, char* argv[] ){
  // Pseudo-code 
  ProgOptions po("ray_urchin: fire rays from a point outward through a geometry along many rays");
  std::string dagversion;
  moab::DagMC::version( &dagversion );
  po.setVersion( dagversion );

  //   cadfile - filename/path of DAGMC geometry file
  std::string input_file;
  po.addRequiredArg<std::string>( "input_file", "Path to input file for preprocessing", &input_file );

  //   rayfile - filename/path file with ray-directions
  std::string ray_file;
  po.addRequiredArg<std::string>( "ray_file", "Path to input file containing ray directions", &ray_file);

  //   (x,y,z) - starting point
  std::string position;
  // moab does not allow vector of doubles or floats, only ints
  // need to read as string then decompose into vector 
  po.addOpt<std::string>( "starting_point", "(x,y,z) location of ray origins",&position);

  // todo: optional volume ID as input
  int start_vol_id = -1;
  po.addOpt<int>("vols,V", "Volume ID of starting volume", &start_vol_id);

  std::string run = "";
  po.addOpt<std::string>( "output", "Type of output hzetrn, sectors",&run); 
  po.parseCommandLine( argc, argv );

  // get the argument
  po.getOpt("starting_point", &position);
  std::vector<double> starting_point = string_to_double_vec(position);
  if(starting_point.size() != 3 ) {
    std::cout << "Insufficient number of entries for coordinate" << std::endl;
    return 1;
  }

  // start point
  moab::CartVect start_pt(starting_point[0],
                          starting_point[1],
                          starting_point[2]);
  // create urchin to fire rays
  Ray_Urchin *urch = new Ray_Urchin(input_file,ray_file,
                                    start_pt,start_vol_id);

  // read DAGMC file and setup for transport
  moab::ErrorCode ret;
  ret = urch->Setup();
  MB_CHK_SET_ERR_CONT(ret,"Failed to perform setup for Urchin.");

  // process all the rays
  ret = urch->process_rays();
  MB_CHK_SET_ERR_CONT(ret,"Error During ray processing for Urchin.");

  // write out all the ray data
  if(run == "hzetrn")
    urch->write_hzetrn2015();
  else if ( run == "sector" )
    urch->write_sector();
  else {
    std::cout << "Unknown output type" << std::endl;
  }
  delete urch;
  std::cout << "All Finished." << std::endl;
}
