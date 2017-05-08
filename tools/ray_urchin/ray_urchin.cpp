#include "ray_urchin.hpp"

// constructor
Ray_Urchin::Ray_Urchin(const std::string cad_file,
                       const std::string ray_file,
                       const moab::CartVect start_pt_in,
                       const int start_vol_id) {
  
  dagmc = new DagMC();
  dmd = new dagmcMetaData(dagmc);
  
  dagmc_file = cad_file;
  urchin_file = ray_file;
  start_position = start_pt_in;
}

// destructor
Ray_Urchin::~RayUrchin() {
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
  int ray_nums[3];

  std::vector<moab::CartVect> ray_list;
  ret = rays_load_and_init(urchin_file, ray_nums, ray_list);
  MB_CHK_ERR_CONT(ret);

  // get the handle of the start volume
  if (start_vol_id != -1 ) {
    start_vol_handle = dagmc->handle_by_id(3,start_vol_id);
  } else {
    // given the start point 
    ret = find_start_vol();
    MB_CHK_SET_ERR_CONT(ret,"Can't find starting volume for provided starting point.");
  }
  // find graveyard handle
  moab::EntityHandle graveyard_handle = find_graveyard(dagmc);
}

// function to open a text file and read a set of unit vectors
moab::Errorcode Ray_Urchin::rays_load_and_init( int[3] &ray_nums, std::vector< moab::CartVect > &ray_list)
{
  std::ifstream ray_file_stream;
  ray_file_stream.open(urchin_file);

  if (!ray_file_stream)
    MB_SET_ERR("Unable to open ray file " + urchin_file);

  // read the total number of rays as well as the longitudinal and lattitudinal numbers
  ray_file_stream >> ray_nums[0] >> ray_nums[1] >> ray_nums[2];

  // populate with the up/down unit vectors
  ray_list.push_back(CartVect(0,0,1));
  ray_list.push_back(CartVect(0,0,-1));
  
  // read unit vectors until the file ends
  while !(ray_file_stream.eof())
    {
      double u, v, w;
      ray_file_stream >> u >> v >> w;
      ray_list.push_back(CartVect(u,v,w));
    }

  ray_file_stream.close();

  return MB_SUCCESS;
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
  std::vector<moab::EntityHandle> vols = dagmc->volumes();
  int inside = 0;

  // if vol id already set
  if (start_vol_id > 0 ) {
    ret = dagmc->point_in_vol(start_vol_handle,start_position.array(),inside);
    MB_CHK_ERR(ret);
    return moab::MB_SUCCESS;
  }
  // otherwise search through them all
  if(!inside) MB_SET_ERR_CONT("Starting point is not in provided starting volume.");
  moab::EntityHandle vol_handle;
  int i = 0;
  int num_vols = dagmc->num_entities(3);
  for ( i = 1 ; i < num_entities ; i++) {
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
  for ( i = 1 ; i < num_entities ; i++) {
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
Ray_History get_ray_hist(const moab::CartVect start_pt,
                         const moab::CartVect dir,
                         std::set< EntityHandle >& include_vols){

  Ray_History ray_hist;

  moab::EntityHandle vol = start_vol_handle;
  double dist = 0.0;
  moab::EntityHandle new_vol = vol;
  moab::EntityHandle surf = 0;
  
  //   fire ray to get all intersections to graveyard : {dir => {cell => track length}}
  while (vol != graveyard_handle)){
    ret = dagmc->ray_fire( vol, start_poisition.array(), dir.array(), surf, dist);
    MB_CHK_ERR_CONT(ret);
    ray_hist.push_back(std::pair(vol,dist));
    include_vols.insert(vol);

    ret = dagmc->next_vol(surf,vol,new_vol);
    MB_CHK_ERR_CONT(ret);
    vol = new_vol;
  }
  
  // pop all the void regions off the end of the list
  while (dmd->volume_material_property_data_eh[ray_hist.end()->first].find("Vacuum") != std::string::npos) {
    include_vols.remove(ray_hist.end()->first);
    ray_hist.pop();
  }
  
  return ray_hist;
}

// function to write out all ray data in the HZETRN2015 format
// input:
//    int[3]& ray_nums - total number of rays and lat/long rays
//    EntityHandle start_vol_handle - handle for starting volume
//    std::set< EntityHandle >& include_vols - superset of all volumes touched by alll rays
//    std::vector< CartVect >& ray_list - list of ray unit directions
//    std::map< CartVect, Ray_History >& ray_hist - map of directions to ray histories
void Ray_Urchin::write_hzetrn2015(const int[3]& ray_nums,
                  const std::set< moab::EntityHandle >& include_vols,
                  const std::vector< moab::CartVect >& ray_list,
                  const std::map< moab::CartVect, Ray_History >& ray_hist) {

  // - write header
  std::cout << "***GEOMETRY_DEFN***" << std::endl;
  std::cout << ray_nums[0] << "t" << ray_nums[1] << "\t" << ray_nums[2] << std::endl << std::endl;
  
  // - write 'raytrace' flag
  std::cout << "raytrace" << std::endl << std::endl;
  
  // - write start point & start cell
  std::cout << start_position[0] << "\t" << start_position[1] << "\t" << start_position[2] << "\t" << start_vol << std::endl << std::endl;
  
  // - write cell list
  std::cout << include_vols.size() << std::endl;
  
  for (std::set< int >::iterator vol_id = include_vols.sort().begin();
       vol_id != include_vols.sort().end();
       vol_id++) {
    //   get cell ID, cell name, material
    std::cout << *vol_id << "\t" << get_name(*vol_id) << "\t" << get_material(*vol_id) << std::endl;
  }

  std::cout << std::endl << ray_list.size() << std::endl;
  
  // - write number of rays
  for (std::vector< moab::CartVect >::iterator dir = ray_list.begin();
       dir != ray_list.end();
       dir++) {
    //    - write ray direction
    std::cout << *dir << "\t" << *(dir+1) <<  "\t" << *(dir+2) << "\t" << ray_hist[*dir].size() << std::endl;
    //    - write track-map for that ray
    for (Ray_History_iter slab = ray_hist[*dir].begin();
         slab != ray_hist[*dir].end();
         slab++) {
      std::cout << slab->first << "\t" slab->second << std::endl;
    }
  }  
}

int main(int argc, char* argvp[] ){
  // Pseudo-code 
  ProgOptions po("ray_urchin: fire rays from a point outward through a geometry along many rays");
  std::string dagversion;
  DagMC::version( &dagversion );
  po.setVersion( dagversion );

  //   cadfile - filename/path of DAGMC geometry file
  std::string input_file;
  po.addRequiredArg<std::string>( "input_file", "Path to input file for preprocessing", &input_file );

  //   rayfile - filename/path file with ray-directions
  std::string ray_file;
  po.addRequiredArg<std::string>( "ray_file", "Path to input file containing ray directions", &ray_file);

  //   (x,y,z) - starting point
  std::vector<double> starting_point;
  po.addRequiredArg<std::vector<double>>( "starting_point", "(x,y,z) location of ray origins", &starting_point);

  // todo: optional volume ID as input
  int start_vol_id = -1;
  po.addOpt<int>("vols,V", "Volume ID of starting volume", &start_vol_id)
  
  po.parseCommandLine( argc, argv );

  moab::CartVect start_pt(starting_point);
  // create urchin to fire rays
  Ray_Urchin *urch = Ray_Urchin(input_file,ray_file,
                                start_pt,start_vol_id);

  // read DAGMC file and setup for transport
  moab::ErroCode ret;
  ret = urch->Setup();
  MB_CHK_SET_ERR_CONT(ret,"Failed to perform setup for Urchin.");
  
  int ray_num[3];
  std::vector<moab::CartVect> ray_list;
  // load input needed for rays
  ret = urch->rays_load_and_init(ray_nums,ray_list);
  MB_CHK_SET_ERR_CONT(ret,"Failed to load & initalise rays from file.");

  // start list of touched volumes
  std::set< EntityHandle > include_vols;
  
  // for each ray
  std::map< CartVect, Ray_History > ray_hist;
  for (std::vector< CartVect >::iterator dir = ray_list.begin();
       dir != ray_list.end();
       dir++){

    ray_hist[*dir] = get_ray_hist(*dir, graveyard_handle, include_vols);  
  }

  // write file
  write_hzetrn2015(ray_nums, include_vols, ray_list, ray_hist);
  
}
