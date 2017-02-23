
// function to open a text file and read a set of unit vectors
// input:
//   std::string ray_file - path/filename of ray information
// output:
//   int[3]& ray_nums - total number of rays + lat/long number of rays
//   std::vector< CartVect >& ray_list - vector of ray unit vectors
Errorcode rays_load_and_init( std::string ray_file, int[3] &ray_nums, std::vector< CartVect > &ray_list)
{

  ifstream ray_file_stream;
  ray_file_stream.open(ray_file);

  if (!ray_file_stream)
    MB_SET_ERR("Unable to open ray file " + ray_file);

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
ErrorCode find_start_vol(Dagmc& dagmc, const CartVect start_p, const int start_vol_id,
                         EntityHandle& start_vol_handle) {

  ErrorCode ret;
  std::vector<EntityHandle> vols = dagmc.volumes();
  int inside = 0;
  
  start_vol_handle = dagmc.entity_by_id(3,start_vol_id);
  
  if (start_vol_id >= 0) {
    ret = dagmc.point_in_volume(start_vol_handle,start_p.array(), inside);
    MB_CHK_ERR(ret);
    if (!inside)
      MB_SET_ERR_CONT("Starting point is not in provided starting volume.");
  }
  
  if (!inside) {
    for (std::vector<EntityHandle>::iterator vol = vols.begin(); vol != vols.end && !inside; vol++) {
      start_vol_handle = *vol;
      ret = dagmc.point_in_volume(*vol, start_p.array(), inside);
      MB_CHK_ERR(ret);
    }
  }
  
  if (!inside) {
    return MB_ENTITY_NOT_FOUND;
  }
  
  return MB_SUCCESS;

}

// function to find the graveyard volume
// input:
//    DagMC dagmc - dagmc instance
EntityHandle find_graveyard(DagMC& dagmc) {

  // get group with name "graveyard"
  // get members of that group
  // return entity handle of that group

  
}

// function to repeatedly fire a ray from a point to the boundary and
//    accumulate a list of pairs for each volume traversed,
//    with each pair containing the volume handle and the distance traversed
// input:
//    DagMC& dagmc - the DagMC instance
//    CartVect start_pt - the starting point of the rays
//    EntityHandle start_vol_handle - the handle for the starting volume
//    CartVect dir - the unit vector for the current ray
//    EntityHandle graveyard_handle - handle for the graveyard volume
// output:
//    std::set< EntityHandle >& include_vols - list of volume entities touched so far
//    (return) Ray_History& - ray history from start_pt to boundary
Ray_History& get_ray_hist(DagMC& dagmc,
                                                              const CartVect start_pt,
                                                              const EntityHandle start_vol_handle,
                                                              const CartVect dir,
                                                              const EntityHandle graveyard_handle,
                                                              std::set< EntityHandle >& include_vols){

  Ray_History ray_hist;

  EntityHandle vol = start_vol_handle;
  CartVect pt = start_pt;
  double dist;
  EntityHandle new_vol = vol;
  EntityHandle surf = 0;
  
  //   fire ray to get all intersections to graveyard : {dir => {cell => track length}}
  while (vol != graveyard_handle)){
    ret = dagmc.ray_fire( vol, start_pt.array(), dir.array(), surf, dist);
    MB_CHK_ERR_CONT(ret);
    ray_hist.push_back(std::pair(vol,dist));
    include_vols.insert(vol);

    ret = dagmc.next_vol(surf,vol,new_vol);
    MB_CHK_ERR_CONT(ret);
    vol = new_vol;
  }
  
  // pop all the void regions off the end of the list
  int void_matl = 0;
  while (get_material(ray_hist.end()->first) == void_matl)
    {
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
void write_urchin(const int[3]& ray_nums,
                  const CartVect& start_pt,
                  const EntityHandle start_vol_handle,
                  const std::set< EntityHandle >& include_vols,
                  const std::vector< CartVect >& ray_list,
                  const std::map< CartVect, Ray_History >& ray_hist) {

  // - write header
  std::cout << "***GEOMETRY_DEFN***" << std::endl;
  std::cout << ray_nums[0] << "t" << ray_nums[1] << "\t" << ray_nums[2] << std::endl << std::endl;
  
  // - write 'raytrace' flag
  std::cout << "raytrace" << std::endl << std::endl;
  
  // - write start point & start cell
  std::cout << start_pt[0] << "\t" << start_pt[1] << "\t" << start_pt[2] << "\t" << start_vol << std::endl << std::endl;
  
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
  for (std::vector< CartVect >::iterator dir = ray_list.begin();
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

Ray_Urchin::Ray_Urchin(const std::string cad_file,
                       const std::string ray_file,
                       const CartVect start_pt_in,
                       const int start_vol_id) {
  
  dagmc = new DagMC();
  
  ret = dagmc.load_file( input_file.c_str() );
  MB_CHK_SET_ERR_CONT(ret,"Failed to load file " + input_file);
  
  ret = dagmc.init_OBBTree();
  MB_CHK_SET_ERR_CONT(ret, "Failed to build OBB tree and/or indices.");
  
  // read rayfile
  int ray_nums[3];
  std::vector< CartVect > ray_list;
  ret = rays_load_and_init(ray_file, ray_nums, ray_list);
  MB_CHK_ERR_CONT(ret);

  // find cell of starting point
  CartVect start_p(starting_point);
  EntityHandle start_vol_handle = 0;
  // search for start vol
  ret = find_start_vol(start_p, start_vol_id, &start_vol_handle);
  MB_CHK_SET_ERR_CONT(ret,"Can't find starting volume for provided starting point.");

  // find graveyard handle
  EntityHandle graveyard_handle = find_graveyard(dagmc);

  
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

  // check for valid input ????
  
  // read DAGMC file and setup for transport
  ErroCode ret;
  
  
  // start list of touched volumes
  std::set< EntityHandle > include_vols;
  
  // for each ray
  std::map< CartVect, Ray_History > ray_hist;
  for (std::vector< CartVect >::iterator dir = ray_list.begin();
       dir != ray_list.end();
       dir++){

    ray_hist[*dir] = get_ray_hist(start_pt, start_vol, *dir, graveyard_handle, include_vols);
    
  }

  // write file
  write_urchin(ray_nums, start_pt, start_vol_handle, include_vols, ray_list, ray_hist);
  
}
