void dagmc_load_and_init(DagMC &dagmc, std::string filename)
{

  EntityHandle input_file_set;
  ErroCode ret;
  
  ret = dagmc.load_file( filename.c_str() );
  CHKERR( ret, "Failed to load file ", filename);

  ret = dagmc.init_OBBTree();
  CHKERR( ret, "Failed to build OBB tree and/or indices.", "");

}

void rays_load_and_init( std::string ray_file, std::vector< int > &ray_nums, std::vector< CartVect > &ray_list)
{

  ifstream ray_file_stream;
  ray_file_stream.open(ray_file);

  double x, y, z;

  ray_num.resize(3);
  ray_file_stream >> ray_num[0] >> ray_num[1] >> ray_num[2];
  
  while !(ray_file_stream.eof())
    {
      ray_file_stream >> x >> y >> z;
      ray_list.push_bach(CartVect(x,y,z));
    }
  
}

std::vector< std::pair< int, double> >& get_ray_hist(CartVect start_pt, int start_vol, CartVect& dir, std::set< int >& includ_vols){

  std::vector< std::pair< int, double> > ray_hist;
  
  int vol = start_vol;
  CartVect pt = start_pt;
  double dist;
  EntityHandle surf;
  
  //   fire ray to get all intersections to graveyard : {dir => {cell => track length}}
  while (vol != graveyard_vol){
    ret = dagmc.ray_fire( vol, start_pt, dir, surf, dist);
    ray_hist.push_back(std::pair(vol,dist));
    vol = dagmc.next_vol(vol, surf);
    include_vols.insert(vol);
  }
  
  // pop all the void regions off the end of the list
  while (get_material(ray_hist.end()->first) == 0)
    {
      include_vols.remove(ray_hist.end()->first);
      ray_hist.pop();
    }
  
  return ray_hist;
}


void write_urchin(std::vector<int>& ray_nums,
                  CartVect& start_pt,
                  int start_vol,
                  std::set< int >& includ_vols,
                  std::vector< CartVect >& ray_list,
                  std::map< CartVect, std::vector< std::pair< int, double > > >& ray_hist) {

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
    for (std::vector< std::pair< int, double > >::iterator slab = ray_hist[*dir].begin();
         slab != ray_hist[*dir].end();
         slab++) {
      std::cout << slab->first << "\t" slab->second << std::endl;
    }
  }

  
  
}

int main(int argc, char* argvp[] ){

  // Pseudo-code
  // initialize MOAB & DAGMC
  
  // read input

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
  int start_vol;
  po.addOpt<int>("vols,V", "Volume ID of starting volume", &start_vol)
  
  po.parseCommandLine( argc, argv );

  // check for valid input ????
  
  // read DAGMC file and setup for transport
  dagmc_load_and_init(input_file);
  
  // read rayfile
  std::vector< int > ray_nums;
  std::vector< CartVect > ray_list;
  ret = rays_load_and_init(ray_file, ray_nums, ray_list);

  // find cell of starting point
  CartVect start_p(starting_point);
  if (start_vol) {
    // check that start point is in vol
  } else {
    // search for start vol
    start_vol = find_start_vol(start_p);
  }

  // start list of touched volumes
  std::set< int > include_vols;
  
  // for each ray
  std::map< CartVect, std::vector< std::pair< int, double > > > ray_hist;
  for (std::vector< CartVect >::iterator dir = ray_list.begin();
       dir != ray_list.end();
       dir++){

    ray_hist[*dir] = get_ray_hist(start_pt, start_vol, *dir, include_vols);
    
  }

  // write file
  write_urchin();
  
}
