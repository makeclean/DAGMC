#include "DagMC.hpp"
#include "dagmcmetadata.hpp"
#include "moab/Core.hpp"
#include <utility>

typedef std::vector< std::pair< moab::EntityHandle, double > > Ray_History; ///< collection of volumes that the ray has encountered
typedef std::vector< std::pair< moab::EntityHandle, double> >::iterator Ray_History_iter; ///< iterator for Ray_History type

class Ray_Urchin {

public:
  // Default constructor
  Ray_Urchin(const std::string cad_file,
             const std::string ray_file,
             const moab::CartVect start_pt_in,
             const int start_vol_id);
  // Destructor
  ~Ray_Urchin();

  // function to setup the internal state of urchin to 
  // be ready to fire rays
  moab::ErrorCode Setup();

 
  moab::ErrorCode process_rays();
  
  // function to write output from urchin
  void write_hzetrn2015();


// private functions
private:
  // get the material name for a given volume
  std::string get_material(moab::EntityHandle volume);

  // find the graveyard volume
  void find_graveyard();

  // function to fine the volume that contains the starting point
  moab::ErrorCode find_start_vol();

 // function to open a text file and read a set of unit vectors
  // input:
  //   std::string ray_file - path/filename of ray information
  // output:
  //   int[3]& ray_nums - total number of rays + lat/long number of rays
  //   std::vector< CartVect >& ray_list - vector of ray unit vectors
  moab::ErrorCode rays_load_and_init(int (&ray_nums)[3], 
                                    std::vector<moab::CartVect> &ray_list);
  
  // 
  Ray_History get_ray_hist(const moab::CartVect dir,
                           std::set<moab::EntityHandle> &include_vols );

// private member variables
private:
  moab::DagMC* dagmc;
  dagmcMetaData* dmd;

  std::string dagmc_file;
  std::string urchin_file; 

  int start_vol_id;
  moab::CartVect start_position;
  moab::EntityHandle start_vol_handle;
  moab::EntityHandle graveyard_handle;

  int ray_nums[3]; // number of rays for each direction?
  std::vector<moab::CartVect> ray_list; // list of rays fired
  std::map<int,Ray_History> ray_hist; /// collection of histories for each direction
  std::set<moab::EntityHandle> include_vols; /// list of touched volumes



};


