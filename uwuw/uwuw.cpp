#include <unistd.h>
#include <iostream>
#include "uwuw.hpp"

// Empty Constructor
UWUW::UWUW(bool verbose) {
  num_tallies = 0;
  num_materials = 0;
  verbosity = verbose;
};

// Default constructor
UWUW::UWUW(char* file, bool verbose) {
  std::string filename(file);

  // turn the filename into a full filepath
  full_filepath = get_full_filepath(filename);

  if(!check_file_exists(full_filepath)) {
    std::cerr << "The file " << full_filepath << " does not exist or is read protected" << std::endl;
    exit(1);
  }

  // load materials
  material_library = load_pyne_materials(full_filepath);

  // load tallies
  tally_library = load_pyne_tallies(full_filepath);

  verbosity = verbose;
};

// Default constructor
UWUW::UWUW(std::string filename, bool verbose) {
  // turn the filename into a full filepath
  full_filepath = get_full_filepath(filename);

  // check for file existence
  if(!check_file_exists(full_filepath)) {
    std::cerr << "The file " << full_filepath << " does not exist or is read protected" << std::endl;
    exit(1);
  }

  // load materials
  material_library = load_pyne_materials(full_filepath);

  // load tallies
  tally_library = load_pyne_tallies(full_filepath);

  verbosity = verbose;
};

// Destructor
UWUW::~UWUW()
{
};

// convert convert a filename into path+filename (for pyne)
std::string UWUW::get_full_filepath(char *filename)
{
  std::string file(filename);
  return UWUW::get_full_filepath(file);
}

// convert convert a filename into path+filename (for pyne)
std::string UWUW::get_full_filepath(std::string filename)
{
  // remove all extra whitespace
  filename.erase(std::remove(filename.begin(), filename.end(),' '), filename.end());
  // use stdlib call
  const char* full_filepath = realpath(filename.c_str(),NULL);
  return std::string(full_filepath);
}

// see if file exists
bool UWUW::check_file_exists(std::string filename)
{
  // from http://stackoverflow.com/questions/12774207/
  // fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
  std::ifstream infile(filename.c_str());
  return infile.good();
}

// loads all materials into map
std::map<std::string, pyne::Material> UWUW::load_pyne_materials(std::string filename, std::string datapath)
{
  std::map<std::string, pyne::Material> library; // material library

  const char* data_path = datapath.c_str();

  if(!hdf5_path_exists(filename,data_path))
    return library;

  num_materials = get_length_of_table(filename,datapath);

  for ( int i = 0 ; i < num_materials ; i++ ) {
    pyne::Material mat; // from file
    mat.from_hdf5(filename,datapath,i);
    // renumber material number by position in the library
    mat.metadata["mat_number"]=i+1;
    library[mat.metadata["name"].asString()]=mat;
  }

  return library;
}

// loads all tallies into map
std::list<pyne::Tally> UWUW::load_pyne_tallies(std::string filename, std::string datapath)
{
  std::list<pyne::Tally> library; // material library

  if(!hdf5_path_exists(filename,datapath))
    return library;

  num_tallies = get_length_of_table(filename,datapath);

  for ( int i = 0 ; i < num_tallies ; i++) {
    pyne::Tally tally; // from file
    tally.from_hdf5(filename,datapath,i);
    library.push_back(tally);
  }
  return library;
}

// see if path exists before we go on
bool UWUW::hdf5_path_exists(std::string filename, std::string datapath)
{
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  //Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  bool datapath_exists = h5wrap::path_exists(db, datapath.c_str());

  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return datapath_exists;
}

// get the length of the datapath found in datapath
int UWUW::get_length_of_table(std::string filename, std::string datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  //Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);

  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  hid_t ds = H5Dopen2(db, datapath.c_str(), H5P_DEFAULT);

  // Initilize to dataspace, to find the indices we are looping over
  hid_t arr_space = H5Dget_space(ds);

  hsize_t arr_dims[1];
  int arr_ndim = H5Sget_simple_extent_dims(arr_space, arr_dims, NULL);

  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);

  return arr_dims[0];
}

// writes all uwuw materials and tallies to file
void UWUW::write_uwuw(std::string filename){
  write_materials(filename);
  write_tallies(filename);
}

//
void UWUW::write_uwuw(char *filename) {
  std::string filename_str(filename);
  write_uwuw(filename_str);
}

// write the new material library
void UWUW::write_materials(std::string output_filename, std::string datapath) {
  std::map<std::string,pyne::Material> :: iterator it;

  // loop over the processed material library and write each one to the file
  for ( it = material_library.begin() ; it != material_library.end() ; ++it ) {
    // the current material
    pyne::Material mat = it->second;
    // write the material to the file
    if(verbosity) {
      std::cout << "writing material, " << mat.metadata["name"].asString();
      std::cout << "writing material, " << mat.metadata["fluka_name"].asString();
      std::cout << " to file " << output_filename << std::endl;
    }
    // write the uwuw materials to the geometry
    mat.write_hdf5(output_filename,datapath);
  }

  // write the nucid table
  write_nucids(output_filename);
  return;
}

// write the nucids to the file and datapath specified
void UWUW::write_nucids(std::string filename, std::string datapath) {
  // Turn off annoying HDF5 errors
  herr_t status;
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  
  //Set file access properties so it closes cleanly
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);

  // open the file
  hid_t db = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl);

  // actual nucid's to write
  std::vector<int> nuc_data = get_nucid_vector();
  // get vector into array
  int* nucids = &nuc_data[0];
  // number of nuclides to write
  hsize_t nuc_dims[1];
  nuc_dims[0] = nuc_data.size();
  
  // crate space to write data 
  hid_t nuc_space = H5Screate_simple(1, nuc_dims, NULL);
  // crate dataspace to write data to
  hid_t nuc_set = H5Dcreate2(db, datapath.c_str(), H5T_NATIVE_INT, nuc_space,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // write the data to the dataset
  H5Dwrite(nuc_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nucids);
  // flush the file
  H5Fflush(db, H5F_SCOPE_GLOBAL);
  //
  status = H5Eclear(H5E_DEFAULT);

  // Close the database
  status = H5Fclose(db);
}

// return sorted array of nucids
std::vector<int> UWUW::get_nucid_vector() {
  // loop over the materials in the material library
  // make a set of nucids

  // set of nucids
  std::set<int> nucids;
  std::map<std::string, pyne::Material>::iterator it;
  for ( it = material_library.begin() ;
	it != material_library.end() ;
	it++ ) {
    // a material object
    pyne::Material mat = it->second;
    //loop over the comp and extract the nucids
    pyne::comp_map comp = mat.comp;
    pyne::comp_map::iterator comp_it;
    for ( comp_it = comp.begin() ;
	  comp_it != comp.end() ;
	  comp_it++ ) {
      // insert the nucs into the list
      nucids.insert(comp_it->first);
    }
  }
  
  // having looped through all the materials and their
  // compositions we now have a set of nucids, insert
  // into an array
  std::vector<int> nucarr(nucids.begin(), nucids.end());
  // std::set is already sorted, can just insert into array
  // return the array
  return nucarr;
}


// write the new material library
void UWUW::write_tallies(std::string output_filename, std::string tally_destination) {
  std::list<pyne::Tally> :: iterator it;

  // loop over the processed material library and write each one to the file
  for ( it = tally_library.begin() ; it != tally_library.end() ; ++it ) {
    // the current tally
    pyne::Tally tally = *it;
    if(verbosity) {
      std::cout << "Writing tally " << tally.tally_name;
      std::cout << " to file " << output_filename << std::endl;
    }

    // write the material to the file
    tally.write_hdf5(output_filename,tally_destination);
  }

  return;
}
