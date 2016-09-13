//----------------------------------*-C++, Fortran-*----------------------------------//
/*!
 * \file   ~/DAGMC/FluDAG/src/cpp/mainFluDAG.cpp
 * \author Julie Zachman
 * \date   Apr 5 2013
 * \brief  Functions called by fluka
 * \note   Unittested
 */
//---------------------------------------------------------------------------//
#include "fluka_funcs.h"

#include "DagMC.hpp"

#include <cstring>
#include <fstream>
#include <time.h>       // for timing the routine

#define flukam flukam_

moab::DagMC *DAG = new moab::DagMC();

#ifdef __cplusplus
extern "C" {
#endif
  void flukam(const int &GeoFlag);
#ifdef __cplusplus
}
#endif

// This has been modified so that all runs are now fluka runs.
// Getting the geometry out of the file is done by dagmc_define
int main(int argc, char* argv[]) {
  
  bool flukarun = false;
  moab::ErrorCode error;
  time_t time_before,time_after;

  // Default h5m filename is for fluka runs
  std::string infile = "dagmc.h5m";

  if ( argc == 1 ) { // then its a fluka run
    // fluka creates a run dir one lvl higher
    infile = "../"+infile;
    flukarun = true;
  } else if ( argc > 2 ) {
    std::cout << "run as main_fludag <facet_file>  to produce"
              << " material assignments" << std::endl;
    std::cout << "too many arguments provided" << std::endl;
    exit(1);
  } else { // its a pre process run
    infile = argv[1]; // must be the 2nd argument
  }

  // check for file existence
  std::ifstream h5mfile (infile.c_str()); // filestream for mesh geom
  if ( !h5mfile.good() ) {
    std::cout << "h5m file does not exist" << std::endl;
    exit(1);
  }

  int max_pbl = 1;

  // get the current time
  time(&time_before);  /* get current time; same as: timer = time(NULL)  */
  // DAG call to load the file
  std::cout << "Loading the faceted geometry file " << infile << "..." << std::endl;
  error = DAG->load_file(infile.c_str()); // load the dag file takeing the faceting from h5m
  if ( error != moab::MB_SUCCESS ) {
    std::cerr << "DAGMC failed to read input file: " << infile << std::endl;
    exit(EXIT_FAILURE);
  }

  // get time to load file 
  time(&time_after);
  double seconds = difftime(time_after,time_before); //get the time in seconds to load file
  time_before = time_after; // reset time to now for the next call
  std::cout << "Time to load the h5m file = " << seconds << " seconds" << std::endl;

  
  // if we are doing a full run, init obb's etc
  if ( flukarun) {
    // DAG call to initialize geometry
    error = DAG->init_OBBTree();
    if ( error != moab::MB_SUCCESS ) {
      std::cerr << "DAGMC failed to initialize geometry and create OBB tree" <<  std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    // if we aren't doing a full run, dont need obbs - much quicker to start
    error = DAG->setup_impl_compl();
    if( error != moab::MB_SUCCESS ) {
      std::cerr << "DAGMC failed to setup implicit compliment" <<  std::endl;
      exit(EXIT_FAILURE);
    }
    error = DAG->setup_indices();
    if( error != moab::MB_SUCCESS ) {
      std::cerr << "DAGMC failed to setup problem indices" <<  std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // figure out how long we took to load
  time(&time_after);
  seconds = difftime(time_after,time_before);
  std::cout << "Time to initialise the geometry" << seconds << std::endl;

  // if fluka preprocess run then create mat file to paste into input deck
  if (!flukarun) {
    std::string lcad = "mat.inp";
    std::cout << "Producing material snippets" << std::endl;
    std::cout << "please paste these into your input deck" << std::endl;
    fludag_write(infile, lcad);

    std::string vol_id = "vol_id_idx.txt";
    std::cout << "Producing volume index & id correspondences" << std::endl;
    fludag_write_ididx(vol_id);
  } else {
    // call fluka run
    // flugg mode is flag = 1
    const int flag = 1;
    flukam(flag);
  }

  return 0;
}

