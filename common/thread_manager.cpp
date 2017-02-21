#include <iostream>
#include "thread_manager.hpp"

// constructor
DagThreadManager::DagThreadManager(int num_threads, moab::Interface *MBI){
  // if the moab pointer is not null point to it
  if(MBI != NULL){
    MOAB->MBI;
  } else {
    // make new moab instance to share all DAGMC data
    MOAB = new moab::Core();
  }
  // number of threads
  num_threads = thread_count;

  // initalise the Dagmc instance vector
  dagmc_instances.reserve(num_threads);
  
  // loop over the number of threads and make a new DAGMC instance for each
  for ( int i = 0 ; i < num_threads ; i++ ) {
    dagmc_instances.push_back(new DagMC(MOAB));
  }
}

// destructor
DagThreadManager::~DagThreadManager(){
  for ( int i = 0 ; i < num_threads ; i++ ) {
    // delete each dagmc instance
    delete dagmc_instances[i];
  }
}

// initalise threads with the loaded data
void DagThreadManager::initalise_threads(bool init_master) {
  for ( int i = 1 ; i < num_threads ; i++ ) {
    moab::ErrorCode rval;
    rval = get_dagmc_instance(i)->load_existing_contents();
    rval = get_dagmc_instance(i)->init_OBBTree();
  }
}
