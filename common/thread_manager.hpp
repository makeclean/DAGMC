#include "DagMC.hpp"
#include "Core.hpp"

class DagThreadManager {
  public: 
  DagThreadManager();
 ~DagThreadManager();

  // setup the dagmc state for the threads
  void initialise_threads();
    
  inline DagMC* get_dagmc_instance(int thread_id) {
    return dagmc_instances[thread_id];
  }

  private:
  int num_threads; ///< Number of threads to be held by the manager
  std::vector<DagMC*> dagmc_instances; ///< vector of dagmc instances
  moab::Interface *MOAB; ///< moab pointer
}
