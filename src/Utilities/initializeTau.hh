// This is a simple little method we use to initialize the Tau profiling package.
#include "TAU.h"
#include "Utilities/Communicator.hh"

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace Spheral {
inline
void
initializeTau() {
#ifdef PROFILING_ON
  TAU_PROFILE("initializeTau", "", TAU_USER);
  int myid = 0;
#ifdef USE_MPI
  MPI_Comm_rank(Communicator::communicator(),&myid);  
#endif
  TAU_PROFILE_SET_NODE(myid);
#endif
}

}
