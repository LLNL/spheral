//---------------------------------Spheral++----------------------------------//
// waitAllWithDeadlockDetection
// 
// Replacement for MPI_Waitall, which automatically spits out after a 
// predetermined time interval any processes that are stuck at the wait.
// Useful for detecting, you know, deadlocks.
//
// Directly cribbed from Jeff Johnsons original code, and cast here as a 
// standalone dohicky.
//----------------------------------------------------------------------------//
#ifndef __Spheral_waitAllWithDeadlockDetection__
#define __Spheral_waitAllWithDeadlockDetection__

#include <string>
#include "mpi.h"

namespace Spheral {
void 
waitallWithDeadlockDetection(const std::string label,
                             const std::vector<int>& sendProcIDs,
                             const std::vector<int>& recvProcIDs, 
                             std::vector<MPI_Request>& sendRequests, 
                             std::vector<MPI_Request>& recvRequests, 
                             std::vector<MPI_Status>& sendStatuses,
                             std::vector<MPI_Status>& recvStatuses,
                             MPI_Comm comm);
}

#endif
