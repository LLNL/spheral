//---------------------------------Spheral++----------------------------------//
// mpiUtilities
//
// Useful methods for doing MPI stuff.  These methods are dummied to do the
// right thing if we're actually building serial as well.
//
// Created by JMO, Thu Aug 21 17:17:57 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_mpiUtilities__
#define __Spheral_mpiUtilities__

#ifdef USE_MPI 
//==============================================================================
// MPI versions.
//==============================================================================
#include "mpi.h"
#include "Distributed/Communicator.hh"

namespace Spheral {

inline
double
safeAllReduceMax(double val) {
  double result;
  MPI_Allreduce(&val, &result, 1, MPI_DOUBLE, MPI_MAX, Communicator::communicator());
  return result;
}

}

#else

//==============================================================================
// Dummy serial versions.
//==============================================================================
namespace Spheral {
inline
double
safeAllReduceMax(double val) {
  std::cerr << "Not here!" << std::endl;
  return val;
}
}

#endif
#endif
