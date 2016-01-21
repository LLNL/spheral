//---------------------------------Spheral++----------------------------------//
// allReduce
//
// Hide (some) of the details about doing MPI all reduces.
//
// Created by JMO, Wed Feb 10 14:38:05 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_allReduce__
#define __Spheral_allReduce__

#include "Utilities/DataTypeTraits.hh"

#ifdef USE_MPI
//------------------------------------------------------------------------------
// MPI version
//------------------------------------------------------------------------------

#include "mpi.h"

namespace Spheral {

template<typename Value>
Value
allReduce(const Value& value, const MPI_Op op, const MPI_Comm comm) {
  Value tmp = value;
  Value result;
  MPI_Allreduce(&tmp, &result, 1, DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

}

#else
//------------------------------------------------------------------------------
// Non-MPI version
//------------------------------------------------------------------------------

namespace Spheral {

#define MPI_MIN 1
#define MPI_MAX 2
#define MPI_SUM 3
#define MPI_COMM_WORLD 0

template<typename Value>
Value
allReduce(const Value& value, const int op, const int comm) {
  return value;
}

}

#endif
#endif

