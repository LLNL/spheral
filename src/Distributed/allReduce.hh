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
#include "Communicator.hh"

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace Spheral {
#ifdef USE_MPI
//------------------------------------------------------------------------------
// MPI version
//------------------------------------------------------------------------------

#define SPHERAL_OP_MIN MPI_MIN
#define SPHERAL_OP_MAX MPI_MAX
#define SPHERAL_OP_SUM MPI_SUM
#define SPHERAL_OP_PROD MPI_PROD
#define SPHERAL_OP_LAND MPI_LAND
#define SPHERAL_OP_LOR MPI_LOR

template<typename Value>
constexpr Value
allReduce(const Value& value, const MPI_Op op,
          const MPI_Comm comm = Communicator::communicator()) {
  Value tmp = value;
  Value result;
  MPI_Allreduce(&tmp, &result, 1,
                DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

template<typename Value>
constexpr Value
distScan(const Value& value, const MPI_Op op,
     const MPI_Comm comm = Communicator::communicator()) {
  Value tmp = value;
  Value result;
  MPI_Scan(&tmp, &result, 1, DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

#else
//------------------------------------------------------------------------------
// Non-MPI version
//------------------------------------------------------------------------------

#define SPHERAL_OP_MIN 1
#define SPHERAL_OP_MAX 2
#define SPHERAL_OP_SUM 3
#define SPHERAL_OP_PROD 4
#define SPHERAL_OP_LAND 5
#define SPHERAL_OP_LOR 6

template<typename Value>
constexpr Value
allReduce(const Value& value, const int /*op*/, const int = 0) {
  return value;
}

template<typename Value>
constexpr Value
distScan(const Value& value, const int /*op*/, const int = 0) {
  return value;
}

inline void
MPI_Barrier(const int = 0) {
  return;
}
#endif
}
#endif

