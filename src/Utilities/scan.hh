//---------------------------------Spheral++----------------------------------//
// scan
//
// Hide (some) of the details about doing MPI scan
//----------------------------------------------------------------------------//
#ifndef __Spheral_scan__
#define __Spheral_scan__

#include "Utilities/DataTypeTraits.hh"

#ifdef USE_MPI
//------------------------------------------------------------------------------
// MPI version
//------------------------------------------------------------------------------

#include "mpi.h"

namespace Spheral {

template<typename Value>
Value
scan(const Value& value, const MPI_Op op, const MPI_Comm comm) {
  Value tmp = value;
  Value result;
  MPI_Scan(&tmp, &result, 1, DataTypeTraits<Value>::MpiDataType(), op, comm);
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
#define MPI_PROD 4
#define MPI_LAND 5
#define MPI_LOR 6
#define MPI_COMM_WORLD 0

template<typename Value>
Value
scan(const Value& value, const int op, const int comm) {
  return value;
}

}

#endif
#endif

