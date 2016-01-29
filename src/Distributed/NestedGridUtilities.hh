//---------------------------------Spheral++----------------------------------//
// NestedGridUtilities -- A set of useful common helper methods for working
// with NestedGridNeighbor structures in parallel.
//
// Created by JMO, Wed Nov 24 14:00:34 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_Distributed_NestedGridUtilities_hh
#define Spheral_Distributed_NestedGridUtilities_hh

#include "DataBase/DataBase.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using NeighborSpace::NestedGridNeighbor;

//------------------------------------------------------------------------------
// Get the NestedGridNeighbor associated with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NeighborSpace::NestedGridNeighbor<Dimension>&
getNestedGridNeighbor(const NodeSpace::NodeList<Dimension>* nodeListPtr) {
  NeighborSpace::NestedGridNeighbor<Dimension>& neighbor = dynamic_cast<NeighborSpace::NestedGridNeighbor<Dimension>&>(nodeListPtr->neighbor());
  return neighbor;
}

//------------------------------------------------------------------------------
// Determine the max number of occupied grid levels for all NodeLists in a
// DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
maxNumGridLevels(const DataBaseSpace::DataBase<Dimension>& dataBase,
                 MPI_Comm& communicator) {

  int numGridLevels = 0;

  // Loop over the NodeLists in the DataBase.
  for (typename DataBaseSpace::DataBase<Dimension>::ConstNodeListIterator nodeListItr = dataBase.nodeListBegin();
       nodeListItr != dataBase.nodeListEnd();
       ++nodeListItr) {
    NeighborSpace::NestedGridNeighbor<Dimension>& neighbor = getNestedGridNeighbor(*nodeListItr);
    numGridLevels = std::max(numGridLevels, neighbor.numGridLevels());
  }

  // Find the global maximum of the number of grid levels.
  int globalNumGridLevels;
  MPI_Allreduce(&numGridLevels, &globalNumGridLevels, 1, MPI_INT, MPI_MAX, communicator);

  CHECK(globalNumGridLevels > 0);
  return globalNumGridLevels;
}

}

#endif
