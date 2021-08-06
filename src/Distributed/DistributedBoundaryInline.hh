#include "NodeList/NeighborNodeList.hh"
#include "Field/Field.hh"
#include "DataBase/DataBase.hh"
#include "Communicator.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Do not use the ghost nodes from the parallel boundary for mesh generation.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
DistributedBoundary<Dimension>::
meshGhostNodes() const {
  return false;
}

#ifdef USE_MPI
//------------------------------------------------------------------------------
// Get the domain ID.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
DistributedBoundary<Dimension>::domainID() const {
  return mDomainID;
}

//------------------------------------------------------------------------------
// Get the number of domains.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
DistributedBoundary<Dimension>::numDomains() const {
  int nProcs;
  MPI_Comm_size(Communicator::communicator(), &nProcs);
  return nProcs;
}

//------------------------------------------------------------------------------
// Get the current map of NodeList <-> (domain, DomainBoundaryNodes) pairs.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename DistributedBoundary<Dimension>::NodeListDomainBoundaryNodeMap&
DistributedBoundary<Dimension>::nodeListDomainBoundaryNodeMap() const {
  return mNodeListDomainBoundaryNodeMap;
}

//------------------------------------------------------------------------------
// Descendent classes can get read/write access to the 
// NodeListDomainBoundaryNodeMap.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename DistributedBoundary<Dimension>::NodeListDomainBoundaryNodeMap&
DistributedBoundary<Dimension>::accessNodeListDomainBoundaryNodeMap() {
  return mNodeListDomainBoundaryNodeMap;
}

#endif // USE_MPI

}
