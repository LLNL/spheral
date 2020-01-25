#include "Utilities/DomainNode.hh"
#include "Communicator.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

#ifdef USE_MPI
//------------------------------------------------------------------------------
// Get the domain ID.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RedistributeNodes<Dimension>::domainID() const {
  int domainID;
  MPI_Comm_rank(Communicator::communicator(), &domainID);
  ENSURE(domainID >= 0 && domainID < numDomains());
  return domainID;
}

//------------------------------------------------------------------------------
// Get the number of domains.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RedistributeNodes<Dimension>::numDomains() const {
  int nProcs;
  MPI_Comm_size(Communicator::communicator(), &nProcs);
  return nProcs;
}
#endif // USE_MPI

//------------------------------------------------------------------------------
// Flag controlling how we compute the work.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RedistributeNodes<Dimension>::computeWork() const {
  return mComputeWork;
}

template<typename Dimension>
inline
void
RedistributeNodes<Dimension>::computeWork(bool x) {
  mComputeWork = x;
}

}
