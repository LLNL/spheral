#include "DomainNode.hh"
#include "DBC.hh"
#include "cdebug.hh"

namespace Spheral {
namespace PartitionSpace {

#ifdef USE_MPI
//------------------------------------------------------------------------------
// Get the domain ID.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
RedistributeNodes<Dimension>::domainID() const {
//   cdebug << "RedistributeNodes::domainID()" << endl;
  int domainID;
  MPI_Comm_rank(mCommunicator, &domainID);
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
//   cdebug << "RedistributeNodes::numDomains()" << endl;
  int nProcs;
  MPI_Comm_size(mCommunicator, &nProcs);
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
RedistributeNodes<Dimension>::computeWork(const bool x) {
  mComputeWork = x;
}


}
}
