//---------------------------------Spheral++----------------------------------//
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"

#include <vector>
#include <tuple>

#ifdef USE_MPI
#include <mpi.h>
#include "Distributed/Communicator.hh"
#endif

namespace Spheral {


//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for all nodes in the given set of
// NodeLists, returning the result as a FieldList<int>.
//------------------------------------------------------------------------------
template<typename Dimension>
// inline
void
setUniqueNodeIDs(FieldList<Dimension,size_t>& uniqueIndex) {

  // Prepare the result.
  const auto numFields = uniqueIndex.numFields();
  const auto maxUnique = uniqueIndex.max();

#ifdef USE_MPI

  // This processors domain id.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // Count up how many nodes have uninitialized unique indices
  int numDomainNodes = 0;
  for (auto fieldi = 0u; fieldi < numFields; ++fieldi){
      const auto& nodeList = uniqueIndex[fieldi]->nodeList();
      const auto ni = nodeList.numInternalNodes();
    for (auto i = 0u; i < ni; ++i) {
      if(uniqueIndex(fieldi,i)==0) numDomainNodes += 1;
    }  
  }
  
  // loop processors and send cumulative number to the next
  int beginID = maxUnique;
  for (int sendProc = 0; sendProc < numProcs - 1; ++sendProc) {
    int sendProcDomainNodes = numDomainNodes;
    MPI_Bcast(&sendProcDomainNodes, 1, MPI_INT, sendProc, Communicator::communicator());
    if (procID > sendProc) beginID += sendProcDomainNodes;
  }

  // initialize the previously uninitialized
  int k = 1;
  for (auto fieldi = 0u; fieldi < numFields; ++fieldi) {
    const auto& nodeList = uniqueIndex[fieldi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    for (auto i = 0u; i != ni; ++i){
      if(uniqueIndex(fieldi,i)==0){
        uniqueIndex(fieldi,i) = beginID + k;
        k++;
      }
    }
  }

#else

  int k = 1;
  for (auto fieldi = 0u; fieldi < numFields; ++fieldi) {
    const auto& nodeList = uniqueIndex[fieldi]->nodeList();
    const auto ni = nodeList.numInternalNodes();

    for (auto i = 0u; i != ni; ++i){
      if(uniqueIndex(fieldi,i)==0){
        uniqueIndex(fieldi,i) = maxUnique + k;
        k++;
      }
    }
  }


#endif

}

}
