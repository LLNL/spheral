//---------------------------------Spheral++----------------------------------//
// globalNodeIDs (contains numGlobalNodes and globalNodeIDs)
//
// Helper methods to assign unique global node IDs to all node in a NodeList,
// serial or parallel.  *Should* compute the same unique ID for each node no
// matter how the parallel domains are carved up.  
//
// Note this is accomplished by sorting based on position, so nodes occupying
// the same position are considered degenerate and may get different IDs for
// different decompositions.
//
// Created by JMO, Fri Oct 29 13:08:33 2004
//----------------------------------------------------------------------------//
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Utilities/peanoHilbertOrderIndices.hh"
#include "Utilities/KeyTraits.hh"
#include "Utilities/DBC.hh"

#include <vector>
#include <tuple>

#ifdef USE_MPI
#include <mpi.h>
#include "Distributed/Communicator.hh"
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Return the total global number of nodes in the NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
numGlobalNodes(const NodeList<Dimension>& nodeList) {
  int localResult = nodeList.numInternalNodes();
#ifdef USE_MPI
  int result;
  MPI_Allreduce(&localResult, &result, 1, MPI_INT, MPI_SUM, Communicator::communicator());
  CHECK(result >= localResult);
#else
  const int result = localResult;
#endif
  return result;
}
  
//------------------------------------------------------------------------------
// Return the total global number of nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
numGlobalNodes(const DataBase<Dimension>& dataBase) {
  return numGlobalNodes<Dimension, typename DataBase<Dimension>::ConstNodeListIterator>(dataBase.nodeListBegin(), dataBase.nodeListEnd());
}

//------------------------------------------------------------------------------
// Return the total global number of nodes in the set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator>
inline
int
numGlobalNodes(const NodeListIterator& begin,
               const NodeListIterator& end) {
  int result = 0;
  for (NodeListIterator nodeListItr = begin; nodeListItr != end; ++nodeListItr) {
    result += numGlobalNodes(**nodeListItr);
  }
  return result;
}
  
//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for the given NodeList, and return
// the set of them on this process.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Field<Dimension, int>
globalNodeIDs(const NodeList<Dimension>& nodeList) {

  typedef typename KeyTraits::Key Key;

  // Get the local domain ID and number of processors.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // Build keys to sort the nodes by.
  DataBase<Dimension> db;
  db.appendNodeList(const_cast<NodeList<Dimension>&>(nodeList));
  FieldList<Dimension, Key> keys = peanoHilbertOrderIndices(db);

  // Build the local list of node info.
  typedef std::vector<std::tuple<Key, int, int> > InfoType;
  InfoType nodeInfo;
  int numLocalNodes = nodeList.numInternalNodes();
  for (int i = 0; i != numLocalNodes; ++i) {
    nodeInfo.push_back(std::tuple<Key, int, int>(keys(0, i), i, procID));
  }

  // Reduce the list of node info to processor 0.
#ifdef USE_MPI
  int numGlobalNodes = numLocalNodes;
  if (procID == 0) {

    // Process 0 receives and builds the global info.
    for (int recvDomain = 1; recvDomain != numProcs; ++recvDomain) {
      MPI_Status status;
      int numRecvNodes;
      MPI_Recv(&numRecvNodes, 1, MPI_INT, recvDomain, 10, Communicator::communicator(), &status);
      CHECK(numRecvNodes >= 0);
      numGlobalNodes += numRecvNodes;
      std::vector<Key> packedKeys(numRecvNodes);
      std::vector<int> packedLocalIDs(numRecvNodes);
      if (numRecvNodes > 0) {
        MPI_Recv(&(*packedKeys.begin()), numRecvNodes, DataTypeTraits<Key>::MpiDataType(),
                 recvDomain, 11, Communicator::communicator(), &status);
        MPI_Recv(&(*packedLocalIDs.begin()), numRecvNodes, MPI_INT,
                 recvDomain, 12, Communicator::communicator(), &status);
      }
      for (int i = 0; i != numRecvNodes; ++i) {
        nodeInfo.push_back(std::tuple<Key, int, int>(packedKeys[i], packedLocalIDs[i], recvDomain));
      }
    }

  } else {
             
    // Send our local info to processor 0.
    std::vector<Key> packedKeys;
    std::vector<int> packedLocalIDs;
    for (typename InfoType::const_iterator itr = nodeInfo.begin();
         itr != nodeInfo.end();
         ++itr) {
      packedKeys.push_back(std::get<0>(*itr));
      packedLocalIDs.push_back(std::get<1>(*itr));
    }
    CHECK(packedKeys.size() == nodeInfo.size());
    CHECK(packedLocalIDs.size() == nodeInfo.size());

    MPI_Send(&numLocalNodes, 1, MPI_INT, 0, 10, Communicator::communicator());
    if (numLocalNodes > 0) {
      MPI_Send(&(*packedKeys.begin()), numLocalNodes, DataTypeTraits<Key>::MpiDataType(),
               0, 11, Communicator::communicator());
      MPI_Send(&(*packedLocalIDs.begin()), numLocalNodes,
               MPI_INT, 0, 12, Communicator::communicator());
    }
  }
  CHECK((int)nodeInfo.size() == numGlobalNodes);
#endif

  // Sort the node info.
  if (nodeInfo.size() > 0) {
    sort(nodeInfo.begin(), nodeInfo.end());
    BEGIN_CONTRACT_SCOPE
    for (int i = 0; i < (int)nodeInfo.size() - 1; ++i) {
      CHECK(std::get<0>(nodeInfo[i]) <= std::get<0>(nodeInfo[i + 1]));
    }
    END_CONTRACT_SCOPE
  }

  // Now we can assign consecutive global IDs based on the sorted list.
  std::vector< std::vector<int> > globalIDs(numProcs);
  for (auto i = 0u; i != nodeInfo.size(); ++i) {
    const int recvProc = std::get<2>(nodeInfo[i]);
    const unsigned int localID = std::get<1>(nodeInfo[i]);
    CHECK(recvProc < (int)globalIDs.size());
    if (localID + 1 > globalIDs[recvProc].size()) globalIDs[recvProc].resize(localID + 1);
    globalIDs[recvProc][localID] = i;
  }

  // Assign process 0's Ids.
  Field<Dimension, int> result("global IDs", nodeList);
  if (procID == 0) {
    for (int i = 0; i != numLocalNodes; ++i) result(i) = globalIDs[0][i];
  }

#ifdef USE_MPI
  // Farm the ID info back to all processors.
  if (procID == 0) {

    // Process 0 sends the info.
    for (int recvProc = 1; recvProc != numProcs; ++recvProc) {
      int numRecvNodes = globalIDs[recvProc].size();
      MPI_Send(&numRecvNodes, 1, MPI_INT, recvProc, 20, Communicator::communicator());
      if (numRecvNodes > 0) MPI_Send(&(*globalIDs[recvProc].begin()), numRecvNodes, MPI_INT,
                                     recvProc, 21, Communicator::communicator());
    }

  } else {

    // Get our ids from process 0.
    MPI_Status status;
    int numRecvNodes;
    MPI_Recv(&numRecvNodes, 1, MPI_INT, 0, 20, Communicator::communicator(), &status);
    CHECK(numRecvNodes == numLocalNodes);
    if (numRecvNodes > 0) MPI_Recv(&(*result.begin()), numRecvNodes, MPI_INT, 0, 21,
                                   Communicator::communicator(), &status);

  }
#endif

  // We're done.
  return result;
}

//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for all nodes across all NodeLists in
// a DataBase, returning the result as a FieldList<int>.
//------------------------------------------------------------------------------
template<typename Dimension>
// inline
FieldList<Dimension, int>
globalNodeIDs(const DataBase<Dimension>& dataBase) {
  return globalNodeIDs<Dimension, typename DataBase<Dimension>::ConstNodeListIterator>(dataBase.nodeListBegin(), dataBase.nodeListEnd());
}

//------------------------------------------------------------------------------
// Compute a unique set of global node IDs for all nodes in the given set of
// NodeLists, returning the result as a FieldList<int>.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator>
// inline
FieldList<Dimension, int>
globalNodeIDs(const NodeListIterator& begin,
              const NodeListIterator& end) {

  // Prepare the result.
  const size_t numNodeLists = std::distance(begin, end);
  FieldList<Dimension, int> result(FieldStorageType::CopyFields);
  for (NodeListIterator itr = begin; itr != end; ++itr) {
    result.appendField(Field<Dimension, int>("global IDs", **itr));
  }
  CONTRACT_VAR(numNodeLists);
  CHECK(result.numFields() == numNodeLists);

#ifdef USE_MPI
  // This processors domain id.
  const int procID = Process::getRank();
  const int numProcs = Process::getTotalNumberOfProcesses();

  // Count up how many nodes are on this domain.
  int numDomainNodes = 0;
  for (NodeListIterator nodeListItr = begin; nodeListItr != end; ++nodeListItr)
    numDomainNodes += (*nodeListItr)->numInternalNodes();

  // Now loop over each processor, and have each send the cumulative number
  // on to the next.
  int beginID = 0;
  for (int sendProc = 0; sendProc < numProcs - 1; ++sendProc) {
    int sendProcDomainNodes = numDomainNodes;
    MPI_Bcast(&sendProcDomainNodes, 1, MPI_INT, sendProc, Communicator::communicator());
    if (procID > sendProc) beginID += sendProcDomainNodes;
  }
  const int endID = beginID + numDomainNodes;

  // Now assign global IDs to each node on this domain in the DataBase.
  // Loop over each NodeList in the DataBase.
  int nodeListID = 0;
  for (NodeListIterator nodeListItr = begin; nodeListItr != end; ++nodeListItr, ++nodeListID) {
    const NodeList<Dimension>& nodeList = **nodeListItr;
    Field<Dimension, int>& globalIDs = **result.fieldForNodeList(nodeList);

    // Construct the set of global IDs for this NodeList on this domain.
    CONTRACT_VAR(endID);
    CHECK(endID - beginID >= (int)nodeList.numInternalNodes());
    for (auto i = 0u; i != nodeList.numInternalNodes(); ++i) globalIDs(i) = beginID + i;
    beginID += nodeList.numInternalNodes();
  }
  CHECK(beginID == endID);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    ENSURE(result.numFields() == numNodeLists);

    // Make sure each global ID is unique.
    const int nGlobal = numGlobalNodes<Dimension, NodeListIterator>(begin, end);
    for (int checkProc = 0; checkProc != numProcs; ++checkProc) {
      for (typename FieldList<Dimension, int>::const_iterator fieldItr = result.begin();
           fieldItr != result.end();
           ++fieldItr) {
        int n = (*fieldItr)->nodeListPtr()->numInternalNodes();
        typename Field<Dimension, int>::const_iterator fieldBegin = (*fieldItr)->begin();
        typename Field<Dimension, int>::const_iterator fieldEnd = fieldBegin + n;
        MPI_Bcast(&n, 1, MPI_INT, checkProc, Communicator::communicator());
        for (int i = 0; i != n; ++i) {
          int id;
          if (procID == checkProc) id = (**fieldItr)(i);
          MPI_Bcast(&id, 1, MPI_INT, checkProc, Communicator::communicator());
          CONTRACT_VAR(nGlobal);
          ENSURE(i >= 0 && i < nGlobal);
          if (procID != checkProc){
            CONTRACT_VAR(fieldEnd);
            ENSURE(find(fieldBegin, fieldEnd, id) == fieldEnd);
          }
        }
      }
    }
  }
  END_CONTRACT_SCOPE

#else

  // Provide a serial version (pretty simple mindedly).
  int numCumulativeNodes = 0;
  for (NodeListIterator nodeListItr = begin; nodeListItr != end; ++nodeListItr) {
    const NodeList<Dimension>& nodeList = **nodeListItr;
    Field<Dimension, int>& globalIDs = **result.fieldForNodeList(nodeList);
    globalIDs = globalNodeIDs(nodeList);
    for (auto i = 0u; i != globalIDs.numElements(); ++i) globalIDs(i) += numCumulativeNodes;
    numCumulativeNodes += globalIDs.numElements();
  }

#endif

  // That's it.
  return result;
}

}
