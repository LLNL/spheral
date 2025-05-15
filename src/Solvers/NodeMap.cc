//---------------------------------Spheral++----------------------------------//
// SpheralNodeMap
//
// Creates a matrix map based on the connectivity of the Spheral nodes.
// Assumes that the global indexing is sequential on one process and
// contiguous between processors. See Utilities/globalNodeIDs.hh. 
//----------------------------------------------------------------------------//
#include "NodeMap.hh"

#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "Registrar/CombinedRegistrar.hh"
#include "Solvers/containsConstantNodes.hh"
#include "Utilities/DBC.hh"
#include "Utilities/globalNodeIDs.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension> 
NodeMap<Dimension>::
NodeMap(const FlatConnectivity<Dimension>& flatConnectivity):
  mConnectivity(flatConnectivity) {
}

//------------------------------------------------------------------------------
// Check that the indices make sense.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeMap<Dimension>::
checkClassInvariants() const {
}

//------------------------------------------------------------------------------
// Get first global index belonging to this processor.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeMap<Dimension>::
firstGlobalIndex() const {
  return mConnectivity.firstGlobalIndex();
}
  
//------------------------------------------------------------------------------
// Get last global index belonging to this processor.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeMap<Dimension>::
lastGlobalIndex() const {
  return mConnectivity.lastGlobalIndex();
}

//------------------------------------------------------------------------------
// Get local number of nodes for this problem.
// Performs a summation over the NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeMap<Dimension>::
numLocalElements() const {
  return mConnectivity.numInternalNodes();
}

//------------------------------------------------------------------------------
// Get global number of nodes for this problem.
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeMap<Dimension>::
numGlobalElements() const {
  return mConnectivity.numGlobalNodes();
}

//------------------------------------------------------------------------------
// Get number of nonzero columns for this row
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeMap<Dimension>::
numElementsPerRow(int localRowIndex) const {
  return mConnectivity.numNonConstNeighbors(localRowIndex);
}

template<typename Dimension>
int
NodeMap<Dimension>::
numElementsPerRow(int nodeListIndex,
                  int nodeIndex) const {
  const auto localIndex = mConnectivity.nodeToLocal(nodeListIndex, nodeIndex);
  return numElementsPerRow(localIndex);
}

//------------------------------------------------------------------------------
// Get nonzero columns corresponding to this row
//------------------------------------------------------------------------------
template<typename Dimension>
void
NodeMap<Dimension>::
getColumnIndices(const int localRowIndex,
                 std::vector<int>& globalColumnIndices) const {
  mConnectivity.globalNeighborIndices(localRowIndex,
                                      globalColumnIndices);
  return;
}

//------------------------------------------------------------------------------
// Get the local matrix index for a given node index. 
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeMap<Dimension>::
getLocalIndex(int nodeListIndex,
              int nodeIndex) const {
  return mConnectivity.nodeToLocal(nodeListIndex, nodeIndex);
}

//------------------------------------------------------------------------------
// Get the node index for a given local index
//------------------------------------------------------------------------------
template<typename Dimension>
std::pair<int, int>
NodeMap<Dimension>::
getNodeIndex(int localIndex) const {
  return mConnectivity.localToNode(localIndex);
}

//------------------------------------------------------------------------------
// Convert local indices to global indices
//------------------------------------------------------------------------------
template<typename Dimension>
int
NodeMap<Dimension>::
getGlobalIndex(int localRowIndex) const {
  return mConnectivity.localToGlobal(localRowIndex);
}

template<typename Dimension>
int
NodeMap<Dimension>::
getGlobalIndex(int nodeListIndex,
               int nodeIndex) const {
  const auto localIndex = mConnectivity.nodeToLocal(nodeListIndex, nodeIndex);
  return mConnectivity.localToGlobal(localIndex);
}
//------------------------------------------------------------------------------
// Return whether node is a constant boundary
//------------------------------------------------------------------------------
template<typename Dimension>
bool 
NodeMap<Dimension>::
isConstantBoundaryNode(int nodeListIndex,
                       int nodeIndex) const {
  const auto localIndex = mConnectivity.nodeToLocal(nodeListIndex, nodeIndex);
  return mConnectivity.isConstantBoundaryNode(localIndex);
}

} // end namespace Spheral
