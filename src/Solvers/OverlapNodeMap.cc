//---------------------------------Spheral++----------------------------------//
// OverlapNodeMap
//
// Node map for overlap connectivity
//----------------------------------------------------------------------------//
#include "OverlapNodeMap.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
OverlapNodeMap<Dimension>::
OverlapNodeMap(const FlatConnectivity<Dimension>& connectivity):
  NodeMap<Dimension>(connectivity) {
}

//------------------------------------------------------------------------------
// Get number of nonzero columns for this row
//------------------------------------------------------------------------------
template<typename Dimension>
int
OverlapNodeMap<Dimension>::
numElementsPerRow(int localRowIndex) const {
  return this->mConnectivity.numNonConstOverlapNeighbors(localRowIndex);
}

template<typename Dimension>
int
OverlapNodeMap<Dimension>::
numElementsPerRow(int nodeListIndex,
                  int nodeIndex) const {
  const auto localIndex = this->mConnectivity.nodeToLocal(nodeListIndex, nodeIndex);
  return numElementsPerRow(localIndex);
}

//------------------------------------------------------------------------------
// Get nonzero columns corresponding to this row
//------------------------------------------------------------------------------
template<typename Dimension>
void
OverlapNodeMap<Dimension>::
getColumnIndices(const int localRowIndex,
                 std::vector<int>& globalColumnIndices) const {
  this->mConnectivity.globalOverlapNeighborIndices(localRowIndex,
                                                   globalColumnIndices);
  return;
}

} // end namespace Spheral
