#include "NodeList/FluidNodeList.hh"
#include "NodeList/NodeListRegistrar.hh"

using Spheral::FieldSpace::FieldList;

namespace Spheral {
namespace NeighborSpace {

//------------------------------------------------------------------------------
// Constructor.
// The input iterators must dereference to const NodeList*, and must be 
// const FluidNodeLists underneath.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIterator>
inline
ConnectivityMap<Dimension>::
ConnectivityMap(const NodeListIterator& begin,
                const NodeListIterator& end,
                const bool buildGhostConnectivity):
  mNodeLists(),
  mBuildGhostConnectivity(buildGhostConnectivity),
  mOffsets(),
  mConnectivity(),
  mNodeTraversalIndices(),
  mKeys(FieldSpace::FieldStorageType::CopyFields) {

  // The private method does the grunt work of filling in the connectivity once we have
  // established the set of NodeLists.
  this->rebuild(begin, end, buildGhostConnectivity);

  // We'd better be valid after the constructor is finished!
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Rebuild for a given set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIterator>
inline
void
ConnectivityMap<Dimension>::
rebuild(const NodeListIterator& begin,
        const NodeListIterator& end, 
        const bool buildGhostConnectivity) {
  mBuildGhostConnectivity = buildGhostConnectivity;

  // Copy the set of NodeLists in the order prescribed by the NodeListRegistrar.
  NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
  const bool domainDecompIndependent = registrar.domainDecompositionIndependent();
  const unsigned numNodeLists = std::distance(begin, end);
  mNodeLists.clear();
  mOffsets.resize(numNodeLists);
  std::vector<unsigned> numNodes(numNodeLists);
  for (NodeListIterator itr = begin; itr != end; ++itr) {
    typename std::vector<const NodeSpace::NodeList<Dimension>*>::iterator posItr = 
      registrar.findInsertionPoint(*itr, mNodeLists.begin(), mNodeLists.end());
    const unsigned i = std::distance(mNodeLists.begin(), posItr);
    CHECK(i < numNodeLists);
    mNodeLists.insert(posItr, *itr);
    numNodes[i] = ((domainDecompIndependent or mBuildGhostConnectivity) ?
                   (*itr)->numNodes() :
                   (*itr)->numInternalNodes());
  }

  // Construct the offsets.
  mOffsets[0] = 0;
  for (unsigned i = 1; i != numNodeLists; ++i) mOffsets[i] = mOffsets[i - 1] + numNodes[i - 1];
  CHECK(mOffsets.size() == mNodeLists.size());

  this->computeConnectivity();
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Are we computing ghost connectivity?
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
buildGhostConnectivity() const {
  return mBuildGhostConnectivity;
}

//------------------------------------------------------------------------------
// Get the set of NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<const NodeSpace::NodeList<Dimension>*>&
ConnectivityMap<Dimension>::
nodeLists() const {
  return mNodeLists;
}

//------------------------------------------------------------------------------
// Get the set of neighbors for the given node in the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector<int> >&
ConnectivityMap<Dimension>::
connectivityForNode(const NodeSpace::NodeList<Dimension>* nodeListPtr,
                    const int nodeID) const {
  const bool ghostValid = (mBuildGhostConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  REQUIRE(nodeID >= 0 and 
          (nodeID < nodeListPtr->numInternalNodes()) or
          (ghostValid and nodeID < nodeListPtr->numNodes()));
  const int nodeListID = std::distance(mNodeLists.begin(),
                                       std::find(mNodeLists.begin(), mNodeLists.end(), nodeListPtr));
  REQUIRE(nodeListID < mConnectivity.size() and nodeListID < mOffsets.size());
  REQUIRE(mOffsets[nodeListID] + nodeID < mConnectivity.size());
  return mConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Same as above, indicating the NodeList as an integer index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector< std::vector<int> >&
ConnectivityMap<Dimension>::
connectivityForNode(const int nodeListID,
                    const int nodeID) const {
  const bool ghostValid = (mBuildGhostConnectivity or
                           NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
  REQUIRE(nodeListID >= 0 and nodeListID < mConnectivity.size());
  REQUIRE(nodeID >= 0 and 
          (nodeID < mNodeLists[nodeListID]->numInternalNodes()) or
          (ghostValid and nodeID < mNodeLists[nodeListID]->numNodes()));
  REQUIRE(mOffsets[nodeListID] + nodeID < mConnectivity.size());
  return mConnectivity[mOffsets[nodeListID] + nodeID];
}

//------------------------------------------------------------------------------
// Compute the number of neighbors for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
ConnectivityMap<Dimension>::
numNeighborsForNode(const NodeSpace::NodeList<Dimension>* nodeListPtr,
                    const int nodeID) const {
  const std::vector< std::vector<int> >& neighbors = connectivityForNode(nodeListPtr, nodeID);
  int result = 0;
  for (std::vector< std::vector<int> >::const_iterator itr = neighbors.begin();
       itr != neighbors.end();
       ++itr) result += itr->size();
  return result;
}

template<typename Dimension>
inline
int
ConnectivityMap<Dimension>::
numNeighborsForNode(const int nodeListID,
                    const int nodeID) const {
  REQUIRE(nodeListID < mNodeLists.size());
  return this->numNeighborsForNode(mNodeLists[nodeListID], nodeID);
}

//------------------------------------------------------------------------------
// A single point to determine if in looping over nodes and neighbors the given
// pair should be calculated or not when we are doing pairs simultaneously.
//------------------------------------------------------------------------------

//#pragma omp declare target
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
calculatePairInteraction(const int nodeListi, const int i,
                         const int nodeListj, const int j,
                         const int firstGhostNodej) const {
//  const bool domainDecompIndependent = NodeListRegistrar<Dimension>::getInstance().domainDecompositionIndependent();
  const bool domainDecompIndependent = true;
  if (domainDecompIndependent) {
    return ((nodeListj > nodeListi) or
            (nodeListj == nodeListi and 
             (mKeys(nodeListj, j) == mKeys(nodeListi, i) ? j > 1 : mKeys(nodeListj, j) > mKeys(nodeListi, i))));
  } else {
    return ((nodeListj > nodeListi) or
            (nodeListj == nodeListi and j > i) or
            (nodeListj < nodeListi and j >= firstGhostNodej));
  }
}
//#pragma omp end declare target

//------------------------------------------------------------------------------
// Iterators for walking nodes in a prescribed order (domain decomposition 
// independent when needed).
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename ConnectivityMap<Dimension>::const_iterator
ConnectivityMap<Dimension>::
begin(const int nodeList) const {
  REQUIRE(nodeList >= 0 and nodeList < mNodeTraversalIndices.size());
  return mNodeTraversalIndices[nodeList].begin();
}

template<typename Dimension>
inline
typename ConnectivityMap<Dimension>::const_iterator
ConnectivityMap<Dimension>::
end(const int nodeList) const {
  REQUIRE(nodeList >= 0 and nodeList < mNodeTraversalIndices.size());
  return mNodeTraversalIndices[nodeList].end();
}

//------------------------------------------------------------------------------
// Get the ith NodeList/FluidNodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::NodeList<Dimension>&
ConnectivityMap<Dimension>::
nodeList(const int index) const {
  REQUIRE(index >= 0 and index < mNodeLists.size());
  return *(mNodeLists[index]);
}

}
}
