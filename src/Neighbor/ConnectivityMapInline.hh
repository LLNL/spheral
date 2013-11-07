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
                const NodeListIterator& end):
  mNodeLists(),
  mConnectivity(FieldSpace::Copy),
  mNodeTraversalIndices(),
  mDomainDecompIndependent(false),
  mKeys(FieldSpace::Copy) {

  // Copy the set of NodeLists in the order prescribed by the NodeListRegistrar.
  NodeListRegistrar<Dimension>& registrar = NodeListRegistrar<Dimension>::instance();
  for (NodeListIterator itr = begin; itr != end; ++itr) {
    typename std::vector<const NodeSpace::NodeList<Dimension>*>::iterator posItr = registrar.findInsertionPoint(*itr,
                                                                                                                mNodeLists.begin(),
                                                                                                                mNodeLists.end());
    mNodeLists.insert(posItr, *itr);
  }

  // Copy the flag about domain decomposition independence.
  mDomainDecompIndependent = registrar.domainDecompositionIndependent();

  // The private method does the grunt work of filling in the connectivity once we have
  // established the set of NodeLists.
  computeConnectivity();

  // We'd better be valid after the constructor is finished!
  ENSURE(valid());
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
  REQUIRE(nodeID >= 0 and 
          (nodeID < nodeListPtr->numInternalNodes()) or
          (mDomainDecompIndependent and nodeID < nodeListPtr->numNodes()));
  const int nodeListID = std::distance(mNodeLists.begin(),
                                       std::find(mNodeLists.begin(), mNodeLists.end(), nodeListPtr));
  REQUIRE(nodeListID < mNodeLists.size());
  return mConnectivity(nodeListID, nodeID);
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
  REQUIRE(nodeListID >= 0 and nodeListID < mNodeLists.size());
  REQUIRE(nodeID >= 0 and 
          (nodeID < mNodeLists[nodeListID]->numInternalNodes()) or
          (mDomainDecompIndependent and nodeID < mNodeLists[nodeListID]->numNodes()));
  return mConnectivity(nodeListID, nodeID);
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

//------------------------------------------------------------------------------
// A single point to determine if in looping over nodes and neighbors the given
// pair should be calculated or not when we are doing pairs simultaneously.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
ConnectivityMap<Dimension>::
calculatePairInteraction(const int nodeListi, const int i,
                         const int nodeListj, const int j,
                         const int firstGhostNodej) const {
  if (mDomainDecompIndependent) {
    return ((nodeListj > nodeListi) or
            (nodeListj == nodeListi and 
             (mKeys(nodeListj, j) == mKeys(nodeListi, i) ? j > 1 : mKeys(nodeListj, j) > mKeys(nodeListi, i))));
  } else {
    return ((nodeListj > nodeListi) or
            (nodeListj == nodeListi and j > i) or
            (nodeListj < nodeListi and j >= firstGhostNodej));
  }
}

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
