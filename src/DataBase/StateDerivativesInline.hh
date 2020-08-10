namespace Spheral {

//------------------------------------------------------------------------------
// Test if the given pair of nodes is flagged as being calculated or not.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
StateDerivatives<Dimension>::
nodePairCalculated(const NodeIteratorBase<Dimension>& node1,
                   const NodeIteratorBase<Dimension>& node2) const {
  if (node1.ghostNode() or node2.ghostNode()) return false;
  typename CalculatedPairType::const_iterator itr1 = mCalculatedNodePairs.find(node1);
  CHECK((itr1 == mCalculatedNodePairs.end() && itr1 == mCalculatedNodePairs.end()) or
        (itr1 != mCalculatedNodePairs.end() && itr1 != mCalculatedNodePairs.end()));
  if (itr1 == mCalculatedNodePairs.end()) {
    return false;
  } else {
    CHECK(itr1 != mCalculatedNodePairs.end());
    return std::find(itr1->second.begin(), itr1->second.end(), node2) != itr1->second.end();
  }
}

//------------------------------------------------------------------------------
// Flag the given set of nodes as being calculated.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
StateDerivatives<Dimension>::
flagNodePairCalculated(const NodeIteratorBase<Dimension>& node1,
                       const NodeIteratorBase<Dimension>& node2) {
  REQUIRE(node1.internalNode() && node2.internalNode());

  // Have these nodes already been added to the maps?
  // If these are new nodes, create entries for them.
  typename CalculatedPairType::const_iterator itr1 = mCalculatedNodePairs.find(const_cast<NodeIteratorBase<Dimension>&>(node1));
  typename CalculatedPairType::const_iterator itr2 = mCalculatedNodePairs.find(const_cast<NodeIteratorBase<Dimension>&>(node2));
  if (itr1 == mCalculatedNodePairs.end()) mCalculatedNodePairs[node1] = std::vector<NodeIteratorBase<Dimension> >();
  if (itr2 == mCalculatedNodePairs.end()) mCalculatedNodePairs[node2] = std::vector<NodeIteratorBase<Dimension> >();

  // Now put them in each others interaction list.
  mCalculatedNodePairs[node1].push_back(node2);
  if (node1 != node2) mCalculatedNodePairs[node2].push_back(node1);

  ENSURE(std::count(mCalculatedNodePairs[node1].begin(),
                    mCalculatedNodePairs[node1].end(),
                    node2) == 1);
  ENSURE(std::count(mCalculatedNodePairs[node2].begin(),
                    mCalculatedNodePairs[node2].end(),
                    node1) == 1);
}

//------------------------------------------------------------------------------
// Return the current number of significant neighbors for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
StateDerivatives<Dimension>::
numSignificantNeighbors(const NodeIteratorBase<Dimension>& node) const {
  REQUIRE(mNumSignificantNeighbors.find(node) != mNumSignificantNeighbors.end());
  return mNumSignificantNeighbors.find(node)->second;
}

//------------------------------------------------------------------------------
// Increment the number of significant neighbors for the given node.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
StateDerivatives<Dimension>::
incrementSignificantNeighbors(const NodeIteratorBase<Dimension>& node) {
  // Have we seen this node before?
  if (mNumSignificantNeighbors.find(node) == mNumSignificantNeighbors.end()) 
    mNumSignificantNeighbors[node] = 0;
  ++(mNumSignificantNeighbors[node]);
}

}
