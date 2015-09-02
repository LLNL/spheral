namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// Return the NodeList we're using.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeSpace::NodeList<Dimension>&
ConstantVelocityBoundary<Dimension>::
nodeList() const {
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// Return the set of node IDs we're controlling.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>
ConstantVelocityBoundary<Dimension>::
nodeIndices() const {
  std::vector<int> result;
  for (int i = 0; i != mNodeListPtr->numInternalNodes(); ++i) {
    if (mNodes(i) == 1) result.push_back(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the set of velocities being enforced for the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<typename Dimension::Vector>
ConstantVelocityBoundary<Dimension>::
velocityCondition() const {
  std::vector<Vector> result;
  for (int i = 0; i != mNodeListPtr->numInternalNodes(); ++i) {
    if (mNodes(i) == 1) result.push_back(mVelocity(i));
  }
  return result;
}

}
}
