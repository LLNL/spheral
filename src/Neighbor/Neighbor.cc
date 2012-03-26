//---------------------------------Spheral++----------------------------------//
// Neighbor -- Abstract interface base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
//----------------------------------------------------------------------------//
#include "Neighbor.hh"

#include "DBC.hh"
#include "Field/Field.hh"
#include "NodeList/NodeList.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/testBoxIntersection.hh"

namespace Spheral {
namespace NeighborSpace {

using namespace std;

using FieldSpace::Field;
using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Construct with the given NodeList and search type.
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>::
Neighbor(NodeList<Dimension>& nodeList,
         const NeighborSearchType searchType,
         const double kernelExtent):
  mSearchType(searchType),
  mKernelExtent(kernelExtent),
  mMasterListPtr(new vector<int>()),
  mCoarseNeighborListPtr(new vector<int>()),
  mRefineNeighborListPtr(new vector<int>()), // (mCoarseNeighborListPtr),
  mNodeListPtr(&nodeList),
  mNodeExtent("Node Extent", nodeList) {
  mNodeListPtr->registerNeighbor(*this);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Neighbor<Dimension>::
~Neighbor() {
  if (mNodeListPtr) mNodeListPtr->unregisterNeighbor();
  delete mMasterListPtr;
  delete mCoarseNeighborListPtr;
  delete mRefineNeighborListPtr;
}

//------------------------------------------------------------------------------
// Get the type of search.
//------------------------------------------------------------------------------
template<typename Dimension>
NeighborSearchType
Neighbor<Dimension>::
neighborSearchType() const {
  return mSearchType;
}

//------------------------------------------------------------------------------
// Set the type of search.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
neighborSearchType(NeighborSearchType searchType) {
  mSearchType = searchType;
}

//------------------------------------------------------------------------------
// Calculate the maximum radial extent of a given smoothing tensor.
//------------------------------------------------------------------------------
// Make the general case the SPH approximation: assume the argument is a 
// scalar representing the inverse smoothing scale.
template<typename Dimension>
typename Dimension::Vector
Neighbor<Dimension>::
HExtent(const typename Dimension::Scalar& H, 
        const double kernelExtent) {
  CHECK(H > 0.0);
  const double r = kernelExtent/H;
  return Vector(r);
}

// Specializations for symmetric tensor smoothing transformation (ASPH).
// The ASPH tensor has units of inverse length.
template<>
Dim<1>::Vector
Neighbor< Dim<1> >::
HExtent(const Dim<1>::SymTensor& H,
        const double kernelExtent) {
  CHECK(H.Determinant() > 0.0);
  const double r = kernelExtent/H.xx();
  return Vector(r);
}

template<>
Dim<2>::Vector
Neighbor< Dim<2> >::
HExtent(const Dim<2>::SymTensor& H,
        const double kernelExtent) {
  const double Hdet = H.Determinant();
  const SymTensor M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Vector(sqrt(M.yy()), sqrt(M.xx()));
}
  
template<>
Dim<3>::Vector
Neighbor< Dim<3> >::
HExtent(const Dim<3>::SymTensor& H,
        const double kernelExtent) {
  const double Hdet = H.Determinant();
  const SymTensor M = H.square();
  CHECK(Hdet > 0.0);
  return kernelExtent/Hdet*Vector(sqrt(M.yy()*M.zz() - M.yz()*M.zy()),
                        sqrt(M.xx()*M.zz() - M.xz()*M.zx()),
                        sqrt(M.xx()*M.yy() - M.xy()*M.yx()));
}

//------------------------------------------------------------------------------
// Return the field of node extents.
//------------------------------------------------------------------------------
template<typename Dimension>
const Field<Dimension, typename Dimension::Vector>&
Neighbor<Dimension>::
nodeExtentField() const {
  return mNodeExtent;
}

// Allow read/write access to the node extent Field for descendent classes.
template<typename Dimension>
Field<Dimension, typename Dimension::Vector>&
Neighbor<Dimension>::
accessNodeExtentField() {
  return mNodeExtent;
}

//------------------------------------------------------------------------------
// Access the node list.
//------------------------------------------------------------------------------
template<typename Dimension>
const NodeList<Dimension>&
Neighbor<Dimension>::
nodeList() const {
  CHECK(mNodeListPtr);
  return *mNodeListPtr;
}

template<typename Dimension>
const NodeList<Dimension>*
Neighbor<Dimension>::
nodeListPtr() const {
  return mNodeListPtr;
}

//------------------------------------------------------------------------------
// Set the node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
nodeList(NodeList<Dimension>& nodeList) {
  CHECK(&nodeList);
  mNodeListPtr = &nodeList;
  mNodeExtent.setNodeList(nodeList);
}

template<typename Dimension>
void
Neighbor<Dimension>::
nodeListPtr(NodeList<Dimension>* nodeListPtr) {
  CHECK(nodeListPtr);
  mNodeListPtr = nodeListPtr;
  mNodeExtent.setNodeList(*nodeListPtr);
}

//------------------------------------------------------------------------------
// Unregister the node list.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
unregisterNodeList() {
  mNodeListPtr = 0;
}

//------------------------------------------------------------------------------
// Determine the extent of the given node's smoothing tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Vector
Neighbor<Dimension>::
nodeExtent(int nodeID) const {
  CHECK(nodeID >= 0 and nodeID < nodeList().numNodes());
  return HExtent(nodeList().Hfield()(nodeID), kernelExtent());
}

//------------------------------------------------------------------------------
// Force the node extent field to be computed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
setNodeExtents() {
  for (int nodeID = 0; nodeID < nodeList().numNodes(); ++nodeID) {
    mNodeExtent(nodeID) = nodeExtent(nodeID);
  }
}

template<typename Dimension>
void
Neighbor<Dimension>::
setNodeExtents(const vector<int>& nodeIDs) {
  for (typename vector<int>::const_iterator nodeIDItr = nodeIDs.begin();
       nodeIDItr < nodeIDs.end();
       ++nodeIDItr) {
    CHECK(*nodeIDItr >= 0 and *nodeIDItr < nodeList().numNodes());
    mNodeExtent(*nodeIDItr) = nodeExtent(*nodeIDItr);
  }
}

template<typename Dimension>
void
Neighbor<Dimension>::
setInternalNodeExtents() {
  for (int nodeID = 0; nodeID < nodeList().numInternalNodes(); ++nodeID) {
    mNodeExtent(nodeID) = nodeExtent(nodeID);
  }
}

template<typename Dimension>
void
Neighbor<Dimension>::
setGhostNodeExtents() {
  for (int nodeID = nodeList().firstGhostNode(); nodeID < nodeList().numNodes(); ++nodeID) {
    mNodeExtent(nodeID) = nodeExtent(nodeID);
  }
}

//------------------------------------------------------------------------------
// Set the master list of nodes based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
setMasterList(int nodeID) {
  CHECK(valid());
  CHECK(nodeID >= 0 and nodeID < nodeList().numInternalNodes());
  vector<int>& masterList = accessMasterList();
  masterList = vector<int>();
  const Field<Dimension, Vector>& positions = nodeList().positions();
  const Field<Dimension, SymTensor>& Hfield = nodeList().Hfield();
  setMasterList(positions(nodeID), Hfield(nodeID));
}

//------------------------------------------------------------------------------
// Set the refined list of potential neighbors based on a particular node index.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
setRefineNeighborList(int nodeID) {
  CHECK(valid());
  CHECK(nodeID >= 0 and nodeID < nodeList().numInternalNodes());
  CHECK(find(accessMasterList().begin(), accessMasterList().end(), nodeID) !=
         accessMasterList().end());
  const Field<Dimension, Vector>& positions = nodeList().positions();
  const Field<Dimension, SymTensor>& Hfield = nodeList().Hfield();
  setRefineNeighborList(positions(nodeID), Hfield(nodeID));
}

//------------------------------------------------------------------------------
// General precull routine.  Select nodes based on a range of positions and
// extents.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<int>
Neighbor<Dimension>::
precullList(const Vector& minMasterPosition, const Vector& maxMasterPosition,
            const Vector& minMasterExtent, const Vector& maxMasterExtent,
            const vector<int>& coarseList) const {

  // Empty vector to accumulate the result in.
  vector<int> cullList;
  cullList.reserve(coarseList.size());

  // Get a reference to the node positions and node extent field.
  const Field<Dimension, Vector>& positions = nodeList().positions();
  const Field<Dimension, Vector>& nodeExtents = nodeExtentField();

  // What kind of preculling we're doing determines the applied test.
  if (neighborSearchType() == GatherScatter) {
    
    // Gather-Scatter.
    for (typename vector<int>::const_iterator coarseItr = coarseList.begin();
         coarseItr != coarseList.end();
         ++coarseItr) {
      const int j = *coarseItr;
      const Vector& nodePosition = positions(j);
      const Vector minNodeExtent = nodePosition - nodeExtents(j);
      const Vector maxNodeExtent = nodePosition + nodeExtents(j);
      const bool gatherTest = testBoxIntersection(nodePosition,
                                                  nodePosition,
                                                  minMasterExtent,
                                                  maxMasterExtent);
      const bool scatterTest = testBoxIntersection(minMasterPosition,
                                                   maxMasterPosition,
                                                   minNodeExtent,
                                                   maxNodeExtent);
      if (gatherTest or scatterTest) cullList.push_back(j);
    }

  } else if (neighborSearchType() == Gather) {

    // Gather.
    for (typename vector<int>::const_iterator coarseItr = coarseList.begin();
         coarseItr != coarseList.end();
         ++coarseItr) {
      const int j = *coarseItr;
      const Vector& nodePosition = positions(j);
      const bool gatherTest = testBoxIntersection(nodePosition,
                                                  nodePosition,
                                                  minMasterExtent,
                                                  maxMasterExtent);
      if (gatherTest) cullList.push_back(j);
    }

  } else {

    // Scatter.
    CHECK(neighborSearchType() == Scatter);
    for (typename vector<int>::const_iterator coarseItr = coarseList.begin();
         coarseItr != coarseList.end();
         ++coarseItr) {
      const int j = *coarseItr;
      const Vector& nodePosition = positions(j);
      const Vector minNodeExtent = nodePosition - nodeExtents(j);
      const Vector maxNodeExtent = nodePosition + nodeExtents(j);
      const bool scatterTest = testBoxIntersection(minMasterPosition,
                                                   maxMasterPosition,
                                                   minNodeExtent,
                                                   maxNodeExtent);
      if (scatterTest) cullList.push_back(j);
    }

  }

  ENSURE(cullList.size() <= coarseList.size());
  return cullList;
}

//------------------------------------------------------------------------------
// Cull the local (to this NodeList) neighbor info based on the current master
// state.
// *NOTE* -- this is not safe to do when you want to use this neighbor info 
// with different NodeLists!
//------------------------------------------------------------------------------
template<typename Dimension>
void
Neighbor<Dimension>::
precullForLocalNodeList() {

  // Grab the state.
  const Field<Dimension, Vector>& r = this->nodeList().positions();
  const Field<Dimension, Vector>& extent = nodeExtentField();

  // Find the min/max master node positions and extents.
  Vector minMasterPosition = DBL_MAX;
  Vector maxMasterPosition = -DBL_MAX;
  Vector minMasterExtent = DBL_MAX;
  Vector maxMasterExtent = -DBL_MAX;
  for (const_iterator masterItr = masterBegin();
       masterItr != masterEnd();
       ++masterItr) {
    const int i = *masterItr;
    const Vector& ri = r(i);
    const Vector minExtent = ri - extent(i);
    const Vector maxExtent = ri + extent(i);
    for (int j = 0; j != Dimension::nDim; ++j) {
      minMasterPosition(j) = min(minMasterPosition(j), ri(j));
      maxMasterPosition(j) = max(maxMasterPosition(j), ri(j));
      minMasterExtent(j) = min(minMasterExtent(j), minExtent(j));
      maxMasterExtent(j) = max(maxMasterExtent(j), maxExtent(j));
    }
  }

  // Now use this info to precull the coarse neighbor list.
  *mCoarseNeighborListPtr = precullList(minMasterPosition, maxMasterPosition,
                                        minMasterExtent, maxMasterExtent,
                                        *mCoarseNeighborListPtr);
}

//------------------------------------------------------------------------------
// Provide a basic test to determine whether the neighbor base is in a valid
// (i.e., ready to use) state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Neighbor<Dimension>::
valid() const {
  return (kernelExtent() > 0.0 and
          neighborSearchType() != None);
}
}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  namespace NeighborSpace {
    template class Neighbor< Dim<1> >;
    template class Neighbor< Dim<2> >;
    template class Neighbor< Dim<3> >;
  }
}
