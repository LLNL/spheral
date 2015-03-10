#include "Utilities/DBC.hh"

namespace Spheral {
namespace NeighborSpace {

//------------------------------------------------------------------------------
// Access the sampling kernel extent.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
Neighbor<Dimension>::kernelExtent() const {
  return mKernelExtent;
}

template<typename Dimension>
inline
void
Neighbor<Dimension>::kernelExtent(double kernelExtent) {
  CHECK(kernelExtent > 0.0);
  mKernelExtent = kernelExtent;
}

//------------------------------------------------------------------------------
// Return the current sizes of the neighboring lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
Neighbor<Dimension>::numMaster() const {
  return mMasterListPtr->size();
}

template<typename Dimension>
inline
int
Neighbor<Dimension>::numCoarse() const {
  return mCoarseNeighborListPtr->size();
}

template<typename Dimension>
inline
int
Neighbor<Dimension>::numRefine() const {
  return mRefineNeighborListPtr->size();
}

//------------------------------------------------------------------------------
// Access the current master list of nodes we have neighbors for.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<int>&
Neighbor<Dimension>::masterList() const {
  return *mMasterListPtr;
}

//------------------------------------------------------------------------------
// Access the current coarse list of potential neighbors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<int>&
Neighbor<Dimension>::coarseNeighborList() const {
  return *mCoarseNeighborListPtr;
}

//------------------------------------------------------------------------------
// Access the current refined list of potential neighbors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<int>&
Neighbor<Dimension>::refineNeighborList() const {
  return *mRefineNeighborListPtr;
}

//------------------------------------------------------------------------------
// Provide iterators over the master, coarse, and fine neighbor lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Neighbor<Dimension>::const_iterator
Neighbor<Dimension>::masterBegin() const {
  return mMasterListPtr->begin();
}

template<typename Dimension>
inline
typename Neighbor<Dimension>::const_iterator
Neighbor<Dimension>::masterEnd() const {
  return mMasterListPtr->end();
}

template<typename Dimension>
inline
typename Neighbor<Dimension>::const_iterator
Neighbor<Dimension>::coarseNeighborBegin() const {
  return mCoarseNeighborListPtr->begin();
}

template<typename Dimension>
inline
typename Neighbor<Dimension>::const_iterator
Neighbor<Dimension>::coarseNeighborEnd() const {
  return mCoarseNeighborListPtr->end();
}

template<typename Dimension>
inline
typename Neighbor<Dimension>::const_iterator
Neighbor<Dimension>::refineNeighborBegin() const {
  return mRefineNeighborListPtr->begin();
}

template<typename Dimension>
inline
typename Neighbor<Dimension>::const_iterator
Neighbor<Dimension>::refineNeighborEnd() const {
  return mRefineNeighborListPtr->end();
}

//------------------------------------------------------------------------------
// Provide read/write access to the master list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>&
Neighbor<Dimension>::accessMasterList() {
  return *mMasterListPtr;
}

//------------------------------------------------------------------------------
// Provide read/write access to the coarse neighbor node index list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>&
Neighbor<Dimension>::accessCoarseNeighborList() {
  return *mCoarseNeighborListPtr;
}

//------------------------------------------------------------------------------
// Provide read/write access to the refined neighbor node index list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<int>&
Neighbor<Dimension>::accessRefineNeighborList() {
  return *mRefineNeighborListPtr;
}

//------------------------------------------------------------------------------
// Static method for setting master/coarse on multiple NodeLists/Neighbors.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIteratorType>
void
Neighbor<Dimension>::
setMasterNeighborGroup(const typename Dimension::Vector& position,
                       const typename Dimension::SymTensor& H,
                       const NodeListIteratorType& nodeListBegin,
                       const NodeListIteratorType& nodeListEnd,
                       double kernelExtent) {

  // Loop over all the NodeLists.
  typedef typename Dimension::Vector Vector;
  Vector minMasterPosition(FLT_MAX);
  Vector maxMasterPosition(-FLT_MAX);
  Vector minMasterExtent, maxMasterExtent;
  for (NodeListIteratorType nodeListItr = nodeListBegin;
       nodeListItr != nodeListEnd;
       ++nodeListItr) {

    // Begin by setting the basic master/coarse neighbor info for each individual
    // NodeList.
    Neighbor<Dimension>& neighbor = (*nodeListItr)->neighbor();
    neighbor.setMasterList(position, H);

    // Accumulate the min/max master node positions and extents.
    const FieldSpace::Field<Dimension, Vector>& positions = (*nodeListItr)->positions();
    const FieldSpace::Field<Dimension, Vector>& extents = neighbor.nodeExtentField();
    for (typename Neighbor<Dimension>::const_iterator masterItr = 
           neighbor.masterBegin();
         masterItr < neighbor.masterEnd();
         ++masterItr) {
      const Vector& ri = positions(*masterItr);
      const Vector minExtenti = ri - extents(*masterItr);
      const Vector maxExtenti = ri + extents(*masterItr);
      for (int i = 0; i < Dimension::nDim; ++i) {
        minMasterPosition(i) = std::min(minMasterPosition(i), ri(i));
        maxMasterPosition(i) = std::max(maxMasterPosition(i), ri(i));
        minMasterExtent(i) = std::min(minMasterExtent(i), minExtenti(i));
        maxMasterExtent(i) = std::max(maxMasterExtent(i), maxExtenti(i));
      }
    }
  }

  // Don't forget to compare these positions and extents with the given 
  // position/extent. :-P
  const Vector extent = HExtent(H, kernelExtent);
  const Vector minExtenti = position - extent;
  const Vector maxExtenti = position + extent;
  for (int i = 0; i < Dimension::nDim; ++i) {
    minMasterPosition(i) = std::min(minMasterPosition(i), position(i));
    maxMasterPosition(i) = std::max(maxMasterPosition(i), position(i));
    minMasterExtent(i) = std::min(minMasterExtent(i), minExtenti(i));
    maxMasterExtent(i) = std::max(maxMasterExtent(i), maxExtenti(i));
  }

  // Loop over the nodes again, and cull the coarse neighbor lists according
  // to the overall min/max master distribution we just calcuated.
  // WARNING!  After this step the neighbor information in the NodeLists is only
  // guaranteed complete for the set of NodeLists passed to this method!
  for (NodeListIteratorType nodeListItr = nodeListBegin;
       nodeListItr != nodeListEnd;
       ++nodeListItr) {

    // Get the current set of coarse neighbors for this NodeList.
    Neighbor<Dimension>& neighbor = (*nodeListItr)->neighbor();
    std::vector<int>& coarseNeighbors = neighbor.accessCoarseNeighborList();

    // Now cull the set of coarse neighbors.
    coarseNeighbors = neighbor.precullList(minMasterPosition, maxMasterPosition,
                                           minMasterExtent, maxMasterExtent,
                                           coarseNeighbors);

//     // Set the per field coarse data caches for this NodeList.
//     (*nodeListItr)->notifyFieldsCacheCoarseValues();
  }
}

}
}
