#include "Neighbor/NestedGridNeighbor.hh"
#include "Geometry/Dimension.hh"

using namespace Spheral::NeighborSpace;

//------------------------------------------------------------------------------
// Get the grid level for the given node ID.
//------------------------------------------------------------------------------
inline
int
gridLevelForNode1d(const NestedGridNeighbor<Spheral::Dim<1> >* self,
                   const int nodeID) {
  return self->gridLevel(nodeID);
}

inline
int
gridLevelForNode2d(const NestedGridNeighbor<Spheral::Dim<2> >* self,
                   const int nodeID) {
  return self->gridLevel(nodeID);
}

inline
int
gridLevelForNode3d(const NestedGridNeighbor<Spheral::Dim<3> >* self,
                   const int nodeID) {
  return self->gridLevel(nodeID);
}

