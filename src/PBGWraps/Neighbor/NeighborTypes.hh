#ifndef __PBGWRAPS_NEIGHBORTYPES__
#define __PBGWRAPS_NEIGHBORTYPES__

#include "Geometry/Dimension.hh"
#include "Neighbor/GridCellIndex.hh"
#include "Neighbor/GridCellPlane.hh"
#include "Neighbor/Neighbor.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef GridCellIndex<Dim<1> > GridCellIndex1d;
typedef GridCellIndex<Dim<2> > GridCellIndex2d;
typedef GridCellIndex<Dim<3> > GridCellIndex3d;

typedef GridCellPlane<Dim<1> > GridCellPlane1d;
typedef GridCellPlane<Dim<2> > GridCellPlane2d;
typedef GridCellPlane<Dim<3> > GridCellPlane3d;

typedef Neighbor<Dim<1> > Neighbor1d;
typedef Neighbor<Dim<2> > Neighbor2d;
typedef Neighbor<Dim<3> > Neighbor3d;

typedef NestedGridNeighbor<Dim<1> > NestedGridNeighbor1d;
typedef NestedGridNeighbor<Dim<2> > NestedGridNeighbor2d;
typedef NestedGridNeighbor<Dim<3> > NestedGridNeighbor3d;

typedef TreeNeighbor<Dim<1> > TreeNeighbor1d;
typedef TreeNeighbor<Dim<2> > TreeNeighbor2d;
typedef TreeNeighbor<Dim<3> > TreeNeighbor3d;

typedef ConnectivityMap<Dim<1> > ConnectivityMap1d;
typedef ConnectivityMap<Dim<2> > ConnectivityMap2d;
typedef ConnectivityMap<Dim<3> > ConnectivityMap3d;

}

typedef std::vector<Spheral::GridCellIndex1d> vector_of_GridCellIndex1d;
typedef std::vector<Spheral::GridCellIndex2d> vector_of_GridCellIndex2d;
typedef std::vector<Spheral::GridCellIndex3d> vector_of_GridCellIndex3d;

typedef std::vector<std::vector<Spheral::GridCellIndex1d> > vector_of_vector_of_GridCellIndex1d;
typedef std::vector<std::vector<Spheral::GridCellIndex2d> > vector_of_vector_of_GridCellIndex2d;
typedef std::vector<std::vector<Spheral::GridCellIndex3d> > vector_of_vector_of_GridCellIndex3d;

namespace Spheral {

//------------------------------------------------------------------------------
// Get a NodeList from a ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
NodeList<Dimension>*
nodeListFromConnectivityMap(ConnectivityMap<Dimension>& self, const int index) {
  return const_cast<NodeList<Dimension>*>(&(self.nodeList(index)));
}

//------------------------------------------------------------------------------
// Get a FluidNodeList from a ConnectivityMap.
//------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FluidNodeList<Dimension>*
// fluidNodeListFromConnectivityMap(ConnectivityMap<Dimension>& self, const int index) {
//   return const_cast<FluidNodeList<Dimension>*>(&(self.fluidNodeList(index)));
// }

//------------------------------------------------------------------------------
// Get the number of NodeLists from a ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
numNodeListsInConnectivityMap(ConnectivityMap<Dimension>& self) {
  return self.nodeLists().size();
}

}

#endif
