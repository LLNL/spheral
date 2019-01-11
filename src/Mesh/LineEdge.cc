//---------------------------------Spheral++----------------------------------//
// LineEdge -- 1-D mesh node class.
//
// Created by JMO, Tue Oct 12 23:07:22 PDT 2010
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include <algorithm>

namespace Spheral {

using std::abs;

//------------------------------------------------------------------------------
// LineMesh::Edge::Edge(mesh, node1ID, node2ID)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<1> >::Edge::
Edge(const Mesh<Dim<1> >& mesh,
     const unsigned ID,
     const unsigned node1ID,
     const unsigned node2ID):
  mMeshPtr(&mesh),
  mID(ID),
  mNode1ID(node1ID),
  mNode2ID(node2ID) {
  REQUIRE2(mID <= mMeshPtr->numEdges(), "Bad Edge ID:  " << mID << " " << mMeshPtr->numEdges());
  REQUIRE(node1ID < mMeshPtr->numNodes() or node1ID == UNSETID);
  REQUIRE(node2ID < mMeshPtr->numNodes() or node2ID == UNSETID);
  REQUIRE(node1ID == node2ID);
}

}
