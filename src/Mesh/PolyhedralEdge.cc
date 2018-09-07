//---------------------------------Spheral++----------------------------------//
// PolyhedralEdge -- 2-D mesh node class.
//
// Created by JMO, Wed Jan  5 21:09:11 PST 2011
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// PolyhedralMesh::Edge::Edge(mesh, node1ID, node2ID)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<3> >::Edge::
Edge(const Mesh<Dim<3> >& mesh,
     const unsigned ID,
     const unsigned node1ID,
     const unsigned node2ID):
  mMeshPtr(&mesh),
  mID(ID),
  mNode1ID(node1ID),
  mNode2ID(node2ID) {
  REQUIRE2(mID <= mMeshPtr->numEdges(), "Bad Edge ID:  " << mID << " " << mMeshPtr->numEdges());
  REQUIRE2(node1ID < mMeshPtr->numNodes() or node1ID == UNSETID, "Bad node ID:  " << node1ID << " " << mMeshPtr->numNodes());
  REQUIRE2(node2ID < mMeshPtr->numNodes() or node2ID == UNSETID, "Bad node ID:  " << node2ID << " " << mMeshPtr->numNodes());
  REQUIRE2(node1ID != node2ID, "Edge can't have same node twice:  " << node1ID << " " << node2ID);
}

}
