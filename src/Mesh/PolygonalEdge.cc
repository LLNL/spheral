//---------------------------------Spheral++----------------------------------//
// PolygonalEdge -- 2-D mesh node class.
//
// Created by JMO, Tue Nov 16 15:00:14 PST 2010
//----------------------------------------------------------------------------//
#include "Mesh.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// PolygonalMesh::Edge::Edge(mesh, node1ID, node2ID)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<2> >::Edge::
Edge(const Mesh<Dim<2> >& mesh,
     const unsigned ID,
     const unsigned node1ID,
     const unsigned node2ID):
  mMeshPtr(&mesh),
  mID(ID),
  mNode1ID(node1ID),
  mNode2ID(node2ID) {
  REQUIRE2(mID <= mMeshPtr->numEdges(), "Bad Edge ID:  " << mID << " " << mMeshPtr->numEdges());
  REQUIRE2(node1ID < mMeshPtr->mNodePositions.size() or node1ID == UNSETID, "Bad node ID:  " << node1ID << " " << mMeshPtr->mNodePositions.size());
  REQUIRE2(node2ID < mMeshPtr->mNodePositions.size() or node2ID == UNSETID, "Bad node ID:  " << node2ID << " " << mMeshPtr->mNodePositions.size());
  REQUIRE2(node1ID != node2ID, "Ack!  " << node1ID << " " << node2ID);
}

}
