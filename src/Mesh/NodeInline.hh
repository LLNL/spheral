//---------------------------------Spheral++----------------------------------//
// Mesh::Node
//
// Created by JMO, Fri Oct 15 23:24:11 PDT 2010
//----------------------------------------------------------------------------//
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Mesh::Node::Node(mesh, nodeID)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Mesh<Dimension>::Node::
Node(const Mesh<Dimension>& mesh,
     const unsigned ID,
     const std::vector<unsigned>& zoneIDs):
  mMeshPtr(&mesh),
  mID(ID),
  mZoneIDs(zoneIDs) {
  REQUIRE(mID < mMeshPtr->mNodePositions.size());
}

//------------------------------------------------------------------------------
// Mesh::Node::ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Node::
ID() const {
  return mID;
}

//------------------------------------------------------------------------------
// Mesh::Node::position
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::Vector
Mesh<Dimension>::Node::
position() const {
  REQUIRE(mID < mMeshPtr->numNodes());
  return mMeshPtr->mNodePositions[mID];
}

//------------------------------------------------------------------------------
// Mesh::Node::zoneIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<unsigned>&
Mesh<Dimension>::Node::
zoneIDs() const {
  return mZoneIDs;
}

}
