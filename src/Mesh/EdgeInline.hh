//---------------------------------Spheral++----------------------------------//
// Mesh::Edge
//
// Created by JMO, Fri Oct 15 23:24:11 PDT 2010
//----------------------------------------------------------------------------//
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Mesh::Edge::ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Edge::
ID() const {
  return mID;
}

//------------------------------------------------------------------------------
// Mesh::Edge::node1ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Edge::
node1ID() const {
  return mNode1ID;
}

//------------------------------------------------------------------------------
// Mesh::Edge::node2ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Edge::
node2ID() const {
  return mNode2ID;
}

//------------------------------------------------------------------------------
// Mesh::Edge::node1
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Mesh<Dimension>::Node&
Mesh<Dimension>::Edge::
node1() const {
  return mMeshPtr->node(mNode1ID);
}

//------------------------------------------------------------------------------
// Mesh::Edge::node2
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Mesh<Dimension>::Node&
Mesh<Dimension>::Edge::
node2() const {
  return mMeshPtr->node(mNode2ID);
}

//------------------------------------------------------------------------------
// Mesh::Edge::position
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::Vector
Mesh<Dimension>::Edge::
position() const {
  REQUIRE(mNode1ID < mMeshPtr->numNodes());
  REQUIRE(mNode2ID < mMeshPtr->numNodes());
  return 0.5*(mMeshPtr->mNodePositions[mNode1ID] +
              mMeshPtr->mNodePositions[mNode2ID]);
}

//------------------------------------------------------------------------------
// Mesh::Edge::length
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
Mesh<Dimension>::Edge::
length() const {
  REQUIRE(mNode1ID < mMeshPtr->numNodes());
  REQUIRE(mNode2ID < mMeshPtr->numNodes());
  return (mMeshPtr->mNodePositions[mNode2ID] -
          mMeshPtr->mNodePositions[mNode1ID]).magnitude();
}

}
