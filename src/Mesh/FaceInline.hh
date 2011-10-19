//---------------------------------Spheral++----------------------------------//
// Mesh::Face
//
// Created by JMO, Fri Oct 15 23:24:11 PDT 2010
//----------------------------------------------------------------------------//
#include "Utilities/DBC.hh"

namespace Spheral {
namespace MeshSpace {

//------------------------------------------------------------------------------
// Mesh::Face::ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Face::
ID() const {
  return mID;
}

//------------------------------------------------------------------------------
// Mesh::Face::numNodes
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Face::
numNodes() const {
  return mNodeIDs.size();
}

//------------------------------------------------------------------------------
// Mesh::Face::numEdges
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Face::
numEdges() const {
  return mEdgeIDs.size();
}

//------------------------------------------------------------------------------
// Face::nodeIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<unsigned>&
Mesh<Dimension>::Face::
nodeIDs() const {
  return mNodeIDs;
}

//------------------------------------------------------------------------------
// Face::edgeIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<unsigned>&
Mesh<Dimension>::Face::
edgeIDs() const {
  return mEdgeIDs;
}

//------------------------------------------------------------------------------
// Mesh::Face::zone1ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Face::
zone1ID() const {
  return mZone1ID;
}

//------------------------------------------------------------------------------
// Mesh::Face::zone2ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Face::
zone2ID() const {
  return mZone2ID;
}

//------------------------------------------------------------------------------
// Face::oppositeZoneID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Face::
oppositeZoneID(const unsigned zoneID) const {
  REQUIRE(zoneID == mZone1ID or zoneID == mZone2ID);
  return zoneID == mZone1ID ? mZone2ID : mZone1ID;
}

}
}
