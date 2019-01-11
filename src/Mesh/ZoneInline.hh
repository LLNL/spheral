//---------------------------------Spheral++----------------------------------//
// Mesh::Zone
//
// Created by JMO, Fri Oct 15 23:24:11 PDT 2010
//----------------------------------------------------------------------------//
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Mesh::Zone::ID
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Zone::
ID() const {
  return mID;
}

//------------------------------------------------------------------------------
// Mesh::Zone::numNodes
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Zone::
numNodes() const {
  return mNodeIDs.size();
}

//------------------------------------------------------------------------------
// Mesh::Zone::numEdges
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Zone::
numEdges() const {
  return mEdgeIDs.size();
}

//------------------------------------------------------------------------------
// Mesh::Zone::numFaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::Zone::
numFaces() const {
  return mFaceIDs.size();
}

//------------------------------------------------------------------------------
// LineZone::nodeIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<unsigned>&
Mesh<Dimension>::Zone::
nodeIDs() const {
  return mNodeIDs;
}

//------------------------------------------------------------------------------
// LineZone::edgeIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<unsigned>&
Mesh<Dimension>::Zone::
edgeIDs() const {
  return mEdgeIDs;
}

//------------------------------------------------------------------------------
// LineZone::faceIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<int>&
Mesh<Dimension>::Zone::
faceIDs() const {
  return mFaceIDs;
}

//------------------------------------------------------------------------------
// LineZone::convexHull
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::ConvexHull
Mesh<Dimension>::Zone::
convexHull() const {
  std::vector<Vector> points;
  for (std::vector<unsigned>::const_iterator itr = mNodeIDs.begin();
       itr != mNodeIDs.end();
       ++itr) points.push_back(mMeshPtr->mNodePositions[*itr]);
  return ConvexHull(points);
}

}
