//---------------------------------Spheral++----------------------------------//
// LineFace -- 1-D face class.
//
// Created by JMO, Thu Oct 14 21:21:02 PDT 2010
//----------------------------------------------------------------------------//
#include <algorithm>

#include "Mesh.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;

//------------------------------------------------------------------------------
// Mesh::Face(...)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<1> >::Face::
Face(const Mesh<Dim<1> >& mesh,
     const unsigned ID,
     const unsigned zone1ID,
     const unsigned zone2ID,
     const vector<unsigned>& edgeIDs):
  mMeshPtr(&mesh),
  mID(ID),
  mZone1ID(zone1ID),
  mZone2ID(zone2ID),
  mNodeIDs(std::vector<unsigned>(1, mesh.mEdges[edgeIDs[0]].node1ID())),
  mEdgeIDs(edgeIDs) {
  REQUIRE(mEdgeIDs.size() == 1);
  REQUIRE(mNodeIDs.size() == 1);
}

//------------------------------------------------------------------------------
// LineFace::position
//------------------------------------------------------------------------------
template<>
Mesh<Dim<1> >::Vector
Mesh<Dim<1> >::Face::
position() const {
  REQUIRE(mNodeIDs.size() == 1);
  return mMeshPtr->mNodePositions[mNodeIDs[0]];
}

//------------------------------------------------------------------------------
// LineFace::area
//------------------------------------------------------------------------------
template<>
double
Mesh<Dim<1> >::Face::
area() const {
  return 1.0;
}

//------------------------------------------------------------------------------
// LineFace::unitNormal
//------------------------------------------------------------------------------
template<>
Dim<1>::Vector
Mesh<Dim<1> >::Face::
unitNormal() const {
  return Vector(1.0);
}

}
}
