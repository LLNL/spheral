//---------------------------------Spheral++----------------------------------//
// PolyhedralZone -- 3-D zone class.
//
// Created by JMO, Thu Jan  6 09:12:50 PST 2011
//----------------------------------------------------------------------------//
#include <vector>
#include <algorithm>

#include "Mesh.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Comparator to help sorting the faces counter-clockwise about a position.
//------------------------------------------------------------------------------
template<typename Element, typename Vector>
struct CounterClockwiseCompareElements {
  CounterClockwiseCompareElements(const vector<Element>& elements,
                                  const unsigned originID):
    mElements(elements),
    mCentroid(elements[originID].position()) {}
  bool operator()(const unsigned id1, const unsigned id2) {
    REQUIRE(id1 < mElements.size());
    REQUIRE(id2 < mElements.size());
    return ((mElements[id1].position() - mCentroid).cross(mElements[id2].position() - mCentroid).z() > 0.0);
  }
  const vector<Element>& mElements;
  Vector mCentroid;
};

//------------------------------------------------------------------------------
// Mesh::Zone(...)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<3> >::Zone::
Zone(const Mesh<Dim<3> >& mesh,
     const unsigned ID,
     const vector<unsigned>& faceIDs):
  mMeshPtr(&mesh),
  mID(ID),
  mNodeIDs(),
  mEdgeIDs(),
  mFaceIDs(faceIDs) {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(mFaceIDs.size() > 3);
    for (vector<unsigned>::const_iterator itr = mFaceIDs.begin();
         itr != mFaceIDs.end();
         ++itr) REQUIRE(*itr < mMeshPtr->mFaces.size());
  }
  END_CONTRACT_SCOPE;
  
  // Construct the edge and node IDs.
  for (vector<unsigned>::const_iterator faceItr = mFaceIDs.begin();
       faceItr != mFaceIDs.end();
       ++faceItr) {
    const Face& face = mMeshPtr->mFaces[*faceItr];
    const vector<unsigned>& edgeIDs = face.edgeIDs();
    const vector<unsigned>& nodeIDs = face.nodeIDs();
    CHECK(edgeIDs.size() > 2);
    CHECK(nodeIDs.size() == edgeIDs.size());
    copy(edgeIDs.begin(), edgeIDs.end(), back_inserter(mEdgeIDs));
    copy(nodeIDs.begin(), nodeIDs.end(), back_inserter(mNodeIDs));
  }
  sort(mEdgeIDs.begin(), mEdgeIDs.end());
  mEdgeIDs.erase(unique(mEdgeIDs.begin(), mEdgeIDs.end()), mEdgeIDs.end());
  sort(mNodeIDs.begin(), mNodeIDs.end());
  mNodeIDs.erase(unique(mNodeIDs.begin(), mNodeIDs.end()), mNodeIDs.end());

  // Post-conditions.
  ENSURE(mNodeIDs.size() > 3);  // At least a tet.
  ENSURE(mEdgeIDs.size() > 5);  // At least a tet.
  ENSURE(mFaceIDs.size() > 3);  // At least a tet.
}

//------------------------------------------------------------------------------
// PolyhedralZone::position
//------------------------------------------------------------------------------
template<>
Dim<3>::Vector
Mesh<Dim<3> >::Zone::
position() const {
  Vector result;
  for (unsigned i = 0; i != mFaceIDs.size(); ++i) result += mMeshPtr->mFaces[mFaceIDs[i]].position();
  result /= mFaceIDs.size();
  return result;
}

//------------------------------------------------------------------------------
// PolyhedralZone::volume
// This calculation relies on the faces being planar!
//------------------------------------------------------------------------------
template<>
double
Mesh<Dim<3> >::Zone::
volume() const {
  double result = 0.0;
  const Vector xzone = this->position();
  for (vector<unsigned>::const_iterator faceItr = mFaceIDs.begin();
       faceItr != mFaceIDs.end();
       ++faceItr) {
    const Face& face = mMeshPtr->mFaces[*faceItr];
    const Vector xface = face.position();
    const Vector faceHat = face.unitNormal();
    const double faceArea = face.area();
    result += abs((xface - xzone).dot(faceHat)) * faceArea;
  }
  result /= 3.0;
  return result;
}

}
}
