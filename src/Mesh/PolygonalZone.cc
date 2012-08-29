//---------------------------------Spheral++----------------------------------//
// PolygonalZone -- 2-D zone class.
//
// Created by JMO, Tue Nov 16 15:30:39 PST 2010
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
//     return isgn0((mElements[id1].position() - mCentroid).cross(mElements[id2].position() - mCentroid).z());
  }
  const vector<Element>& mElements;
  Vector mCentroid;
};

//------------------------------------------------------------------------------
// Mesh::Zone(...)
//------------------------------------------------------------------------------
template<>
Mesh<Dim<2> >::Zone::
Zone(const Mesh<Dim<2> >& mesh,
     const unsigned ID,
     const vector<int>& faceIDs):
  mMeshPtr(&mesh),
  mID(ID),
  mNodeIDs(),
  mEdgeIDs(),
  mFaceIDs(faceIDs) {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(mFaceIDs.size() > 2);
    BOOST_FOREACH(int i, mFaceIDs) {
      int j = (i < 0 ? ~i : i);
      REQUIRE(j < mMeshPtr->mFaces.size());
      REQUIRE(mMeshPtr->mFaces[j].mEdgeIDs.size() == 1);
      REQUIRE(mMeshPtr->mFaces[j].mEdgeIDs[0] == j);
    }
  }
  END_CONTRACT_SCOPE;
  
  // Copy the face IDs as the edge IDs (they are degenerate after all!).
  int faceID;
  BOOST_FOREACH(faceID, mFaceIDs) mEdgeIDs.push_back(faceID < 0 ? ~faceID : faceID);

  // We need the nodes sorted counter-clockwise around the zone.  The faces
  // are already sorted, so we can leverage that info.  Note we are using the
  // fact that faces and edges are degenerate here.
  const unsigned numEdges = mEdgeIDs.size();
  int i;
  unsigned j, n1i, n2i, n1j, n2j;
  BOOST_FOREACH(faceID, mFaceIDs) {
    i = (faceID < 0 ? ~faceID : faceID);
    j = (i + 1) % numEdges;
    n1i = mMeshPtr->mEdges[mEdgeIDs[i]].node1ID();
    n2i = mMeshPtr->mEdges[mEdgeIDs[i]].node2ID();
    n1j = mMeshPtr->mEdges[mEdgeIDs[j]].node1ID();
    n2j = mMeshPtr->mEdges[mEdgeIDs[j]].node2ID();
//     if (!(n1i == n1j or n1i == n2j or
//           n2i == n1j or n2i == n2j)) {
//       cerr << "Bad nodes/edges to zone:  " << endl
//            << "  Nodes for edges:" << endl;
//       for (unsigned ii = 0; ii != numEdges; ++ii) {
//         cerr << "  --> " 
//              << mMeshPtr->mEdges[mEdgeIDs[ii]].node1ID() << " "
//              << mMeshPtr->mEdges[mEdgeIDs[ii]].node2ID() << " : "
//              << mMeshPtr->mNodePositions[mMeshPtr->mEdges[mEdgeIDs[ii]].node1ID()] << " "
//              << mMeshPtr->mNodePositions[mMeshPtr->mEdges[mEdgeIDs[ii]].node2ID()] << endl;
//       }
//     }
    CHECK2(n1i == n1j or n1i == n2j or
           n2i == n1j or n2i == n2j, "Bad node IDs:  " << n1i << " " << n2i << " " << n1j << " " << n2j);
    mNodeIDs.push_back((n1i == n1j) or (n1i == n2j) ? n1i : n2i);
  }
  CHECK(mNodeIDs.size() == mEdgeIDs.size());

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    ENSURE(mNodeIDs.size() == mFaceIDs.size());
    ENSURE(mEdgeIDs.size() == mFaceIDs.size());
    ENSURE(mFaceIDs.size() > 2);
    for (unsigned i = 0; i != mFaceIDs.size(); ++i) {
      ENSURE(mNodeIDs[i] < mMeshPtr->mNodes.size());
      ENSURE(mEdgeIDs[i] < mMeshPtr->mEdges.size());
      ENSURE((mFaceIDs[i] < 0 and ~mFaceIDs[i] < mMeshPtr->mFaces.size()) or
             mFaceIDs[i] < mMeshPtr->mFaces.size());
    }

    // Make sure the elements are unique!
    vector<unsigned> nodeIDs(mNodeIDs);
    sort(nodeIDs.begin(), nodeIDs.end());
    ENSURE(unique(nodeIDs.begin(), nodeIDs.end()) == nodeIDs.end());
    vector<unsigned> edgeIDs(mEdgeIDs);
    sort(edgeIDs.begin(), edgeIDs.end());
    ENSURE(unique(edgeIDs.begin(), edgeIDs.end()) == edgeIDs.end());
    vector<unsigned> faceIDs;
    BOOST_FOREACH(i, mFaceIDs) faceIDs.push_back(i < 0 ? ~i : i);
    sort(faceIDs.begin(), faceIDs.end());
    ENSURE(unique(faceIDs.begin(), faceIDs.end()) == faceIDs.end());

    // Make sure elements are listed counter-clockwise.
    CounterClockwiseCompareElements<Node, Vector> nodeComparator(mMeshPtr->mNodes, mNodeIDs[0]);
    CounterClockwiseCompareElements<Edge, Vector> edgeComparator(mMeshPtr->mEdges, mEdgeIDs[0]);
    CounterClockwiseCompareElements<Face, Vector> faceComparator(mMeshPtr->mFaces, mFaceIDs[0]);
    for (unsigned i = 0; i < mFaceIDs.size() - 1; ++i) {
      ENSURE2(nodeComparator(mNodeIDs[i], mNodeIDs[i + 1]) >= 0,
              nodeComparator(mNodeIDs[i], mNodeIDs[i + 1]) << " "
              << mMeshPtr->mNodes[mNodeIDs[0]].position() << " "
              << mMeshPtr->mNodes[mNodeIDs[i]].position() << " "
              << mMeshPtr->mNodes[mNodeIDs[i + 1]].position());
      ENSURE2(edgeComparator(mEdgeIDs[i], mEdgeIDs[i + 1]) >= 0,
              edgeComparator(mEdgeIDs[i], mEdgeIDs[i + 1]) << " "
              << mMeshPtr->mEdges[mEdgeIDs[0]].position() << " "
              << mMeshPtr->mEdges[mEdgeIDs[i]].position() << " "
              << mMeshPtr->mEdges[mEdgeIDs[i + 1]].position());
      ENSURE2(faceComparator(mFaceIDs[i], mFaceIDs[i + 1]) >= 0,
              faceComparator(mFaceIDs[i], mFaceIDs[i + 1]) << " "
              << mMeshPtr->mFaces[mFaceIDs[0]].position() << " "
              << mMeshPtr->mFaces[mFaceIDs[i]].position() << " "
              << mMeshPtr->mFaces[mFaceIDs[i + 1]].position());
    }
  }
  END_CONTRACT_SCOPE;
}

//------------------------------------------------------------------------------
// PolygonalZone::position
//------------------------------------------------------------------------------
template<>
Dim<2>::Vector
Mesh<Dim<2> >::Zone::
position() const {
  unsigned i, j;
  double weight, weightSum = 0.0;
  Vector result;
  for (i = 0; i != mNodeIDs.size(); ++i) {
    j = (i + 1) % mNodeIDs.size();
    weight = (mMeshPtr->mNodePositions[mNodeIDs[j]] - mMeshPtr->mNodePositions[mNodeIDs[i]]).magnitude();
    weightSum += weight;
    result += weight*(mMeshPtr->mNodePositions[mNodeIDs[j]] + mMeshPtr->mNodePositions[mNodeIDs[i]]);
  }
  CHECK(weightSum > 0.0);
  result /= 2.0*weightSum;
  return result;
}

//------------------------------------------------------------------------------
// PolygonalZone::volume
//------------------------------------------------------------------------------
template<>
double
Mesh<Dim<2> >::Zone::
volume() const {
  double result = 0.0;
  const Vector centroid = this->position();
  unsigned i, j;
  for (i = 0; i != mNodeIDs.size(); ++i) {
    j = (i + 1) % mNodeIDs.size();
    const Vector x1 = mMeshPtr->mNodePositions[mNodeIDs[i]];
    const Vector x2 = mMeshPtr->mNodePositions[mNodeIDs[j]];
    const double dv = (x1 - centroid).cross(x2 - centroid).z();
    CHECK2(fuzzyGreaterThanOrEqual(dv, 0.0, 1.0e-8),
           "Negative triangle!  " << dv << " " << centroid << " " << this->convexHull().contains(centroid) << " " << this->convexHull().distance(centroid));
    result += dv;
  }
  result *= 0.5;
  ENSURE2(result > 0.0, "Nonpositive volume!  " << result);
  return result;
}

}
}
