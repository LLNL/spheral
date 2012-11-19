//---------------------------------Spheral++----------------------------------//
// LineMesh -- 1-D mesh class.
//
// Created by JMO, Tue Oct 12 23:07:22 PDT 2010
//----------------------------------------------------------------------------//
#include <algorithm>
#include <limits>

#include "Mesh.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/DBC.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;

//------------------------------------------------------------------------------
// Comparator function to sort indices.
//------------------------------------------------------------------------------
struct CompareIndicesByPositions {
  CompareIndicesByPositions(const vector<Dim<1>::Vector>& positions):
    mPositions(positions) {}
  bool operator()(const unsigned lhs, const unsigned rhs) {
    REQUIRE(lhs < mPositions.size());
    REQUIRE(rhs < mPositions.size());
    return mPositions[lhs].x() < mPositions[rhs].x();
  }
  const vector<Dim<1>::Vector>& mPositions;
};

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, xmin, xmax)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<1> >::
reconstructInternal(const vector<Mesh<Dim<1> >::Vector>& generators,
                    const Mesh<Dim<1> >::Vector& xmin,
                    const Mesh<Dim<1> >::Vector& xmax) {

  // Is there anything to do?
  if (generators.size() == 0) return;

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(xmin < xmax);
    for (vector<Vector>::const_iterator itr = generators.begin();
         itr != generators.end();
         ++itr) {
      REQUIRE2(xmin <= *itr and *itr <= xmax, "Node out of bounds:  " << *itr << " not in [" << xmin << " " << xmax << "]");
    }
  }
  END_CONTRACT_SCOPE;

  // Find the sorted order for the zones by index.
  vector<unsigned> zoneOrder;
  for (unsigned i = 0; i != generators.size(); ++i) zoneOrder.push_back(i);
  CompareIndicesByPositions zoneComparator(generators);
  sort(zoneOrder.begin(), zoneOrder.end(), zoneComparator);
  CHECK(zoneOrder.size() == generators.size());
  BEGIN_CONTRACT_SCOPE;
  {
    // This check not only checks that the positions are sorted but
    // also that there are no duplicates.
    for (vector<unsigned>::const_iterator itr = zoneOrder.begin();
         itr < zoneOrder.end() - 1;
         ++itr) {
      CHECK2(generators[*itr].x() < generators[*(itr + 1)].x(), 
             "Bad generator positions:  " << generators[*itr] << " " << generators[*(itr + 1)]);
    }
  }
  END_CONTRACT_SCOPE;

  // It's also useful to have the reverse mapping.
  vector<unsigned> sortedZoneIDs(zoneOrder.size(), UNSETID);
  for (unsigned i = 0; i != zoneOrder.size(); ++i) sortedZoneIDs[zoneOrder[i]] = i;
  CHECK(find(sortedZoneIDs.begin(), sortedZoneIDs.end(), UNSETID) == sortedZoneIDs.end());

  // Compute the node positions from the zones.
  mNodePositions.push_back(xmin);
  for (vector<unsigned>::const_iterator itr = zoneOrder.begin();
       itr < zoneOrder.end() - 1;
       ++itr) mNodePositions.push_back(0.5*(generators[*itr] + generators[*(itr + 1)]));
  mNodePositions.push_back(xmax);
  // cerr << "Sorted generators:  ";
  // for (unsigned i = 0; i != zoneOrder.size(); ++i) cerr << " " << generators[zoneOrder[i]];
  // cerr << endl << "Node positions : ";
  // for (unsigned i = 0; i != mNodePositions.size(); ++i) cerr << " " << mNodePositions[i];
  // cerr << endl;
  CHECK2(mNodePositions.size() == generators.size() + 1, "Blago!  " << mNodePositions.size() << " " << generators.size());

  // Construct the nodes from the positions.
  mNodes.reserve(mNodePositions.size());
  for (unsigned i = 0; i != mNodePositions.size(); ++i) {
    vector<unsigned> zones;
    if (i > 0)                         zones.push_back(zoneOrder[i - 1]);
    if (i < mNodePositions.size() - 1) zones.push_back(zoneOrder[i]);
    mNodes.push_back(Node(*this, i, zones));
    // cerr << "Node -> zones:  " << i << " " << zones[0] << " " << zones[1] << " : " << mNodePositions[i] << " "
    //      << (zones[0] != UNSETID ? generators[zones[0]] : -1.0) << " "
    //      << (zones[1] != UNSETID ? generators[zones[1]] : -1.0) << " "
    //      << endl;
  }
  CHECK(mNodes.size() == mNodePositions.size());

  // Construct the edges.  For a 1-D LineMesh edges are degenerate -- they only
  // have one node.
  mEdges.reserve(mNodePositions.size());
  for (unsigned i = 0; i != mNodePositions.size(); ++i) mEdges.push_back(Edge(*this, i, i, i));
  CHECK(mEdges.size() == mNodePositions.size());

  // Construct the faces.  Same as the edges -- one per node.
  mFaces.reserve(mNodePositions.size());
  mFaces.push_back(Face(*this,
                        0,
                        UNSETID,
                        ~zoneOrder[0],
                        vector<unsigned>(1, 0)));
  for (unsigned i = 1; i < mNodePositions.size() - 1; ++i) mFaces.push_back(Face(*this,
                                                                                 i,
                                                                                 zoneOrder[i - 1],
                                                                                 ~zoneOrder[i],
                                                                                 vector<unsigned>(1, i)));
  mFaces.push_back(Face(*this,
                        mNodePositions.size() - 1,
                        zoneOrder[mNodePositions.size() - 2],
                        ~UNSETID,
                        vector<unsigned>(1, mNodePositions.size() - 1)));
  CHECK(mFaces.size() == mNodePositions.size());
  // for (unsigned i = 0; i != mFaces.size(); ++i) {
  //   cerr << "Face:  " << i << " " << mFaces[i].zone1ID() << " " << mFaces[i].zone2ID() << " : "
  //        << mFaces[i].position() << " " 
  //        << (mFaces[i].zone1ID() < generators.size() ? generators[mFaces[i].zone1ID()] : -1.0) << " "
  //        << (mFaces[i].zone2ID() < generators.size() ? generators[mFaces[i].zone2ID()] : -1.0) << " "
  //        << endl;
  // }

  // Finally construct the zones.
  for (unsigned igen = 0; igen != generators.size(); ++igen) {
    const int node1 = sortedZoneIDs[igen];
    const int node2 = node1 + 1;
    CHECK2(generators[igen].x() > mNodePositions[node1].x() and
           generators[igen].x() < mNodePositions[node2].x(),
           "Generator outside node boundaries!  "
           << igen << " " << node1 << " " << node2 << " : "
           << generators[igen] << " "
           << mNodePositions[node1] << " "
           << mNodePositions[node2]);

    // Now we can build the zone.
    // We use the fact that there is a one to one mapping of nodes->faces.
    vector<int> faceIDs;
    faceIDs.push_back(node1);
    faceIDs.push_back(~node2);
    mZones.push_back(Zone(*this, igen, faceIDs));
  }

  // Build the parallel info.
  this->generateDomainInfo();
  // MPI_Barrier(MPI_COMM_WORLD);
  // for (unsigned irank = 0; irank != Process::getTotalNumberOfProcesses(); ++irank) {
  //   if (Process::getRank() == irank) {
  //     cerr << "LineMesh neighborDomains : " << mNeighborDomains.size() << " : ";
  //     copy(mNeighborDomains.begin(), mNeighborDomains.end(), ostream_iterator<unsigned>(cerr, " "));
  //     cerr << endl;
  //   }
  //   MPI_Barrier(MPI_COMM_WORLD);
  // }

  // That's it.
  ENSURE(mNodePositions.size() == generators.size() + 1);
  ENSURE(mNodes.size() == generators.size() + 1);
  ENSURE(mEdges.size() == generators.size() + 1);
  ENSURE(mFaces.size() == generators.size() + 1);
  ENSURE(mZones.size() == generators.size());
  ENSURE(mSharedNodes.size() == mNeighborDomains.size());
  ENSURE(mSharedFaces.size() == mNeighborDomains.size());
  BEGIN_CONTRACT_SCOPE;
  {
    // In 1D we know that every shared node should correspond to a shared face!
    for (unsigned i = 0; i != mSharedNodes.size(); ++i) ENSURE(mSharedFaces[i].size() == mSharedNodes[i].size());
  }
  END_CONTRACT_SCOPE;
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, polygon)
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<1> >::
reconstruct(const vector<Dim<1>::Vector>& generators,
            const Dim<1>::FacetedVolume& boundary) {
  this->reconstruct(generators, boundary.xmin(), boundary.xmax());
}

//------------------------------------------------------------------------------
// Compute the bounding surface of the mesh.
//------------------------------------------------------------------------------
template<>
Dim<1>::FacetedVolume
Mesh<Dim<1> >::
boundingSurface() const {
  // The best we can do is the bounding vertex positions.
  double xmin =  std::numeric_limits<double>::max(), 
         xmax = -std::numeric_limits<double>::max();
  for (unsigned i = 0; i != mNodePositions.size(); ++i) {
    xmin = std::min(xmin, mNodePositions[i].x());
    xmax = std::max(xmax, mNodePositions[i].x());
  }
  xmin = allReduce(xmin, MPI_MIN, MPI_COMM_WORLD);
  xmax = allReduce(xmax, MPI_MAX, MPI_COMM_WORLD);
  return FacetedVolume(Vector(0.5*(xmin + xmax)), 0.5*(xmax - xmin));
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<> const unsigned Mesh<Dim<1> >::minFacesPerZone = 2;
template<> const unsigned Mesh<Dim<1> >::minEdgesPerZone = 2;
template<> const unsigned Mesh<Dim<1> >::minNodesPerZone = 2;
template<> const unsigned Mesh<Dim<1> >::minEdgesPerFace = 1;
template<> const unsigned Mesh<Dim<1> >::minNodesPerFace = 1;

}
}

//------------------------------------------------------------------------------
// Instantiate the generic mesh non-inlined methods.
//------------------------------------------------------------------------------
#include "Mesh.cc"
template class Spheral::MeshSpace::Mesh<Spheral::Dim<1> >;
