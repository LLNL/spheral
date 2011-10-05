//---------------------------------Spheral++----------------------------------//
// PolygonalMesh -- 2-D mesh class.
//
// Created by JMO, Tue Nov 16 14:18:20 PST 2010
//----------------------------------------------------------------------------//
#include <limits>
#include <algorithm>
#include <set>
#include <sstream>

#include "boost/foreach.hpp"
#include "boost/random.hpp"
#include "boost/random/uniform_01.hpp"

#include "Mesh.hh"
#include "MeshConstructionUtilities.hh"
#include "findMatchingVertex.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/boundPointWithinBox.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/hashes.hh"
#include "Utilities/DBC.hh"

#include "VoroPP.hh"

#include "Utilities/timingUtilities.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------
// Mesh::reconstructInternal
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<2> >::
reconstructInternal(const vector<Dim<2>::Vector>& generators,
                    const Dim<2>::Vector& xmin0,
                    const Dim<2>::Vector& xmax0) {

  // Some useful typedefs.
  typedef Dim<2> Dimension;

  // Is there anything to do?
  if (generators.size() == 0) return;

  // The tolerance on which we consider positions to be degenerate.
  const double xtol = 1.0e-12*(xmax0 - xmin0).maxElement();
  const double xtol2 = xtol*xtol;

  // Pre-conditions.
  unsigned i, j, k, igen, numGens = generators.size();
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(xmin0.x() < xmax0.x() and
            xmin0.y() < xmax0.y());
    for (igen = 0; igen != numGens; ++igen)
      REQUIRE2(testPointInBox(generators[igen], xmin0, xmax0), 
               "Generator out of bounds:  " << generators[igen] << " not in [" << xmin0 << " " << xmax0 << "]");
    for (i = 0; i < numGens - 1; ++i) {
      for (j = i + 1; j < numGens; ++j) {
        REQUIRE2((generators[i] - generators[j]).magnitude2() > xtol2, 
                 "Degenerate generator positions:  " << i << " " << j << " " << generators[i] << " " << generators[j]);
      }
    }
  }
  END_CONTRACT_SCOPE;

  // This version is currently a hack 'cause the 2D implementation of Voro++ 
  // doesn't really support walls.  We're just cutting down the bounding box.
  Vector xmin = xmin0, xmax = xmax0;
//   for (vector<MeshWallPtr>::iterator wallItr = mWallPtrs.begin();
//        wallItr != mWallPtrs.end();
//        ++wallItr) {
//     xmin = elementWiseMax(xmin, (**wallItr).xmin());
//     xmax = elementWiseMin(xmax, (**wallItr).xmax());
//   }
//   cerr << "PolygonalMesh choosing (xmin, xmax) range to be : " << xmin << " " << xmax << endl;

  // The inverse box scale.
  const Vector box = xmax - xmin;
  const Vector boxInv(safeInv(box.x()),
                      safeInv(box.y()));

  // Construct a random number generator.
  typedef boost::mt19937 base_generator_type;
  base_generator_type basegen(4181989193);
  boost::uniform_01<base_generator_type> rangen(basegen);

  // Use Voro++ to generate our tessellation information.
//   Timing::Time t0 = Timing::currentTime();
  vector<unsigned> badGenerators(1U, 1U);
  vector<vector<Vector> > cellVertices;
  vector<vector<vector<unsigned> > > cellFaceVertices;
  vector<double> cellVolumes;
  vector<Vector> cellCentroids;
  vector<vector<unsigned> > neighborCells;
  unsigned iter = 0;
  double costheta, sintheta;
  vector<Vector>& modgenerators = const_cast<vector<Vector>&>(generators);  // This is deliberately bad to remind us this whole loop is a hack!
  while (badGenerators.size() > 0 and iter < 10) {
    ++iter;

    // Build the Voro++ object.
    VoroPP<Dimension> voro(generators, xmin, xmax,
                           20, 20, 1, 1.0e-8);

//     // Add walls to the container.
//     for (vector<MeshWallPtr>::iterator wallItr = mWallPtrs.begin();
//          wallItr != mWallPtrs.end();
//          ++wallItr) voro.addBoundary(**wallItr);

    // Read out the Voro++ data.
    badGenerators = voro.allCells(cellVertices, cellFaceVertices, cellVolumes, cellCentroids, neighborCells);

    // Did Voro++ miss any generators?
    if (badGenerators.size() > 0) {
      cerr << "  --> Tweaking " << badGenerators.size() << " generators..." << endl;
      for (k = 0; k != badGenerators.size(); ++k) {
        igen = badGenerators[k];
//         cerr << "PolygonalMesh: Tweaking bad generator " << igen << " " << modgenerators[igen];
        costheta = 2.0*(rangen() - 0.5);
        CHECK(abs(costheta) <= 1.0);
        sintheta = sqrt(1.0 - costheta*costheta) * sgn(rangen());
        modgenerators[igen] += Vector(1.0e-5*box.x() * costheta,
                                      1.0e-5*box.y() * sintheta);
//         modgenerators[igen] = boundPointWithinBox(modgenerators[igen] + Vector(1.0e-5*box.x() * costheta,
//                                                                                1.0e-5*box.y() * sintheta),
//                                                   xmin, xmax);
//         cerr << " --> " << modgenerators[igen] << endl;
      }
    }
  }
  VERIFY2(badGenerators.size() == 0, "Voro++ unable to handle all generators: " << badGenerators.size() << " bad of " << numGens);
  CHECK(cellVertices.size() == numGens);
  CHECK(cellFaceVertices.size() == numGens);
  CHECK(cellVolumes.size() == numGens);
  CHECK(cellCentroids.size() == numGens);
  CHECK(neighborCells.size() == numGens);
//   if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to construct Voro++ data." << endl;

  // Prepare storage for the mesh face IDs for each cell, and the set of shared nodes
  // between cells.
//   t0 = Timing::currentTime();
  const int UNSETFACEID = numeric_limits<int>::max();
  vector<vector<int> > cellFaceIDs;
  vector<vector<set<pair<unsigned, unsigned> > > > sharedLocalNodes;
  cellFaceIDs.reserve(numGens);
  sharedLocalNodes.reserve(numGens);
  unsigned nfacesi, nfacesj;
  for (igen = 0; igen != numGens; ++igen) {
    const vector<Vector>& verticesi = cellVertices[igen];
    nfacesi = verticesi.size();
    cellFaceIDs.push_back(vector<int>(nfacesi, UNSETFACEID));
    sharedLocalNodes.push_back(vector<set<pair<unsigned, unsigned> > >(nfacesi));
    for (i = 0; i != nfacesi; ++i) {
      sharedLocalNodes[igen][i].insert(make_pair(igen, i));
    }
  }
  CHECK(cellFaceIDs.size() == numGens);
  CHECK(sharedLocalNodes.size() == numGens);

  // Find the common faces and nodes between zones.  We can assign
  // mesh face IDs at this time as well.
  unsigned numFaces = 0, jgen, iface, jface;
  for (igen = 0; igen != numGens; ++igen) {
    const vector<unsigned>& neighborsi = neighborCells[igen];
    const vector<Vector>& verticesi = cellVertices[igen];
    nfacesi = neighborsi.size();
    CHECK(verticesi.size() == nfacesi);
    CHECK(cellFaceIDs[igen].size() == nfacesi);
    CHECK(sharedLocalNodes[igen].size() == nfacesi);
    for (iface = 0; iface != nfacesi; ++iface) {
      jgen = neighborsi[iface];
//       cerr << "Generator neighbors:  " << igen << " " << jgen << endl;
      if (jgen == UNSETID) {
        // This face is on a boundary, so just assign it a new ID and we're done.
        cellFaceIDs[igen][iface] = numFaces++;
      } else if (jgen > igen) {
        // Find which local face this corresponds to in the neighbor cell.
        const vector<Vector>& verticesj = cellVertices[jgen];
        nfacesj = verticesj.size();
        k = findMatchingVertex(iface, verticesi, verticesj);
        jface = (k == 0 ? nfacesj - 1 : k - 1);
        // if (!(((verticesi[iface] - verticesj[k]).magnitude2() < 1.0e-5) and
        //       ((verticesi[(iface + 1) % nfacesi] - verticesj[jface]).magnitude2() < 1.0e-5))) {
        //   cerr << "Generators:  " << igen << " " << jgen << endl;
        //   cerr << "neighborsi:  ";
        //   for (unsigned ii = 0; ii != neighborsi.size(); ++ii) cerr << neighborsi[ii] << " ";
        //   cerr << endl;
        //   cerr << "verticesi : ";
        //   for (unsigned t = 0; t != verticesi.size(); ++t) cerr << verticesi[t] << " ";
        //   cerr << endl
        //        << "verticesj : ";
        //   for (unsigned t = 0; t != verticesj.size(); ++t) cerr << verticesj[t] << " ";
        //   cerr << endl
        //        << iface << " <-> " << k << endl
        //        << (iface + 1) % nfacesi << " <-> " << jface << endl;
        // }
        CHECK2((verticesi[iface] - verticesj[k]).magnitude2() < 1.0e-5,
               "Possible bad node match:  " << verticesi[iface] << " " << verticesj[k] << " " << (verticesi[iface] - verticesj[k]).magnitude2());
        CHECK2((verticesi[(iface + 1) % nfacesi] - verticesj[jface]).magnitude2() < 1.0e-5,
               "Possible bad node match:  " << verticesi[(iface + 1) % nfacesi] << " " << verticesj[jface] << " " << (verticesi[(iface + 1) % nfacesi] - verticesj[jface]).magnitude2());

        // Assign the new face to each cell.
        cellFaceIDs[igen][iface] = numFaces++;
        cellFaceIDs[jgen][jface] = ~(cellFaceIDs[igen][iface]);

        // Store the matching node info.
        sharedLocalNodes[igen][iface].insert(make_pair(jgen, k));
        sharedLocalNodes[jgen][k].insert(make_pair(igen, iface));

        sharedLocalNodes[igen][(iface + 1) % nfacesi].insert(make_pair(jgen, jface));
        sharedLocalNodes[jgen][jface].insert(make_pair(igen, (iface + 1) % nfacesi));
      }
    }
  }
//   if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to assign face IDs." << endl;

//   // Blago!
//   for (igen = 0; igen != numGens; ++igen) {
//     cerr << igen << " sharedLocalNodes: " << endl;
//     for (i = 0; i != cellVertices[igen].size(); ++i) {
//       cerr << "  --> " << i << " " << cellVertices[igen][i] << " : ";
//       for (set<pair<unsigned, unsigned> >::const_iterator itr = sharedLocalNodes[igen][i].begin();
//            itr != sharedLocalNodes[igen][i].end();
//            ++itr) cerr << "(" << itr->first << " " << itr->second << ") ";
//       cerr << endl;
//     }
//   }
//   // Blago!

  // Prepare storage for the mesh node IDs in each cell.
//   t0 = Timing::currentTime();
  vector<vector<unsigned> > cellMeshNodeIDs;
  vector<vector<unsigned> > nodeZoneIDs;
  cellMeshNodeIDs.reserve(numGens);
  for (igen = 0; igen != numGens; ++igen) cellMeshNodeIDs.push_back(vector<unsigned>(cellVertices[igen].size(), UNSETID));
  CHECK(cellMeshNodeIDs.size() == numGens);

  // Now we know enough to build the mesh node information.
  unsigned numNodes = 0, n;
  for (igen = 0; igen != numGens; ++igen) {
    n = sharedLocalNodes[igen].size();
    for (i = 0; i != n; ++i) {
      if (cellMeshNodeIDs[igen][i] == UNSETID) {
        set<pair<unsigned, unsigned> > fullNodeSet, currentNodeSet, newNodeSet = sharedLocalNodes[igen][i];
        while (not newNodeSet.empty()) {
          currentNodeSet = newNodeSet;
          newNodeSet = set<pair<unsigned, unsigned> >();

          // Create the union of the full set with the current.
          set<pair<unsigned, unsigned> > tmp;
          set_union(fullNodeSet.begin(), fullNodeSet.end(),
                    currentNodeSet.begin(), currentNodeSet.end(),
                    inserter(tmp, tmp.begin()));
          fullNodeSet.swap(tmp);

          // Look for any new shared node info we don't have yet.
          for (set<pair< unsigned, unsigned> >::const_iterator citr = currentNodeSet.begin();
               citr != currentNodeSet.end();
               ++citr) {
            const set<pair<unsigned, unsigned> >& otherNodeSet = sharedLocalNodes[citr->first][citr->second];
            for (set<pair<unsigned, unsigned> >::const_iterator oitr = otherNodeSet.begin();
                 oitr != otherNodeSet.end();
                 ++oitr) {
              if (fullNodeSet.find(*oitr) == fullNodeSet.end()) newNodeSet.insert(*oitr);
            }
          }
        }
        CHECK(fullNodeSet.size() > 0);

        // Now we should know all cells that share this node in common.  Assign them all this as
        // the next mesh node ID.
        mNodePositions.push_back(Vector());
        nodeZoneIDs.push_back(vector<unsigned>());
        nodeZoneIDs.back().reserve(fullNodeSet.size());
        for (set<pair<unsigned, unsigned> >::const_iterator itr = fullNodeSet.begin();
             itr != fullNodeSet.end();
             ++itr) {
          CHECK(itr->first < cellMeshNodeIDs.size() and itr->second < cellMeshNodeIDs[itr->first].size());
          CHECK(cellMeshNodeIDs[itr->first][itr->second] == UNSETID);
          cellMeshNodeIDs[itr->first][itr->second] = numNodes;
          mNodePositions.back() += cellVertices[itr->first][itr->second];
          nodeZoneIDs.back().push_back(itr->first);
        }
        ++numNodes;
        mNodePositions.back() /= fullNodeSet.size();
        boundPointWithinBox(mNodePositions.back(), xmin, xmax);
        CHECK(nodeZoneIDs.back().size() == fullNodeSet.size());
      }
    }
  }
  CHECK(mNodePositions.size() == numNodes);
  CHECK(nodeZoneIDs.size() == numNodes);
//   if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to construct unique Node IDs." << endl;

  BEGIN_CONTRACT_SCOPE;
  {
    // Make sure the unique node positions are consistent with the cell by cell 
    // vertices.
    const double tol = 1.0e-8*(xmax - xmin).maxElement();
    for (igen = 0; igen != numGens; ++igen) {
      const vector<Vector>& vertices = cellVertices[igen];
      const vector<unsigned>& nodeIDs = cellMeshNodeIDs[igen];
      n = vertices.size();
      CHECK2(nodeIDs.size() == n, nodeIDs.size() << " " << n);
      for (i = 0; i != n; ++i) {
        CHECK(nodeIDs[i] < mNodePositions.size());
        CHECK2((vertices[i] - mNodePositions[nodeIDs[i]]).magnitude2() < tol,
               "Blago!  " << igen << " " << i << " " << nodeIDs[i] << " "
               << vertices[i] << " " << mNodePositions[i]);
      }
    }
  }
  END_CONTRACT_SCOPE;

  // Construct the nodes.
//   t0 = Timing::currentTime();
  mNodes.reserve(numNodes);
  unsigned inode;
  for (inode = 0; inode != numNodes; ++inode) {
    CHECK(inode < nodeZoneIDs.size());
    mNodes.push_back(Node(*this, inode, nodeZoneIDs[inode]));
    BEGIN_CONTRACT_SCOPE;
    {
      BOOST_FOREACH(igen, mNodes.back().mZoneIDs) { CHECK(igen < numGens); }
    }
    END_CONTRACT_SCOPE;
  }
  CHECK(mNodes.size() == mNodePositions.size());
//   if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to construct nodes." << endl;

  // Build the edges and faces simultaneously.
//   t0 = Timing::currentTime();
  mEdges.reserve(numFaces);
  mFaces.reserve(numFaces);
  unsigned iedge = 0;
  for (igen = 0; igen != numGens; ++igen) {
    const vector<int>& faceIDs = cellFaceIDs[igen];
    const vector<unsigned>& nodeIDs = cellMeshNodeIDs[igen];
    const vector<unsigned>& neighborsi = neighborCells[igen];
    n = faceIDs.size();
    CHECK(nodeIDs.size() == n);
    CHECK(neighborsi.size() == n);
    for (i = 0; i != n; ++i) {
      if (faceIDs[i] >= 0) {
        CHECK(faceIDs[i] == iedge);
        j = (i + 1) % n;
        CHECK2(nodeIDs[i] != nodeIDs[j], nodeIDs[i] << " " << nodeIDs[j] << " : " << i << " " << j << " " << n);
        mEdges.push_back(Edge(*this, iedge, nodeIDs[i], nodeIDs[j]));
        mFaces.push_back(Face(*this, iedge, igen, neighborsi[i], vector<unsigned>(1, iedge)));
        ++iedge;
      }
    }
  }
  CHECK(mEdges.size() == numFaces);
  CHECK(mFaces.size() == numFaces);
//   if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to construct edges." << endl;

  // Build the zones.
//   t0 = Timing::currentTime();
  mZones.reserve(numGens);
  for (igen = 0; igen != numGens; ++igen) {
    n = cellFaceIDs[igen].size();
    CHECK2(n >= 3, "Bad number of faces for zone:  " << igen << " " << n);
    vector<unsigned> faceIDs;
    faceIDs.reserve(n);
    for (i = 0; i != n; ++i) faceIDs.push_back(cellFaceIDs[igen][i] >= 0 ?
                                               cellFaceIDs[igen][i] :
                                               ~cellFaceIDs[igen][i]);
    CHECK(faceIDs.size() == n);
    mZones.push_back(Zone(*this, igen, faceIDs));
  }
//   if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to construct zones." << endl;
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<>
const unsigned Mesh<Dim<2> >::UNSETID = numeric_limits<unsigned>::max();

}
}

//------------------------------------------------------------------------------
// Instantiate the generic mesh non-inlined methods.
//------------------------------------------------------------------------------
#include "Mesh.cc"
template class Spheral::MeshSpace::Mesh<Spheral::Dim<2> >;
