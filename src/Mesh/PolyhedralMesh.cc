//---------------------------------Spheral++----------------------------------//
// PolyhedralMesh -- 3-D mesh class.
//
// Created by JMO, Wed Jan  5 21:09:11 PST 2011
//----------------------------------------------------------------------------//
#include <limits>
#include <algorithm>
#include <vector>
#include <set>
#include <sstream>
#include <functional>

#include <ctime>

#include "boost/foreach.hpp"
#include "boost/bimap.hpp"
#include "boost/random.hpp"
#include "boost/random/uniform_01.hpp"

#include "Mesh.hh"
#include "VoroPP.hh"
#include "Cell.hh"
#include "findMatchingVertex.hh"
#include "Geometry/Dimension.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/boundPointWithinBox.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/CounterClockwiseComparator.hh"
#include "Utilities/DBC.hh"
#include "MeshConstructionUtilities.hh"

#include "Utilities/timingUtilities.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;
using namespace boost;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Find the active node corresponding to the given node ID.
//------------------------------------------------------------------------------
unsigned
findActiveNode(unsigned i, const vector<unsigned>& nodeMap) {
  const unsigned n = nodeMap.size();
  REQUIRE(i < n);
  while (nodeMap[i] != i) i = nodeMap[i];
  return i;
}

//------------------------------------------------------------------------------
// Find the unique and ordered set of IDs from a list.
//------------------------------------------------------------------------------
vector<unsigned>
uniqueIDset(const vector<unsigned>& cellNodeIDs,
            const vector<unsigned>& indices,
            const map<unsigned, unsigned>& mapping) {
  unsigned i;
  const unsigned nids = cellNodeIDs.size();
  const unsigned n = indices.size();
  BEGIN_CONTRACT_SCOPE;
  {
    BOOST_FOREACH(i, indices) { REQUIRE(i < nids); }
  }
  END_CONTRACT_SCOPE;
  vector<unsigned> result;
  if (n == 0) return result;
  map<unsigned, unsigned>::const_iterator itr = mapping.find(cellNodeIDs[indices[0]]);
  CHECK(itr != mapping.end());
  result.push_back(itr->second);
  i = 1;
  while (i < n) {
    itr = mapping.find(cellNodeIDs[indices[i]]);
    CHECK(itr != mapping.end());
    if (itr->second != result.back()) result.push_back(itr->second);
    ++i;
  }
  if (result.back() == result.front()) result.pop_back();
  BEGIN_CONTRACT_SCOPE;
  {
    BOOST_FOREACH(i, result) { ENSURE(count(result.begin(), result.end(), i) == 1); }
  }
  END_CONTRACT_SCOPE;
  return result;
}

// //------------------------------------------------------------------------------
// // Remove degenerate edges around each face of each zone.
// //------------------------------------------------------------------------------
// void
// removeDegenerateCellVertices(vector<vector<Dim<3>::Vector> >& cellVertices,
//                              vector<vector<vector<unsigned> > >& cellFaceVertices,
//                              vector<vector<unsigned> >& neighborCells,
//                              const double& edgeTolFrac) {

//   typedef Dim<3>::Vector Vector;
//   const unsigned numGens = cellVertices.size();
//   REQUIRE(cellFaceVertices.size() == numGens);
//   REQUIRE(neighborCells.size() == numGens);

//   // Walk the generators.
//   bool change;
//   unsigned igen, i, j, k, ki, kj, nf, nv, iface, numUnique;
//   double edgeTol;
//   for (igen = 0; igen != numGens; ++igen) {
//     nf = cellFaceVertices[igen].size();
//     CHECK(neighborCells[igen].size() == nf);

//     // Scan for the maximum edge and set our removal threshold.
//     for (iface = 0; iface != nf; ++iface) {
//       nv = cellFaceVertices[igen][iface].size();
//       edgeTol = 0.0;
//       for (i = 0; i != nv; ++i) {
//         j = (i + 1) % nv;
//         ki = cellFaceVertices[igen][iface][i];
//         kj = cellFaceVertices[igen][iface][j];
//         edgeTol = max(edgeTol, (cellVertices[igen][ki] - cellVertices[igen][kj]).magnitude2());
//       }
//     }
//     CHECK(edgeTol > 0.0);
//     edgeTol = edgeTolFrac*sqrt(edgeTol);

//     // Prepare the mapping of nodes to their single surviving value.
//     vector<unsigned> nodeMap;
//     for (i = 0; i != cellVertices.size(); ++i) nodeMap.push_back(i);
//     change = false;

//     // Look for any nodes we need to remove.
//     nv = cellVertices[igen].size();
//     for (i = 0; i != nv; ++i) {
//       for (j = i + 1; j < nv; ++j) {
//         if ((cellVertices[igen][i] - cellVertices[igen][j]).magnitude() < edgeTol) {
//           change = true;
//           k = min(nodeMap[i], nodeMap[j]);
//           nodeMap[i] = k;
//           nodeMap[j] = k;
//         }
//       }
//     }

//     // If we need to remove nodes, do it!
//     if (change) {

//       // Force the node map to be consistent, with all redundant nodes pointing
//       // to a single result.
//       nv = cellVertices[igen].size();
//       while (change) {
//         change = false;
//         numUnique = 0;
//         for (i = 0; i != nv; ++i) {
//           if (nodeMap[i] != nodeMap[nodeMap[i]]) {
//             change = true;
//             nodeMap[i] = nodeMap[nodeMap[i]];
//           }
//           if (nodeMap[i] == i) ++numUnique;
//         }
//       }

//       // Remove the postions we are eliminating.
//       vector<unsigned> nodes2kill;
//       numUnique = 0;
//       for (i = 0; i != nv; ++i) {
//         if (nodeMap[i] == i) {
//           nodeMap[i] = numUnique++;
//         } else {
//           nodeMap[i] = nv + 1;
//           nodes2kill.push_back(i);
//         }
//       }
//       CHECK(numUnique + nodes2kill.size() == nv);
//       removeElements(cellVertices[igen], nodes2kill);

//       // Renumber the face vertex indexing.
//       vector<unsigned> faces2kill;
//       for (iface = 0; iface != nf; ++iface) {
//         nv = cellFaceVertices[igen][iface].size();
//         vector<unsigned> newresult;
//         for (i = 0; i != nv; ++i) {
//           j = nodeMap[cellFaceVertices[igen][iface][i]];
//           if (j < nv) newresult.push_back(j);
//         }
//         if (newresult.size() >= 3) {
//           cellFaceVertices[igen][iface] = newresult;
//         } else {
//           faces2kill.push_back(iface);
//         }
//       }

//       // Are we killing any faces?
//       if (faces2kill.size() > 0) {
//         removeElements(cellFaceVertices[igen], faces2kill);
//         removeElements(neighborCells, faces2kill);
//       }
//     }
//   }
// }

//------------------------------------------------------------------------------
// Mesh::reconstructInternal
//------------------------------------------------------------------------------
template<>
void
Mesh<Dim<3> >::
reconstructInternal(const vector<Dim<3>::Vector>& generators,
                    const Dim<3>::Vector& xmin0,
                    const Dim<3>::Vector& xmax0) {

  // Some useful typedefs.
  typedef Dim<3> Dimension;
  typedef ConvexHull::Facet Facet;
  typedef pair<unsigned, unsigned> EdgeHash;
  typedef vector<int> FaceHash;

  // Is there anything to do?
  if (generators.size() == 0) return;

  // The tolerance on which we consider positions to be degenerate.
  const double xtol = 1.0e-12*(xmax0 - xmin0).maxElement();
  const double xtol2 = xtol*xtol;

  // Pre-conditions.
  unsigned i, j, k, igen, jgen, numGens = generators.size();
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(xmin0.x() < xmax0.x() and
            xmin0.y() < xmax0.y() and
            xmin0.z() < xmax0.z());
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

  // Enforce the boundary condition limits on volume.
  Vector xmin = xmin0, xmax = xmax0;
  for (vector<MeshWallPtr>::iterator wallItr = mWallPtrs.begin();
       wallItr != mWallPtrs.end();
       ++wallItr) {
    xmin = elementWiseMax(xmin, (**wallItr).xmin());
    xmax = elementWiseMin(xmax, (**wallItr).xmax());
  }
  // cerr << "PolyhedralMesh choosing (xmin, xmax) range to be : " << xmin << " " << xmax << endl;

  // The inverse box scale.
  const Vector box = xmax - xmin;
  const Vector boxInv(safeInv(box.x()),
                      safeInv(box.y()),
                      safeInv(box.z()));

  // Construct a random number generator.
  typedef boost::mt19937 base_generator_type;
  base_generator_type basegen(4181989193);
  boost::uniform_01<base_generator_type> rangen(basegen);

  // Use Voro++ to generate our tessellation information.
  Timing::Time t0 = Timing::currentTime();
  const double edgeTol = 1.0e-8;
  vector<unsigned> badGenerators(1U, 1U);
  vector<Cell<Dimension> > cells;
  unsigned iter = 0;
  double costheta, sintheta, cosphi, sinphi;
  vector<Vector>& modgenerators = const_cast<vector<Vector>&>(generators);  // This is deliberately bad to remind us this whole loop is a hack!
  while (badGenerators.size() > 0 and iter < 10) {
    ++iter;

    // Build the Voro++ object.
    VoroPP<Dimension> voro(generators, xmin, xmax,
                           20, 20, 20, edgeTol);

//      // Add walls to the container.
//      for (vector<MeshWallPtr>::iterator wallItr = mWallPtrs.begin();
//           wallItr != mWallPtrs.end();
//           ++wallItr) voro.addBoundary(**wallItr);

    // Read out the Voro++ data.
    badGenerators = voro.allCells(cells);

    // Did Voro++ miss any generators?
    if (badGenerators.size() > 0) {
      cerr << "  --> Tweaking " << badGenerators.size() << " generators..." << endl;
      for (k = 0; k != badGenerators.size(); ++k) {
        igen = badGenerators[k];
//         cerr << "PolyhedralMesh: Tweaking bad generator " << igen << " " << modgenerators[igen];
        costheta = 2.0*(rangen() - 0.5);
        cosphi = 2.0*(rangen() - 0.5);
        CHECK(abs(costheta) <= 1.0);
        CHECK(abs(cosphi) <= 1.0);
        sintheta = sqrt(1.0 - costheta*costheta) * sgn(rangen());
        sinphi = sqrt(1.0 - cosphi*cosphi) * sgn(rangen());
        modgenerators[igen] += Vector(1.0e-5*box.x() * costheta*sinphi,
                                      1.0e-5*box.y() * sintheta*sinphi,
                                      1.0e-5*box.z() * cosphi);
//         modgenerators[igen] = boundPointWithinBox(modgenerators[igen] + Vector(1.0e-5*box.x() * costheta,
//                                                                                1.0e-5*box.y() * sintheta),
//                                                   xmin, xmax);
//         cerr << " --> " << modgenerators[igen] << endl;
      }
    }
  }
  VERIFY2(badGenerators.size() == 0, "Voro++ unable to handle all generators: " << badGenerators.size() << " bad of " << numGens);
  CHECK(cells.size() == numGens);
  if (Process::getRank() == 0) cerr << "PolyhedralMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct Voro++ data." << endl;

  // Now we need to take the (probably inconsistent) cell by cell information
  // from Voro++ and stitch it together into a fully consistent topology.
  t0 = Timing::currentTime();

  // Enforce consistency for the cells neighboring across faces.
  for (igen = 0; igen != numGens; ++igen) cells[igen].cullDegenerateNeighbors(cells);

  // Match up face info between cells.
  unsigned iface, nf;
  for (igen = 0; igen != numGens; ++igen) {
    const vector<unsigned>& neighbors = cells[igen].newNeighbors();
    nf = cells[igen].numOldFaces();
    CHECK(neighbors.size() == nf);
    for (iface = 0; iface != nf; ++iface) {
      jgen = neighbors[iface];
      CHECK(jgen == Cell<Dimension>::UNSETID or
            jgen == Cell<Dimension>::DELETED or
            jgen < numGens);
      if (jgen < numGens and
          jgen > igen) cells[igen].matchFace(iface, cells[jgen]);
    }
  }

  // Propagate the info of which cells share vertices until it stops 
  // changing.
  bool done = false;
  while (not done) {
    done = true;
    for (igen = 0; igen != numGens; ++igen) {
      done = done and cells[igen].findMinCellsForVertices(cells);
    }
  }

  // Lock the cells -- they should have consistent info now.
  for (igen = 0; igen != numGens; ++igen) cells[igen].lock(cells);
  Cell<Dimension>::lockMinCellsForVertices(cells);
  if (Process::getRank() == 0) cerr << "PolyhedralMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to stitch together cell topology." << endl;

  // Find the unique node positions, and build up the IDs of the cells that
  // share these nodes.
  t0 = Timing::currentTime();
  unsigned nv;
  vector<vector<unsigned> > nodeZoneIDs;
  for (igen = 0; igen != numGens; ++igen) {
    const vector<Vector>& vertices = cells[igen].newVertices();
    nv = vertices.size();
    for (i = 0; i != nv; ++i) {
      jgen = cells[igen].minCellForVertex(i);
      if (jgen == igen) {
        cells[igen].realNodeID(i, mNodePositions.size());
        // cerr << igen << " ASSIGNING real node ID " << cells[igen].realNodeID(i) << " @ " << vertices[i] << endl;
        mNodePositions.push_back(vertices[i]);
        nodeZoneIDs.push_back(vector<unsigned>(1, igen));
      } else {
        CHECK(jgen < igen);
        j = cells[igen].localVertexForMinCell(i);
        k = cells[jgen].realNodeID(j);
        CHECK2(k < mNodePositions.size(), k << " " << igen << " " << jgen << " " << mNodePositions.size());
        cells[igen].realNodeID(i, k);
        // cerr << igen << " COPYING real node ID " << cells[igen].realNodeID(i) << " @ " << vertices[i] << endl;
        nodeZoneIDs[k].push_back(igen);
      }
    }
  }
  const unsigned numNodes = mNodePositions.size();
  CHECK(nodeZoneIDs.size() == numNodes);
  
  // Create the nodes.
  for (i = 0; i != numNodes; ++i) mNodes.push_back(Node(*this, i, nodeZoneIDs[i]));
  CHECK(mNodes.size() == numNodes);

  // Create the edges.
  unsigned numEdges = 0, inode, jnode;
  EdgeHash ehash;
  map<EdgeHash, unsigned> edgeHash2ID;
  for (igen = 0; igen != numGens; ++igen) {
    const vector<vector<unsigned> >& indices = cells[igen].newFaceVertices();
    nf = indices.size();
    for (iface = 0; iface != nf; ++iface) {
      nv = indices[iface].size();
      for (i = 0; i != nv; ++i) {
        j = (i + 1) % nv;
        inode = cells[igen].realNodeID(indices[iface][i]);
        jnode = cells[igen].realNodeID(indices[iface][j]);
        CHECK(inode < numNodes and jnode < numNodes);
        ehash = hashEdge(inode, jnode);
        if (edgeHash2ID.find(ehash) == edgeHash2ID.end()) {
          mEdges.push_back(Edge(*this, numEdges, inode, jnode));
          edgeHash2ID[ehash] = numEdges++;
        }
      }
    }
  }
  CHECK(mEdges.size() == numEdges);
  CHECK(edgeHash2ID.size() == numEdges);
  
  // Create the faces.
  unsigned numFaces = 0;
  vector<vector<unsigned> > zoneFaceIDs(numGens);
  for (igen = 0; igen != numGens; ++igen) {
    const vector<vector<unsigned> >& indices = cells[igen].newFaceVertices();
    const vector<unsigned>& neighbors = cells[igen].newNeighbors();
    nf = indices.size();
    CHECK(neighbors.size() == nf);
    for (iface = 0; iface != nf; ++iface) {
      jgen = neighbors[iface];
      if (jgen == UNSETID or jgen > igen) {
        vector<unsigned> edgeIDs;
        nv = indices[iface].size();
        for (i = 0; i != nv; ++i) {
          j = (i + 1) % nv;
          inode = cells[igen].realNodeID(indices[iface][i]);
          jnode = cells[igen].realNodeID(indices[iface][j]);
          CHECK(inode < numNodes and jnode < numNodes);
          ehash = hashEdge(inode, jnode);
          CHECK(edgeHash2ID.find(ehash) != edgeHash2ID.end());
          edgeIDs.push_back(edgeHash2ID[ehash]);
        }
        CHECK(edgeIDs.size() >= Cell<Dimension>::minVerticesPerFace);
        mFaces.push_back(Face(*this, numFaces, igen, jgen, edgeIDs));
        zoneFaceIDs[igen].push_back(numFaces);
        if (jgen != UNSETID) zoneFaceIDs[jgen].push_back(numFaces);
        ++numFaces;
      }
    }
  }
  CHECK(mFaces.size() == numFaces);

  // Create the zones.
  for (igen = 0; igen != numGens; ++igen) {
    CHECK(zoneFaceIDs[igen].size() == cells[igen].numNewFaces());
    mZones.push_back(Zone(*this, igen, zoneFaceIDs[igen]));
  }
  if (Process::getRank() == 0) cerr << "PolyhedralMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct mesh elements." << endl;
}

}
}

//------------------------------------------------------------------------------
// Instantiate the generic mesh non-inlined methods.
//------------------------------------------------------------------------------
#include "Mesh.cc"
template class Spheral::MeshSpace::Mesh<Spheral::Dim<3> >;
