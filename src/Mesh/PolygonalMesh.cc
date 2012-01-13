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
  typedef pair<unsigned, unsigned> EdgeHash;

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

  // Enforce the boundary condition limits on volume.
  Vector xmin = xmin0, xmax = xmax0;
  for (vector<MeshWallPtr>::iterator wallItr = mWallPtrs.begin();
       wallItr != mWallPtrs.end();
       ++wallItr) {
    xmin = elementWiseMax(xmin, (**wallItr).xmin());
    xmax = elementWiseMin(xmax, (**wallItr).xmax());
  }
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

//     if (badGenerators.size() > 0) {
//       stringstream fname;
//       fname << "BadGenerators_" << Process::getRank() << ".txt";
//       ofstream os(fname.str().c_str());
//       for (i = 0; i != numGens; ++i) {
//         os << generators[i].x() << " " << generators[i].y() << " " << generators[i].z() << endl;
//       }
//       os.close();
//       VERIFY(false);
//     }

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
      }
    }
  }
  VERIFY2(badGenerators.size() == 0, "Voro++ unable to handle all generators: " << badGenerators.size() << " bad of " << numGens);
  CHECK(cells.size() == numGens);
  if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
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
  if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
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
  
//   // BLAGO!
//   cerr << "------------------------------ LOCKED CELLS ------------------------------" << endl;
//   for (igen = 0; igen != numGens; ++igen) {
//     cerr << cells[igen].dumpCell() << endl;
//   }
//   // BLAGO!

  // Create the nodes.
  for (i = 0; i != numNodes; ++i) mNodes.push_back(Node(*this, i, nodeZoneIDs[i]));
  CHECK(mNodes.size() == numNodes);

  // Create the edges, faces, and zones.
  // This logic relies on vertices being order counter-clockwise around the cell.
  unsigned numEdges = 0, inode, jnode;
  EdgeHash ehash;
  map<EdgeHash, unsigned> edgeHash2ID;
  for (igen = 0; igen != numGens; ++igen) {
    nv = cells[igen].newVertices().size();
    const vector<unsigned>& neighbors = cells[igen].newNeighbors();
    CHECK(neighbors.size() == nv);
    vector<unsigned> zoneFaceIDs;
    for (i = 0; i != nv; ++i) {
      j = (i + 1) % nv;
      inode = cells[igen].realNodeID(i);
      jnode = cells[igen].realNodeID(j);
      CHECK(inode < numNodes and jnode < numNodes);
      CHECK2(inode != jnode,
             "Can't have an edge with the same node twice: " << inode << " " << jnode << endl
             << igen << " " << iface << " " << i << " " << j << endl
             << mNodePositions[inode] << endl
             << cells[igen].dumpCell());
      ehash = hashEdge(inode, jnode);
      if (edgeHash2ID.find(ehash) == edgeHash2ID.end()) {
        jgen = neighbors[i];
        vector<unsigned> edgeIDs(1, numEdges);
        mEdges.push_back(Edge(*this, numEdges, inode, jnode));
        mFaces.push_back(Face(*this, numEdges, igen, jgen, edgeIDs));
        edgeHash2ID[ehash] = numEdges++;
      }
      zoneFaceIDs.push_back(edgeHash2ID[ehash]);
    }
    CHECK(zoneFaceIDs.size() == nv);
    mZones.push_back(Zone(*this, igen, zoneFaceIDs));
  }
  CHECK(mEdges.size() == numEdges);
  CHECK(mFaces.size() == numEdges);
  CHECK(edgeHash2ID.size() == numEdges);
  CHECK(mZones.size() == numGens);

//     const vector<vector<unsigned> >& indices = cells[igen].newFaceVertices();
//     nf = indices.size();
//     for (iface = 0; iface != nf; ++iface) {
//       nv = indices[iface].size();
//       for (i = 0; i != nv; ++i) {
//         j = (i + 1) % nv;
//         inode = cells[igen].realNodeID(indices[iface][i]);
//         jnode = cells[igen].realNodeID(indices[iface][j]);
//         CHECK(inode < numNodes and jnode < numNodes);
//         CHECK2(inode != jnode,
//                "Can't have an edge with the same node twice: " << inode << " " << jnode << endl
//                << igen << " " << iface << " " << i << " " << j << endl
//                << mNodePositions[inode] << endl
//                << cells[igen].dumpCell());
//         ehash = hashEdge(inode, jnode);
//         if (edgeHash2ID.find(ehash) == edgeHash2ID.end()) {
//           mEdges.push_back(Edge(*this, numEdges, inode, jnode));
//           edgeHash2ID[ehash] = numEdges++;
//         }
//       }
//     }
//   }
  
//   // Create the faces.
//   for (igen = 0; igen != numGens; ++igen) {
//     const vector<vector<unsigned> >& indices = cells[igen].newFaceVertices();
//     const vector<unsigned>& neighbors = cells[igen].newNeighbors();
//     nf = indices.size();
//     CHECK(neighbors.size() == nf);
//     for (iface = 0; iface != nf; ++iface) {
//       jgen = neighbors[iface];
//       if (jgen == UNSETID or jgen > igen) {
//         vector<unsigned> edgeIDs;
//         nv = indices[iface].size();
//         for (i = 0; i != nv - (3 - Dimension::nDim); ++i) {
//           j = (i + 1) % nv;
//           inode = cells[igen].realNodeID(indices[iface][i]);
//           jnode = cells[igen].realNodeID(indices[iface][j]);
//           CHECK(inode < numNodes and jnode < numNodes);
//           ehash = hashEdge(inode, jnode);
//           CHECK(edgeHash2ID.find(ehash) != edgeHash2ID.end());
//           edgeIDs.push_back(edgeHash2ID[ehash]);
//         }
//         CHECK(edgeIDs.size() >= Cell<Dimension>::minVerticesPerFace);
//         mFaces.push_back(Face(*this, numFaces, igen, jgen, edgeIDs));
//         zoneFaceIDs[igen].push_back(numFaces);
//         if (jgen != UNSETID) zoneFaceIDs[jgen].push_back(numFaces);
//         ++numFaces;
//       }
//     }
//   }

//   // Create the zones.
//   for (igen = 0; igen != numGens; ++igen) {
//     CHECK(zoneFaceIDs[igen].size() == cells[igen].numNewFaces());
//     mZones.push_back(Zone(*this, igen, zoneFaceIDs[igen]));
//   }

  if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct mesh elements." << endl;
}

//------------------------------------------------------------------------------
// Compute the bounding surface of the mesh.
//------------------------------------------------------------------------------
template<>
Dim<2>::FacetedVolume
Mesh<Dim<2> >::
boundingSurface() const {
  // Flatten the set of communicated nodes into a set.
  set<unsigned> sharedNodes;
  unsigned domainID;
  for (domainID = 0; domainID != mSharedNodes.size(); ++domainID) {
    std::copy(mSharedNodes[domainID].begin(), mSharedNodes[domainID].end(),
              std::inserter(sharedNodes, sharedNodes.begin()));
  }

  // Look for the faces that bound the mesh.
  set<unsigned> nodeIDs;
  vector<vector<unsigned> > facetIndices;
  BOOST_FOREACH(const Face& face, mFaces) {
    CHECK(face.mNodeIDs.size() == 2);
    if ((face.zone1ID() == UNSETID or face.zone2ID() == UNSETID) and
        (sharedNodes.find(face.mNodeIDs[0]) == sharedNodes.end() and
         sharedNodes.find(face.mNodeIDs[1]) == sharedNodes.end())) {
      nodeIDs.insert(face.mNodeIDs[0]);
      nodeIDs.insert(face.mNodeIDs[1]);
      vector<unsigned> ids;
      ids.push_back(face.mNodeIDs[0]);
      ids.push_back(face.mNodeIDs[1]);
      facetIndices.push_back(ids);
    }
  }

  // Extract the subset of vertices we actually need, and renumber the indices.
  unsigned i;
  vector<Vector> vertices;
  map<unsigned, unsigned> old2new;
  vertices.reserve(nodeIDs.size());
  for (set<unsigned>::const_iterator itr = nodeIDs.begin();
       itr != nodeIDs.end();
       ++itr) {
    old2new[*itr] = vertices.size();
    vertices.push_back(mNodePositions[*itr]);
  }
  for (i = 0; i != facetIndices.size(); ++i) {
    CHECK(facetIndices[i].size() == 2);
    facetIndices[i][0] = old2new[facetIndices[i][0]];
    facetIndices[i][1] = old2new[facetIndices[i][1]];
  }

#ifdef USE_MPI
  // In the parallel case we have to construct the total surface and distribute
  // it to everyone.
  const unsigned rank = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  if (numDomains > 1) {
    unsigned j, bufSize, nvertices, nfacets;
    vector<char> localBuffer, buffer;
    vector<char>::const_iterator bufItr;
    vector<Vector> otherVertices;
    vector<unsigned> otherIndices;
    packElement(vertices, localBuffer);
    packElement(facetIndices.size(), localBuffer);
    for (i = 0; i != facetIndices.size(); ++i) packElement(facetIndices[i], localBuffer);
    for (domainID = 0; domainID != numDomains; ++domainID) {
      buffer = localBuffer;
      bufSize = localBuffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, domainID, MPI_COMM_WORLD);
      if (bufSize > 0) {
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, domainID, MPI_COMM_WORLD);
        if (rank != domainID) {
          nvertices = vertices.size();
          otherVertices = vector<Vector>();
          bufItr = buffer.begin();
          unpackElement(otherVertices, bufItr, buffer.end());
          copy(otherVertices.begin(), otherVertices.end(), back_inserter(vertices));
          unpackElement(nfacets, bufItr, buffer.end());
          for (i = 0; i != nfacets; ++i) {
            otherIndices = vector<unsigned>();
            unpackElement(otherIndices, bufItr, buffer.end());
            CHECK2(otherIndices.size() == 2, "Bad size: " << otherIndices.size());
            facetIndices.push_back(vector<unsigned>());
            facetIndices.back().reserve(otherIndices.size());
            for (j = 0; j != otherIndices.size(); ++j) facetIndices.back().push_back(otherIndices[j] + nvertices);
            CHECK(facetIndices.back().size() == otherIndices.size());
          }
        }
      }
    }
  }
#endif

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    cerr << "Blago!" << endl;
    BOOST_FOREACH(const vector<unsigned>& indices, facetIndices) {
      ENSURE(indices.size() == 2);
      cerr << "  indices.size() = " << indices.size() << endl;
    }
  }
  END_CONTRACT_SCOPE;

  // That's it.
  return FacetedVolume(vertices, facetIndices);
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<> const unsigned Mesh<Dim<2> >::minFacesPerZone = 3;
template<> const unsigned Mesh<Dim<2> >::minEdgesPerZone = 3;
template<> const unsigned Mesh<Dim<2> >::minNodesPerZone = 3;
template<> const unsigned Mesh<Dim<2> >::minEdgesPerFace = 1;
template<> const unsigned Mesh<Dim<2> >::minNodesPerFace = 2;

}
}

//------------------------------------------------------------------------------
// Instantiate the generic mesh non-inlined methods.
//------------------------------------------------------------------------------
#include "Mesh.cc"
template class Spheral::MeshSpace::Mesh<Spheral::Dim<2> >;
