//---------------------------------Spheral++----------------------------------//
// PolyhedralMesh -- 3-D mesh class.
//
// Created by JMO, Wed Jan  5 21:09:11 PST 2011
//----------------------------------------------------------------------------//
#include <limits>
#include <set>
#include <sstream>
#include "boost/foreach.hpp"

#include "polytope.hh"

#include "Mesh.hh"
#include "MeshConstructionUtilities.hh"
#include "Utilities/DBC.hh"
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
Mesh<Dim<3> >::
reconstructInternal(const vector<Dim<3>::Vector>& generators,
                    const Dim<3>::Vector& xmin,
                    const Dim<3>::Vector& xmax) {

  // Some useful typedefs.
  typedef Dim<3> Dimension;
  typedef pair<unsigned, unsigned> EdgeHash;

  // Is there anything to do?
  if (generators.size() == 0) return;

  // The tolerance on which we consider positions to be degenerate.
  const double xtol = 1.0e-12*(xmax - xmin).maxElement();
  const double xtol2 = xtol*xtol;

  // Pre-conditions.
  unsigned i, j, k, igen, jgen, numGens = generators.size();
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(xmin.x() < xmax.x() and
            xmin.y() < xmax.y() and
            xmin.z() < xmax.z());
    for (igen = 0; igen != numGens; ++igen)
      REQUIRE2(testPointInBox(generators[igen], xmin, xmax), 
               "Generator out of bounds:  " << generators[igen] << " not in [" << xmin << " " << xmax << "]");
    for (i = 0; i < numGens - 1; ++i) {
      for (j = i + 1; j < numGens; ++j) {
        REQUIRE2((generators[i] - generators[j]).magnitude2() > xtol2, 
                 "Degenerate generator positions:  " << i << " " << j << " " << generators[i] << " " << generators[j]);
      }
    }
  }
  END_CONTRACT_SCOPE;

  // The inverse box scale.
  const Vector box = xmax - xmin;
  const Vector boxInv(safeInv(box.x()),
                      safeInv(box.y()),
                      safeInv(box.z()));

  // Build the normalized generator positions.
  vector<double> normGens;
  normGens.reserve(3*generators.size());
  for (igen = 0; igen != numGens; ++igen) {
    normGens.push_back((generators[igen].x() - xmin.x())*boxInv.x());
    normGens.push_back((generators[igen].y() - xmin.y())*boxInv.y());
    normGens.push_back((generators[igen].z() - xmin.z())*boxInv.z());
  }
  CHECK(normGens.size() == 3*numGens);

  // Do the polytope tessellation.  We use the TetGen based tessellator for now.
  Timing::Time t0 = Timing::currentTime();
  polytope::Tessellation<3, double> tessellation;
  {
#ifdef USE_MPI
    polytope::DistributedTessellator<3, double> tessellator(new polytope::TetgenTessellator<double>(),
                                                            true,     // Manage memory for serial tessellator
                                                            true);    // Build parallel connectivity
    tessellator.tessellate(normGens, const_cast<double*>(xmin.begin()), const_cast<double*>(xmax.begin()), tessellation);
#else
    polytope::TetGenTessellator<double> tessellator;
    tessellator.tessellate(normGens, const_cast<double*>(xmin.begin()), const_cast<double*>(xmax.begin()), tessellation);
#endif
  }
  CHECK(tessellation.cells.size() == numGens);
  if (Process::getRank() == 0) cerr << "PolyhedralMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct polytope tessellation." << endl;

  // Read out and construct the Mesh node positions.
  t0 = Timing::currentTime();
  const unsigned numNodes = tessellation.nodes.size()/3;
  mNodePositions.reserve(numNodes);
  for (i = 0; i != numNodes; ++i) {
    mNodePositions.push_back(Vector(max(xmin.x(), min(xmax.x(), xmin.x() + tessellation.nodes[3*i]*box.x())),
                                    max(xmin.y(), min(xmax.y(), xmin.y() + tessellation.nodes[3*i+1]*box.y())),
                                    max(xmin.z(), min(xmax.z(), xmin.z() + tessellation.nodes[3*i+2]*box.z()))));
  }

  // Build the edges and faces.
  const unsigned numFaces = tessellation.faces.size();
  unsigned inode, jnode, iedge, n, numEdges = 0;
  EdgeHash ehash;
  map<EdgeHash, unsigned> edgeHash2ID;
  map<EdgeHash, unsigned>::iterator emapItr;
  map<unsigned, set<unsigned> > nodeZones;
  for (i = 0; i != numFaces; ++i) {
    n = tessellation.faces[i].size();
    CHECK(n >= 3);
    CHECK(tessellation.faceCells[i].size() == 1 or
          tessellation.faceCells[i].size() == 2);
    igen = tessellation.faceCells[i][0];
    jgen = (tessellation.faceCells[i].size() == 2 ? 
            tessellation.faceCells[i][1] :
            UNSETID);
    vector<unsigned> faceEdges;
    for (j = 0; j != n; ++j) {
      inode = tessellation.faces[i][j];
      jnode = tessellation.faces[i][(j + 1) % n];
      ehash = hashEdge(inode, jnode);
      emapItr = edgeHash2ID.find(ehash);
      if (emapItr == edgeHash2ID.end()) {
        iedge = numEdges;
        edgeHash2ID[ehash] = numEdges++;
        mEdges.push_back(Edge(*this, iedge, inode, jnode));
      } else {
        iedge = emapItr->second;
      }
      faceEdges.push_back(iedge);
    }
    mFaces.push_back(Face(*this, i, igen, jgen, faceEdges));
    BOOST_FOREACH(j, tessellation.faceCells[i]) {
      nodeZones[inode].insert(j);
      nodeZones[jnode].insert(j);
    }
  }
  CHECK(mEdges.size() == numEdges);
  CHECK(mFaces.size() == numFaces);

  // Construct the nodes.
  mNodes.reserve(numNodes);
  for (i = 0; i != numNodes; ++i) mNodes.push_back(Node(*this, i, vector<unsigned>(nodeZones[i].begin(),
                                                                                   nodeZones[i].end())));
  CHECK(mNodes.size() == numNodes);

  // Construct the zones.
  for (i = 0; i != numGens; ++i) {
    vector<unsigned> faceIDs;
    BOOST_FOREACH(int f, tessellation.cells[i]) {
      faceIDs.push_back(f >= 0 ? f : ~f);
    }
    mZones.push_back(Zone(*this, i, faceIDs));
  }
  CHECK(mZones.size() == numGens);

  // Copy the parallel info.
  mNeighborDomains = tessellation.neighborDomains;
  mSharedNodes = tessellation.sharedNodes;

  // Report our final timing and we're done.
  if (Process::getRank() == 0) cerr << "PolyhedralMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct mesh elements." << endl;
}

//------------------------------------------------------------------------------
// Compute the bounding surface of the mesh.
//------------------------------------------------------------------------------
template<>
Dim<3>::FacetedVolume
Mesh<Dim<3> >::
boundingSurface() const {

  // Flatten the set of communicated nodes into a set.
  set<unsigned> sharedNodes;
  unsigned domainID;
  for (domainID = 0; domainID != mSharedNodes.size(); ++domainID) {
    std::copy(mSharedNodes[domainID].begin(), mSharedNodes[domainID].end(),
              std::inserter(sharedNodes, sharedNodes.begin()));
  }

  // Build the global IDs for the mesh nodes.
  const vector<unsigned> local2globalIDs = this->globalMeshNodeIDs();

  // Look for the faces that bound the mesh.  We build up the global
  // vertex indices, and the associated positions.
  bool useFace;
  unsigned i, j, iglobal, izone;
  map<unsigned, Vector> globalVertexPositions;
  vector<vector<unsigned> > facetIndices;
  BOOST_FOREACH(const Face& face, mFaces) {
    const vector<unsigned>& nodeIDs = face.nodeIDs();
    useFace = (face.zone1ID() == UNSETID or face.zone2ID() == UNSETID);
    i = 0;
    while (useFace and i != nodeIDs.size()) {
      useFace = (sharedNodes.find(nodeIDs[i]) == sharedNodes.end());
      ++i;
    }
    if (useFace) {
      vector<unsigned> ids;
      BOOST_FOREACH(i, nodeIDs) {
        CHECK(i < local2globalIDs.size());
        iglobal = local2globalIDs[i];
        CHECK(globalVertexPositions.find(iglobal) == globalVertexPositions.end() or
              fuzzyEqual((globalVertexPositions[iglobal] - mNodePositions[i]).magnitude2(), 0.0));
        globalVertexPositions[iglobal] = mNodePositions[i];
        ids.push_back(iglobal);
      }

      // Do we need to reverse the face orientation?
      izone = std::min(face.zone1ID(), face.zone2ID());
      CHECK(izone != UNSETID);
      if (face.compare(mZones[izone].position()) == 1) std::reverse(ids.begin(), ids.end());

      facetIndices.push_back(ids);
    }
  }

#ifdef USE_MPI
  // In the parallel case we have to construct the total surface and distribute
  // it to everyone.
  const unsigned rank = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  if (numDomains > 1) {
    unsigned bufSize, nfacets;
    vector<char> localBuffer, buffer;
    vector<char>::const_iterator bufItr;
    vector<Vector> otherVertices;
    vector<unsigned> otherIndices;

    // Pack our local data, and then erase our local copy to be rebuilt 
    // consistently for everyone.
    packElement(globalVertexPositions, localBuffer);
    packElement(unsigned(facetIndices.size()), localBuffer);
    for (i = 0; i != facetIndices.size(); ++i) packElement(facetIndices[i], localBuffer);
    globalVertexPositions = map<unsigned, Vector>();
    facetIndices = vector<vector<unsigned> >();

    // Distribute the complete data to everyone.
    for (domainID = 0; domainID != numDomains; ++domainID) {
      buffer = localBuffer;
      bufSize = localBuffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, domainID, MPI_COMM_WORLD);
      if (bufSize > 0) {
        buffer.resize(bufSize);
        MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, domainID, MPI_COMM_WORLD);
        bufItr = buffer.begin();
        unpackElement(globalVertexPositions, bufItr, buffer.end());
        unpackElement(nfacets, bufItr, buffer.end());
        for (i = 0; i != nfacets; ++i) {
          otherIndices = vector<unsigned>();
          unpackElement(otherIndices, bufItr, buffer.end());
          CHECK2(otherIndices.size() >= 3, "Bad size: " << otherIndices.size());

          facetIndices.push_back(otherIndices);
        }
      }
    }
  }
#endif

  // Extract the vertex positions as an array, and map the global IDs 
  // to index in this array.
  map<unsigned, unsigned> global2vertexID;
  vector<Vector> vertices;
  vertices.reserve(globalVertexPositions.size());
  i = 0;
  for (map<unsigned, Vector>::const_iterator itr = globalVertexPositions.begin();
       itr != globalVertexPositions.end();
       ++itr, ++i) {
    global2vertexID[itr->first] = i;
    vertices.push_back(itr->second);
  }
  CHECK(i == globalVertexPositions.size());
  CHECK(global2vertexID.size() == globalVertexPositions.size());
  CHECK(vertices.size() == globalVertexPositions.size());

  // Transform the facet node indices to the vertex array numbering.
  BOOST_FOREACH(vector<unsigned>& indices, facetIndices) {
    CHECK(indices.size() >= 3);
    for (j = 0; j != indices.size(); ++j) {
      CHECK(global2vertexID.find(indices[j]) != global2vertexID.end());
      indices[j] = global2vertexID[indices[j]];
      CHECK(indices[j] < vertices.size());
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    BOOST_FOREACH(const vector<unsigned>& indices, facetIndices) {
      ENSURE(indices.size() >= 3);
      ENSURE(*max_element(indices.begin(), indices.end()) < vertices.size());
    }
  }
  END_CONTRACT_SCOPE;

  // That's it.
  return FacetedVolume(vertices, facetIndices);
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<> const unsigned Mesh<Dim<3> >::minFacesPerZone = 4;
template<> const unsigned Mesh<Dim<3> >::minEdgesPerZone = 6;
template<> const unsigned Mesh<Dim<3> >::minNodesPerZone = 4;
template<> const unsigned Mesh<Dim<3> >::minEdgesPerFace = 3;
template<> const unsigned Mesh<Dim<3> >::minNodesPerFace = 3;

}
}

//------------------------------------------------------------------------------
// Instantiate the generic mesh non-inlined methods.
//------------------------------------------------------------------------------
#include "Mesh.cc"
template class Spheral::MeshSpace::Mesh<Spheral::Dim<3> >;
