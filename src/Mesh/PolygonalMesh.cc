//---------------------------------Spheral++----------------------------------//
// PolygonalMesh -- 2-D mesh class.
//
// Created by JMO, Tue Nov 16 14:18:20 PST 2010
//----------------------------------------------------------------------------//
#include <limits>
#include <set>
#include <sstream>
#include "boost/foreach.hpp"

#include "polytope.hh"

#include "Mesh.hh"
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
Mesh<Dim<2> >::
reconstructInternal(const vector<Dim<2>::Vector>& generators,
                    const Dim<2>::Vector& xmin,
                    const Dim<2>::Vector& xmax) {

  // Some useful typedefs.
  typedef Dim<2> Dimension;

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
            xmin.y() < xmax.y());
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
                      safeInv(box.y()));

  // Build the normalized generator positions.
  vector<double> normGens;
  normGens.reserve(2*generators.size());
  for (igen = 0; igen != numGens; ++igen) {
    normGens.push_back((generators[igen].x() - xmin.x())*boxInv.x());
    normGens.push_back((generators[igen].y() - xmin.y())*boxInv.y());
  }
  CHECK(normGens.size() == 2*numGens);

  // Do the polytope tessellation.  We use the Triangle based tessellator for now.
  Timing::Time t0 = Timing::currentTime();
  polytope::Tessellation<2, double> tessellation;
  {
#ifdef USE_MPI
    polytope::DistributedTessellator<2, double> tessellator(new polytope::TriangleTessellator<double>(),
                                                            true,     // Manage memory for serial tessellator
                                                            true);    // Build parallel connectivity
    tessellator.tessellate(normGens, const_cast<double*>(xmin.begin()), const_cast<double*>(xmax.begin()), tessellation);
#else
    polytope::TriangleTessellator<double> tessellator;
    tessellator.tessellate(normGens, const_cast<double*>(xmin.begin()), const_cast<double*>(xmax.begin()), tessellation);
#endif
  }
  CHECK(tessellation.cells.size() == numGens);
  if (Process::getRank() == 0) cerr << "PolygonalMesh:: required " 
                                    << Timing::difference(t0, Timing::currentTime())
                                    << " seconds to construct polytope tessellation." << endl;

  // Read out and construct the Mesh node positions.
  t0 = Timing::currentTime();
  const unsigned numNodes = tessellation.nodes.size()/2;
  mNodePositions.reserve(numNodes);
  for (i = 0; i != numNodes; ++i) {
    mNodePositions.push_back(Vector(max(xmin.x(), min(xmax.x(), xmin.x() + tessellation.nodes[2*i]*box.x())),
                                    max(xmin.y(), min(xmax.y(), xmin.y() + tessellation.nodes[2*i+1]*box.y()))));
  }

  // Build the edges and faces, which are degnerate in 2D.
  // Simultaneously build up the node->zone connectivity.
  const unsigned numEdges = tessellation.faces.size();
  unsigned inode, jnode;
  map<unsigned, set<unsigned> > nodeZones;
  for (i = 0; i != numEdges; ++i) {
    CHECK(tessellation.faces[i].size() == 2);
    CHECK(tessellation.faceCells[i].size() == 1 or
          tessellation.faceCells[i].size() == 2);
    inode = tessellation.faces[i][0];
    jnode = tessellation.faces[i][1];
    igen = tessellation.faceCells[i][0];
    jgen = (tessellation.faceCells[i].size() == 2 ? 
            tessellation.faceCells[i][1] :
            UNSETID);
    mEdges.push_back(Edge(*this, i, inode, jnode));
    mFaces.push_back(Face(*this, i, igen, jgen, vector<unsigned>(1, i)));
    BOOST_FOREACH(j, tessellation.faceCells[i]) {
      nodeZones[inode].insert(j);
      nodeZones[jnode].insert(j);
    }
  }
  CHECK(mEdges.size() == numEdges);
  CHECK(mFaces.size() == numEdges);
  CHECK(nodeZones.size() == numNodes);

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

  // Build the global IDs for the mesh nodes.
  const vector<unsigned> local2globalIDs = this->globalMeshNodeIDs();

  // Look for the faces that bound the mesh.  We build up the global
  // vertex indices, and the associated positions.
  unsigned i, j, iglobal, jglobal;
  map<unsigned, Vector> globalVertexPositions;
  vector<vector<unsigned> > facetIndices;
  BOOST_FOREACH(const Face& face, mFaces) {
    CHECK(face.mNodeIDs.size() == 2);
    i = face.mNodeIDs[0];
    j = face.mNodeIDs[1];
    if ((face.zone1ID() == UNSETID or face.zone2ID() == UNSETID) and
        (sharedNodes.find(i) == sharedNodes.end() or
         sharedNodes.find(j) == sharedNodes.end())) {
      iglobal = local2globalIDs[i];
      jglobal = local2globalIDs[j];
      globalVertexPositions[iglobal] = mNodePositions[i];
      globalVertexPositions[jglobal] = mNodePositions[j];
      vector<unsigned> ids;
      ids.push_back(iglobal);
      ids.push_back(jglobal);
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
          CHECK2(otherIndices.size() == 2, "Bad size: " << otherIndices.size());
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
    CHECK(indices.size() == 2);
    CHECK(global2vertexID.find(indices[0]) != global2vertexID.end());
    CHECK(global2vertexID.find(indices[1]) != global2vertexID.end());
    indices[0] = global2vertexID[indices[0]];
    indices[1] = global2vertexID[indices[1]];
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    BOOST_FOREACH(const vector<unsigned>& indices, facetIndices) {
      ENSURE(indices.size() == 2);
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
