//---------------------------------Spheral++----------------------------------//
// Mesh
//
// Created by JMO, Tue Oct 12 23:07:22 PDT 2010
//----------------------------------------------------------------------------//
#include <limits>
#include <numeric>
#include "boost/unordered_map.hpp"
#include "boost/tuple/tuple_comparison.hpp"

#include "MeshConstructionUtilities.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/DBC.hh"
#include "Utilities/allReduce.hh"
#include "NodeList/NodeList.hh"

#ifdef USE_MPI
extern "C" {
#include "mpi.h"
}
#include "Utilities/packElement.hh"
#endif

namespace Spheral {
namespace MeshSpace {

using namespace std;
using namespace boost;
using ::boost::unordered_map;
using std::min;
using std::max;
using std::abs;

using NodeSpace::NodeList;

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh():
  mNodePositions(),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNeighborDomains(),
  mSharedNodes(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets(),
  mWallPtrs() {
}

//------------------------------------------------------------------------------
// Mesh::Mesh(generators, xmin, xmax)
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh(const vector<typename Dimension::Vector>& generators,
     const typename Dimension::Vector& xmin,
     const typename Dimension::Vector& xmax):
  mNodePositions(),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNeighborDomains(),
  mSharedNodes(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets(),
  mWallPtrs() {
  this->reconstruct(generators, xmin, xmax);
}

//------------------------------------------------------------------------------
// Mesh::Mesh(generators, polygon)
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh(const vector<typename Dimension::Vector>& generators,
     const typename Dimension::FacetedVolume& boundary):
  mNodePositions(),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNeighborDomains(),
  mSharedNodes(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets(),
  mWallPtrs() {
  this->reconstruct(generators, boundary);
}

//------------------------------------------------------------------------------
// The explict constructor where the user directly passes in all the 
// information.
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>::
Mesh(const vector<Vector>& nodePositions,
     const vector<vector<unsigned> >& edgeNodes,
     const vector<vector<unsigned> >& faceEdges,
     const vector<vector<unsigned> >& zoneFaces):
  mNodePositions(nodePositions),
  mNodes(),
  mEdges(),
  mFaces(),
  mZones(),
  mNodeListNameOffsets(),
  mNodeListIndexOffsets(),
  mWallPtrs() {

  // Reverse the zoneFace structure.
  vector<vector<unsigned> > faceZones(faceEdges.size());
  for (unsigned izone = 0; izone != zoneFaces.size(); ++izone) {
    const vector<unsigned>& zf = zoneFaces[izone];
    for (unsigned j = 0; j != zf.size(); ++j) {
      VERIFY(zf[j] < faceZones.size());
      faceZones[zf[j]].push_back(izone);
    }
  }

  // Back out the set of zones for each node.
  vector<set<unsigned> > nodeZones(nodePositions.size());
  for (unsigned izone = 0; izone != zoneFaces.size(); ++izone) {
    const vector<unsigned>& zf = zoneFaces[izone];
    for (unsigned i = 0; i != zf.size(); ++i) {
      const unsigned iface = zf[i];
      CHECK(iface < faceEdges.size());
      const vector<unsigned>& fe = faceEdges[iface];
      for (unsigned j = 0; j != fe.size(); ++j) {
        const unsigned iedge = fe[j];
        CHECK(iedge < edgeNodes.size());
        const vector<unsigned>& en = edgeNodes[iedge];
        CHECK(en.size() == 2);
        CHECK(en[0] < nodeZones.size());
        CHECK(en[1] < nodeZones.size());
        nodeZones[en[0]].insert(izone);
        nodeZones[en[1]].insert(izone);
      }
    }
  }

  // Construct the nodes.
  for (unsigned inode = 0; inode != nodePositions.size(); ++inode) {
    mNodes.push_back(Node(*this, inode, vector<unsigned>(nodeZones[inode].begin(), nodeZones[inode].end())));
  }

  // Construct the edges.
  for (unsigned iedge = 0; iedge != edgeNodes.size(); ++iedge) {
    VERIFY(edgeNodes[iedge].size() == 2);
    mEdges.push_back(Edge(*this, iedge, edgeNodes[iedge][0], edgeNodes[iedge][1]));
  }

  // Construct the faces.
  for (unsigned iface = 0; iface != faceEdges.size(); ++iface) {
    VERIFY(faceZones[iface].size() == 1 or faceZones[iface].size() == 2);
    if (faceZones[iface].size() == 1) {
      mFaces.push_back(Face(*this, iface, faceZones[iface][0], UNSETID, faceEdges[iface]));
    } else {
      mFaces.push_back(Face(*this, iface, faceZones[iface][0], faceZones[iface][1], faceEdges[iface]));
    }
  }

  // Construct the zones.
  for (unsigned izone = 0; izone != zoneFaces.size(); ++izone) {
    mZones.push_back(Zone(*this, izone, zoneFaces[izone]));
  }

  // Post-conditions.
  ENSURE(mNodePositions.size() == nodePositions.size());
  ENSURE(mNodes.size() == nodePositions.size());
  ENSURE(mEdges.size() == edgeNodes.size());
  ENSURE(mFaces.size() == faceEdges.size());
  ENSURE(mZones.size() == zoneFaces.size());
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>&
Mesh<Dimension>::
operator=(const Mesh<Dimension>& rhs) {
  if (this != &rhs) {
    mNodePositions = rhs.mNodePositions;
    mNodes = rhs.mNodes;
    mEdges = rhs.mEdges;
    mFaces = rhs.mFaces;
    mZones = rhs.mZones;
    mNeighborDomains = rhs.mNeighborDomains;
    mSharedNodes = rhs.mSharedNodes;
    mNodeListNameOffsets = rhs.mNodeListNameOffsets;
    mNodeListIndexOffsets = rhs.mNodeListIndexOffsets;
    mWallPtrs = rhs.mWallPtrs;

    // Set the mesh pointers appropriately.
    for (typename vector<Node>::iterator itr = mNodes.begin();
         itr != mNodes.end();
         ++itr) itr->mMeshPtr = this;
    for (typename vector<Edge>::iterator itr = mEdges.begin();
         itr != mEdges.end();
         ++itr) itr->mMeshPtr = this;
    for (typename vector<Face>::iterator itr = mFaces.begin();
         itr != mFaces.end();
         ++itr) itr->mMeshPtr = this;
    for (typename vector<Zone>::iterator itr = mZones.begin();
         itr != mZones.end();
         ++itr) itr->mMeshPtr = this;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Mesh::clear
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
clear() {
  mNodePositions =        vector<Vector>();
  mNodes =                NodeContainer();
  mEdges =                EdgeContainer();
  mFaces =                FaceContainer();
  mZones =                ZoneContainer();
  mNeighborDomains =      vector<unsigned>();
  mSharedNodes =          vector<vector<unsigned> >();
  mNodeListNameOffsets =  map<string, unsigned>();
  mNodeListIndexOffsets = vector<unsigned>();
  mWallPtrs =             vector<MeshWallPtr>();
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, xmin, xmax)
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
reconstruct(const vector<typename Dimension::Vector>& generators,
            const typename Dimension::Vector& xmin,
            const typename Dimension::Vector& xmax) {
  this->clear();
  this->reconstructInternal(generators, xmin, xmax);
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, polygon)
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
reconstruct(const vector<typename Dimension::Vector>& generators,
            const typename Dimension::FacetedVolume& boundary) {
  this->clear();
  this->addWall(MeshWallPtr(new FacetedMeshWall<Dimension>(boundary)));
  this->reconstructInternal(generators, boundary.xmin(), boundary.xmax());
}

//------------------------------------------------------------------------------
// Mesh::removeZonesByMask
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
removeZonesByMask(const vector<unsigned>& zoneMask) {

  // Pre-conditions.
  VERIFY2(zoneMask.size() == mZones.size(), "Mask wrong size:  " << zoneMask.size() << " != " << mZones.size());

  // Construct masks indicating which nodes, edges, and faces we need to keep
  // or remove.
  vector<unsigned> nodeMask(mNodes.size(), 0);
  vector<unsigned> edgeMask(mEdges.size(), 0);
  vector<unsigned> faceMask(mFaces.size(), 0);
  for (size_t izone = 0; izone != zoneMask.size(); ++izone) {
    if (zoneMask[izone] == 1) {
      const Zone& zone = mZones[izone];
      const vector<unsigned>& nodeIDs = zone.nodeIDs();
      const vector<unsigned>& edgeIDs = zone.edgeIDs();
      const vector<unsigned>& faceIDs = zone.faceIDs();
      for (vector<unsigned>::const_iterator itr = nodeIDs.begin();
           itr != nodeIDs.end();
           ++itr) {
        CHECK(*itr < mNodes.size());
        nodeMask[*itr] = 1;
      }
      for (vector<unsigned>::const_iterator itr = edgeIDs.begin();
           itr != edgeIDs.end();
           ++itr) {
        CHECK(*itr < mEdges.size());
        edgeMask[*itr] = 1;
      }
      for (vector<unsigned>::const_iterator itr = faceIDs.begin();
           itr != faceIDs.end();
           ++itr) {
        CHECK(*itr < mFaces.size());
        faceMask[*itr] = 1;
      }
    }
  }
  
  // We need the new IDs for all elements.
  const vector<unsigned> newNodeIDs = this->recomputeIDs(nodeMask);
  const vector<unsigned> newEdgeIDs = this->recomputeIDs(edgeMask);
  const vector<unsigned> newFaceIDs = this->recomputeIDs(faceMask);
  const vector<unsigned> newZoneIDs = this->recomputeIDs(zoneMask);
  CHECK2(newNodeIDs.size() == mNodes.size(), newNodeIDs.size() << " " << mNodes.size());
  CHECK2(newEdgeIDs.size() == mEdges.size(), newEdgeIDs.size() << " " << mEdges.size());
  CHECK2(newFaceIDs.size() == mFaces.size(), newFaceIDs.size() << " " << mFaces.size());
  CHECK2(newZoneIDs.size() == mZones.size(), newZoneIDs.size() << " " << mZones.size());

  // Update the IDs of nodes and their internal data.
  for (unsigned i = 0; i != mNodes.size(); ++i) {
    Node& node = mNodes[i];
    node.mID = newNodeIDs[i];
    reassignIDs(node.mZoneIDs, newZoneIDs);
  }
  
  // Update the IDs of edges and their internal data.
  for (unsigned i = 0; i != mEdges.size(); ++i) {
    Edge& edge = mEdges[i];
    edge.mID = newEdgeIDs[i];
    edge.mNode1ID = newNodeIDs[edge.mNode1ID];
    edge.mNode2ID = newNodeIDs[edge.mNode2ID];
  }

  // Update the IDs of faces and their internal data.
  for (unsigned i = 0; i != mFaces.size(); ++i) {
    Face& face = mFaces[i];
    face.mID = newFaceIDs[i];
    if (face.mZone1ID != UNSETID) face.mZone1ID = newZoneIDs[face.mZone1ID];
    if (face.mZone2ID != UNSETID) face.mZone2ID = newZoneIDs[face.mZone2ID];
    reassignIDs(face.mNodeIDs, newNodeIDs);
    reassignIDs(face.mEdgeIDs, newEdgeIDs);
  }

  // Update the IDs of the zones and their internal data.
  for (unsigned i = 0; i != mZones.size(); ++i) {
    Zone& zone = mZones[i];
    zone.mID = newZoneIDs[i];
    reassignIDs(zone.mNodeIDs, newNodeIDs);
    reassignIDs(zone.mEdgeIDs, newEdgeIDs);
    reassignIDs(zone.mFaceIDs, newFaceIDs);
  }

  // Erase the nodes.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != nodeMask.size(); ++i) {
      if (nodeMask[i] == 0) kill.push_back(i);
    }
    removeElements(mNodes, kill);
    removeElements(mNodePositions, kill);
  }

  // Erase the edges.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != edgeMask.size(); ++i) {
      if (edgeMask[i] == 0) kill.push_back(i);
    }
    removeElements(mEdges, kill);
  }

  // Erase the faces.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != faceMask.size(); ++i) {
      if (faceMask[i] == 0) kill.push_back(i);
    }
    removeElements(mFaces, kill);
  }

  // Erase the zones.
  {
    vector<unsigned> kill;
    for (unsigned i = 0; i != zoneMask.size(); ++i) {
      if (zoneMask[i] == 0) kill.push_back(i);
    }
    removeElements(mZones, kill);
  }

  // Any pre-existing parallel info is now invalid, so just clear out the old
  // data.
  mNeighborDomains = vector<unsigned>();
  mSharedNodes = vector<vector<unsigned> >();
}

//------------------------------------------------------------------------------
// Mesh::generateDomainInfo
// Create the parallel domain info: i.e., which nodes on each domain are linked.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
generateDomainInfo() {
  REQUIRE(mNodePositions.size() == numNodes());

  // This method is empty and a no-op unless we're building a parallel code!
#ifdef USE_MPI

  // Start out by determining the global extent of the mesh.
  Vector xmin(numeric_limits<double>::max(),
              numeric_limits<double>::max(),
              numeric_limits<double>::max());
  Vector xmax(-numeric_limits<double>::max(),
              -numeric_limits<double>::max(),
              -numeric_limits<double>::max());
  for (unsigned i = 0; i != numNodes(); ++i) {
    xmin = elementWiseMin(xmin, mNodePositions[i]);
    xmax = elementWiseMax(xmax, mNodePositions[i]);
  }
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    xmin(i) = allReduce(xmin(i), MPI_MIN, MPI_COMM_WORLD);
    xmax(i) = allReduce(xmax(i), MPI_MAX, MPI_COMM_WORLD);
  }

  // Define the hashing scale.
  const double dxhash = (xmax - xmin).maxElement() / numeric_limits<KeyElement>::max();

  // Puff out the bounds a bit.  We do the all reduce just to ensure
  // bit perfect consistency across processors.
  Vector boxInv;
  for (unsigned i = 0; i != Dimension::nDim; ++i) {
    xmin(i) = allReduce(xmin(i) - dxhash, MPI_MIN, MPI_COMM_WORLD);
    xmax(i) = allReduce(xmax(i) + dxhash, MPI_MAX, MPI_COMM_WORLD);
    boxInv(i) = safeInv(xmax(i) - xmin(i));
  }

  // Hash the node positions.  We want these sorted by key as well
  // to make testing if a key is present fast.
  vector<Key> nodeHashes;
  unordered_map<Key, unsigned> key2nodeID;
  nodeHashes.reserve(numNodes());
  for (unsigned i = 0; i != numNodes(); ++i) {
    nodeHashes.push_back(hashPosition(mNodePositions[i], xmin, xmax, boxInv));
    key2nodeID[nodeHashes.back()] = i;
  }
  sort(nodeHashes.begin(), nodeHashes.end());
  CHECK2(nodeHashes.size() == numNodes(), "Bad sizes:  " << nodeHashes.size() << " " << numNodes());
  CHECK2(key2nodeID.size() == numNodes(), "Bad sizes:  " << key2nodeID.size() << " " << numNodes());

  // Puff out our domain positions a bit from their centroid to try and
  // ensure our intersection tests don't miss something.
  vector<Vector> hullPoints;
  {
    Vector centroid;
    for (unsigned i = 0; i != mNodePositions.size(); ++i) centroid += mNodePositions[i];
    CHECK(mNodePositions.size() > 0);
    centroid /= mNodePositions.size();
    hullPoints.reserve(mNodePositions.size());
    for (unsigned i = 0; i != mNodePositions.size(); ++i) 
      hullPoints.push_back(1.05*(mNodePositions[i] - centroid) + centroid);
  }

  // Generate convex hulls enclosing each domain.
  const unsigned numDomains = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
  CHECK(rank < numDomains);
  vector<ConvexHull> domainHulls(numDomains);
  domainHulls[rank] = ConvexHull(hullPoints);

  // Globally exchange the convex hulls.  This might be a bottle neck!
  {
    vector<char> localBuffer;
    packElement(domainHulls[rank], localBuffer);
    for (int sendProc = 0; sendProc != numDomains; ++sendProc) {
      unsigned bufSize = localBuffer.size();
      MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      vector<char> buffer = localBuffer;
      buffer.resize(bufSize);
      MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, MPI_COMM_WORLD);
      vector<char>::const_iterator itr = buffer.begin();
      unpackElement(domainHulls[sendProc], itr, buffer.end());
      CHECK(itr == buffer.end());
    }
  }

  // Check which bounding volume hulls intersect our own.  This represents the
  // set we'll potentially have to communicate with.
  vector<unsigned> potentialNeighborDomains;
  for (unsigned idomain = 0; idomain != numDomains; ++idomain) {
    if (idomain != rank and domainHulls[idomain].intersect(domainHulls[rank])) potentialNeighborDomains.push_back(idomain);
  }

  // Ensure consistency in the potential exchange pattern!
  BEGIN_CONTRACT_SCOPE;
  {
    for (int sendProc = 0; sendProc != numDomains; ++sendProc) {
      unsigned num = potentialNeighborDomains.size();
      MPI_Bcast(&num, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      vector<unsigned> otherNeighbors = potentialNeighborDomains;
      otherNeighbors.resize(num);
      MPI_Bcast(&otherNeighbors.front(), num, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      CHECK((binary_search(potentialNeighborDomains.begin(), potentialNeighborDomains.end(), sendProc) == true and
             binary_search(otherNeighbors.begin(), otherNeighbors.end(), rank) == true) or
            (binary_search(potentialNeighborDomains.begin(), potentialNeighborDomains.end(), sendProc) == false and
             binary_search(otherNeighbors.begin(), otherNeighbors.end(), rank) == false));
    }
  }
  END_CONTRACT_SCOPE;

  // Exchange the keys.
  vector<vector<Key> > neighborHashes;
  exchangeTuples(nodeHashes, potentialNeighborDomains, neighborHashes);

  // Determine which of our nodes are shared with each neighbor domain.
  const unsigned numNeighborDomains = potentialNeighborDomains.size();
  mNeighborDomains.reserve(numNeighborDomains);
  mSharedNodes.reserve(numNeighborDomains);
  for (unsigned k = 0; k != numNeighborDomains; ++k) {
    const unsigned otherProc = potentialNeighborDomains[k];
    const vector<Key>& otherHashes = neighborHashes[k];
    vector<unsigned> sharedIDs;
    sharedIDs.reserve(otherHashes.size());
    for (vector<Key>::const_iterator keyItr = nodeHashes.begin();
         keyItr != nodeHashes.end();
         ++keyItr) {
      if (binary_search(otherHashes.begin(), otherHashes.end(), *keyItr) == true)
        sharedIDs.push_back(key2nodeID[*keyItr]);
    }
    if (sharedIDs.size() > 0) {
      mNeighborDomains.push_back(otherProc);
      mSharedNodes.push_back(sharedIDs);
    }
  }
  CHECK(mNeighborDomains.size() == mSharedNodes.size());

  // Check that the communicated info is consistent.
  CHECK2(validDomainInfo(xmin, xmax, false) == "", validDomainInfo(xmin, xmax, false));
        
  // Determine which process owns each shared node.  We use the convention that the process
  // with the lowest rank wins.
  unordered_map<unsigned, unsigned> nodeOwner;
  for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
    for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
      CHECK(k < mSharedNodes.size() and j < mSharedNodes[k].size());
      nodeOwner[mSharedNodes[k][j]] = numDomains;
    }
  }
  for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
    const unsigned potentialOwner = min(rank, mNeighborDomains[k]);
    for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
      CHECK(k < mSharedNodes.size() and j < mSharedNodes[k].size());
      const unsigned i = mSharedNodes[k][j];
      CHECK(nodeOwner.find(i) != nodeOwner.end());
      nodeOwner[i] = min(nodeOwner[i], potentialOwner);
    }
  }

  // Turns out we actually want to keep the duplicates nodes!
//   // Cull the lists of shared nodes so that we only talk with the owner about any individual
//   // node.
//   {
//     vector<unsigned> procsToKill;
//     for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
//       const unsigned otherProc = mNeighborDomains[k];
//       vector<unsigned> nodesToKill;
//       CHECK(k < mSharedNodes.size());
//       for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
//         CHECK(k < mSharedNodes.size() and j < mSharedNodes[k].size());
//         const unsigned nodeID = mSharedNodes[k][j];
//         CHECK(nodeOwner.find(nodeID) != nodeOwner.end());
//         const unsigned owner = nodeOwner[nodeID];
//         CHECK(owner < numDomains);
//         if (owner != rank and owner != otherProc) nodesToKill.push_back(j);
//       }
//       if (nodesToKill.size() == mSharedNodes[k].size()) {
//         procsToKill.push_back(k);
//       } else {
//         removeElements(mSharedNodes[k], nodesToKill);
//       }
//     }
//     removeElements(mNeighborDomains, procsToKill);
//     removeElements(mSharedNodes, procsToKill);
//   }
//   CHECK(mNeighborDomains.size() == mSharedNodes.size());

//   // Check that the communicated info is consistent.
//   CHECK2(validDomainInfo(xmin, xmax, true) == "", validDomainInfo(xmin, xmax, true));

  // For the sake of bit-perfectness, exchange node positions such that every shared
  // node has the position it takes on the owner domain.
  {
    // Start by sending our info for all nodes we own.
    list<vector<char> > sendBuffers;
    list<unsigned> sendBufferSizes;
    vector<MPI_Request> sendRequests;
    sendRequests.reserve(2*mNeighborDomains.size());
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      CHECK(k < mSharedNodes.size());
      const unsigned otherProc = mNeighborDomains[k];
      if (otherProc > rank) {
        vector<char> buf;
        for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
          const unsigned i = mSharedNodes[k][j];
          CHECK(i < mNodePositions.size());
          CHECK(nodeOwner.find(i) != nodeOwner.end());
          if (nodeOwner[i] == rank) packElement(mNodePositions[i], buf);
        }
        sendBufferSizes.push_back(buf.size());
        sendBuffers.push_back(buf);
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&sendBufferSizes.back(), 1, MPI_UNSIGNED, otherProc, 10, MPI_COMM_WORLD, &sendRequests.back());
        if (buf.size() > 0) {
          sendRequests.push_back(MPI_Request());
          MPI_Isend(&(sendBuffers.back().front()), buf.size(), MPI_CHAR, otherProc, 11, MPI_COMM_WORLD, &sendRequests.back());
        }
      }
    }
    CHECK(sendRequests.size() <= 2*mNeighborDomains.size());
    CHECK(sendBufferSizes.size() == sendBuffers.size());

    // Now get the postions from all domains sending to us.
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      const unsigned otherProc = mNeighborDomains[k];
      if (otherProc < rank) {
        unsigned bufSize;
        MPI_Status stat1, stat2;
        MPI_Recv(&bufSize, 1, MPI_UNSIGNED, otherProc, 10, MPI_COMM_WORLD, &stat1);
        if (bufSize > 0) {
          vector<char> buffer(bufSize);
          MPI_Recv(&buffer.front(), bufSize, MPI_CHAR, otherProc, 11, MPI_COMM_WORLD, &stat2);
          vector<char>::const_iterator itr = buffer.begin();
          for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
            const unsigned i = mSharedNodes[k][j];
            CHECK(i < mNodePositions.size());
            CHECK(nodeOwner.find(i) != nodeOwner.end());
            if (nodeOwner[i] == otherProc) unpackElement(mNodePositions[i], itr, buffer.end());
          }
          CHECK(itr == buffer.end());
        }
      }
    }

    // Wait for our sends to complete.
    vector<MPI_Status> status(sendRequests.size());
    MPI_Waitall(sendRequests.size(), &sendRequests.front(), &status.front());
  }

  // That's it.
  ENSURE2(validDomainInfo(xmin, xmax, false) == "", validDomainInfo(xmin, xmax, false));
#endif
}

//------------------------------------------------------------------------------
// Check that the communicated information is consistent.
// Optionally we can also check that each send node is owned by only one 
// process.
//------------------------------------------------------------------------------
template<typename Dimension>
string
Mesh<Dimension>::
validDomainInfo(const typename Dimension::Vector& xmin,
                const typename Dimension::Vector& xmax,
                const bool checkUniqueSendProc) const {
  string result = "";

#ifdef USE_MPI
  const unsigned rank = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();

  // The inverse box scale.
  Vector boxInv;
  for (unsigned i = 0; i != Dimension::nDim; ++i) boxInv(i) = safeInv(xmax(i) - xmin(i));

  // First check that all processors agree about who is talking to whom.
  {
    for (unsigned sendProc = 0; sendProc != numDomains; ++sendProc) {
      unsigned numNeighbors = mNeighborDomains.size();
      vector<unsigned> checkNeighbors = mNeighborDomains;
      MPI_Bcast(&numNeighbors, 1, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
      if (numNeighbors > 0) {
        checkNeighbors.resize(numNeighbors);
        MPI_Bcast(&checkNeighbors.front(), numNeighbors, MPI_UNSIGNED, sendProc, MPI_COMM_WORLD);
        if (not (binary_search(checkNeighbors.begin(), checkNeighbors.end(), rank) == 
                 binary_search(mNeighborDomains.begin(), mNeighborDomains.end(), sendProc))) {
          result = "Processors don't agree about who is talking to whom!";
        }
      }
    }
  }

  // Check that the processors that are talking to agree about the number of shared nodes.
  if (result == "") {
    vector<MPI_Request> requests(2*mNeighborDomains.size());
    vector<unsigned> numLocalSharedNodes, numOtherSharedNodes(mNeighborDomains.size());
    for (unsigned k = 0; k != mSharedNodes.size(); ++k) numLocalSharedNodes.push_back(mSharedNodes[k].size());
    CHECK(numLocalSharedNodes.size() == mNeighborDomains.size());
    CHECK(numOtherSharedNodes.size() == mNeighborDomains.size());
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      const unsigned otherProc = mNeighborDomains[k];
      MPI_Isend(&numLocalSharedNodes[k], 1, MPI_UNSIGNED, otherProc,
                (rank + 1)*numDomains + otherProc, MPI_COMM_WORLD, &requests[k]);
      MPI_Irecv(&numOtherSharedNodes[k], 1, MPI_UNSIGNED, otherProc,
                (otherProc + 1)*numDomains + rank, MPI_COMM_WORLD, &requests[mNeighborDomains.size() + k]);
    }
    vector<MPI_Status> status(requests.size());
    MPI_Waitall(requests.size(), &requests.front(), &status.front());
    if (not (numLocalSharedNodes == numOtherSharedNodes)) {
      result = "Processors don't agree about number of shared nodes.";
    }
  }
      
  // Check that we agree with all our neighbors about which nodes we share.
  if (result == "") {

    vector<Key> allLocalHashes;
    vector<vector<Key> > otherHashes;
    for (unsigned i = 0; i != mNodePositions.size(); ++i) 
      allLocalHashes.push_back(hashPosition(mNodePositions[i], xmin, xmax, boxInv));
    CHECK(allLocalHashes.size() == mNodePositions.size());
    exchangeTuples(allLocalHashes, mNeighborDomains, mSharedNodes, otherHashes);

    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      for (unsigned i = 0; i != mSharedNodes[k].size(); ++i) {
        if (not (allLocalHashes[mSharedNodes[k][i]] == otherHashes[k][i]))
          result = "Hash neighbor comparisons fail!";
      }
    }
  }

  // Check the uniqueness of shared nodes (only shared with one other domain unless owned).
  if (result == "" and checkUniqueSendProc) {
    unordered_map<unsigned, unsigned> shareCount;
    unordered_map<unsigned, vector<unsigned> > nodeSharedWithDomains;
    for (unsigned k = 0; k != mNeighborDomains.size(); ++k) {
      const unsigned otherProc = mNeighborDomains[k];
      CHECK(k < mSharedNodes.size());
      for (unsigned j = 0; j != mSharedNodes[k].size(); ++j) {
        CHECK(j < mSharedNodes[k].size());
        const unsigned i = mSharedNodes[k][j];
        nodeSharedWithDomains[i].push_back(otherProc);
      }
    }
    for (unordered_map<unsigned, vector<unsigned> >::const_iterator itr = nodeSharedWithDomains.begin();
         itr != nodeSharedWithDomains.end();
         ++itr) {
      const unsigned i = itr->first;
      const vector<unsigned>& otherDomains = itr->second;
      if (otherDomains.size() > 1 and
          *min_element(otherDomains.begin(), otherDomains.end()) < rank) {
        stringstream thpt;
        thpt << " Node " << i << " on domain " << rank << " shared with domains [";
        for (unsigned k = 0; k != otherDomains.size(); ++k) thpt << " " << otherDomains[k];
        thpt << "]" << endl;
        result += thpt.str();
      }
    }
  }
#endif

  return result;
}

//------------------------------------------------------------------------------
// Fill in the NodeList offsets.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Mesh<Dimension>::
storeNodeListOffsets(const vector<NodeList<Dimension>*>& nodeListPtrs,
                     const vector<unsigned>& offsets) {
  VERIFY(nodeListPtrs.size() == offsets.size());
  this->storeNodeListOffsets(nodeListPtrs.begin(), nodeListPtrs.end(), offsets);
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename Dimension>
const unsigned Mesh<Dimension>::UNSETID = numeric_limits<unsigned>::max();
}
}
