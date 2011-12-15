//---------------------------------Spheral++----------------------------------//
// Mesh
//
// Created by JMO, Tue Oct 12 23:07:22 PDT 2010
//----------------------------------------------------------------------------//
#include <numeric>
#include <limits>
#include <algorithm>

#include "Utilities/removeElements.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/allReduce.hh"
#include "Utilities/DBC.hh"
#include "MeshConstructionUtilities.hh"

namespace Spheral {
namespace MeshSpace {

//------------------------------------------------------------------------------
// Mesh::~Mesh()
//------------------------------------------------------------------------------
template<typename Dimension>
inline
Mesh<Dimension>::
~Mesh() {
}

//------------------------------------------------------------------------------
// Add a wall for mesh construction.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Mesh<Dimension>::
addWall(typename Mesh<Dimension>::MeshWallPtr wallPtr) {
  mWallPtrs.push_back(wallPtr);
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, xmin, xmax, boundaryBegin, boundaryEnd)
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename BoundaryIterator>
inline
void
Mesh<Dimension>::
reconstruct(const std::vector<typename Dimension::Vector>& generators,
            const typename Dimension::Vector& xmin,
            const typename Dimension::Vector& xmax,
            const BoundaryIterator boundaryBegin,
            const BoundaryIterator boundaryEnd) {
  this->clear();

  // Add the boundary conditions.
  for (BoundaryIterator itr = boundaryBegin; itr != boundaryEnd; ++itr) {
    this->addWall((**itr).meshWall());
  }

  // Dispatch the build.
  this->reconstructInternal(generators, xmin, xmax);
}

//------------------------------------------------------------------------------
// Mesh::reconstruct(generators, FacetedVolume, boundaryBegin, boundaryEnd)
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename BoundaryIterator>
inline
void
Mesh<Dimension>::
reconstruct(const std::vector<typename Dimension::Vector>& generators,
            const typename Dimension::FacetedVolume& boundary,
            const BoundaryIterator boundaryBegin,
            const BoundaryIterator boundaryEnd) {
  this->clear();

  // Add the facted volume as a wall.
  this->addWall(MeshWallPtr(new FacetedMeshWall<Dimension>(boundary)));

  // Add the boundary conditions.
  for (BoundaryIterator itr = boundaryBegin; itr != boundaryEnd; ++itr) {
    this->addWall((**itr).meshWall());
  }

  // Dispatch the build.
  this->reconstructInternal(generators, boundary.xmin(), boundary.xmax());
}

//------------------------------------------------------------------------------
// Mesh::numNodes
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::
numNodes() const {
  return mNodes.size();
}

//------------------------------------------------------------------------------
// Mesh::numEdges
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::
numEdges() const {
  return mEdges.size();
}

//------------------------------------------------------------------------------
// Mesh::numFaces
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::
numFaces() const {
  return mFaces.size();
}

//------------------------------------------------------------------------------
// Mesh::numZones
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::
numZones() const {
  return mZones.size();
}

//------------------------------------------------------------------------------
// Mesh::node(i)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Mesh<Dimension>::Node&
Mesh<Dimension>::
node(const unsigned i) const {
  REQUIRE2(i < mNodes.size(), i << " " << mNodes.size());
  return mNodes[i];
}

//------------------------------------------------------------------------------
// Mesh::nodeBegin
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::NodeIterator
Mesh<Dimension>::
nodeBegin() const {
  return mNodes.begin();
}

//------------------------------------------------------------------------------
// Mesh::nodeEnd
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::NodeIterator
Mesh<Dimension>::
nodeEnd() const {
  return mNodes.end();
}

//------------------------------------------------------------------------------
// Mesh::edge(i)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Mesh<Dimension>::Edge&
Mesh<Dimension>::
edge(const unsigned i) const {
  REQUIRE2(i < mEdges.size(), i << " " << mEdges.size());
  return mEdges[i];
}

//------------------------------------------------------------------------------
// Mesh::edgeBegin
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::EdgeIterator
Mesh<Dimension>::
edgeBegin() const {
  return mEdges.begin();
}

//------------------------------------------------------------------------------
// Mesh::edgeEnd
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::EdgeIterator
Mesh<Dimension>::
edgeEnd() const {
  return mEdges.end();
}

//------------------------------------------------------------------------------
// Mesh::face(i)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Mesh<Dimension>::Face&
Mesh<Dimension>::
face(const unsigned i) const {
  REQUIRE2(i < mFaces.size(), i << " " << mFaces.size());
  return mFaces[i];
}

//------------------------------------------------------------------------------
// Mesh::faceBegin
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::FaceIterator
Mesh<Dimension>::
faceBegin() const {
  return mFaces.begin();
}

//------------------------------------------------------------------------------
// Mesh::faceEnd
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::FaceIterator
Mesh<Dimension>::
faceEnd() const {
  return mFaces.end();
}

//------------------------------------------------------------------------------
// Mesh::zone(i)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Mesh<Dimension>::Zone&
Mesh<Dimension>::
zone(const unsigned i) const {
  REQUIRE2(i < mZones.size(), i << " " << mZones.size());
  return mZones[i];
}

//------------------------------------------------------------------------------
// Mesh::zoneBegin
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::ZoneIterator
Mesh<Dimension>::
zoneBegin() const {
  return mZones.begin();
}

//------------------------------------------------------------------------------
// Mesh::zoneEnd
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Mesh<Dimension>::ZoneIterator
Mesh<Dimension>::
zoneEnd() const {
  return mZones.end();
}

//------------------------------------------------------------------------------
// Mesh::zone(nodelist, i)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Mesh<Dimension>::Zone&
Mesh<Dimension>::
zone(const NodeSpace::NodeList<Dimension>& nodeList, const unsigned i) const {
  std::map<std::string, unsigned>::const_iterator itr = mNodeListNameOffsets.find(nodeList.name());
  REQUIRE(itr != mNodeListNameOffsets.end());
  REQUIRE(itr->second     < mZones.size());
  REQUIRE(itr->second + i < mZones.size());
  return mZones[itr->second + i];
}

template<typename Dimension>
inline
const typename Mesh<Dimension>::Zone&
Mesh<Dimension>::
zone(const unsigned nodeListi, const unsigned i) const {
  REQUIRE(nodeListi < mNodeListIndexOffsets.size());
  REQUIRE(mNodeListIndexOffsets[nodeListi]     < mZones.size());
  REQUIRE(mNodeListIndexOffsets[nodeListi] + i < mZones.size());
  return mZones[mNodeListIndexOffsets[nodeListi] + i];
}

//------------------------------------------------------------------------------
// Mesh::offset(nodelist)
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned
Mesh<Dimension>::
offset(const NodeSpace::NodeList<Dimension>& nodeList) const {
  std::map<std::string, unsigned>::const_iterator itr = mNodeListNameOffsets.find(nodeList.name());
  REQUIRE(itr != mNodeListNameOffsets.end());
  REQUIRE(itr->second < mZones.size());
  return itr->second;
}

template<typename Dimension>
inline
unsigned
Mesh<Dimension>::
offset(const unsigned nodeListi) const {
  REQUIRE(nodeListi < mNodeListIndexOffsets.size());
  return mNodeListIndexOffsets[nodeListi];
}

//------------------------------------------------------------------------------
// Mesh::neighborDomains
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<unsigned>&
Mesh<Dimension>::
neighborDomains() const {
  return mNeighborDomains;
}

//------------------------------------------------------------------------------
// Mesh::sharedNodes
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<std::vector<unsigned> >&
Mesh<Dimension>::
sharedNodes() const {
  return mSharedNodes;
}

//------------------------------------------------------------------------------
// Mesh::minimumScale
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
Mesh<Dimension>::
minimumScale() const {
  
  // Minimum edge length.
  double result = std::numeric_limits<double>::max();
  for (unsigned i = 0; i != numEdges(); ++i) {
    result = std::min(result, (mNodePositions[mEdges[i].node1ID()] - 
                               mNodePositions[mEdges[i].node2ID()]).magnitude2());
  }

  // Minimum face->zone distance.
  for (unsigned izone = 0; izone != numZones(); ++izone) {
    const Zone& zone = this->zone(izone);
    const Vector zonePosition = zone.position();
    const std::vector<unsigned>& faces = zone.faceIDs();
    for (std::vector<unsigned>::const_iterator itr = faces.begin();
         itr != faces.end();
         ++itr) {
      const Face& face = this->face(*itr);
      result = std::min(result, (face.position() - zonePosition).magnitude2());
    }
  }
  result = allReduce(0.5*sqrt(result), MPI_MIN, MPI_COMM_WORLD);

  // That's it.
  ENSURE(result > 0.0);
  return result;
}

//------------------------------------------------------------------------------
// 1-D specialization.
//------------------------------------------------------------------------------
template<>
inline
double
Mesh<Dim<1> >::
minimumScale() const {
  double result = std::numeric_limits<double>::max();
  for (unsigned i = 0; i != numZones(); ++i) {
    result = std::min(result, std::abs(mNodePositions[mZones[i].mNodeIDs[0]].x() - 
                                       mNodePositions[mZones[i].mNodeIDs[1]].x()));
  }
  result = allReduce(0.5*result, MPI_MIN, MPI_COMM_WORLD);

  // That's it.
  ENSURE(result > 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Mesh::recomputeIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::vector<unsigned>
Mesh<Dimension>::
recomputeIDs(const std::vector<unsigned>& mask) const {
  std::vector<unsigned> result(mask.size(), UNSETID);
  unsigned newID = 0;
  for (unsigned i = 0; i != mask.size(); ++i) {
    if (mask[i] == 1) {
      result[i] = newID;
      ++newID;
    }
  }
  ENSURE(std::accumulate(mask.begin(), mask.end(), 0) == newID);
  return result;
}

//------------------------------------------------------------------------------
// Mesh::reassignIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Mesh<Dimension>::
reassignIDs(std::vector<unsigned>& ids,
            const std::vector<unsigned>& old2new) const {
  for (size_t k = 0; k != ids.size(); ++k) {
    if (ids[k] != UNSETID) {
      CHECK(ids[k] < old2new.size());
      ids[k] = old2new[ids[k]];
    }
  }
}

//------------------------------------------------------------------------------
// Mesh::removeUNSETIDs
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
Mesh<Dimension>::
removeUNSETIDs(std::vector<unsigned>& ids) const {
  std::vector<unsigned> kill;
  for (size_t k = 0; k != ids.size(); ++k) {
    if (ids[k] == UNSETID) kill.push_back(k);
  }
  removeElements(ids, kill);
}

//------------------------------------------------------------------------------
// Fill in the NodeList offsets.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename NodeListIterator>
inline
void
Mesh<Dimension>::
storeNodeListOffsets(NodeListIterator begin,
                     NodeListIterator end,
                     const std::vector<unsigned>& offsets) {
  VERIFY2(offsets.size() > 1, "offsets size:  " << offsets.size());
  VERIFY2(std::distance(begin, end) + 1 == offsets.size(),
          "Bad sizes:  " << std::distance(begin, end) << " " << offsets.size());
//   for (typename std::vector<unsigned>::const_iterator itr = offsets.begin();
//        itr < offsets.end() - 1;
//        ++itr) VERIFY2(*itr < mZones.size() or (*itr == mZones.size() and *itr == *(itr + 1)), "Bad offset:  " << *itr << " " << mZones.size());
//   VERIFY2(offsets.back() == mZones.size(), "Bad last offset: " << offsets.back() << " " << mZones.size());
  mNodeListNameOffsets = std::map<std::string, unsigned>();
  mNodeListIndexOffsets = std::vector<unsigned>();
  NodeListIterator itr = begin;
  for (unsigned i = 0; i != offsets.size() - 1; ++i, ++itr) {
    mNodeListNameOffsets[(**itr).name()] = offsets[i];
    mNodeListIndexOffsets.push_back(offsets[i]);
  }
}

}
}
