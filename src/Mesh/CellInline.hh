#include "Utilities/DBC.hh"

namespace Spheral {
namespace MeshSpace {

//------------------------------------------------------------------------------
// ID, volume, and centroid.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
Cell<Dimension>::
ID() const {
  return mID;
}

template<typename Dimension>
double
Cell<Dimension>::
volume() const {
  return mVolume;
}

template<typename Dimension>
const typename Dimension::Vector&
Cell<Dimension>::
centroid() const {
  return mCentroid;
}

//------------------------------------------------------------------------------
// Access original cell info.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
Cell<Dimension>::
numOldVertices() const {
  return mOldVertices.size();
}

template<typename Dimension>
unsigned
Cell<Dimension>::
numOldFaces() const {
  return mOldFaceVertices.size();
}

template<typename Dimension>
const std::vector<typename Dimension::Vector>& 
Cell<Dimension>::
oldVertices() const {
  return mOldVertices;
}

template<typename Dimension>
const std::vector<std::vector<unsigned> >& 
Cell<Dimension>::
oldFaceVertices() const {
  return mOldFaceVertices;
}

template<typename Dimension>
const std::vector<unsigned>& 
Cell<Dimension>::
oldNeighbors() const {
  return mOldNeighbors;
}

//------------------------------------------------------------------------------
// Access new cell info.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
Cell<Dimension>::
numNewVertices() const {
  return mNewVertices.size();
}

template<typename Dimension>
unsigned
Cell<Dimension>::
numNewFaces() const {
  return mNewFaceVertices.size();
}

template<typename Dimension>
const std::vector<typename Dimension::Vector>& 
Cell<Dimension>::
newVertices() const {
  return mNewVertices;
}

template<typename Dimension>
const std::vector<std::vector<unsigned> >& 
Cell<Dimension>::
newFaceVertices() const {
  return mNewFaceVertices;
}

template<typename Dimension>
const std::vector<unsigned>& 
Cell<Dimension>::
newNeighbors() const {
  return mNewNeighbors;
}

//------------------------------------------------------------------------------
// The map of old -> vertices.
//------------------------------------------------------------------------------
template<typename Dimension>
const std::vector<unsigned>&
Cell<Dimension>::
vertexMap() const {
  return mVertexMap;
}

//------------------------------------------------------------------------------
// The set of cells and vertices that share one of our vertices.
//------------------------------------------------------------------------------
template<typename Dimension>
const std::map<unsigned, unsigned>&
Cell<Dimension>::
cellsForVertex(const unsigned i) const {
  REQUIRE(i < mCellsForVertex.size());
  return mCellsForVertex[i];
}

//------------------------------------------------------------------------------
// The mesh based (or "real") node ID associated with one of our vertices.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned
Cell<Dimension>::
realNodeID(const unsigned ivertex) const {
  REQUIRE(ivertex < mRealNodeIDs.size());
  return mRealNodeIDs[ivertex];
}

}
}
