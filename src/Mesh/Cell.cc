//---------------------------------Spheral++----------------------------------//
// Cell
// 
// This class is designed to take (possibly) inconsistent cell by cell 
// information from tools like Voro++ and help tie them together into 
// consistent connectivity between cells.
//
// Created by JMO, Sun Sep 11 18:32:22 PDT 2011
//----------------------------------------------------------------------------//
#include <limits>
#include <string>
#include <sstream>
#include "boost/foreach.hpp"

#include "Cell.hh"
#include "findMatchingVertex.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/removeElements.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Default constructor -- only provided for use with creating std::vector!
//------------------------------------------------------------------------------
template<typename Dimension>
Cell<Dimension>::
Cell():
  mID(UNSETID),
  mVolume(0.0),
  mMaxEdge(0.0),
  mCentroid(),
  mOldVertices(),
  mNewVertices(),
  mOldFaceVertices(),
  mNewFaceVertices(),
  mOldNeighbors(),
  mNewNeighbors(),
  mVertexMap(),
  mMinCellForVertex(),
  mRealNodeIDs() {
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Cell<Dimension>::
Cell(const unsigned ID,
     const double vol,
     const typename Dimension::Vector& centroid,
     const vector<typename Dimension::Vector>& vertices,
     const vector<vector<unsigned> >& faceVertices,
     const vector<unsigned>& neighbors,
     const double edgeTol):
  mID(ID),
  mVolume(vol),
  mMaxEdge(0.0),
  mCentroid(centroid),
  mOldVertices(vertices),
  mNewVertices(),
  mOldFaceVertices(faceVertices),
  mNewFaceVertices(faceVertices.size()),
  mOldNeighbors(neighbors),
  mNewNeighbors(neighbors),
  mVertexMap(vertices.size(), UNSETID),
  mMinCellForVertex(),
  mRealNodeIDs() {

  // Pre-conditions.
  const unsigned nf = faceVertices.size();
  unsigned iface, i;
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(neighbors.size() == nf);
    REQUIRE(count(neighbors.begin(), neighbors.end(), ID) == 0);
    for (iface = 0; iface != nf; ++iface) {
      REQUIRE(faceVertices[iface].size() >= minVerticesPerFace);
      BOOST_FOREACH(i, faceVertices[iface]) { REQUIRE(i < vertices.size()); }
    }
  }
  END_CONTRACT_SCOPE;

  // Identify the maximum edge length.
  unsigned j, nv;
  for (iface = 0; iface != nf; ++iface) {
    nv = mOldFaceVertices[iface].size();
    for (i = 0; i != nv; ++i) {
      j = (i + 1) % nv;
      mMaxEdge = max(mMaxEdge, (mOldVertices[mOldFaceVertices[iface][i]] -
                                mOldVertices[mOldFaceVertices[iface][j]]).magnitude2());
    }
  }
  CHECK(mMaxEdge > 0.0);
  mMaxEdge = sqrt(mMaxEdge);

  // Initialize all vertices to be unique.
  nv = mOldVertices.size();
  for (i = 0; i != nv; ++i) mVertexMap[i] = i;

  // Look for any vertices that we need to remove because they are degenerate.
  const double threshold = FastMath::square(edgeTol*mMaxEdge);
  for (i = 0; i != nv; ++i) {
    for (j = i + 1; j < nv; ++j) {
      if ((mOldVertices[i] - mOldVertices[j]).magnitude2() < threshold) mVertexMap[j] = mVertexMap[i];
    }
  }

  // Reduce the vertex map to the unique set.
  updateVertexMap();

  // Create the new (non-degenerate) vertices for this cell.
  // In the process update the mVertexMap to point from the old degenerate
  // vertex ordering to the new non-degenerate set.
  unsigned nv1 = 0, k;
  for (i = 0; i != nv; ++i) {
    if (mVertexMap[i] == i) {
      mNewVertices.push_back(mOldVertices[i]);
      mVertexMap[i] = nv1++;
      CHECK(mNewVertices.size() == nv1);
    } else {
      k = mVertexMap[mVertexMap[i]];
      CHECK(k < nv1);
      CHECK(mVertexMap[k] <= k);
      mVertexMap[i] = k;
    }
  }
  CHECK(mNewVertices.size() == nv1);

  // Initialize the new face vertices to reflect the new non-degenerate 
  // vertex ordering.  This *may* result in degenerate (line or point)
  // faces, which we keep for now.
  unsigned ki, kj;
  for (iface = 0; iface != nf; ++iface) {
    const unsigned nv0 = mOldFaceVertices[iface].size();
    for (i = 0; i != nv0; ++i) {
      j = (i + 1) % nv0;
      ki = mVertexMap[mOldFaceVertices[iface][i]];
      kj = mVertexMap[mOldFaceVertices[iface][j]];
      if (ki != kj) {
        mNewFaceVertices[iface].push_back(ki);
      }
    }
  }

  // Initialize the min cell for each of our vertices to be ourself.
  for (i = 0; i != nv1; ++i) {
    mMinCellForVertex.push_back(make_pair(mID, i));
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    ENSURE(mVolume > 0.0);
    ENSURE(mMaxEdge > 0.0);
    ENSURE(mOldVertices == vertices);
    ENSURE(mNewVertices.size() > 0 and mNewVertices.size() <= mOldVertices.size());
    ENSURE(mOldFaceVertices == faceVertices);
    ENSURE(mNewFaceVertices.size() == mOldFaceVertices.size());
    ENSURE(mOldNeighbors == neighbors);
    ENSURE(mNewNeighbors == mOldNeighbors);
    ENSURE(count(mOldNeighbors.begin(), mOldNeighbors.end(), mID) == 0);
    ENSURE(count(mNewNeighbors.begin(), mNewNeighbors.end(), mID) == 0);
    for (i = 0; i != mNewFaceVertices.size(); ++i) {
      for (j = 0; j != mNewFaceVertices[i].size(); ++j) {
        ENSURE(mNewFaceVertices[i][j] < mNewVertices.size());
        ENSURE2(count(mNewFaceVertices[i].begin(), mNewFaceVertices[i].end(), mNewFaceVertices[i][j]) == 1, 
                "Bad node count : " << i << " " << j << endl << dumpCell());
      }
    }
    ENSURE(mVertexMap.size() == mOldVertices.size());
    ENSURE(mMinCellForVertex.size() == mNewVertices.size());
    for (i = 0; i != mMinCellForVertex.size(); ++i) {
      ENSURE(mMinCellForVertex[i].first == mID);
      ENSURE(mMinCellForVertex[i].second == i);
    }
    ENSURE(mRealNodeIDs.size() == 0);
  }
  END_CONTRACT_SCOPE;
}

//------------------------------------------------------------------------------
// Cull any inconsistent neighbor connectivity.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Cell<Dimension>::
cullDegenerateNeighbors(vector<Cell<Dimension> >& cells) {

  // Pre-conditions.
  const unsigned ncells = cells.size();
  REQUIRE(mID < ncells);
  BEGIN_CONTRACT_SCOPE;
  {
    unsigned i;
    BOOST_FOREACH(i, mNewNeighbors) { REQUIRE(i == UNSETID or
                                              i == DELETED or
                                              i < ncells); }
  }
  END_CONTRACT_SCOPE;

  const unsigned nfi = mNewNeighbors.size();
  unsigned i, j, jcell, iface, jface, nfj, nvi, nvj;

  // Walk our neighbors.
  CHECK(mOldFaceVertices.size() == nfi);
  for (iface = 0; iface != nfi; ++iface) {
    jcell = mNewNeighbors[iface];
    CHECK(jcell != mID);

    // Check if both cells have valid info for this neighbor relation.
    if (jcell < ncells) {
      CHECK(cells[jcell].mID == jcell);
      nfj = cells[jcell].mNewNeighbors.size();
      jface = distance(cells[jcell].mNewNeighbors.begin(),
                       find(cells[jcell].mNewNeighbors.begin(), cells[jcell].mNewNeighbors.end(), mID));

      // If the other cell doesn't have us listed as a neighbor, we have to remove
      // this face and neighbor association.
      if (jface == nfj) {
        mNewNeighbors[iface] = DELETED;
        mNewFaceVertices[iface] = vector<unsigned>();

      } else {
        // If we got here both cells acknowledge the cross-face relationship,
        // so ensure we have enough vertices in the face to define a 
        // valid area.
        CHECK(jface < nfj);
        nvi =              mNewFaceVertices[iface].size();
        nvj = cells[jcell].mNewFaceVertices[jface].size();
        if (min(nvi, nvj) < minVerticesPerFace) {
          // A least one of these cells doesn't have enough unique vertices on this face
          // to form a valid area, so delete the face from both cells.
          mNewNeighbors[iface] = DELETED;
          cells[jcell].mNewNeighbors[jface] = DELETED;
          mNewFaceVertices[iface] = vector<unsigned>();
          cells[jcell].mNewFaceVertices[jface] = vector<unsigned>();
        }
      }
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    // We should have a valid set of neighbors when we're done.
    for (iface = 0; iface != nfi; ++iface) {
      jcell = mNewNeighbors[iface];
      ENSURE(jcell == DELETED or jcell == UNSETID or jcell < ncells);
      if (jcell < ncells) {
        nfj = cells[jcell].mNewNeighbors.size();
        jface = distance(cells[jcell].mNewNeighbors.begin(),
                         find(cells[jcell].mNewNeighbors.begin(), cells[jcell].mNewNeighbors.end(), mID));
        ENSURE(jface < nfj);
        ENSURE(cells[jcell].mNewNeighbors[jface] == mID);
        nvi =              mNewFaceVertices[iface].size();
        nvj = cells[jcell].mNewFaceVertices[jface].size();
        ENSURE(nvi >= minVerticesPerFace);
        ENSURE(nvj >= minVerticesPerFace);
      }
    }
  }
  END_CONTRACT_SCOPE;

}

//------------------------------------------------------------------------------
// Match up one of our faces with the corollary in a neighbor cell.  When
// we're done both cells will fill in their new vertex arrays for this face 
// and be consistent with one and other.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Cell<Dimension>::
matchFace(const unsigned iface,
          Cell<Dimension>& otherCell) {

  // Pre-conditions.
  const unsigned nfi = mNewNeighbors.size();
  REQUIRE(iface < nfi);
  REQUIRE(mNewFaceVertices[iface].size() >= minVerticesPerFace);
  REQUIRE2(count(otherCell.mNewNeighbors.begin(), otherCell.mNewNeighbors.end(),mID) == 1, 
           '\n' << "This cell:  " << this->dumpCell() << '\n' << "Other cell:  " << otherCell.dumpCell());

  // Locate the matching face in the other cell.
  const unsigned nfj = otherCell.mNewNeighbors.size();
  const unsigned jface = distance(otherCell.mNewNeighbors.begin(),
                                  find(otherCell.mNewNeighbors.begin(),
                                       otherCell.mNewNeighbors.end(),
                                       mID));
  CHECK(jface < nfj);
  CHECK(otherCell.mNewFaceVertices[jface].size() >= minVerticesPerFace);
  const unsigned nvi =           mNewFaceVertices[iface].size();
  const unsigned nvj = otherCell.mNewFaceVertices[jface].size();

  // Determine which cell is going to be in control of this face.
  unsigned i, j, k, nv, masterFace, slaveFace;
  Cell<Dimension>* masterCellPtr;
  Cell<Dimension>* slaveCellPtr;
  if ((nvi < nvj) or (nvi == nvj and mID < otherCell.mID)) {
    masterCellPtr = this;
    slaveCellPtr = &otherCell;
    nv = nvi;
    masterFace = iface;
    slaveFace = jface;
  } else {
    masterCellPtr = &otherCell;
    slaveCellPtr = this;
    nv = nvj;
    masterFace = jface;
    slaveFace = iface;
  }

  // Walk the vertices of the control face, and match to the closest one in
  // the slave face.
  slaveCellPtr->mNewFaceVertices[slaveFace] = vector<unsigned>();
  for (k = 0; k != nv; ++k) {
    i = masterCellPtr->mNewFaceVertices[masterFace][k];
    j = findMatchingVertex(masterCellPtr->mNewVertices[i], slaveCellPtr->mNewVertices);
    slaveCellPtr->mNewFaceVertices[slaveFace].push_back(j);
    if (masterCellPtr->mMinCellForVertex[i].first < slaveCellPtr->mMinCellForVertex[j].first) {
      slaveCellPtr->mMinCellForVertex[j] = masterCellPtr->mMinCellForVertex[i];
    } else {
      masterCellPtr->mMinCellForVertex[i] = slaveCellPtr->mMinCellForVertex[j];
    }
  }
  CHECK(slaveCellPtr->mNewFaceVertices[slaveFace].size() == masterCellPtr->mNewFaceVertices[masterFace].size());

  // Reverse the order of the face vertices on the slave generator.
  reverse(slaveCellPtr->mNewFaceVertices[slaveFace].begin(),
          slaveCellPtr->mNewFaceVertices[slaveFace].end());
}

//------------------------------------------------------------------------------
// Enforce consistency in the minimum cell for each of our vertices.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Cell<Dimension>::
findMinCellsForVertices(vector<Cell<Dimension> >& cells) {
  bool result = true;
  const unsigned ncells = cells.size();
  const unsigned nv = mNewVertices.size();
  REQUIRE(mMinCellForVertex.size() == nv);
  unsigned i, j, jgen, ni, nj;
  for (i = 0; i != nv; ++i) {
    jgen = mMinCellForVertex[i].first;
    j = mMinCellForVertex[i].second;
    CHECK(jgen < ncells);
    if (cells[jgen].mMinCellForVertex[j] != mMinCellForVertex[i]) {
      result = false;
      if (cells[jgen].mMinCellForVertex[j].first < jgen) {
        mMinCellForVertex[i] = cells[jgen].mMinCellForVertex[j];
      } else {
        cells[jgen].mMinCellForVertex[j] = mMinCellForVertex[i];
      }
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Finish off our new vertex and face information.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Cell<Dimension>::
lock(vector<Cell<Dimension> >& cells) {

  // Pre-conditions.
  unsigned i, j, k;
  const unsigned nv0 = mNewVertices.size();
  const unsigned nf0 = mNewFaceVertices.size();
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(mNewVertices.size() > minVerticesPerFace);
    REQUIRE(mNewFaceVertices.size() == mOldFaceVertices.size());
    REQUIRE(mNewNeighbors.size() == nf0);
    REQUIRE(mMinCellForVertex.size() == nv0);
  }
  END_CONTRACT_SCOPE;

  const unsigned ncells = cells.size();

  // Remove any deleted faces.
  vector<unsigned> faces2kill;
  for (i = 0; i != nf0; ++i) {
    if (mNewNeighbors[i] == DELETED) {
      CHECK(mNewFaceVertices[i].size() == 0);
      faces2kill.push_back(i);
    }
  }
  removeElements(mNewFaceVertices, faces2kill);
  removeElements(mNewNeighbors, faces2kill);
  CHECK(mNewFaceVertices.size() == mNewNeighbors.size());

  // Initialize the real node ID info, though it is not set yet.
  mRealNodeIDs = vector<unsigned>(nv0, UNSETID);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    ENSURE(mNewVertices.size() <= mOldVertices.size());
    ENSURE(mNewFaceVertices.size() <= mOldFaceVertices.size());
    for (i = 0; i != mNewFaceVertices.size(); ++i) {
      BOOST_FOREACH(j, mNewFaceVertices[i]) { ENSURE(j < mNewVertices.size()); }
    }
    ENSURE(mNewNeighbors.size() == mNewFaceVertices.size());
    BOOST_FOREACH(i, mNewNeighbors) { ENSURE(i != DELETED); }
    ENSURE(mMinCellForVertex.size() == mNewVertices.size());
    BOOST_FOREACH(i, mRealNodeIDs) { ENSURE(i == UNSETID); }
  }
  END_CONTRACT_SCOPE;
}

//------------------------------------------------------------------------------
// Set the "real" (mesh based) ID associated with one of our vertices.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Cell<Dimension>::
realNodeID(const unsigned ivertex, const unsigned ID) {
  REQUIRE(ivertex < mRealNodeIDs.size());
  REQUIRE(mRealNodeIDs.size() == mNewVertices.size());
  mRealNodeIDs[ivertex] = ID;
}

//------------------------------------------------------------------------------
// Dump this cells info to a string.
//------------------------------------------------------------------------------
template<typename Dimension>
string
Cell<Dimension>::
dumpCell() const {
  stringstream result;
  const unsigned nv0 = mOldVertices.size();
  const unsigned nf0 = mOldFaceVertices.size();
  const unsigned nv1 = mNewVertices.size();
  const unsigned nf1 = mNewFaceVertices.size();
  unsigned i, j;
  result << "ID: " << mID << '\n'
         << "Volume: " << mVolume << '\n'
         << "maxEdge: " << mMaxEdge << '\n';
  for (i = 0; i != nv0; ++i) {
    result << (i == 0 ? "Old vertices: " : "              ")
           << i << " : " << mOldVertices[i] << '\n';
  }
  for (i = 0; i != nv1; ++i) {
    result << (i == 0 ? "New vertices: " : "              ")
           << i << " : " << mNewVertices[i] << '\n';
  }
  for (i = 0; i != nf0; ++i) {
    result << (i == 0 ? "Old face vertices: " : "                   ")
           << i << " : ";
    for (j = 0; j != mOldFaceVertices[i].size(); ++j) result << mOldFaceVertices[i][j] << " ";
    result << '\n';
  }
  for (i = 0; i != nf1; ++i) {
    result << (i == 0 ? "New face vertices: " : "                   ")
           << i << " : ";
    for (j = 0; j != mNewFaceVertices[i].size(); ++j) result << mNewFaceVertices[i][j] << " ";
    result << '\n';
  }
  result << "Old neighbors: ";
  for (i = 0; i != nf0; ++i) result << mOldNeighbors[i] << " ";
  result << '\n'
         << "New neighbors: ";
  for (i = 0; i != nf1; ++i) result << mNewNeighbors[i] << " ";
  result << '\n' << "mVertexMap: ";
  for (i = 0; i != nv0; ++i) result << "(" << i << " " << mVertexMap[i] << ") ";
  result << '\n' << "Real Node IDs: ";
  for (i = 0; i != mRealNodeIDs.size(); ++i) result << mRealNodeIDs[i] << " ";
  result << '\n';
  return result.str();
}

//------------------------------------------------------------------------------
// Reduce the new face indices to point to the lowest common denominators.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Cell<Dimension>::
updateVertexMap() {
  const unsigned nv = mOldVertices.size();
  REQUIRE(mVertexMap.size() == nv);
  unsigned i, j;
  bool done = false;
  while (not done) {
    done = true;
    for (i = 0; i != nv; ++i) {
      j = mVertexMap[i];
      if (mVertexMap[j] != j) {
        done = false;
        mVertexMap[i] = mVertexMap[j];
      }
    }
  }

  // Post-conditions.
  for (i = 0; i != nv; ++i) {
    ENSURE(mVertexMap[i] == mVertexMap[mVertexMap[i]]);
  }
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename Dimension> const unsigned Cell<Dimension>::UNSETID = numeric_limits<unsigned>::max();
template<typename Dimension> const unsigned Cell<Dimension>::DELETED = numeric_limits<unsigned>::max() - 1U;

template<> const unsigned Cell<Dim<2> >::minVerticesPerFace = 2U;
template<> const unsigned Cell<Dim<3> >::minVerticesPerFace = 3U;

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class Cell<Dim<2> >;
template class Cell<Dim<3> >;
}
}
