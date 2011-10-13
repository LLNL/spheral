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
  mCellsForVertex(vertices.size()),
  mRealNodeIDs() {

  // Pre-conditions.
  const unsigned nf = faceVertices.size();
  unsigned iface, i;
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(neighbors.size() == nf);
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

  // Intitialize the cells per vertex to include ourself.
  for (i = 0; i != nv; ++i) {
    mCellsForVertex[i][mID] = i;
  }

  // Any faces on the boundary of the tesselation should go ahead and assign the
  // original face vertices as the new ones.
  for (iface = 0; iface != nf; ++iface) {
    if (mNewNeighbors[iface] == UNSETID) mNewFaceVertices[iface] = mOldFaceVertices[iface];
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    ENSURE(mVolume > 0.0);
    ENSURE(mMaxEdge > 0.0);
    ENSURE(mOldVertices == vertices);
    ENSURE(mNewVertices.size() == 0);
    ENSURE(mOldFaceVertices == faceVertices);
    ENSURE(mNewFaceVertices.size() == mOldFaceVertices.size());
    ENSURE(mOldNeighbors == neighbors);
    ENSURE(mNewNeighbors == mOldNeighbors);
    for (i = 0; i != mNewFaceVertices.size(); ++i) {
      ENSURE((mNewNeighbors[i] == UNSETID and mNewFaceVertices[i] == mOldFaceVertices[i]) or
             mNewFaceVertices[i].size() == 0);
    }
    ENSURE(mVertexMap.size() == vertices.size());
    ENSURE(mCellsForVertex.size() == vertices.size());
    for (i = 0; i != mCellsForVertex.size(); ++i) {
      ENSURE(mCellsForVertex[i].size() == 1);
      ENSURE(mCellsForVertex[i].find(mID) != mCellsForVertex[i].end());
      ENSURE(mCellsForVertex[i][mID] == i);
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
  CHECK(mNewNeighbors.size() == nfi);
  CHECK(mOldFaceVertices.size() == nfi);
  for (iface = 0; iface != nfi; ++iface) {
    jcell = mNewNeighbors[iface];

    // Check if both cells have valid info for this neighbor relation.
    if (jcell < ncells) {
      nfj = cells[jcell].mNewNeighbors.size();
      jface = distance(cells[jcell].mNewNeighbors.begin(),
                       find(cells[jcell].mNewNeighbors.begin(), cells[jcell].mNewNeighbors.end(), mID));

      // If the other cell doesn't have us listed as a neighbor, we have to remove
      // this face and neighbor association.
      if (jface == nfj) {
        mNewNeighbors[iface] = DELETED;

      } else {
        // If we got here both cells acknowledge the cross-face relationship,
        // so ensure we have enough vertices in the face to define a 
        // valid area.
        CHECK(jface < nfj);
        set<unsigned> uniquei, uniquej;
        nvi = mOldFaceVertices[iface].size();
        nvj = cells[jcell].mOldFaceVertices[jface].size();
        for (i = 0; i != nvi; ++i) uniquei.insert(mVertexMap[mOldFaceVertices[iface][i]]);
        for (j = 0; j != nvj; ++j) uniquej.insert(cells[jcell].mVertexMap[cells[jcell].mOldFaceVertices[jface][j]]);
        if (min(nvi, nvj) < minVerticesPerFace) {
          // A least one of these cells doesn't have enough unique vertices on this face
          // to form a valid area, so delete the face from both cells.
          CHECK(max(nvi, nvj) < minVerticesPerFace);  // Hopefully we don't have to change the vertices!
          mNewNeighbors[iface] = DELETED;
          cells[jcell].mNewNeighbors[iface] = DELETED;
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
        set<unsigned> uniquei, uniquej;
        nvi = mOldFaceVertices[iface].size();
        nvj = cells[jcell].mOldFaceVertices[jface].size();
        for (i = 0; i != nvi; ++i) uniquei.insert(mVertexMap[mOldFaceVertices[iface][i]]);
        for (j = 0; j != nvj; ++j) uniquej.insert(cells[jcell].mVertexMap[cells[jcell].mOldFaceVertices[jface][j]]);
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
  REQUIRE(mNewFaceVertices[iface].size() == 0);
  REQUIRE2(count(otherCell.mNewNeighbors.begin(), otherCell.mNewNeighbors.end(),mID) == 1, 
           '\n' << "This cell:  " << this->dumpCell() << '\n' << "Other cell:  " << otherCell.dumpCell());

  // Locate the matching face in the other cell.
  const unsigned nfj = otherCell.mNewNeighbors.size();
  const unsigned jface = distance(otherCell.mNewNeighbors.begin(),
                                  find(otherCell.mNewNeighbors.begin(),
                                       otherCell.mNewNeighbors.end(),
                                       mID));
  CHECK(jface < nfj);
  CHECK(otherCell.mNewFaceVertices[jface].size() == 0);
  const unsigned nvi = mOldFaceVertices[iface].size();
  const unsigned nvj = otherCell.mOldFaceVertices[jface].size();

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
  for (k = 0; k != nv; ++k) {
    i = masterCellPtr->mOldFaceVertices[masterFace][k];
    j = findMatchingVertex(masterCellPtr->mOldVertices[i], 
                           slaveCellPtr->mOldVertices,
                           slaveCellPtr->mOldFaceVertices[slaveFace]);
    masterCellPtr->mNewFaceVertices[masterFace].push_back(i);
    slaveCellPtr->mNewFaceVertices[slaveFace].push_back(j);
    masterCellPtr->mCellsForVertex[i][slaveCellPtr->mID] = j;
    slaveCellPtr->mCellsForVertex[j][masterCellPtr->mID] = i;
  }
  CHECK(masterCellPtr->mNewFaceVertices[masterFace].size() == nv);
  CHECK(slaveCellPtr->mNewFaceVertices[slaveFace].size() == nv);

  // Reverse the order of the face vertices on the slave generator.
  reverse(slaveCellPtr->mNewFaceVertices[slaveFace].begin(),
          slaveCellPtr->mNewFaceVertices[slaveFace].end());
}

//------------------------------------------------------------------------------
// Enforce symmetry in the set of known cells that share our vertices.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Cell<Dimension>::
distributeSharedVertices(vector<Cell<Dimension> >& cells) {
  bool result = true;
  const unsigned ncells = cells.size();
  const unsigned nv = mOldVertices.size();
  REQUIRE(mCellsForVertex.size() == nv);
  unsigned i, j, jgen, ni, nj;
  for (i = 0; i != nv; ++i) {
    ni = mCellsForVertex[i].size();
    for (map<unsigned, unsigned>::const_iterator itr = mCellsForVertex[i].begin();
         itr != mCellsForVertex[i].end();
         ++itr) {
      jgen = itr->first;
      j = itr->second;
      CHECK(jgen < ncells and j < cells[jgen].mCellsForVertex.size());
      nj = cells[jgen].mCellsForVertex[j].size();
      mCellsForVertex[i].insert(cells[jgen].mCellsForVertex[j].begin(),
                                cells[jgen].mCellsForVertex[j].end());
      cells[jgen].mCellsForVertex[j].insert(mCellsForVertex[i].begin(),
                                            mCellsForVertex[i].end());
      if (cells[jgen].mCellsForVertex[j].size() != nj) result = false;
    }
    if (mCellsForVertex[i].size() != ni) result = false;
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
  const unsigned nv0 = mOldVertices.size();
  const unsigned nf0 = mOldFaceVertices.size();
  BEGIN_CONTRACT_SCOPE;
  {
    REQUIRE(mNewVertices.size() == 0);
    REQUIRE(mNewFaceVertices.size() <= nf0);
    REQUIRE(mNewNeighbors.size() == nf0);
    REQUIRE(mVertexMap.size() == nv0);
    REQUIRE(mCellsForVertex.size() == nv0);
    for (i = 0; i != nv0; ++i) {
      REQUIRE(mVertexMap[i] < nv0);
      REQUIRE(mVertexMap[i] <= i);
      REQUIRE(mVertexMap[i] == i or mVertexMap[mVertexMap[i]] == mVertexMap[i]);
    }
  }
  END_CONTRACT_SCOPE;

  const unsigned ncells = cells.size();

  // // Blago!
  // cerr << "Dumping INITIAL cells for vertices for cell " << mID << endl;
  // for (i = 0; i != nv0; ++i) {
  //   cerr << "  Vertex " << i << " @ " << mOldVertices[i] << endl;
  //   for (map<unsigned, unsigned>::const_iterator itr = mCellsForVertex[i].begin();
  //        itr != mCellsForVertex[i].end();
  //        ++itr) {
  //     cerr << "       " << itr->first << " " << itr->second << endl;
  //   }
  // }
  // // Blago!

  // Create the new vertices and update the vertex map to point from the 
  // old vertex numbering to the new.
  unsigned nv1 = 0;
  for (i = 0; i != nv0; ++i) {
    if (mVertexMap[i] == i) {
      mNewVertices.push_back(mOldVertices[i]);
      mVertexMap[i] = nv1++;
      CHECK(mNewVertices.size() == nv1);
    } else {
      k = mVertexMap[mVertexMap[i]];
      CHECK(k < nv1);
      CHECK(mVertexMap[k] = k);
      mVertexMap[i] = k;
      mCellsForVertex[k].insert(mCellsForVertex[i].begin(), mCellsForVertex[i].end());
    }
  }
  mCellsForVertex.resize(nv1);
  BEGIN_CONTRACT_SCOPE;
  {
    CHECK(mNewVertices.size() == nv1);
    CHECK(mCellsForVertex.size() == nv1);
    BOOST_FOREACH(i, mVertexMap) { CHECK(i < nv1); }
    for (i = 0; i != nv1; ++i) {
      CHECK(mCellsForVertex[i].size() > 0);
    }
  }
  END_CONTRACT_SCOPE;

  // Update the new face info.
  vector<unsigned> faces2kill;
  for (i = 0; i != nf0; ++i) {
    if (mNewNeighbors[i] == DELETED) {
      CHECK(mNewFaceVertices[i].size() == 0);
      faces2kill.push_back(i);
    } else {
      CHECK(mNewFaceVertices[i].size() >= minVerticesPerFace);
      for (j = 0; j != mNewFaceVertices[i].size(); ++j) {
        mNewFaceVertices[i][j] = mVertexMap[mNewFaceVertices[i][j]];
      }
    }
  }
  removeElements(mNewFaceVertices, faces2kill);
  removeElements(mNewNeighbors, faces2kill);
  CHECK(mNewFaceVertices.size() == mNewNeighbors.size());

  // // Blago!
  // cerr << "Dumping FINAL cells for vertices for cell " << mID << endl;
  // for (i = 0; i != nv1; ++i) {
  //   cerr << "  Vertex " << i << " @ " << mNewVertices[i] << endl;
  //   for (typename map<unsigned, unsigned>::const_iterator itr = mCellsForVertex[i].begin();
  //        itr != mCellsForVertex[i].end();
  //        ++itr) {
  //     cerr << "       " << itr->first << " " << itr->second << endl;
  //   }
  // }
  // // Blago!

  // Walk all neighbor cells and update their vertex info for us.
  unsigned nvj;
  for (i = 0; i != nv1; ++i) {
    for (typename map<unsigned, unsigned>::const_iterator itr = mCellsForVertex[i].begin();
         itr != mCellsForVertex[i].end();
         ++itr) {
      CHECK(itr->first < ncells);
      if (itr->first != mID) {
        Cell<Dimension>& otherCell = cells[itr->first];
        nvj = otherCell.mCellsForVertex.size();
        for (j = 0; j != nvj; ++j) {
          map<unsigned, unsigned> newCells;
          for (map<unsigned, unsigned>::const_iterator otherItr = otherCell.mCellsForVertex[j].begin();
               otherItr != otherCell.mCellsForVertex[j].end();
               ++otherItr) {
            if (otherItr->first == mID) {
              newCells[mID] = mVertexMap[otherItr->second];
            } else {
              newCells.insert(*otherItr);
            }
          }
          otherCell.mCellsForVertex[j] = newCells;
        }
      }
    }
  }

  // Initialize the real node ID info, though it is not set yet.
  mRealNodeIDs = vector<unsigned>(nv1, UNSETID);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    ENSURE(mNewVertices.size() <= mOldVertices.size());
    ENSURE(mNewFaceVertices.size() <= mOldFaceVertices.size());
    for (i = 0; i != mNewFaceVertices.size(); ++i) {
      BOOST_FOREACH(j, mNewFaceVertices[i]) { j < mNewVertices.size(); }
    }
    ENSURE(mNewNeighbors.size() == mNewFaceVertices.size());
    BOOST_FOREACH(i, mNewNeighbors) { ENSURE(i != DELETED); }
    ENSURE(mCellsForVertex.size() == mNewVertices.size());
    for (i = 0; i != mCellsForVertex.size(); ++i) {
      ENSURE(mCellsForVertex[i].size() >= 1);
    }
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
// Same as above, but using the vertex numbering of a neighboring cell.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Cell<Dimension>::
realNodeID(const unsigned jvertex, 
           const Cell<Dimension>& jcell,
           const unsigned ID) {
  REQUIRE(mRealNodeIDs.size() == mNewVertices.size());
  REQUIRE(jvertex < jcell.mCellsForVertex.size());

  // Find our vertex in the other cells set.
  const unsigned nj = jcell.mCellsForVertex[jvertex].size();
  const map<unsigned, unsigned>::const_iterator itr = jcell.mCellsForVertex[jvertex].find(mID);
  CHECK(itr != jcell.mCellsForVertex[jvertex].end());
  const unsigned ivertex = itr->second;
  this->realNodeID(ivertex, ID);
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
    result << (i == 0 ? "Old face vertices: " : "                   ");
    for (j = 0; j != mOldFaceVertices[i].size(); ++j) result << mOldFaceVertices[i][j] << " ";
    result << '\n';
  }
  for (i = 0; i != nf1; ++i) {
    result << (i == 0 ? "New face vertices: " : "                   ");
    for (j = 0; j != mNewFaceVertices[i].size(); ++j) result << mNewFaceVertices[i][j] << " ";
    result << '\n';
  }
  result << "Old neighbors: ";
  for (i = 0; i != nf0; ++i) result << mOldNeighbors[i] << " ";
  result << '\n'
         << "New neighbors: ";
  for (i = 0; i != nf1; ++i) result << mNewNeighbors[i] << " ";
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
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<> const unsigned Cell<Dim<3> >::UNSETID = numeric_limits<unsigned>::max();
template<> const unsigned Cell<Dim<3> >::DELETED = numeric_limits<unsigned>::max() - 1U;
template<> const unsigned Cell<Dim<3> >::minVerticesPerFace = 3U;

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
template class Cell<Dim<3> >;
}
}
