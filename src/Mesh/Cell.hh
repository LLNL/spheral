//---------------------------------Spheral++----------------------------------//
// Cell
// 
// This class is designed to take (possibly) inconsistent cell by cell 
// information from tools like Voro++ and help tie them together into 
// consistent connectivity between cells.
//
// Created by JMO, Sun Sep 11 18:32:22 PDT 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_Cell__
#define __Spheral_Cell__

#include <vector>
#include <map>
#include <string>

namespace Spheral {
namespace MeshSpace {

template<typename Dimension>
class Cell {
  //--------------------------- Public Interface ---------------------------//
public:
  typedef typename Dimension::Vector Vector;
  static const unsigned UNSETID;
  static const unsigned DELETED;
  static const unsigned minVerticesPerFace;

  // Constructors.
  Cell();
  Cell(const unsigned ID,
       const double vol,
       const Vector& centroid,
       const std::vector<Vector>& vertices,
       const std::vector<std::vector<unsigned> >& faceVertices,
       const std::vector<unsigned>& neighbors,
       const double edgeTol);

  // Our ID, volume, & centroid.
  unsigned ID() const;
  double volume() const;
  const Vector& centroid() const;

  // Access original cell info.
  unsigned numOldVertices() const;
  unsigned numOldFaces() const;
  const std::vector<Vector>& oldVertices() const;
  const std::vector<std::vector<unsigned> >& oldFaceVertices() const;
  const std::vector<unsigned>& oldNeighbors() const;

  // Access new cell info.
  unsigned numNewVertices() const;
  unsigned numNewFaces() const;
  const std::vector<Vector>& newVertices() const;
  const std::vector<std::vector<unsigned> >& newFaceVertices() const;
  const std::vector<unsigned>& newNeighbors() const;

  // The map of old -> new vertices.
  const std::vector<unsigned>& vertexMap() const;

  // The minimum cell ID which shares the given local vertex.
  unsigned minCellForVertex(const unsigned i) const;
  unsigned localVertexForMinCell(const unsigned i) const;

  // Make cell connectivity across faces consistent.
  void cullDegenerateNeighbors(std::vector<Cell<Dimension> >& cells);

  // Match up one of our faces with another cell, making each cell
  // consistent for this face.
  void matchFace(const unsigned iface, 
                 Cell<Dimension>& otherCell);

  // Enforce consistency in the minimum cell for each of our vertices.
  bool findMinCellsForVertices(std::vector<Cell<Dimension> >& cells);

  // Finish off our internal info,so all new cell info is set.
  void lock(std::vector<Cell<Dimension> >& cells);

  // Access the "real" (mesh based) ID for one of our vertices.
  unsigned realNodeID(const unsigned ivertex) const;

  // set the "real" (mesh based) ID associated with one of our vertices.
  void realNodeID(const unsigned ivertex, const unsigned ID);

  // Dump out the basic state of the cell to a string.
  std::string dumpCell() const;

  //--------------------------- Private Interface ---------------------------//
private:
  unsigned mID;
  double mVolume, mMaxEdge;
  Vector mCentroid;
  std::vector<Vector> mOldVertices, mNewVertices;
  std::vector<std::vector<unsigned> > mOldFaceVertices, mNewFaceVertices;
  std::vector<unsigned> mOldNeighbors, mNewNeighbors;
  std::vector<unsigned> mVertexMap;
  std::vector<std::pair<unsigned, unsigned> > mMinCellForVertex; // (jcell, jvertex)
  std::vector<unsigned> mRealNodeIDs;

  // Update the vertex map so everything points to the lowest common set.
  void updateVertexMap();
};

}
}

#include "CellInline.hh"

#endif
