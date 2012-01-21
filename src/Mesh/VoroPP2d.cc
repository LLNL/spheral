//---------------------------------Spheral++----------------------------------//
// VoroPP2d
//----------------------------------------------------------------------------//
#include <iostream>
#include <sstream>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/tokenizer.hpp"

#include "VoroPP2d.hh"
#include "MeshConstructionUtilities.hh"
#include "Utilities/testBoxIntersection.hh"

#include "voro++_2d/voro++_2d.hh"

#include "Utilities/timingUtilities.hh"

namespace Spheral {
namespace MeshSpace {

using namespace std;
using std::min;
using std::max;
using std::abs;
using namespace boost;

//------------------------------------------------------------------------------
// Provide a few specialized classes to help reading Voro++ formatted data.
//------------------------------------------------------------------------------
struct VoroVector2d {
  double x, y;
  void fillVector(Dim<2>::Vector& vec) { vec.x(x); vec.y(y); }
};

inline
std::istream&
operator>>(std::istream& is, VoroVector2d& val) {
  char c;
  is >> c;
  CHECK(c == '(');
  std::string dummy = "";
  double* valptr = &val.x;
  while (c != ')') {
    is >> c;
    if (c == ',' or c == ')') {
      *valptr = boost::lexical_cast<double>(dummy);
      dummy = "";
      ++valptr;
    } else {
      dummy += c;
    }
  }
  return is;
}

//------------------------------------------------------------------------------
// VorPP::VoroPP(...)  2-D
//------------------------------------------------------------------------------
VoroPP2d::
VoroPP2d(const vector<Dim<2>::Vector>& generators,
         const Dim<2>::Vector& xmin,
         const Dim<2>::Vector& xmax,
         const unsigned nx,
         const unsigned ny):
  mNumGenerators(generators.size()),
  mNx(nx),
  mNy(ny),
  mGeneratorsPtr(&generators),
  mXmin(xmin),
  mXmax(xmax),
  mContainerPtr(new container_2d(xmin.x(), xmax.x(),
                                 xmin.y(), xmax.y(),
                                 mNx, mNy,           // The number of sub-regions in each dimension.
                                 false, false,       // Periodic
                                 16)),                // initial guess to reserve points per subblock?
  mGen2CellGen(generators.size(), 0),
  mVoro2GenOrder() {
  this->construct(generators, mNx, mNy);
}

//------------------------------------------------------------------------------
// VorPP::construct -- private method called by all constructors.
//------------------------------------------------------------------------------
void
VoroPP2d::
construct(const vector<Dim<2>::Vector>& generators,
          const unsigned nx,
          const unsigned ny) {

  // Add the generators to the container.
  unsigned i, j, ij, igen, ncells = nx*ny;
  vector<unsigned> numGensInCell(ncells, 0);
  unordered_map<unsigned, vector<unsigned> > gensInCell;
  for (igen = 0; igen != mNumGenerators; ++igen) {
    subRegion(generators[igen], i, j);
    ij = i + nx*j;
    CHECK(ij < ncells);
    mGen2CellGen[igen] = numGensInCell[ij];
    numGensInCell[ij] += 1;
    gensInCell[ij].push_back(igen);
    mContainerPtr->put(igen, generators[igen].x(), generators[igen].y());
  }

  // Construct the map telling us how to go from Voro++ ordering of the generators
  // back to our own.
  mVoro2GenOrder.reserve(generators.size());
  vector<unsigned>::const_iterator itr;
  for (j = 0; j != ny; ++j) {
    for (i = 0; i != nx; ++i) {
      ij = i + nx*j;
      CHECK(ij < ncells);
      if (numGensInCell[ij] > 0) {
        CHECK(gensInCell[ij].size() == numGensInCell[ij]);
        for (itr = gensInCell[ij].begin(); itr != gensInCell[ij].end(); ++itr) {
          mVoro2GenOrder.push_back(*itr);
        }
      }
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE;
  {
    ENSURE(mVoro2GenOrder.size() == generators.size());
    vector<unsigned> gen2voro(generators.size(), generators.size() + 1);
    for (i = 0; i != generators.size(); ++i) gen2voro[mVoro2GenOrder[i]] = i;
    ENSURE(generators.size() == 0 or
           *max_element(gen2voro.begin(), gen2voro.end()) == generators.size() - 1);
  }
  END_CONTRACT_SCOPE;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
VoroPP2d::
~VoroPP2d() {
  ENSURE(mContainerPtr.use_count() == 1);
}

//------------------------------------------------------------------------------
// Return the vertices for all cells.
// Vertex positions are listed in counter-clockwise order for each cell.
// Note we assign a unique ID to each vertex, but the vertices themselves are
// degenerate!
//------------------------------------------------------------------------------
void
VoroPP2d::
allCells(vector<Dim<2>::Vector>& vertices,
         vector<vector<unsigned> >& cellVertexIndices) const {

  // Clear the results.
  vertices = vector<Vector>();
  cellVertexIndices = vector<vector<unsigned> >(mNumGenerators);
  
  typedef boost::tuple<uint64_t, uint64_t, uint64_t> Key;
  Vector boxInv(1.0/(mXmax.x() - mXmin.x()),
                1.0/(mXmax.y() - mXmin.y()));

  // Walk the Voro++ generator order, converting back to our own.
  double xgen, ygen;
  unsigned ivorogen, igen, i, j, ij, k, nunique, nedges;
  Key hashi, hashj;
  voronoicell_2d cell;
//   Timing::Time t0 = Timing::currentTime();
  for (igen = 0; igen != mNumGenerators; ++igen) {
    CHECK(igen < cellVertexIndices.size());
    xgen = (*mGeneratorsPtr)[igen].x();
    ygen = (*mGeneratorsPtr)[igen].y();
    cell.init(mXmin.x(), mXmax.x(),
              mXmin.y(), mXmax.y());

    // Compute the current cell.
    subRegion((*mGeneratorsPtr)[igen], i, j);
    ij = i + mNx*j;
    // cerr << "Generator:  " << igen << " " << (*mGeneratorsPtr)[igen] << " " << i << " " << j << " " << ij << " " << mGen2CellGen[igen] << endl;
    mContainerPtr->compute_cell_sphere(cell, i, j, ij, mGen2CellGen[igen], xgen, ygen);
    // cerr << "   num nodes : " << cell.p << endl;
    CHECK2(cell.p >= 2,
           "Bad number of nodes for generator: " << igen << " @ " << (*mGeneratorsPtr)[igen] << " : " << cell.p);
    
    // If Voro++ has inflicted us with a degenerate 2 node situation, fill in an extra point.
    if (cell.p == 2) {
      const Vector n1 = boundPointWithinBox(Vector(xgen + 0.5*cell.pts[0],
                                                   ygen + 0.5*cell.pts[1]), mXmin, mXmax);
      const Vector n2 = boundPointWithinBox(Vector(xgen + 0.5*cell.pts[2],
                                                   ygen + 0.5*cell.pts[3]), mXmin, mXmax);
      const Vector n12 = n2 - n1;
      const Vector n3 = boundPointWithinBox(n1 + 1.0*Vector(n12.y(), -(n12.x())), mXmin, mXmax);
      vertices.push_back(n1);
      vertices.push_back(n2);
      vertices.push_back(n3);
      cellVertexIndices[igen].push_back(vertices.size() - 3);
      cellVertexIndices[igen].push_back(vertices.size() - 2);
      cellVertexIndices[igen].push_back(vertices.size() - 1);

    } else {

      // Extract this cell's vertex info.
      // cerr << "   --> ";
      // for (i = 0; i != cell.p; ++i) {
      //   cellVertexIndices[igen].push_back(boundPointWithinBox(Vector(xgen + 0.5*cell.pts[2*i],
      //                                                           ygen + 0.5*cell.pts[2*i + 1]), mXmin, mXmax));
      //   cerr << " " << cellVertexIndices[igen].back();
      // }
      k = 0;
      do {
        vertices.push_back(boundPointWithinBox(Vector(xgen + 0.5*cell.pts[2*k],
                                                      ygen + 0.5*cell.pts[2*k + 1]), mXmin, mXmax));
        cellVertexIndices[igen].push_back(vertices.size() - 1);
        k = cell.ed[2*k];
        // cerr << " " << cellVertexIndices[igen].back();
      } while (k != 0);
      // cerr << endl;
      CHECK(cellVertexIndices[igen].size() == cell.p);
    }

    // Check that this cell has enough unique vertices.
    nunique = 0U;
    nedges = cellVertexIndices[igen].size();
    for (i = 0; i != nedges; ++i) {
      j = (i + 1) % nedges;
      hashi = hashPosition(vertices[cellVertexIndices[igen][i]], mXmin, mXmax, boxInv);
      hashj = hashPosition(vertices[cellVertexIndices[igen][j]], mXmin, mXmax, boxInv);
      if (compare(hashi, hashj) != 0) {
        ++nunique;
      } else if (nedges - i - 1 + nunique < 3) {
        incrementTuple(hashj, Key(1U, 1U, 0U));
        vertices[cellVertexIndices[igen][j]] = quantizedPosition(hashj, mXmin, mXmax);
        ++nunique;
      }
    }
    CHECK2(nunique >= 3, nunique << " " << cellVertexIndices[igen].size());
  }
//   if (Process::getRank() == 0) cerr << "VoroPP2d:: required " 
//                                     << Timing::difference(t0, Timing::currentTime())
//                                     << " seconds to read all vertex info for "
//                                     << mNumGenerators << " generators." << endl;
}

//------------------------------------------------------------------------------
// Find the indices describing the sub-region containing point p.
//------------------------------------------------------------------------------
void
VoroPP2d::
subRegion(const Dim<2>::Vector& p, unsigned& i, unsigned& j) const {
  REQUIRE2(testPointInBox(p, mXmin, mXmax),
           "Point not in box:  " << p << " not in [" << mXmin << ", " << mXmax << "]");
  const Vector delta((mXmax.x() - mXmin.x())/mNx,
                     (mXmax.y() - mXmin.y())/mNy);
  i = max(0U, min(mNx - 1, unsigned((p.x() - mXmin.x())/delta.x())));
  j = max(0U, min(mNy - 1, unsigned((p.y() - mXmin.y())/delta.y())));
}

}
}
