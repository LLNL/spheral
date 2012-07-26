//---------------------------------Spheral++----------------------------------//
// GridCellIndex -- a simple little class to hold the gridcell indicies for
// a gridcell in each dimension.  Basically an integer version of the
// GeomVector<Dimension> class.
//
// Created by:  JMO, Mon Dec 27 10:47:34 PST 1999
//----------------------------------------------------------------------------//
#include <limits.h>
#include <math.h>

#include "GridCellIndex.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
namespace NeighborSpace {

using namespace std;

//------------------------------------------------------------------------------
// Define a global function which returns a vector of GridCellIndex objects
// within a given range of grid cells.
//------------------------------------------------------------------------------
vector<GridCellIndex<Dim<1> > >
GridCellIndexRange(const GridCellIndex<Dim<1> >& gridCellMin,
                   const GridCellIndex<Dim<1> >& gridCellMax) {
  int numCells = gridCellMax.xIndex() - gridCellMin.xIndex() + 1;
  CHECK(numCells > 0);

  vector<GridCellIndex<Dim<1> > > result;
  result.reserve(numCells);
  for (int ix = gridCellMin.xIndex(); ix <= gridCellMax.xIndex(); ++ix) {
    result.push_back(GridCellIndex<Dim<1> >(ix));
    CHECK(result.size() <= numCells);
  }

  return result;  
}

vector<GridCellIndex<Dim<2> > >
GridCellIndexRange(const GridCellIndex<Dim<2> >& gridCellMin,
                   const GridCellIndex<Dim<2> >& gridCellMax) {
  int dx = gridCellMax.xIndex() - gridCellMin.xIndex();
  int dy = gridCellMax.yIndex() - gridCellMin.yIndex();
  CHECK(dx >= 0 && dy >= 0);
  int numCells = (dx + 1)*(dy + 1);
  CHECK(numCells >= 0);

  vector<GridCellIndex<Dim<2> > > result;
  result.reserve(numCells);
  for (int ix = gridCellMin.xIndex(); ix <= gridCellMax.xIndex(); ++ix) {
    for (int iy = gridCellMin.yIndex(); iy <= gridCellMax.yIndex(); ++iy) {
      result.push_back(GridCellIndex<Dim<2> >(ix, iy));
      CHECK(result.size() <= numCells);
    }
  }

  return result;  
}

vector<GridCellIndex<Dim<3> > >
GridCellIndexRange(const GridCellIndex<Dim<3> >& gridCellMin,
                   const GridCellIndex<Dim<3> >& gridCellMax) {
  int dx = gridCellMax.xIndex() - gridCellMin.xIndex();
  int dy = gridCellMax.yIndex() - gridCellMin.yIndex();
  int dz = gridCellMax.zIndex() - gridCellMin.zIndex();
  CHECK(dx >= 0 && dy >= 0 && dz >= 0);
  int numCells = (dx + 1)*(dy + 1)*(dz + 1);
  CHECK(numCells >= 0);

  vector<GridCellIndex<Dim<3> > > result;
  result.reserve(numCells);
  for (int ix = gridCellMin.xIndex(); ix <= gridCellMax.xIndex(); ++ix) {
    for (int iy = gridCellMin.yIndex(); iy <= gridCellMax.yIndex(); ++iy) {
      for (int iz = gridCellMin.zIndex(); iz <= gridCellMax.zIndex(); ++iz) {
        result.push_back(GridCellIndex<Dim<3> >(ix, iy, iz));
        CHECK(result.size() <= numCells);
      }
    }
  }

  return result;  
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
// Instantiate the static variables.
// Set the maximum index value to 2^21/2, so that (mIndexMax - mIndexMin)^3 will
// be in the range of an unsigned long long (2^64 on the machines I'm using).

// 2^20
#define INDEXMAX 1048576 - 1
#define INDEXMIN -INDEXMAX - 1
#define XMULTIPLIER (INDEXMAX - INDEXMIN)

template<> const int GridCellIndex< Dim<1> >::mIndexMax = INDEXMAX; // 2^20 - 1
template<> const int GridCellIndex< Dim<1> >::mIndexMin = INDEXMIN;

template<> const int GridCellIndex< Dim<2> >::mIndexMax = INDEXMAX; // 2^20 - 1
template<> const int GridCellIndex< Dim<2> >::mIndexMin = INDEXMIN;

template<> const int GridCellIndex< Dim<3> >::mIndexMax = INDEXMAX; // 2^20 - 1
template<> const int GridCellIndex< Dim<3> >::mIndexMin = INDEXMIN;

}
}
