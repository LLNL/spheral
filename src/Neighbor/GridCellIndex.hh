//---------------------------------Spheral++----------------------------------//
// GridCellIndex -- a simple little class to hold the gridcell indicies for
// a gridcell in each dimension.  Basically an integer version of the
// GeomVector<Dimension> class.
//
// Created by:  JMO, Mon Dec 27 10:47:34 PST 1999
//----------------------------------------------------------------------------//

#ifndef __Spheral_GridCellIndex_hh__
#define __Spheral_GridCellIndex_hh__

#ifndef __GCCXML__
#include <vector>
#else
#include "fakestl.hh"
#endif

#include "GridCellIndexBase.hh"

namespace Spheral {
namespace NeighborSpace {

template<typename Dimension> class GridCellPlane;

template<typename Dimension>
class GridCellIndex: public GridCellIndexBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  GridCellIndex();
  GridCellIndex(int xIndex);
  GridCellIndex(int xIndex, int yIndex);
  GridCellIndex(int xIndex, int yIndex, int zIndex);
  GridCellIndex(const GridCellIndex& rhs);
  ~GridCellIndex();

  int xIndex() const;
  int yIndex() const;
  int zIndex() const;

  void xIndex(int x);
  void yIndex(int y);
  void zIndex(int z);

  void setIndicies(int xIndex);
  void setIndicies(int xIndex, int yIndex);
  void setIndicies(int xIndex, int yIndex, int zIndex);

  int operator()(int i) const;
  int& operator()(int i);

  GridCellIndex& operator=(const GridCellIndex& rhs);
  GridCellIndex& operator=(int rhs);

  GridCellIndex operator-() const;

  GridCellIndex operator+(const GridCellIndex& rhs) const;
  GridCellIndex operator-(const GridCellIndex& rhs) const;

  GridCellIndex& operator+=(const GridCellIndex& rhs);
  GridCellIndex& operator-=(const GridCellIndex& rhs);

  GridCellIndex operator+(const int rhs) const;
  GridCellIndex operator-(const int rhs) const;
  GridCellIndex operator*(const int rhs) const;
  GridCellIndex operator/(const int rhs) const;

  GridCellIndex& operator+=(const int rhs);
  GridCellIndex& operator-=(const int rhs);

  int dot(const GridCellIndex& rhs) const;

  int compare(const GridCellIndex& rhs) const;
  bool operator==(const GridCellIndex& rhs) const;
  bool operator!=(const GridCellIndex& rhs) const;
  bool operator<(const GridCellIndex& rhs) const;
  bool operator>(const GridCellIndex& rhs) const;
  bool operator<=(const GridCellIndex& rhs) const;
  bool operator>=(const GridCellIndex& rhs) const;

  bool operator<(const GridCellPlane<Dimension>& plane) const;
  bool operator>(const GridCellPlane<Dimension>& plane) const;
  bool operator<=(const GridCellPlane<Dimension>& plane) const;
  bool operator>=(const GridCellPlane<Dimension>& plane) const;

  bool inRange(const GridCellIndex& minGridCell, const GridCellIndex& maxGridCell) const;

  double magnitude() const;
  int magnitude2() const;
  int minElement() const;
  int maxElement() const;
  int sumElements() const;
  int productElements() const;

  int indexMin() const;
  int indexMax() const;

private:
  //--------------------------- Private Interface ---------------------------//
  static const int mIndexMin;
  static const int mIndexMax;

};

// Declare explicit specializations.
template<> GridCellIndex<Dim<1> >::GridCellIndex(int);
template<> GridCellIndex<Dim<2> >::GridCellIndex(int, int);
template<> GridCellIndex<Dim<3> >::GridCellIndex(int, int, int);

template<> int GridCellIndex<Dim<2> >::yIndex() const;
template<> int GridCellIndex<Dim<3> >::yIndex() const;
template<> int GridCellIndex<Dim<3> >::zIndex() const;

template<> void GridCellIndex<Dim<2> >::yIndex(int);
template<> void GridCellIndex<Dim<3> >::yIndex(int);
template<> void GridCellIndex<Dim<3> >::zIndex(int);

template<> void GridCellIndex<Dim<1> >::setIndicies(int);
template<> void GridCellIndex<Dim<2> >::setIndicies(int, int);
template<> void GridCellIndex<Dim<3> >::setIndicies(int, int, int);

template<> GridCellIndex<Dim<1> > GridCellIndex<Dim<1> >::operator-() const;
template<> GridCellIndex<Dim<2> > GridCellIndex<Dim<2> >::operator-() const;
template<> GridCellIndex<Dim<3> > GridCellIndex<Dim<3> >::operator-() const;

template<> GridCellIndex<Dim<1> > GridCellIndex<Dim<1> >::operator+(const GridCellIndex<Dim<1> >&) const;
template<> GridCellIndex<Dim<2> > GridCellIndex<Dim<2> >::operator+(const GridCellIndex<Dim<2> >&) const;
template<> GridCellIndex<Dim<3> > GridCellIndex<Dim<3> >::operator+(const GridCellIndex<Dim<3> >&) const;

template<> GridCellIndex<Dim<1> > GridCellIndex<Dim<1> >::operator-(const GridCellIndex<Dim<1> >&) const;
template<> GridCellIndex<Dim<2> > GridCellIndex<Dim<2> >::operator-(const GridCellIndex<Dim<2> >&) const;
template<> GridCellIndex<Dim<3> > GridCellIndex<Dim<3> >::operator-(const GridCellIndex<Dim<3> >&) const;

template<> GridCellIndex<Dim<1> > GridCellIndex<Dim<1> >::operator+(int) const;
template<> GridCellIndex<Dim<2> > GridCellIndex<Dim<2> >::operator+(int) const;
template<> GridCellIndex<Dim<3> > GridCellIndex<Dim<3> >::operator+(int) const;

template<> GridCellIndex<Dim<1> > GridCellIndex<Dim<1> >::operator-(int) const;
template<> GridCellIndex<Dim<2> > GridCellIndex<Dim<2> >::operator-(int) const;
template<> GridCellIndex<Dim<3> > GridCellIndex<Dim<3> >::operator-(int) const;

template<> GridCellIndex<Dim<1> > GridCellIndex<Dim<1> >::operator*(int) const;
template<> GridCellIndex<Dim<2> > GridCellIndex<Dim<2> >::operator*(int) const;
template<> GridCellIndex<Dim<3> > GridCellIndex<Dim<3> >::operator*(int) const;

template<> GridCellIndex<Dim<1> > GridCellIndex<Dim<1> >::operator/(int) const;
template<> GridCellIndex<Dim<2> > GridCellIndex<Dim<2> >::operator/(int) const;
template<> GridCellIndex<Dim<3> > GridCellIndex<Dim<3> >::operator/(int) const;

template<> int GridCellIndex< Dim<1> >::dot(const GridCellIndex< Dim<1> >&) const;
template<> int GridCellIndex< Dim<2> >::dot(const GridCellIndex< Dim<2> >&) const;
template<> int GridCellIndex< Dim<3> >::dot(const GridCellIndex< Dim<3> >&) const;

template<> int GridCellIndex<Dim<1> >::compare(const GridCellIndex<Dim<1> >&) const;

template<> bool GridCellIndex<Dim<1> >::operator==(const GridCellIndex<Dim<1> >&) const;
template<> bool GridCellIndex<Dim<2> >::operator==(const GridCellIndex<Dim<2> >&) const;
template<> bool GridCellIndex<Dim<3> >::operator==(const GridCellIndex<Dim<3> >&) const;

template<> bool GridCellIndex<Dim<1> >::inRange(const GridCellIndex<Dim<1> >&,
                                                const GridCellIndex<Dim<1> >&) const;
template<> bool GridCellIndex<Dim<2> >::inRange(const GridCellIndex<Dim<2> >&,
                                                const GridCellIndex<Dim<2> >&) const;
template<> bool GridCellIndex<Dim<3> >::inRange(const GridCellIndex<Dim<3> >&,
                                                const GridCellIndex<Dim<3> >&) const;

template<> int GridCellIndex< Dim<1> >::magnitude2() const;
template<> int GridCellIndex< Dim<2> >::magnitude2() const;
template<> int GridCellIndex< Dim<3> >::magnitude2() const;

template<> int GridCellIndex<Dim<1> >::minElement() const;
template<> int GridCellIndex<Dim<2> >::minElement() const;
template<> int GridCellIndex<Dim<3> >::minElement() const;

template<> int GridCellIndex<Dim<1> >::maxElement() const;
template<> int GridCellIndex<Dim<2> >::maxElement() const;
template<> int GridCellIndex<Dim<3> >::maxElement() const;

template<> int GridCellIndex<Dim<1> >::sumElements() const;
template<> int GridCellIndex<Dim<2> >::sumElements() const;
template<> int GridCellIndex<Dim<3> >::sumElements() const;

template<> int GridCellIndex<Dim<1> >::productElements() const;
template<> int GridCellIndex<Dim<2> >::productElements() const;
template<> int GridCellIndex<Dim<3> >::productElements() const;

}
}

#ifndef __GCCXML__
#include "GridCellIndexInline.hh"
#endif

#else

namespace Spheral {
namespace NeighborSpace {

// Forward declaration.
template<typename Dimension> class GridCellIndex;

std::vector<GridCellIndex<Dim<1> > >
GridCellIndexRange(const GridCellIndex<Dim<1> >& gridCellMin,
                   const GridCellIndex<Dim<1> >& gridCellMax);
std::vector<GridCellIndex<Dim<2> > >
GridCellIndexRange(const GridCellIndex<Dim<2> >& gridCellMin,
                   const GridCellIndex<Dim<2> >& gridCellMax);
std::vector<GridCellIndex<Dim<3> > >
GridCellIndexRange(const GridCellIndex<Dim<3> >& gridCellMin,
                   const GridCellIndex<Dim<3> >& gridCellMax);

void
incrementGridCell(GridCellIndex<Dim<1> >& gridCell,
                  const GridCellIndex<Dim<1> >& gridCellMin,
                  const GridCellIndex<Dim<1> >& gridCellMax);
void
incrementGridCell(GridCellIndex<Dim<2> >& gridCell,
                  const GridCellIndex<Dim<2> >& gridCellMin,
                  const GridCellIndex<Dim<2> >& gridCellMax);
void
incrementGridCell(GridCellIndex<Dim<3> >& gridCell,
                  const GridCellIndex<Dim<3> >& gridCellMin,
                  const GridCellIndex<Dim<3> >& gridCellMax);

template<typename Dimension>
GridCellIndex<Dimension> operator+(int lhs, const GridCellIndex<Dimension>& rhs);
template<typename Dimension>
GridCellIndex<Dimension> operator-(int lhs, const GridCellIndex<Dimension>& rhs);
template<typename Dimension>
GridCellIndex<Dimension> operator*(int lhs, const GridCellIndex<Dimension>& rhs);

}
}

#endif
