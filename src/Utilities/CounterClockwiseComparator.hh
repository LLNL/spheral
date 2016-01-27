#ifndef __Spheral_CounterClockwiseComparator__
#define __Spheral_CounterClockwiseComparator__

#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Comparator to help sorting the elements of a container of Vector 
// counter-clockwise about a position.
// Provides two types of comparisons:
//   1.  pos1 -> pos2
//   2.  pos[index1] -> pos[index2]
//------------------------------------------------------------------------------
// General 3-D definition.
template<typename Vector, typename Container>
struct CounterClockwiseComparator {
  CounterClockwiseComparator(const Container& positions,
                             const Vector& origin,
                             const Vector& normal):
    mPositions(positions),
    mOrigin(origin),
    mNormal(normal) {}
  bool operator()(const Vector& pos1, const Vector& pos2) {
    return ((pos1 - mOrigin).cross(pos2 - mOrigin).dot(mNormal) > 0.0);
  }
  bool operator()(const unsigned id1, const unsigned id2) {
    return (*this)(mPositions.at(id1), mPositions.at(id2));
  }
  const Container& mPositions;
  const Vector mOrigin, mNormal;
};

// Specialized 2-D definition (normal implied).
template<typename Container>
struct CounterClockwiseComparator<Dim<2>::Vector, Container> {
  CounterClockwiseComparator(const Container& positions,
                             const Dim<2>::Vector& origin):
    mPositions(positions),
    mOrigin(origin) {}
  bool operator()(const Dim<2>::Vector& pos1, const Dim<2>::Vector& pos2) {
    return ((pos1 - mOrigin).cross(pos2 - mOrigin).z() > 0.0);
  }
  bool operator()(const unsigned id1, const unsigned id2) {
    return (*this)(mPositions.at(id1), mPositions.at(id2));
  }
  const Container& mPositions;
  const Dim<2>::Vector mOrigin;
};

}

#endif
