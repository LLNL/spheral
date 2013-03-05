//---------------------------------Spheral++----------------------------------//
// CounterClockwiseCompareElements.
//
// A comparator functor to help sorting 2-D positions counterclockwise about
// a given origin.
//
// Created by JMO, Thu Nov 18 10:21:51 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral__CounterClockwiseCompareElements__
#define __Spheral__CounterClockwiseCompareElements__

#include "Utilities/DBC.hh"

//------------------------------------------------------------------------------
// Comparator to help sorting the faces counter-clockwise about a position.
//------------------------------------------------------------------------------
template<typename Element, typename Vector>
struct CounterClockwiseCompareElements {
  CounterClockwiseCompareElements(const vector<Element>& elements,
                                  const Vector& centroid):
    mElements(elements),
    mCentroid(centroid) {};
  int operator()(const unsigned id1, const unsigned id2) {
    REQUIRE(id1 < mElements.size());
    REQUIRE(id2 < mElements.size());
    return int(sgn0(((mElements[id1].position() - mCentroid).cross(mElements[id2].position() - mCentroid)).z()));
  }
  const vector<Element>& mElements;
  const Vector& mCentroid;
};

#endif
