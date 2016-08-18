//---------------------------------Spheral++----------------------------------//
// comparisons
// 
// A collection of functors useful for comparing various types in STL operations.
//
// Created by JMO, Sat Mar 13 13:25:08 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_comparisons__
#define __Spheral_comparisons__

namespace Spheral {

template<typename Pair>
struct ComparePairByFirstElement {
  bool operator()(const Pair& lhs, const Pair& rhs) {
    return lhs.first < rhs.first;
  }
};

}

#endif
