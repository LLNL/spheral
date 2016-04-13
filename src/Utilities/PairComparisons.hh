//---------------------------------Spheral++----------------------------------//
// PairComparisons
//
// A set of functors that are handy for comparing std::pair types in STL 
// algorithms.
// 
// Created by J. Michael Owen, Thu Aug 12 16:20:18 PDT 2010
//----------------------------------------------------------------------------//
#ifndef _Spheral_PairComparisons_
#define _Spheral_PairComparisons_

namespace Spheral {

//------------------------------------------------------------------------------
// Compare pairs by their elements, using the first value as most significant.
//------------------------------------------------------------------------------
template<typename Pair>
struct ComparePairs {
  bool operator()(const Pair& lhs, const Pair& rhs) const {
    return (lhs.first < rhs.first   ? true  :
	    lhs.first > rhs.first   ? false :
	    lhs.second < rhs.second ? true  :
	                              false);
  }
};

//------------------------------------------------------------------------------
// Compare pairs by the first element in increasing order.
//------------------------------------------------------------------------------
template<typename Pair>
struct ComparePairsByFirstElement {
  bool operator()(const Pair& lhs, const Pair& rhs) const {
    return lhs.first < rhs.first;
  }
  bool operator()(const Pair& lhs, const typename Pair::first_type& rhs) const {
    return lhs.first < rhs;
  }
};

//------------------------------------------------------------------------------
// Compare pairs by the first element in decreasing order.
//------------------------------------------------------------------------------
template<typename Pair>
struct ComparePairsByFirstElementInDecreasingOrder {
  bool operator()(const Pair& lhs, const Pair& rhs) const {
    return lhs.first > rhs.first;
  }
  bool operator()(const Pair& lhs, const typename Pair::first_type& rhs) const {
    return lhs.first > rhs;
  }
};

//------------------------------------------------------------------------------
// Compare pairs by the second element in increasing order.
//------------------------------------------------------------------------------
template<typename Pair>
struct ComparePairsBySecondElement {
  bool operator()(const Pair& lhs, const Pair& rhs) const {
    return lhs.second < rhs.second;
  }
  bool operator()(const typename Pair::second_type& lhs, const Pair& rhs) const {
    return lhs < rhs.second;
  }
  bool operator()(const Pair& lhs, const typename Pair::second_type& rhs) const {
    return lhs.second < rhs;
  }
};

//------------------------------------------------------------------------------
// Compare pairs by the second element in decreasing order.
//------------------------------------------------------------------------------
template<typename Pair>
struct ComparePairsBySecondElementInDecreasingOrder {
  bool operator()(const Pair& lhs, const Pair& rhs) const {
    return lhs.second > rhs.second;
  }
  bool operator()(const typename Pair::second_type& lhs, const Pair& rhs) const {
    return lhs > rhs.second;
  }
  bool operator()(const Pair& lhs, const typename Pair::second_type& rhs) const {
    return lhs.second > rhs;
  }
};

//------------------------------------------------------------------------------
// Return the result of summing the first element of two pairs.
//------------------------------------------------------------------------------
template<typename Pair>
struct AddPairFirstElements {
  typename Pair::first_type operator()(const typename Pair::first_type& lhs, const Pair& rhs) const {
    return lhs + rhs.first;
  }
};

//------------------------------------------------------------------------------
// Return the result of summing the second element of two pairs.
//------------------------------------------------------------------------------
template<typename Pair>
struct AddPairSecondElements {
  typename Pair::second_type operator()(const typename Pair::second_type& lhs, const Pair& rhs) const {
    return lhs + rhs.second;
  }
};

//------------------------------------------------------------------------------
// Compare pairs for equality.
//------------------------------------------------------------------------------
template<typename Pair>
struct PairsEqualByFirstElement {
  bool operator()(const Pair& lhs, const Pair& rhs) const {
    return (lhs.first == rhs.first);
  }
};

template<typename Pair>
struct PairsEqualBySecondElement {
  bool operator()(const Pair& lhs, const Pair& rhs) const {
    return (lhs.second == rhs.second);
  }
};

template<typename Pair>
struct PairFirstElementEqualTo {
  typename Pair::first_type mVal;
  PairFirstElementEqualTo(const typename Pair::first_type x): mVal(x) {}
  bool operator()(const Pair& lhs) const {
    return (lhs.first == mVal);
  }
  bool operator()(const Pair& lhs, const typename Pair::first_type& rhs) const {
    return lhs.first == rhs;
  }
};

template<typename Pair>
struct PairSecondElementEqualTo {
  typename Pair::second_type mVal;
  PairSecondElementEqualTo(const typename Pair::second_type x): mVal(x) {}
  bool operator()(const Pair& lhs) const {
    return (lhs.second == mVal);
  }
  bool operator()(const Pair& lhs, const typename Pair::second_type& rhs) const {
    return lhs.second == rhs;
  }
};

}

#endif
