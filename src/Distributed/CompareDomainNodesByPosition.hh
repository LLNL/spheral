//---------------------------------Spheral++----------------------------------//
// CompareDomainNodesByPosition.
//
// Comparison functor for use in sorting DomainNodes by position, allowing the
// user to cycle through the significant dimension via a construction 
// parameter.
//
// Created by JMO, Sat Dec  4 23:04:01 PST 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_PartitionSpace_CompareDomainNodesByPosition_hh__
#define __Spheral_PartitionSpace_CompareDomainNodesByPosition_hh__

#include "Utilities/DomainNode.hh"

namespace Spheral {

template<typename Dimension>
class CompareDomainNodesByPosition {
 public:
  CompareDomainNodesByPosition(const int positionIndex):
    mPositionIndex(positionIndex) {}

  bool operator()(const typename Dimension::Vector& lhs,
                  const typename Dimension::Vector& rhs) const {
    for (int i = 0; i != Dimension::nDim; ++i) {
      const int j = (mPositionIndex + i) % Dimension::nDim;
      CHECK(j >= 0 && j < Dimension::nDim);
      if (lhs(j) < rhs(j)) {
        return true;
      } else if (lhs(j) > rhs(j)) {
        return false;
      }
    }
    return false;
  }

  bool operator()(const DomainNode<Dimension>& lhs,
                  const DomainNode<Dimension>& rhs) const {
    return operator()(lhs.position, rhs.position);
  }

 private:
  int mPositionIndex;
};

}

#endif
