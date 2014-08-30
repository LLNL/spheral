#include "TreeNeighbor.cc"

namespace Spheral {
  namespace NeighborSpace {

    //------------------------------------------------------------------------------
    // Define our static members.
    //------------------------------------------------------------------------------
    template<typename Dimension> const unsigned TreeNeighbor<Dimension>::num1dbits = 21U;
    template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::max1dKey = 1U << TreeNeighbor<Dimension>::num1dbits;
    template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::xkeymask = (1U << TreeNeighbor<Dimension>::num1dbits) - 1U;
    template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::ykeymask = TreeNeighbor<Dimension>::xkeymask << TreeNeighbor<Dimension>::num1dbits;
    template<typename Dimension> const uint64_t TreeNeighbor<Dimension>::zkeymask = TreeNeighbor<Dimension>::ykeymask << TreeNeighbor<Dimension>::num1dbits;

    //------------------------------------------------------------------------------
    // Explicit instantiation.
    //------------------------------------------------------------------------------
    template class TreeNeighbor<Dim<1> >;
    template class TreeNeighbor<Dim<2> >;
    template class TreeNeighbor<Dim<3> >;
  }
}

