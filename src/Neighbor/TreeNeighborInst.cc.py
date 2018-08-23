text = """
#include "Neighbor/TreeNeighbor.cc"

namespace Spheral {
  namespace NeighborSpace {

    //------------------------------------------------------------------------------
    // Explicit instantiation.
    //------------------------------------------------------------------------------
    template class TreeNeighbor< Dim< %(ndim)s > >;
  }
}

"""
