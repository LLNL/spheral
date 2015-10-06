text = """
#include "Neighbor.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace NeighborSpace {
    template class Neighbor< Dim< %(ndim)s > >;
  }
}
"""
