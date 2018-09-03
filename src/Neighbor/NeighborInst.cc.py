text = """
#include "Neighbor.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class Neighbor< Dim< %(ndim)s > >;
}
"""
