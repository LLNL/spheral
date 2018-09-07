text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geodyn.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Geodyn<Dim< %(ndim)s > >;
}
"""
