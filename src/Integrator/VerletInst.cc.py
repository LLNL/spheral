text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Verlet.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Verlet< Dim< %(ndim)s > >;
}
"""
