text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SynchronousRK4.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class SynchronousRK4< Dim< %(ndim)s > >;
  }
}
"""
