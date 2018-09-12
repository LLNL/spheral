text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/SynchronousRK4.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SynchronousRK4< Dim< %(ndim)s > >;
}
"""
