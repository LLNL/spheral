text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/SynchronousRK1.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SynchronousRK1< Dim< %(ndim)s > >;
}
"""
