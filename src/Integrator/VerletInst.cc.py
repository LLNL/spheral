text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/Verlet.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Verlet< Dim< %(ndim)s > >;
}
"""
