text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Integrator< Dim< %(ndim)s > >;
}
"""
