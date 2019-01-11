text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/Integrator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Integrator< Dim< %(ndim)s > >;
}
"""
