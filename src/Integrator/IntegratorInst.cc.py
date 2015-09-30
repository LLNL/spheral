text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class Integrator< Dim< %(ndim)s > >;
  }
}
"""
