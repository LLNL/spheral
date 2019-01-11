text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/SynchronousRK2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SynchronousRK2< Dim< %(ndim)s > >;
}
"""
