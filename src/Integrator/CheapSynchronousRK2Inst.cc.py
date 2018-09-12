text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/CheapSynchronousRK2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CheapSynchronousRK2< Dim< %(ndim)s > >;
}
"""
