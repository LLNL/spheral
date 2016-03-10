text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CheapSynchronousRK2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class CheapSynchronousRK2< Dim< %(ndim)s > >;
  }
}
"""
