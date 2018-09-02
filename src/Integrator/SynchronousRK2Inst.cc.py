text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SynchronousRK2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SynchronousRK2< Dim< %(ndim)s > >;
}
"""
