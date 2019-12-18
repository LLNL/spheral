text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Physics.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class Physics<Dim< %(ndim)s > >;
}
"""
