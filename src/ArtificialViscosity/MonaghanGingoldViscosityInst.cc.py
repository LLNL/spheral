text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "MonaghanGingoldViscosity.cc"

namespace Spheral {
  template class MonaghanGingoldViscosity< Dim< %(ndim)s > >;
}
"""
