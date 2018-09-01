text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorMonaghanGingoldViscosity.cc"

namespace Spheral {
  template class TensorMonaghanGingoldViscosity< Dim< %(ndim)s > >;
}
"""
