text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MeshIdealHPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MeshIdealHPolicy<Dim< %(ndim)s > >;
}
"""
