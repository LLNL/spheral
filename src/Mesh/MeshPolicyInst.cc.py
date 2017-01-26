text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MeshPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MeshPolicy<Dim< %(ndim)s > >;
}

"""
