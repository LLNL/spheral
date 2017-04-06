text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ShearModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ShearModulusPolicy<Dim< %(ndim)s > >;
}

"""
