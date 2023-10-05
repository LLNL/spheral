text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/ShearModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ShearModulusPolicy<Dim< %(ndim)s > >;
}
"""
