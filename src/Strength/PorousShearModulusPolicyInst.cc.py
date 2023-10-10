text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PorousShearModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousShearModulusPolicy<Dim< %(ndim)s > >;
}
"""
