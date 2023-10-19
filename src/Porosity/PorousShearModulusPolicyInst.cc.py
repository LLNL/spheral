text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorousShearModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousShearModulusPolicy<Dim< %(ndim)s > >;
}
"""
