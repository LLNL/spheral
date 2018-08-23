text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Material/GammaLawGas.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class GammaLawGas< Dim< %(ndim)s >  >;
  }
}
"""
