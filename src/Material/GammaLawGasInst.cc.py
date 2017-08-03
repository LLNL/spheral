text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GammaLawGas.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class GammaLawGas< Dim< %(ndim)s >  >;
  }
}
"""
