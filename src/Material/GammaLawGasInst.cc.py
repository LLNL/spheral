text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Material/GammaLawGas.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GammaLawGas< Dim< %(ndim)s >  >;
}
"""
