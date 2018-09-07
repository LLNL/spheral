text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GammaLawGas.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GammaLawGas< Dim< %(ndim)s >  >;
}
"""
