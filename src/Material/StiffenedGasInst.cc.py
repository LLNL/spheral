text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Material/StiffenedGas.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class StiffenedGas< Dim< %(ndim)s >  >;
}
"""
