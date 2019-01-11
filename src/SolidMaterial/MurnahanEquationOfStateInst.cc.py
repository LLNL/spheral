text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/MurnahanEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MurnahanEquationOfState<Dim< %(ndim)s > >;
}
"""
