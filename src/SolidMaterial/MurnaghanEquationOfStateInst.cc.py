text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/MurnaghanEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MurnaghanEquationOfState<Dim< %(ndim)s > >;
}
"""
