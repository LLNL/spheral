text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MurnahanEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class MurnahanEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
