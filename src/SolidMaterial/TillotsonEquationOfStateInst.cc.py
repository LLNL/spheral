text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TillotsonEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class TillotsonEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
