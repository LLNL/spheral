text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TillotsonEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class TillotsonEquationOfState<Dim< %(ndim)s > >;
}
"""
