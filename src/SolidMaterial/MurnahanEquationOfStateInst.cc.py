text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MurnahanEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MurnahanEquationOfState<Dim< %(ndim)s > >;
}
"""
