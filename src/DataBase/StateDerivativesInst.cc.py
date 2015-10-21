text = """
//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "StateDerivatives.cc"

namespace Spheral {
  template class StateDerivatives<Dim< %(ndim)s > >;
}
"""
