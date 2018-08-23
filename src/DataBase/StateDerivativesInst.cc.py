text = """
//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/StateDerivatives.cc"

namespace Spheral {
  template class StateDerivatives<Dim< %(ndim)s > >;
}
"""
