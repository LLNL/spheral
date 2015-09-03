text = """
//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "StateBase.cc"

namespace Spheral {
  template class StateBase<Dim< %(ndim)s > >;
}
"""
