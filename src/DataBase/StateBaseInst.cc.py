text = """
//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase/StateBase.cc"

namespace Spheral {
  template class StateBase<Dim< %(ndim)s > >;
}
"""
