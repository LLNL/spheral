text = """
//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "DataBase/State.cc"

namespace Spheral {
  template class State<Dim< %(ndim)s > >;
}
"""
