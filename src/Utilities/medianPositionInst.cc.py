text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "medianPosition.cc"

namespace Spheral {
  template Dim< %(ndim)s >::Vector medianPosition(vector<Dim< %(ndim)s >::Vector>& positions);
}
"""
