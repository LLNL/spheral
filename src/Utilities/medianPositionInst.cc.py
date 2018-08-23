text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/medianPosition.cc"

namespace Spheral {
  template Dim< %(ndim)s >::Vector medianPosition(vector<Dim< %(ndim)s >::Vector>& positions);
}
"""
