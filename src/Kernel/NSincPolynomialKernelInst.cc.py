text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NSincPolynomialKernel.cc"

namespace Spheral {
  template class NSincPolynomialKernel< Dim< %(ndim)s > >;
}
"""
