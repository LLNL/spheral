text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/NSincPolynomialKernel.cc"

namespace Spheral {
  template class NSincPolynomialKernel< Dim< %(ndim)s > >;
}
"""
