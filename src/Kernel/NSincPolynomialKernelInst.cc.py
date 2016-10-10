text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NSincPolynomialKernel.cc"

namespace Spheral {
  namespace KernelSpace {
    template class NSincPolynomialKernel< Dim< %(ndim)s > >;
  }
}
"""
