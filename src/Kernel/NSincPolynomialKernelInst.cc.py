text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/NSincPolynomialKernel.cc"

namespace Spheral {
  namespace KernelSpace {
    template class NSincPolynomialKernel< Dim< %(ndim)s > >;
  }
}
"""
