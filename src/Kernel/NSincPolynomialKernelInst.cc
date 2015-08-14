//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NSincPolynomialKernel.cc"

namespace Spheral {
  namespace KernelSpace {
    template class NSincPolynomialKernel<Dim<1> >;
    template class NSincPolynomialKernel<Dim<2> >;
    template class NSincPolynomialKernel<Dim<3> >;
  }
}
