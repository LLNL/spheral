//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NBSplineKernel.cc"

namespace Spheral {
  namespace KernelSpace {
    template class NBSplineKernel< Dim<1> >;
    template class NBSplineKernel< Dim<2> >;
    template class NBSplineKernel< Dim<3> >;
  }
}


