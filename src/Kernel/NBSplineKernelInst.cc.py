text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NBSplineKernel.cc"

namespace Spheral {
  namespace KernelSpace {
    template class NBSplineKernel<  %(ndim)s  >;
  }
}


"""
