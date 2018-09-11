text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/NBSplineKernel.cc"

namespace Spheral {
  namespace KernelSpace {
    template class NBSplineKernel< Dim< %(ndim)s >  >;
  }
}


"""
