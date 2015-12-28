text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NBSplineKernel.cc"

namespace Spheral {
  namespace KernelSpace {
    template class NBSplineKernel< Dim< %(ndim)s >  >;
  }
}


"""
