text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Kernel/NBSplineKernel.cc"

namespace Spheral {
  template class NBSplineKernel< Dim< %(ndim)s >  >;
}
"""
