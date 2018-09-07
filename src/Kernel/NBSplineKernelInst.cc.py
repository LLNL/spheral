text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NBSplineKernel.cc"

namespace Spheral {
  template class NBSplineKernel< Dim< %(ndim)s >  >;
}
"""
