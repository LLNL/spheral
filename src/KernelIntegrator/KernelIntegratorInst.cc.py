text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "KernelIntegrator.cc"

namespace Spheral {
  template class KernelIntegrator<Dim<%(ndim)s>>;
}
"""
