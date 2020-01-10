text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RK/ReproducingKernelMethods.cc"

namespace Spheral {
template class ReproducingKernelMethods<Dim<%(ndim)s>>;
}
"""
