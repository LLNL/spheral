text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RK/ReproducingKernel.cc"

namespace Spheral {
template class ReproducingKernel<Dim<%(ndim)s>>;
}
"""
