text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "KernelIntegrator/SPHIntegrationKernel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template class SPHIntegrationKernel<Dim<%(ndim)s>>;
}
"""
