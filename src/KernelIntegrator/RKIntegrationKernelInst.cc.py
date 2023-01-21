text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "KernelIntegrator/RKIntegrationKernel.cc"
#include "Geometry/Dimension.hh"
"""

for order in [0, 1, 2, 3, 4, 5, 6, 7]:
    text += """
namespace Spheral {
template class RKIntegrationKernel<Dim<%(ndim)s>, """
    text += """%(order)s>;
}
""" % {"order" : order}
