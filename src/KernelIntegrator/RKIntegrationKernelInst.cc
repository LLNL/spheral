//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "KernelIntegrator/RKIntegrationKernel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class RKIntegrationKernel<Dim<1>, """
    text += """%(order)s>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RKIntegrationKernel<Dim<2>, """
    text += """%(order)s>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RKIntegrationKernel<Dim<3>, """
    text += """%(order)s>;
#endif
}