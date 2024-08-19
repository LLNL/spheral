//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "KernelIntegrator/RKIntegrationKernel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template class RKIntegrationKernel<Dim<1>, 0>;
template class RKIntegrationKernel<Dim<1>, 1>;
template class RKIntegrationKernel<Dim<1>, 2>;
template class RKIntegrationKernel<Dim<1>, 3>;
template class RKIntegrationKernel<Dim<1>, 4>;
template class RKIntegrationKernel<Dim<1>, 5>;
template class RKIntegrationKernel<Dim<1>, 6>;
template class RKIntegrationKernel<Dim<1>, 7>;
#endif

#if defined(SPHERAL_ENABLE_2D)
template class RKIntegrationKernel<Dim<2>, 0>;
template class RKIntegrationKernel<Dim<2>, 1>;
template class RKIntegrationKernel<Dim<2>, 2>;
template class RKIntegrationKernel<Dim<2>, 3>;
template class RKIntegrationKernel<Dim<2>, 4>;
template class RKIntegrationKernel<Dim<2>, 5>;
template class RKIntegrationKernel<Dim<2>, 6>;
template class RKIntegrationKernel<Dim<2>, 7>;
#endif

#if defined(SPHERAL_ENABLE_3D)
template class RKIntegrationKernel<Dim<3>, 0>;
template class RKIntegrationKernel<Dim<3>, 1>;
template class RKIntegrationKernel<Dim<3>, 2>;
template class RKIntegrationKernel<Dim<3>, 3>;
template class RKIntegrationKernel<Dim<3>, 4>;
template class RKIntegrationKernel<Dim<3>, 5>;
template class RKIntegrationKernel<Dim<3>, 6>;
template class RKIntegrationKernel<Dim<3>, 7>;
#endif
}
