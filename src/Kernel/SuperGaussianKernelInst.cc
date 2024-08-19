//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/SuperGaussianKernel.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template<> const double SuperGaussianKernel<Dim<1> >::mKW = 0.5*3.0;
  template<> const double SuperGaussianKernel<Dim<1> >::mKGW = 0.5*1.0;

  template class SuperGaussianKernel<Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
  template<> const double SuperGaussianKernel<Dim<2> >::mKW = 0.5*4.0;
  template<> const double SuperGaussianKernel<Dim<2> >::mKGW = 0.5*2.0;

  template class SuperGaussianKernel<Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
  template<> const double SuperGaussianKernel<Dim<3> >::mKW = 0.5*5.0;
  template<> const double SuperGaussianKernel<Dim<3> >::mKGW = 0.5*3.0;

  template class SuperGaussianKernel<Dim<3> >;
#endif
}
