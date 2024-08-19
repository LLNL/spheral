//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/QuinticSplineKernel.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
template<>
QuinticSplineKernel< Dim<1> >::QuinticSplineKernel():
  Kernel<Dim<1>, QuinticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(FastMath::pow5(3.0)/40.0);
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}

template class QuinticSplineKernel< Dim<1> >;
#endif

#if defined(SPHERAL_ENABLE_2D)
template<>
QuinticSplineKernel< Dim<2> >::QuinticSplineKernel():
  Kernel<Dim<2>, QuinticSplineKernel< Dim<2> > >() {
  setVolumeNormalization(FastMath::pow7(3.0)*7.0/(478.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}

template class QuinticSplineKernel< Dim<2> >;
#endif

#if defined(SPHERAL_ENABLE_3D)
template<>
QuinticSplineKernel< Dim<3> >::QuinticSplineKernel():
  Kernel<Dim<3>, QuinticSplineKernel< Dim<3> > >() {
  setVolumeNormalization(FastMath::pow7(3.0)/(40.0*M_PI));
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
}

template class QuinticSplineKernel< Dim<3> >;
#endif
}
