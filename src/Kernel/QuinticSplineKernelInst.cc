//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/QuinticSplineKernel.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
"""

if ndim == "1":
   text += """
template<>
QuinticSplineKernel< Dim<1> >::QuinticSplineKernel():
  Kernel<Dim<1>, QuinticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(FastMath::pow5(3.0)/40.0);
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
#endif

#if defined(SPHERAL_ENABLE_2D)

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
"""

if ndim == "1":
   text += """
template<>
QuinticSplineKernel< Dim<1> >::QuinticSplineKernel():
  Kernel<Dim<1>, QuinticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(FastMath::pow5(3.0)/40.0);
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
#endif

#if defined(SPHERAL_ENABLE_3D)

//------------------------------------------------------------------------------
// Empty constructor
//------------------------------------------------------------------------------
"""

if ndim == "1":
   text += """
template<>
QuinticSplineKernel< Dim<1> >::QuinticSplineKernel():
  Kernel<Dim<1>, QuinticSplineKernel< Dim<1> > >() {
  setVolumeNormalization(FastMath::pow5(3.0)/40.0);
  setKernelExtent(1.0);
  setInflectionPoint(0.342037); // (2.0/15.0*(7.0 - pow(2.0, 1.0/3.0) - pow(22.0, 2.0/3.0)));
#endif
}