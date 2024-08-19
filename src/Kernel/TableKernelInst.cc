//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/TableKernel.cc"
#include "Kernel/BSplineKernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Kernel/HatKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/ExpInvKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class TableKernel< Dim<1> >;

  template TableKernel<Dim<1>>::TableKernel(const BSplineKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const W4SplineKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const GaussianKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const SuperGaussianKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const PiGaussianKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const SincKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const NSincPolynomialKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const NBSplineKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const HatKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const QuarticKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const QuinticSplineKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const WendlandC2Kernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const WendlandC4Kernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const WendlandC6Kernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const ExpInvKernel<Dim<1>>&, const unsigned);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class TableKernel< Dim<2> >;

  template TableKernel<Dim<2>>::TableKernel(const BSplineKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const W4SplineKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const GaussianKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const SuperGaussianKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const PiGaussianKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const SincKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const NSincPolynomialKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const NBSplineKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const HatKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const QuarticKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const QuinticSplineKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const WendlandC2Kernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const WendlandC4Kernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const WendlandC6Kernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const ExpInvKernel<Dim<2>>&, const unsigned);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class TableKernel< Dim<3> >;

  template TableKernel<Dim<3>>::TableKernel(const BSplineKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const W4SplineKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const GaussianKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const SuperGaussianKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const PiGaussianKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const SincKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const NSincPolynomialKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const NBSplineKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const HatKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const QuarticKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const QuinticSplineKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const WendlandC2Kernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const WendlandC4Kernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const WendlandC6Kernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const ExpInvKernel<Dim<3>>&, const unsigned);
#endif
}
