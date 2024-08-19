//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Kernel/TableKernel.cc"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template class TableKernel< Dim<1> >;

  template TableKernel<Dim<1>>::TableKernel(const WendlandC2Kernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const WendlandC4Kernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const WendlandC6Kernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const GaussianKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const PiGaussianKernel<Dim<1>>&, const unsigned);
  template TableKernel<Dim<1>>::TableKernel(const NBSplineKernel<Dim<1>>&, const unsigned);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template class TableKernel< Dim<2> >;

  template TableKernel<Dim<2>>::TableKernel(const WendlandC2Kernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const WendlandC4Kernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const WendlandC6Kernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const GaussianKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const PiGaussianKernel<Dim<2>>&, const unsigned);
  template TableKernel<Dim<2>>::TableKernel(const NBSplineKernel<Dim<2>>&, const unsigned);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template class TableKernel< Dim<3> >;

  template TableKernel<Dim<3>>::TableKernel(const WendlandC2Kernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const WendlandC4Kernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const WendlandC6Kernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const GaussianKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const PiGaussianKernel<Dim<3>>&, const unsigned);
  template TableKernel<Dim<3>>::TableKernel(const NBSplineKernel<Dim<3>>&, const unsigned);
#endif
}
