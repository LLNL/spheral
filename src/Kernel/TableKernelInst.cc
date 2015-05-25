//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TableKernel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace KernelSpace {
    template class TableKernel< Dim<1> >;
    template class TableKernel< Dim<2> >;
    template class TableKernel< Dim<3> >;

    template TableKernel<Dim<1> >::TableKernel(const BSplineKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const W4SplineKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const GaussianKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const SuperGaussianKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const PiGaussianKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const HatKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const SincKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const NSincPolynomialKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const QuarticSplineKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const QuinticSplineKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const NBSplineKernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const WendlandC4Kernel<Dim<1> >&, const int, const double);
    template TableKernel<Dim<1> >::TableKernel(const WendlandC6Kernel<Dim<1> >&, const int, const double);

    template TableKernel<Dim<2> >::TableKernel(const BSplineKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const W4SplineKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const GaussianKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const SuperGaussianKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const PiGaussianKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const HatKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const SincKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const NSincPolynomialKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const QuarticSplineKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const QuinticSplineKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const NBSplineKernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const WendlandC4Kernel<Dim<2> >&, const int, const double);
    template TableKernel<Dim<2> >::TableKernel(const WendlandC6Kernel<Dim<2> >&, const int, const double);

    template TableKernel<Dim<3> >::TableKernel(const BSplineKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const W4SplineKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const GaussianKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const SuperGaussianKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const PiGaussianKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const HatKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const SincKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const NSincPolynomialKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const QuarticSplineKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const QuinticSplineKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const NBSplineKernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const WendlandC4Kernel<Dim<3> >&, const int, const double);
    template TableKernel<Dim<3> >::TableKernel(const WendlandC6Kernel<Dim<3> >&, const int, const double);
  }
}
