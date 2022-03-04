//---------------------------------Spheral++----------------------------------//
// SphericalKernel
//
// Take a 3D Kernel and build a specialized 1D tabulated version appropriate
// for use with the spherical SPH algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Created by JMO, Wed Dec  2 16:41:20 PST 2020
//----------------------------------------------------------------------------//

#include "SphericalKernel.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/DBC.hh"

#include "Kernel/BSplineKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/HatKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/ExpInvKernel.hh"

namespace Spheral {

namespace {  // anonymous

//------------------------------------------------------------------------------
// Functor for doing our volume integration of the kernel  
//------------------------------------------------------------------------------
template<typename KernelType>
struct W3S1Func {
  const KernelType& mW;
  double metaMax;
  size_t mn;
  W3S1Func(const KernelType& W, const size_t n): mW(W), metaMax(mW.kernelExtent()), mn(n) {}

  // Define a local nested functor we'll use to do the volume integral for
  // an (etaj, etai) lookup.
  struct VolFunc {
    const KernelType& mW;
    VolFunc(const KernelType& W): mW(W) {}
    double operator()(const double q) const {
      CHECK(q >= 0.0);
      return q * mW(q, 1.0);
    }
  };

  // Now the lookup based on (etaj, etai) -- the (rprime/h, r/h) from the paper
  double operator()(const double low, const double high) const {
    if (low >= high) return 0.0;
    return simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                        low,
                                                        high,
                                                        mn);
  }
};

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename KernelType>
SphericalKernel::SphericalKernel(const KernelType& kernel,
                                 const unsigned numIntegral,
                                 const unsigned numKernel):
  mInterp(0.0, kernel.kernelExtent(),
          0.0, kernel.kernelExtent(),
          numKernel, numKernel, 
          W3S1Func<KernelType>(kernel, numIntegral)),
  mBaseKernel3d(kernel, numKernel),
  mBaseKernel1d(kernel, numKernel),
  metamax(kernel.kernelExtent()) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SphericalKernel::SphericalKernel(const SphericalKernel& rhs):
  mInterp(rhs.mInterp),
  mBaseKernel3d(rhs.mBaseKernel3d),
  mBaseKernel1d(rhs.mBaseKernel1d),
  metamax(rhs.metamax) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalKernel::~SphericalKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
SphericalKernel&
SphericalKernel::operator=(const SphericalKernel& rhs) {
  if (this != &rhs) {
    mInterp = rhs.mInterp;
    mBaseKernel3d = rhs.mBaseKernel3d;
    mBaseKernel1d = rhs.mBaseKernel1d;
    metamax = rhs.metamax;
  }
  return *this;
}

// Constructor instantiations
template SphericalKernel::SphericalKernel<TableKernel<Dim<3>>>          (const TableKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<BSplineKernel<Dim<3>>>        (const BSplineKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<NBSplineKernel<Dim<3>>>       (const NBSplineKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<W4SplineKernel<Dim<3>>>       (const W4SplineKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<GaussianKernel<Dim<3>>>       (const GaussianKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<SuperGaussianKernel<Dim<3>>>  (const SuperGaussianKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<PiGaussianKernel<Dim<3>>>     (const PiGaussianKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<HatKernel<Dim<3>>>            (const HatKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<SincKernel<Dim<3>>>           (const SincKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<NSincPolynomialKernel<Dim<3>>>(const NSincPolynomialKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<QuarticSplineKernel<Dim<3>>>  (const QuarticSplineKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<QuinticSplineKernel<Dim<3>>>  (const QuinticSplineKernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<WendlandC2Kernel<Dim<3>>>     (const WendlandC2Kernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<WendlandC4Kernel<Dim<3>>>     (const WendlandC4Kernel<Dim<3>>&, const unsigned, const unsigned);
template SphericalKernel::SphericalKernel<WendlandC6Kernel<Dim<3>>>     (const WendlandC6Kernel<Dim<3>>&, const unsigned, const unsigned);

}

// We need to instantiate the special TableKernel constructor we use
#include "TableKernel.cc"
namespace Spheral {
template TableKernel<Dim<1>>::TableKernel(const TableKernel<Dim<3>>&, const unsigned);
}

