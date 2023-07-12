//---------------------------------Spheral++----------------------------------//
// SphericalRadialKernel
//
// Take a 1D Kernel and build a specialized 1D tabulated version appropriate
// for use with the radial spherical SPH algorithm
//
// Created by JMO, Wed Jul  5 10:46:14 PDT 2023
//----------------------------------------------------------------------------//

#include "SphericalRadialKernel.hh"
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
// Functor for doing our volume integration of the kernel, which is a function
// or R/h
//------------------------------------------------------------------------------
template<typename KernelType>
struct RadialKernelNormalization {
  const KernelType& mW;
  double metamax;
  size_t mn;
  RadialKernelNormalization(const KernelType& W, const size_t n): mW(W), metamax(W.kernelExtent()), mn(n) {}

  // Define a local nested functor we'll use to do the volume integral for
  // an (etaj, etai) lookup.
  struct VolFunc {
    const KernelType& mW;
    VolFunc(const KernelType& W): mW(W) {}
    double operator()(const double q) const {
      return q*q*mW(std::abs(q), 1.0);
    }
  };

  // Integrate the volume function over the radial kernel extent for the given eta0 = R/h
  double operator()(const double eta0) const {
    const auto Ainv = simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                                   eta0 - metamax,
                                                                   eta0 + metamax,
                                                                   mn);
    ENSURE(Ainv > 0.0);
    return 1.0/Ainv;
  }
};

//------------------------------------------------------------------------------
// Functor for doing the volume integration of the kernel normalization gradient,
// which is needed for computing the full spherical gradient.
//------------------------------------------------------------------------------
template<typename KernelType>
struct RadialKernelNormalizationGradient {
  const KernelType& mW;
  double metamax;
  size_t mn;
  RadialKernelNormalizationGradient(const KernelType& W, const size_t n): mW(W), mn(n) {}

  // Define a local nested functor we'll use to do the volume integral for
  // an (etaj, etai) lookup.
  struct VolFunc {
    const KernelType& mW;
    VolFunc(const KernelType& W): mW(W) {}
    double operator()(const double q) const {
      return q*q*mW.grad(std::abs(q), 1.0);
    }
  };

  // Integrate the volume function over the radial kernel extent for the given eta0 = R/h
  double operator()(const double eta0) const {
    const auto gradAInv = simpsonsIntegration<VolFunc, double, double>(VolFunc(mW),
                                                                       eta0 - metamax,
                                                                       eta0 + metamax,
                                                                       mn);
    ENSURE(gradAInv > 0.0);
    return gradAInv;
  }
};

}            // anonymous

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename KernelType>
SphericalRadialKernel::SphericalRadialKernel(const KernelType& kernel,
                                             const unsigned numIntegral,
                                             const unsigned numKernel,
                                             const bool useInterpolation):


  mAInterp(0.0,
           10.0*kernel.kernelExtent(),
           numKernel,
           RadialKernelNormalization<KernelType>(kernel, numIntegral)),
  mGradAInvInterp(0.0,
                  10.0*kernel.kernelExtent(),
                  numKernel,
                  RadialKernelNormalizationGradient<KernelType>(kernel, numIntegral)),
  mBaseKernel(kernel),
  metamax(kernel.kernelExtent()),
  metacutoff(10.0*kernel.kernelExtent()),
  mNumIntegral(numIntegral),
  mUseInterpolation(useInterpolation) {
}

//------------------------------------------------------------------------------
// Copy constructor
//------------------------------------------------------------------------------
SphericalRadialKernel::SphericalRadialKernel(const SphericalRadialKernel& rhs):
  mAInterp(rhs.mAInterp),
  mGradAInvInterp(rhs.mGradAInvInterp),
  mBaseKernel(rhs.mBaseKernel),
  metamax(rhs.metamax),
  metacutoff(rhs.metacutoff),
  mNumIntegral(rhs.mNumIntegral),
  mUseInterpolation(rhs.mUseInterpolation) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalRadialKernel::~SphericalRadialKernel() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
SphericalRadialKernel&
SphericalRadialKernel::operator=(const SphericalRadialKernel& rhs) {
  if (this != &rhs) {
    mAInterp = rhs.mAInterp;
    mGradAInvInterp = rhs.mGradAInvInterp;
    mBaseKernel = rhs.mBaseKernel;
    metamax = rhs.metamax;
    metacutoff = metacutoff;
    mNumIntegral = rhs.mNumIntegral;
    mUseInterpolation = rhs.mUseInterpolation;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Lookup the kernel for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
double
SphericalRadialKernel::operator()(const Dim<1>::Vector& etaj,
                                  const Dim<1>::Vector& etai,
                                  const Dim<1>::Scalar  Hdet) const {
  REQUIRE(Hdet >= 0.0);
  const auto Ai = volumeNormalization(etai[0]);
  return Ai*FastMath::pow3(Hdet)*mBaseKernel(std::abs(etaj[0] - etai[0]), 1.0);
}

//------------------------------------------------------------------------------
// Lookup the grad kernel for (rj/h, ri/h) = (etaj, etai)
// Using the Leibniz integral rule to differentiate this integral.
//------------------------------------------------------------------------------
Dim<1>::Vector
SphericalRadialKernel::grad(const Dim<1>::Vector& etaj,
                            const Dim<1>::Vector& etai,
                            const Dim<1>::SymTensor& H) const {
  const auto Hdet4 = FastMath::pow4(H.Determinant());
  REQUIRE(Hdet4 >= 0.0);
  const auto ei = etai[0];
  const auto ej = etaj[0];
  const auto eji = ej - ei;
  const auto Ai = volumeNormalization(std::abs(ei));
  const auto gradAInvi = gradAInv(std::abs(ei));
  return Vector(sgn(eji)*Ai*Hdet4*(mBaseKernel.grad(eji, 1.0) - Ai*gradAInvi));
}

//------------------------------------------------------------------------------
// Simultaneously lookup (W,  grad W) for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
void
SphericalRadialKernel::kernelAndGrad(const Dim<1>::Vector& etaj,
                                     const Dim<1>::Vector& etai,
                                     const Dim<1>::SymTensor& H,
                                     Dim<1>::Scalar& W,
                                     Dim<1>::Vector& gradW,
                                     Dim<1>::Scalar& deltaWsum) const {
  const auto ei = etai[0];
  const auto ej = etaj[0];
  const auto eji = ej - ei;
  if (eji > metamax) {
    W = 0.0;
    gradW.Zero();
    deltaWsum = 0.0;
  } else {
    const auto Hdet = H.Determinant();
    CHECK(Hdet >= 0.0);
    const auto Ai = volumeNormalization(std::abs(ei));
    const auto gradAInvi = gradAInv(std::abs(ei));
    const auto thpt = Ai*FastMath::pow3(Hdet);
    const auto gW = mBaseKernel.grad(eji, 1.0);
    W = thpt * mBaseKernel(std::abs(eji), 1.0);
    gradW = Vector(sgn(eji)*thpt*Hdet*(gW - Ai*gradAInvi));
    deltaWsum = Hdet*gW;
  }
}

//------------------------------------------------------------------------------
// Return the volume normalization as a function of eta
//------------------------------------------------------------------------------
inline
double
SphericalRadialKernel::volumeNormalization(const double eta) const {
  return (mUseInterpolation ? (eta <= metacutoff ?
                               mAInterp(eta) :
                               mAInterp(eta) * FastMath::pow2(eta/metacutoff)) :
          RadialKernelNormalization<TableKernel<Dim<1>>>(mBaseKernel, mNumIntegral)(eta));
}

//------------------------------------------------------------------------------
// Return gradA(eta)
//------------------------------------------------------------------------------
inline
double
SphericalRadialKernel::gradAInv(const double eta) const {
  return (mUseInterpolation ? (eta <= metacutoff ?
                               mGradAInvInterp(eta) :
                               mGradAInvInterp(eta) * FastMath::pow2(eta/metacutoff)) :
          RadialKernelNormalizationGradient<TableKernel<Dim<1>>>(mBaseKernel, mNumIntegral)(eta));
}

// Constructor instantiations
template SphericalRadialKernel::SphericalRadialKernel(const TableKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const BSplineKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const NBSplineKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const W4SplineKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const GaussianKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const SuperGaussianKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const PiGaussianKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const HatKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const SincKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const NSincPolynomialKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const QuarticSplineKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const QuinticSplineKernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const WendlandC2Kernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const WendlandC4Kernel<Dim<1>>&, const unsigned, const unsigned, const bool);
template SphericalRadialKernel::SphericalRadialKernel(const WendlandC6Kernel<Dim<1>>&, const unsigned, const unsigned, const bool);

}
