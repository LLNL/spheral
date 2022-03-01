#include "Utilities/SpheralFunctions.hh"  // sgn
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
inline
bool
SphericalKernel::operator==(const SphericalKernel& rhs) const {
  return ((mInterp == rhs.mInterp) and
          (metamax == rhs.metamax));
}

//------------------------------------------------------------------------------
// Lookup the kernel for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
inline
double
SphericalKernel::operator()(const Dim<1>::Vector& etaj,
                                 const Dim<1>::Vector& etai,
                                 const Dim<1>::Scalar  Hdet) const {
  REQUIRE(Hdet >= 0.0);
  const auto ei = std::max(1e-8, std::abs(etai[0]));
  const auto ej = std::max(1e-8, std::abs(etaj[0]));
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return 0.0;
  const auto max_bound = std::min(metamax, ei + ej);
  return 2.0*M_PI/(ei*ej)*FastMath::cube(Hdet)*mInterp(min_bound, max_bound);
}

//------------------------------------------------------------------------------
// Lookup the grad kernel for (rj/h, ri/h) = (etaj, etai)
// Using the Leibniz integral rule to differentiate this integral.
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
SphericalKernel::grad(const Dim<1>::Vector& etaj,
                           const Dim<1>::Vector& etai,
                           const Dim<1>::SymTensor& H) const {
  const auto Hdet = H.Determinant();
  REQUIRE(Hdet >= 0.0);
  const auto ei = std::max(1e-8, std::abs(etai[0]));
  const auto ej = std::max(1e-8, std::abs(etaj[0]));
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return Vector::zero;
  const auto max_bound = std::min(metamax, ei + ej);
  const auto B = (ei + ej >= metamax ?
                  0.0 :
                  max_bound*mBaseKernel3d.kernelValue(max_bound, 1.0));
  const auto A = min_bound*mBaseKernel3d.kernelValue(min_bound, 1.0)*(ei > ej ? 1.0 : -1.0);
  return Vector(2.0*M_PI/(ei*ej)*FastMath::pow4(Hdet)*(B - A - 1.0/ei*mInterp(min_bound, max_bound)));
}

//------------------------------------------------------------------------------
// Simultaneously lookup (W,  grad W) for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
inline
void
SphericalKernel::kernelAndGrad(const Dim<1>::Vector& etaj,
                                    const Dim<1>::Vector& etai,
                                    const Dim<1>::SymTensor& H,
                                    Dim<1>::Scalar& W,
                                    Dim<1>::Vector& gradW,
                                    Dim<1>::Scalar& deltaWsum) const {
  const auto Hdet = H.Determinant();
  REQUIRE(Hdet >= 0.0);
  const auto ei = std::max(1e-8, std::abs(etai[0]));
  const auto ej = std::max(1e-8, std::abs(etaj[0]));
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) {
    W = 0.0;
    gradW.Zero();
    deltaWsum = 0.0;
  } else {
    const auto max_bound = std::min(metamax, ei + ej);
    const auto B = (ei + ej >= metamax ?
                    0.0 :
                    max_bound*mBaseKernel3d.kernelValue(max_bound, 1.0));
    const auto A = min_bound*mBaseKernel3d.kernelValue(min_bound, 1.0)*(ei > ej ? 1.0 : -1.0);
    const auto pre = 2.0*M_PI/(ei*ej)*FastMath::cube(Hdet);
    const auto interpVal = mInterp(min_bound, max_bound);
    W = pre*interpVal;
    gradW = Vector(pre*Hdet*(B - A - interpVal/ei));
    deltaWsum = mBaseKernel1d.gradValue((etaj - etai).magnitude(), Hdet);
  }
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
const typename SphericalKernel::InterpolatorType&
SphericalKernel::Winterpolator() const {
  return mInterp;
}

inline
const TableKernel<Dim<3>>&
SphericalKernel::baseKernel3d() const {
  return mBaseKernel3d;
}

inline
const TableKernel<Dim<1>>&
SphericalKernel::baseKernel1d() const {
  return mBaseKernel1d;
}

inline
double
SphericalKernel::etamax() const {
  return metamax;
}

}
