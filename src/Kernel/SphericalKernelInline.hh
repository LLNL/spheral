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
  const auto ei = std::max(1e-10, std::abs(etai[0]));
  const auto ej = std::max(1e-10, std::abs(etaj[0]));
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto a = std::abs(ej - ei);            // Lower integration limit
  if (a > metamax) return 0.0;
  const auto b = std::min(metamax, ei + ej);   // Upper integration limit
  return 2.0*M_PI/(ei*ej)*FastMath::cube(Hdet)*mInterp(a, b);
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
  const auto ei = std::max(1e-10, std::abs(etai[0]));
  const auto ej = std::max(1e-10, std::abs(etaj[0]));
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto a = std::abs(ej - ei);            // Lower integration limit
  if (a > metamax) return Vector::zero;
  const auto b = std::min(metamax, ei + ej);   // Upper integration limit
  const auto A = a*mBaseKernel3d.kernelValue(a, 1.0)*sgn0(ei - ej);
  const auto B = (ei + ej >= metamax ?
                  0.0 :
                  b*mBaseKernel3d.kernelValue(b, 1.0));
  return Vector(2.0*M_PI/(ei*ej)*FastMath::pow4(Hdet)*(B - A - 1.0/ei*mInterp(a, b)));
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
  const auto ei = std::max(1e-10, std::abs(etai[0]));
  const auto ej = std::max(1e-10, std::abs(etaj[0]));
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto a = std::abs(ej - ei);            // Lower integration limit
  if (a > metamax) {
    W = 0.0;
    gradW.Zero();
    deltaWsum = 0.0;
  } else {
    const auto b = std::min(metamax, ei + ej); // Upper integration limit
    const auto B = (ei + ej >= metamax ?
                    0.0 :
                    b*mBaseKernel3d.kernelValue(b, 1.0));
    const auto A = a*mBaseKernel3d.kernelValue(a, 1.0)*sgn0(ei - ej);
    const auto pre = 2.0*M_PI/(ei*ej)*FastMath::cube(Hdet);
    const auto interpVal = mInterp(a, b);
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
