#include "Utilities/SpheralFunctions.hh"  // sgn
#include "Utilities/FastMath.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Lookup the kernel for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
inline
double
SphericalTableKernel::operator()(const Dim<1>::Vector& etaj,
                                 const Dim<1>::Vector& etai,
                                 const Dim<1>::Scalar  Hdeti) const {
  REQUIRE(Hdeti >= 0.0);
  const auto ei = std::max(1e-10, etai[0]);
  const auto ej = std::max(1e-10, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return 0.0;
  const auto max_bound = std::min(metamax, ei + ej);
  return 2.0*M_PI/(ei*ej)*FastMath::cube(Hdeti)*mInterp(min_bound, max_bound);
}

//------------------------------------------------------------------------------
// Lookup the grad kernel for (rj/h, ri/h) = (etaj, etai)
// Using the Leibniz integral rule to differentiate this integral.
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
SphericalTableKernel::grad(const Dim<1>::Vector& etaj,
                           const Dim<1>::Vector& etai,
                           const Dim<1>::Scalar  Hdeti) const {
  REQUIRE(Hdeti >= 0.0);
  const auto ei = std::max(1e-10, etai[0]);
  const auto ej = std::max(1e-10, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return Vector::zero;
  const auto max_bound = std::min(metamax, ei + ej);
  const auto B = (ei + ej >= metamax ?
                  0.0 :
                  max_bound*mKernel.kernelValue(max_bound, 1.0));
  const auto A = min_bound*mKernel.kernelValue(min_bound, 1.0)*(ei > ej ? 1.0 : -1.0);
  return Vector(2.0*M_PI/(ei*ej)*FastMath::pow4(Hdeti)*(B - A - 1.0/ei*mInterp(min_bound, max_bound)));
}

//------------------------------------------------------------------------------
// Simultaneously lookup (W,  grad W) for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
inline
std::pair<double, Dim<1>::Vector>
SphericalTableKernel::kernelAndGradValue(const Dim<1>::Vector& etaj,
                                         const Dim<1>::Vector& etai,
                                         const Dim<1>::Scalar  Hdeti) const {
  REQUIRE(Hdeti >= 0.0);
  const auto ei = std::max(1e-10, etai[0]);
  const auto ej = std::max(1e-10, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return std::make_pair(0.0, Vector::zero);
  const auto max_bound = std::min(metamax, ei + ej);
  const auto B = (ei + ej >= metamax ?
                  0.0 :
                  max_bound*mKernel.kernelValue(max_bound, 1.0));
  const auto A = min_bound*mKernel.kernelValue(min_bound, 1.0)*(ei > ej ? 1.0 : -1.0);
  const auto pre = 2.0*M_PI/(ei*ej)*FastMath::cube(Hdeti);
  const auto interpVal = mInterp(min_bound, max_bound);
  return std::make_pair(pre*interpVal,
                        Vector(pre*Hdeti*(B - A - interpVal/ei)));
}

//------------------------------------------------------------------------------
// kernelValue, gradValue, grad2Value -- pass through to the base TableKernel
//------------------------------------------------------------------------------
inline
double 
SphericalTableKernel::kernelValue(const double etaMagnitude, const double Hdet) const {
  return mKernel.kernelValue(etaMagnitude, Hdet);
}

inline
double 
SphericalTableKernel::gradValue(const double etaMagnitude, const double Hdet) const {
  return mKernel.gradValue(etaMagnitude, Hdet);
}

inline
typename Dim<1>::Scalar
SphericalTableKernel::grad2Value(const double etaMagnitude, const double Hdet) const {
  return mKernel.grad2Value(etaMagnitude, Hdet);
}

//------------------------------------------------------------------------------
// Data accessors
//------------------------------------------------------------------------------
inline
const typename SphericalTableKernel::InterpolatorType&
SphericalTableKernel::Winterpolator() const {
  return mInterp;
}

inline
const TableKernel<Dim<3>>&
SphericalTableKernel::kernel() const {
  return mKernel;
}

inline
double
SphericalTableKernel::etamax() const {
  return metamax;
}

}
