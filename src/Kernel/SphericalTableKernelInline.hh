#include "Utilities/SpheralFunctions.hh"  // sgn

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
  return 2.0*M_PI/(ei*ej)*Hdeti*mInterp(Dim<2>::Vector(min_bound, max_bound));
}

//------------------------------------------------------------------------------
// Simultaneously lookup (W,  grad W) for (rj/h, ri/h) = (etaj, etai)
//------------------------------------------------------------------------------
inline
std::pair<double, double>
SphericalTableKernel::kernelAndGradValue(const Dim<1>::Vector& etaj,
                                         const Dim<1>::Vector& etai,
                                         const Dim<1>::Scalar  Hdeti) const {
  REQUIRE(Hdeti >= 0.0);
  const auto ei = std::max(1e-10, etai[0]);
  const auto ej = std::max(1e-10, etaj[0]);
  CHECK(ei > 0.0);
  CHECK(ej > 0.0);
  const auto min_bound = std::abs(ej - ei);
  if (min_bound > metamax) return std::make_pair(0.0, 0.0);
  const auto max_bound = std::min(metamax, ei + ej);
  const auto Wval = 2.0*M_PI/(ei*ej)*Hdeti*mInterp(Dim<2>::Vector(min_bound, max_bound));
  return std::make_pair(Wval,
                        2.0*M_PI/(ei*ej)*Hdeti*mGradInterp(Dim<2>::Vector(min_bound, max_bound)) + (ei + ej)/(ei*ej)*Wval);
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
const typename SphericalTableKernel::InterpolatorType&
SphericalTableKernel::gradWinterpolator() const {
  return mGradInterp;
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
