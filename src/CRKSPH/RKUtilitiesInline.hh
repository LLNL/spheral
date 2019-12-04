#include <algorithm>
#include <numeric>
#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Get dot product of two vectors with length "polynomialSize", given offsets
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
inline
typename Dimension::Scalar
RKUtilities<Dimension, correctionOrder>::
innerProductRK(const std::vector<double>& x,
               const std::vector<double>& y,
               const int offsetx,
               const int offsety) {
  CHECK(x.size() >= offsetx + polynomialSize);
  CHECK(y.size() >= offsety + polynomialSize);
  return std::inner_product(std::begin(x) + offsetx,
                            std::begin(x) + offsetx + polynomialSize,
                            std::begin(y) + offsety,
                            0.0);
}

//------------------------------------------------------------------------------
// Get flattened index for symmetric matrix
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
flatSymmetricIndex(const int d1, const int d2) {
  const auto dim = Dimension::nDim;
  const auto k1 = std::min(d1, d2);
  const auto k2 = std::max(d1, d2);
  return dim * (dim - 1) / 2 - (dim - k1) * (dim - k1 - 1) / 2 + k2;
}

//------------------------------------------------------------------------------
// Get storage size of a symmetric matrix
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
symmetricMatrixSize(const int d) {
  return d * (d + 1) / 2;
}

//------------------------------------------------------------------------------
// Get expected length of corrections vector
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
correctionsSize(bool needHessian) {
  const auto dim = Dimension::nDim;
  return (needHessian
          ? polynomialSize * (1 + dim + symmetricMatrixSize(dim))
          : polynomialSize * (1 + dim));
}

//------------------------------------------------------------------------------
// Get offsets for coefficients and polynomials to allow inner products
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetGradC(const int d) {
  return polynomialSize * (1 + d);
}

template<typename Dimension, CRKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetHessC(const int d1, const int d2) {
  const auto dim = Dimension::nDim;
  const auto d12 = flatSymmetricIndex(d1, d2);
  return polynomialSize * (1 + dim + d12);
}

template<typename Dimension, CRKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetGradP(const int d) {
  return polynomialSize * d;
}

template<typename Dimension, CRKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetHessP(const int d1, const int d2) {
  const auto dim = Dimension::nDim;
  const auto d12 = flatSymmetricIndex(d1, d2);
  return polynomialSize * d12;
}

//------------------------------------------------------------------------------
// Get the coefficient sizes for each case
//------------------------------------------------------------------------------
// template<> constexpr int RKUtilities<Dim<1>, CRKOrder::ZerothOrder>::polynomialSize = 1;
// template<> constexpr int RKUtilities<Dim<2>, CRKOrder::ZerothOrder>::polynomialSize = 1;
// template<> constexpr int RKUtilities<Dim<3>, CRKOrder::ZerothOrder>::polynomialSize = 1;
// template<> constexpr int RKUtilities<Dim<1>, CRKOrder::LinearOrder>::polynomialSize = 2;
// template<> constexpr int RKUtilities<Dim<2>, CRKOrder::LinearOrder>::polynomialSize = 3;
// template<> constexpr int RKUtilities<Dim<3>, CRKOrder::LinearOrder>::polynomialSize = 4;
// template<> constexpr int RKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::polynomialSize = 3;
// template<> constexpr int RKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::polynomialSize = 6;
// template<> constexpr int RKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::polynomialSize = 10;
// template<> constexpr int RKUtilities<Dim<1>, CRKOrder::CubicOrder>::polynomialSize = 4;
// template<> constexpr int RKUtilities<Dim<2>, CRKOrder::CubicOrder>::polynomialSize = 10;
// template<> constexpr int RKUtilities<Dim<3>, CRKOrder::CubicOrder>::polynomialSize = 20;
// template<> constexpr int RKUtilities<Dim<1>, CRKOrder::QuarticOrder>::polynomialSize = 5;
// template<> constexpr int RKUtilities<Dim<2>, CRKOrder::QuarticOrder>::polynomialSize = 15;
// template<> constexpr int RKUtilities<Dim<3>, CRKOrder::QuarticOrder>::polynomialSize = 35;
// template<> constexpr int RKUtilities<Dim<1>, CRKOrder::QuinticOrder>::polynomialSize = 6;
// template<> constexpr int RKUtilities<Dim<2>, CRKOrder::QuinticOrder>::polynomialSize = 21;
// template<> constexpr int RKUtilities<Dim<3>, CRKOrder::QuinticOrder>::polynomialSize = 56;

//------------------------------------------------------------------------------
// Get the polynomials
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::ZerothOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::ZerothOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::ZerothOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1};
}

// Linear order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::LinearOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1,x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::LinearOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1,x[0],x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::LinearOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1,x[0],x[1],x[2]};
}

// Quadratic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1,x[0],x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2]};
}

// Cubic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::CubicOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::CubicOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::CubicOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2]};
}

// Quartic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuarticOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuarticOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuarticOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]};
}

// Quintic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuinticOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuinticOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuinticOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]};
}

// Sextic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::SexticOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::SexticOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::SexticOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2]};
}

// Septic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::SepticOrder>::
getPolynomials(const Dim<1>::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::SepticOrder>::
getPolynomials(const Dim<2>::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::SepticOrder>::
getPolynomials(const Dim<3>::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2]};
}

//------------------------------------------------------------------------------
// Get the gradients
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::ZerothOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::ZerothOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,0};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::ZerothOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0};
}

// Linear order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::LinearOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0,1};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::LinearOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,1,0,0,0,1};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::LinearOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,1,0,0,0,0,1,0,0,0,0,1};
}

// Quadratic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0,1,2*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,1,0,2*x[0],x[1],0,0,0,1,0,x[0],2*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2]};
}

// Cubic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::CubicOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0,1,2*x[0],3*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::CubicOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::CubicOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2]};
}

// Quartic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuarticOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuarticOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuarticOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2]};
}

// Quintic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuinticOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuinticOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuinticOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],3*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[1]*x[1]*x[1],3*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2]};
}

// Sextic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::SexticOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::SexticOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[1]*x[1],5*x[0]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::SexticOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],3*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[1]*x[1]*x[1],3*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],2*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],0,5*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[1]*x[1]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[2]*x[2],0,x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[2]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2]};
}

// Septic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::SepticOrder>::
getGradPolynomials(const Dim<1>::Vector& x) {
  return {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0],7*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::SepticOrder>::
getGradPolynomials(const Dim<2>::Vector& x) {
  return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],0,7*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[1]*x[1],5*x[0]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],5*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],7*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::SepticOrder>::
getGradPolynomials(const Dim<3>::Vector& x) {
  return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],3*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,7*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[0]*x[0]*x[2],5*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],5*x[0]*x[0]*x[0]*x[0]*x[2]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[1]*x[1]*x[1],3*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],2*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],0,5*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[1]*x[1]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],2*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2],0,5*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],5*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],4*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2],0,7*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],5*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[2]*x[2],0,x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[2]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[0]*x[2]*x[2],0,x[0]*x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2]*x[2],0,x[0]*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[0]*x[2]*x[2]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],7*x[2]*x[2]*x[2]*x[2]*x[2]*x[2]};
}

//------------------------------------------------------------------------------
// Get the hessians
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::ZerothOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::ZerothOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::ZerothOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,0,0};
}

// Linear order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::LinearOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0,0};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::LinearOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0,0,0,0,0,0,0};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::LinearOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
}

// Quadratic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0,0,2};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,2};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2};
}

// Cubic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::CubicOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0,0,2,6*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::CubicOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::CubicOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2]};
}

// Quartic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuarticOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0,0,2,6*x[0],12*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuarticOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuarticOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2]};
}

// Quintic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::QuinticOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0,0,2,6*x[0],12*x[0]*x[0],20*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::QuinticOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],2*x[1]*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,4*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1],0,0,2*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],12*x[0]*x[1]*x[1],20*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::QuinticOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],12*x[0]*x[0]*x[2],6*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],2*x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,6*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],0,6*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,3*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,2*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,12*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,0,20*x[1]*x[1]*x[1],12*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,3*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,0,4*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,0,2*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],12*x[0]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],12*x[1]*x[2]*x[2],20*x[2]*x[2]*x[2]};
}

// Sextic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::SexticOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0,0,2,6*x[0],12*x[0]*x[0],20*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::SexticOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],2*x[1]*x[1]*x[1],0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,4*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,0,5*x[0]*x[0]*x[0]*x[0],8*x[0]*x[0]*x[0]*x[1],9*x[0]*x[0]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1],0,0,2*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],12*x[0]*x[1]*x[1],20*x[1]*x[1]*x[1],0,0,2*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],20*x[0]*x[1]*x[1]*x[1],30*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::SexticOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],12*x[0]*x[0]*x[2],6*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],2*x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],20*x[0]*x[0]*x[0]*x[2],12*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],6*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2],2*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],2*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,6*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],0,6*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,8*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],0,9*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,8*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],4*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,3*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,2*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,4*x[0]*x[0]*x[0]*x[1],8*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],9*x[0]*x[0]*x[2]*x[2],0,2*x[0]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],8*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,12*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,0,20*x[1]*x[1]*x[1],12*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],2*x[0]*x[0]*x[2]*x[2],0,0,20*x[0]*x[1]*x[1]*x[1],12*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,0,30*x[1]*x[1]*x[1]*x[1],20*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,3*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,0,4*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,0,4*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,0,5*x[1]*x[1]*x[1]*x[1],8*x[1]*x[1]*x[1]*x[2],9*x[1]*x[1]*x[2]*x[2],8*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,0,2*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],12*x[0]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],12*x[1]*x[2]*x[2],20*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[2],0,0,2*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],0,0,2*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[2]*x[2],20*x[0]*x[2]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],20*x[1]*x[2]*x[2]*x[2],30*x[2]*x[2]*x[2]*x[2]};
}

// Septic order
template<>
inline
std::vector<double>
RKUtilities<Dim<1>, CRKOrder::SepticOrder>::
getHessPolynomials(const Dim<1>::Vector& x) {
  return {0,0,2,6*x[0],12*x[0]*x[0],20*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0],42*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<2>, CRKOrder::SepticOrder>::
getHessPolynomials(const Dim<2>::Vector& x) {
  return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],2*x[1]*x[1]*x[1],0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1],0,0,42*x[0]*x[0]*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0]*x[1],20*x[0]*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,4*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,0,5*x[0]*x[0]*x[0]*x[0],8*x[0]*x[0]*x[0]*x[1],9*x[0]*x[0]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],10*x[0]*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[1],10*x[0]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1],0,0,2*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],12*x[0]*x[1]*x[1],20*x[1]*x[1]*x[1],0,0,2*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],20*x[0]*x[1]*x[1]*x[1],30*x[1]*x[1]*x[1]*x[1],0,0,2*x[0]*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[0]*x[1]*x[1],20*x[0]*x[0]*x[1]*x[1]*x[1],30*x[0]*x[1]*x[1]*x[1]*x[1],42*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
RKUtilities<Dim<3>, CRKOrder::SepticOrder>::
getHessPolynomials(const Dim<3>::Vector& x) {
  return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],12*x[0]*x[0]*x[2],6*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],2*x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],20*x[0]*x[0]*x[0]*x[2],12*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],6*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2],2*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],2*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,42*x[0]*x[0]*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0]*x[1],30*x[0]*x[0]*x[0]*x[0]*x[2],20*x[0]*x[0]*x[0]*x[1]*x[1],20*x[0]*x[0]*x[0]*x[1]*x[2],20*x[0]*x[0]*x[0]*x[2]*x[2],12*x[0]*x[0]*x[1]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[2],12*x[0]*x[0]*x[1]*x[2]*x[2],12*x[0]*x[0]*x[2]*x[2]*x[2],6*x[0]*x[1]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[1]*x[2]*x[2],6*x[0]*x[1]*x[2]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2]*x[2],2*x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],2*x[1]*x[1]*x[1]*x[2]*x[2],2*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,6*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],0,6*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,8*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],0,9*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,8*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],4*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],0,10*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[2],0,12*x[0]*x[0]*x[0]*x[1]*x[1],8*x[0]*x[0]*x[0]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2],0,12*x[0]*x[0]*x[1]*x[1]*x[1],9*x[0]*x[0]*x[1]*x[1]*x[2],6*x[0]*x[0]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2],0,10*x[0]*x[1]*x[1]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[1]*x[1]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,3*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,2*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,4*x[0]*x[0]*x[0]*x[1],8*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],9*x[0]*x[0]*x[2]*x[2],0,2*x[0]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],8*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],0,5*x[0]*x[0]*x[0]*x[0]*x[1],10*x[0]*x[0]*x[0]*x[0]*x[2],0,4*x[0]*x[0]*x[0]*x[1]*x[1],8*x[0]*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[0]*x[2]*x[2],0,3*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[1]*x[2],9*x[0]*x[0]*x[1]*x[2]*x[2],12*x[0]*x[0]*x[2]*x[2]*x[2],0,2*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[1]*x[2]*x[2],8*x[0]*x[1]*x[2]*x[2]*x[2],10*x[0]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,12*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,0,20*x[1]*x[1]*x[1],12*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],2*x[0]*x[0]*x[2]*x[2],0,0,20*x[0]*x[1]*x[1]*x[1],12*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,0,30*x[1]*x[1]*x[1]*x[1],20*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,0,12*x[0]*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[0]*x[1]*x[2],2*x[0]*x[0]*x[0]*x[2]*x[2],0,0,20*x[0]*x[0]*x[1]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[2],6*x[0]*x[0]*x[1]*x[2]*x[2],2*x[0]*x[0]*x[2]*x[2]*x[2],0,0,30*x[0]*x[1]*x[1]*x[1]*x[1],20*x[0]*x[1]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[1]*x[2]*x[2],6*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],0,0,42*x[1]*x[1]*x[1]*x[1]*x[1],30*x[1]*x[1]*x[1]*x[1]*x[2],20*x[1]*x[1]*x[1]*x[2]*x[2],12*x[1]*x[1]*x[2]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,3*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,0,4*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,0,4*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,0,5*x[1]*x[1]*x[1]*x[1],8*x[1]*x[1]*x[1]*x[2],9*x[1]*x[1]*x[2]*x[2],8*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,0,3*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[2]*x[2],0,0,4*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[1]*x[2],6*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[2]*x[2]*x[2],0,0,5*x[0]*x[1]*x[1]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1]*x[2],9*x[0]*x[1]*x[1]*x[2]*x[2],8*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[2]*x[2]*x[2]*x[2],0,0,6*x[1]*x[1]*x[1]*x[1]*x[1],10*x[1]*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[1]*x[2]*x[2],12*x[1]*x[1]*x[2]*x[2]*x[2],10*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,0,2*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],12*x[0]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],12*x[1]*x[2]*x[2],20*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[2],0,0,2*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],0,0,2*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[2]*x[2],20*x[0]*x[2]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],20*x[1]*x[2]*x[2]*x[2],30*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[0]*x[2],0,0,2*x[0]*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[0]*x[2]*x[2],0,0,2*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[1]*x[2],12*x[0]*x[0]*x[1]*x[2]*x[2],20*x[0]*x[0]*x[2]*x[2]*x[2],0,0,2*x[0]*x[1]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[1]*x[2]*x[2],20*x[0]*x[1]*x[2]*x[2]*x[2],30*x[0]*x[2]*x[2]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[1]*x[2]*x[2],20*x[1]*x[1]*x[2]*x[2]*x[2],30*x[1]*x[2]*x[2]*x[2]*x[2],42*x[2]*x[2]*x[2]*x[2]*x[2]};
}

} // end namespace Spheral
