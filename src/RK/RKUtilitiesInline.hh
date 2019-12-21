#include <algorithm>
#include <numeric>
#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Get dot product of two vectors with length "polynomialSize", given offsets
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
template<typename DataType1, typename DataType2>
inline
typename Dimension::Scalar
RKUtilities<Dimension, correctionOrder>::
innerProductRK(const DataType1& x,
               const DataType2& y,
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
template<typename Dimension, RKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
flatSymmetricIndex(const int d1, const int d2) {
  const auto k1 = std::min(d1, d2);
  const auto k2 = std::max(d1, d2);
  return Dimension::nDim * (Dimension::nDim - 1) / 2 - (Dimension::nDim - k1) * (Dimension::nDim - k1 - 1) / 2 + k2;
}

// //------------------------------------------------------------------------------
// // Get storage size of a symmetric matrix
// //------------------------------------------------------------------------------
// template<typename Dimension, RKOrder correctionOrder>
// inline
// int
// RKUtilities<Dimension, correctionOrder>::
// symmetricMatrixSize(const int d) {
//   return d * (d + 1) / 2;
// }

//------------------------------------------------------------------------------
// Get expected length of corrections vector
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
correctionsSize(bool needHessian) {
  return (needHessian ? hessCorrectionsSize : gradCorrectionsSize);
}
template<typename Dimension, RKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
zerothCorrectionsSize(bool needHessian) {
  return (needHessian
          ? RKUtilities<Dimension, RKOrder::ZerothOrder>::hessCorrectionsSize
          : RKUtilities<Dimension, RKOrder::ZerothOrder>::gradCorrectionsSize);
}

//------------------------------------------------------------------------------
// Get offsets for coefficients and polynomials to allow inner products
//------------------------------------------------------------------------------
template<typename Dimension, RKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetGradC(const int d) {
  return polynomialSize * (1 + d);
}

template<typename Dimension, RKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetHessC(const int d1, const int d2) {
  const auto d12 = flatSymmetricIndex(d1, d2);
  return polynomialSize * (1 + Dimension::nDim + d12);
}

template<typename Dimension, RKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetGradP(const int d) {
  return polynomialSize * d;
}

template<typename Dimension, RKOrder correctionOrder>
inline
int
RKUtilities<Dimension, correctionOrder>::
offsetHessP(const int d1, const int d2) {
  const auto d12 = flatSymmetricIndex(d1, d2);
  return polynomialSize * d12;
}

//------------------------------------------------------------------------------
// Get the polynomials
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::ZerothOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::ZerothOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::ZerothOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1};
}

// Linear order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::LinearOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1,x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::LinearOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::LinearOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[2]};
}

// Quadratic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuadraticOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuadraticOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuadraticOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2]};
}

// Cubic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::CubicOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[0]*x[0],x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::CubicOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::CubicOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2]};
}

// Quartic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuarticOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuarticOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuarticOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]};
}

// Quintic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuinticOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuinticOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuinticOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]};
}

// Sextic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::SexticOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::SexticOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::SexticOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2]};
}

// Septic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::SepticOrder>::
getPolynomials(const Dim<1>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::SepticOrder>::
getPolynomials(const Dim<2>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::SepticOrder>::
getPolynomials(const Dim<3>::Vector& x,
               PolyArray& p) {
  p = {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[0]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[0]*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[0]*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[0]*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2]*x[2]};
}

//------------------------------------------------------------------------------
// Get the gradients
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::ZerothOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::ZerothOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,0};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::ZerothOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,0,0};
}

// Linear order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::LinearOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::LinearOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,0,1};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::LinearOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,0,0,1,0,0,0,0,1};
}

// Quadratic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuadraticOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,2*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuadraticOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,2*x[0],x[1],0,0,0,1,0,x[0],2*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuadraticOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,2*x[0],x[1],x[2],0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2]};
}

// Cubic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::CubicOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,2*x[0],3*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::CubicOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::CubicOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2]};
}

// Quartic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuarticOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuarticOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuarticOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2]};
}

// Quintic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuinticOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuinticOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuinticOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],3*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[1]*x[1]*x[1],3*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2]};
}

// Sextic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::SexticOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::SexticOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[1]*x[1],5*x[0]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::SexticOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],3*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[1]*x[1]*x[1],3*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],2*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],0,5*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[1]*x[1]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[2]*x[2],0,x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[2]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2]};
}

// Septic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::SepticOrder>::
getGradPolynomials(const Dim<1>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0],7*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::SepticOrder>::
getGradPolynomials(const Dim<2>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1],0,7*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[1]*x[1],5*x[0]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1],0,x[0]*x[0]*x[0]*x[0]*x[0]*x[0],2*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],3*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],5*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],7*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::SepticOrder>::
getGradPolynomials(const Dim<3>::Vector& x,
                   GradPolyArray& p) {
  p = {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],3*x[0]*x[0]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],5*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,7*x[0]*x[0]*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[0]*x[0]*x[2],5*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],5*x[0]*x[0]*x[0]*x[0]*x[2]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[1]*x[1]*x[2]*x[2],x[1]*x[1]*x[1]*x[2]*x[2]*x[2],x[1]*x[1]*x[2]*x[2]*x[2]*x[2],x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[1]*x[1]*x[1],3*x[0]*x[1]*x[1]*x[2],2*x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[1]*x[1]*x[2],2*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2],0,5*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[1]*x[1]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[0]*x[0]*x[2]*x[2],0,4*x[0]*x[0]*x[0]*x[1]*x[1]*x[1],3*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],2*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],x[0]*x[0]*x[0]*x[2]*x[2]*x[2],0,5*x[0]*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],2*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],x[0]*x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[0]*x[1]*x[1]*x[1]*x[1]*x[1],5*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],4*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],x[0]*x[2]*x[2]*x[2]*x[2]*x[2],0,7*x[1]*x[1]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],5*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[2]*x[2],0,x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[2]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2],0,0,x[0]*x[0]*x[0]*x[0]*x[0]*x[0],0,x[0]*x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[0]*x[2],0,x[0]*x[0]*x[0]*x[0]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[0]*x[2]*x[2],0,x[0]*x[0]*x[0]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[0]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2]*x[2],0,x[0]*x[0]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[0]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[0]*x[2]*x[2]*x[2]*x[2],0,x[0]*x[1]*x[1]*x[1]*x[1]*x[1],2*x[0]*x[1]*x[1]*x[1]*x[1]*x[2],3*x[0]*x[1]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[1]*x[2]*x[2]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[1]*x[2]*x[2]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2]*x[2]*x[2],7*x[2]*x[2]*x[2]*x[2]*x[2]*x[2]};
}

//------------------------------------------------------------------------------
// Get the hessians
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::ZerothOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::ZerothOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::ZerothOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,0,0};
}

// Linear order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::LinearOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::LinearOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,0,0,0,0,0};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::LinearOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
}

// Quadratic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuadraticOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,2};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuadraticOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,2};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuadraticOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2};
}

// Cubic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::CubicOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,2,6*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::CubicOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,2,0,0,6*x[0],2*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::CubicOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2]};
}

// Quartic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuarticOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,2,6*x[0],12*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuarticOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuarticOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2]};
}

// Quintic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::QuinticOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,2,6*x[0],12*x[0]*x[0],20*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::QuinticOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],2*x[1]*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,4*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1],0,0,2*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],12*x[0]*x[1]*x[1],20*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::QuinticOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],12*x[0]*x[0]*x[2],6*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],2*x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,6*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],0,6*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,3*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,2*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,12*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,0,20*x[1]*x[1]*x[1],12*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,3*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,0,4*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,0,2*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],12*x[0]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],12*x[1]*x[2]*x[2],20*x[2]*x[2]*x[2]};
}

// Sextic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::SexticOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,2,6*x[0],12*x[0]*x[0],20*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::SexticOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],2*x[1]*x[1]*x[1],0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,4*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,0,5*x[0]*x[0]*x[0]*x[0],8*x[0]*x[0]*x[0]*x[1],9*x[0]*x[0]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1],0,0,2*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],12*x[0]*x[1]*x[1],20*x[1]*x[1]*x[1],0,0,2*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],20*x[0]*x[1]*x[1]*x[1],30*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::SexticOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],12*x[0]*x[0]*x[2],6*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],2*x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],20*x[0]*x[0]*x[0]*x[2],12*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],6*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2],2*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],2*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,6*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],0,6*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,8*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],0,9*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,8*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],4*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,3*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,2*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,4*x[0]*x[0]*x[0]*x[1],8*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],9*x[0]*x[0]*x[2]*x[2],0,2*x[0]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],8*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,12*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,0,20*x[1]*x[1]*x[1],12*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],2*x[0]*x[0]*x[2]*x[2],0,0,20*x[0]*x[1]*x[1]*x[1],12*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,0,30*x[1]*x[1]*x[1]*x[1],20*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,3*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,0,4*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,0,4*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,0,5*x[1]*x[1]*x[1]*x[1],8*x[1]*x[1]*x[1]*x[2],9*x[1]*x[1]*x[2]*x[2],8*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,0,2*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],12*x[0]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],12*x[1]*x[2]*x[2],20*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[2],0,0,2*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],0,0,2*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[2]*x[2],20*x[0]*x[2]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],20*x[1]*x[2]*x[2]*x[2],30*x[2]*x[2]*x[2]*x[2]};
}

// Septic order
template<>
inline
void
RKUtilities<Dim<1>, RKOrder::SepticOrder>::
getHessPolynomials(const Dim<1>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,2,6*x[0],12*x[0]*x[0],20*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0],42*x[0]*x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
void
RKUtilities<Dim<2>, RKOrder::SepticOrder>::
getHessPolynomials(const Dim<2>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],2*x[1]*x[1]*x[1],0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1],0,0,42*x[0]*x[0]*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0]*x[1],20*x[0]*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,4*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],6*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,0,5*x[0]*x[0]*x[0]*x[0],8*x[0]*x[0]*x[0]*x[1],9*x[0]*x[0]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1],0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],10*x[0]*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[1],10*x[0]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1],0,0,2*x[0]*x[0]*x[0],6*x[0]*x[0]*x[1],12*x[0]*x[1]*x[1],20*x[1]*x[1]*x[1],0,0,2*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[1]*x[1],20*x[0]*x[1]*x[1]*x[1],30*x[1]*x[1]*x[1]*x[1],0,0,2*x[0]*x[0]*x[0]*x[0]*x[0],6*x[0]*x[0]*x[0]*x[0]*x[1],12*x[0]*x[0]*x[0]*x[1]*x[1],20*x[0]*x[0]*x[1]*x[1]*x[1],30*x[0]*x[1]*x[1]*x[1]*x[1],42*x[1]*x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
void
RKUtilities<Dim<3>, RKOrder::SepticOrder>::
getHessPolynomials(const Dim<3>::Vector& x,
                   HessPolyArray& p) {
  p = {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,20*x[0]*x[0]*x[0],12*x[0]*x[0]*x[1],12*x[0]*x[0]*x[2],6*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],2*x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,30*x[0]*x[0]*x[0]*x[0],20*x[0]*x[0]*x[0]*x[1],20*x[0]*x[0]*x[0]*x[2],12*x[0]*x[0]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],6*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2],2*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],2*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,42*x[0]*x[0]*x[0]*x[0]*x[0],30*x[0]*x[0]*x[0]*x[0]*x[1],30*x[0]*x[0]*x[0]*x[0]*x[2],20*x[0]*x[0]*x[0]*x[1]*x[1],20*x[0]*x[0]*x[0]*x[1]*x[2],20*x[0]*x[0]*x[0]*x[2]*x[2],12*x[0]*x[0]*x[1]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[2],12*x[0]*x[0]*x[1]*x[2]*x[2],12*x[0]*x[0]*x[2]*x[2]*x[2],6*x[0]*x[1]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[1]*x[2]*x[2],6*x[0]*x[1]*x[2]*x[2]*x[2],6*x[0]*x[2]*x[2]*x[2]*x[2],2*x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],2*x[1]*x[1]*x[1]*x[2]*x[2],2*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,6*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],0,6*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,8*x[0]*x[0]*x[0]*x[1],4*x[0]*x[0]*x[0]*x[2],0,9*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,8*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],4*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,5*x[1]*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],0,10*x[0]*x[0]*x[0]*x[0]*x[1],5*x[0]*x[0]*x[0]*x[0]*x[2],0,12*x[0]*x[0]*x[0]*x[1]*x[1],8*x[0]*x[0]*x[0]*x[1]*x[2],4*x[0]*x[0]*x[0]*x[2]*x[2],0,12*x[0]*x[0]*x[1]*x[1]*x[1],9*x[0]*x[0]*x[1]*x[1]*x[2],6*x[0]*x[0]*x[1]*x[2]*x[2],3*x[0]*x[0]*x[2]*x[2]*x[2],0,10*x[0]*x[1]*x[1]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[1]*x[2]*x[2],4*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],0,6*x[1]*x[1]*x[1]*x[1]*x[1],5*x[1]*x[1]*x[1]*x[1]*x[2],4*x[1]*x[1]*x[1]*x[2]*x[2],3*x[1]*x[1]*x[2]*x[2]*x[2],2*x[1]*x[2]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,4*x[0]*x[0]*x[0],0,3*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,2*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],6*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,5*x[0]*x[0]*x[0]*x[0],0,4*x[0]*x[0]*x[0]*x[1],8*x[0]*x[0]*x[0]*x[2],0,3*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],9*x[0]*x[0]*x[2]*x[2],0,2*x[0]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],8*x[0]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[2]*x[2],4*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,6*x[0]*x[0]*x[0]*x[0]*x[0],0,5*x[0]*x[0]*x[0]*x[0]*x[1],10*x[0]*x[0]*x[0]*x[0]*x[2],0,4*x[0]*x[0]*x[0]*x[1]*x[1],8*x[0]*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[0]*x[2]*x[2],0,3*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[1]*x[2],9*x[0]*x[0]*x[1]*x[2]*x[2],12*x[0]*x[0]*x[2]*x[2]*x[2],0,2*x[0]*x[1]*x[1]*x[1]*x[1],4*x[0]*x[1]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[1]*x[2]*x[2],8*x[0]*x[1]*x[2]*x[2]*x[2],10*x[0]*x[2]*x[2]*x[2]*x[2],0,x[1]*x[1]*x[1]*x[1]*x[1],2*x[1]*x[1]*x[1]*x[1]*x[2],3*x[1]*x[1]*x[1]*x[2]*x[2],4*x[1]*x[1]*x[2]*x[2]*x[2],5*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,12*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],0,0,20*x[1]*x[1]*x[1],12*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],2*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,12*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],2*x[0]*x[0]*x[2]*x[2],0,0,20*x[0]*x[1]*x[1]*x[1],12*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2],0,0,30*x[1]*x[1]*x[1]*x[1],20*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0]*x[0],0,0,6*x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,0,12*x[0]*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[0]*x[1]*x[2],2*x[0]*x[0]*x[0]*x[2]*x[2],0,0,20*x[0]*x[0]*x[1]*x[1]*x[1],12*x[0]*x[0]*x[1]*x[1]*x[2],6*x[0]*x[0]*x[1]*x[2]*x[2],2*x[0]*x[0]*x[2]*x[2]*x[2],0,0,30*x[0]*x[1]*x[1]*x[1]*x[1],20*x[0]*x[1]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[1]*x[2]*x[2],6*x[0]*x[1]*x[2]*x[2]*x[2],2*x[0]*x[2]*x[2]*x[2]*x[2],0,0,42*x[1]*x[1]*x[1]*x[1]*x[1],30*x[1]*x[1]*x[1]*x[1]*x[2],20*x[1]*x[1]*x[1]*x[2]*x[2],12*x[1]*x[1]*x[2]*x[2]*x[2],6*x[1]*x[2]*x[2]*x[2]*x[2],2*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,3*x[0]*x[1]*x[1],4*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,0,4*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],6*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[2],0,0,3*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[2]*x[2],0,0,4*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],6*x[0]*x[1]*x[2]*x[2],4*x[0]*x[2]*x[2]*x[2],0,0,5*x[1]*x[1]*x[1]*x[1],8*x[1]*x[1]*x[1]*x[2],9*x[1]*x[1]*x[2]*x[2],8*x[1]*x[2]*x[2]*x[2],5*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[0]*x[1],2*x[0]*x[0]*x[0]*x[0]*x[2],0,0,3*x[0]*x[0]*x[0]*x[1]*x[1],4*x[0]*x[0]*x[0]*x[1]*x[2],3*x[0]*x[0]*x[0]*x[2]*x[2],0,0,4*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[1]*x[2],6*x[0]*x[0]*x[1]*x[2]*x[2],4*x[0]*x[0]*x[2]*x[2]*x[2],0,0,5*x[0]*x[1]*x[1]*x[1]*x[1],8*x[0]*x[1]*x[1]*x[1]*x[2],9*x[0]*x[1]*x[1]*x[2]*x[2],8*x[0]*x[1]*x[2]*x[2]*x[2],5*x[0]*x[2]*x[2]*x[2]*x[2],0,0,6*x[1]*x[1]*x[1]*x[1]*x[1],10*x[1]*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[1]*x[2]*x[2],12*x[1]*x[1]*x[2]*x[2]*x[2],10*x[1]*x[2]*x[2]*x[2]*x[2],6*x[2]*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[1],6*x[0]*x[0]*x[2],0,0,2*x[0]*x[1]*x[1],6*x[0]*x[1]*x[2],12*x[0]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1],6*x[1]*x[1]*x[2],12*x[1]*x[2]*x[2],20*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[2],0,0,2*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[2]*x[2],0,0,2*x[0]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[2]*x[2],20*x[0]*x[2]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[2]*x[2],20*x[1]*x[2]*x[2]*x[2],30*x[2]*x[2]*x[2]*x[2],0,0,0,0,0,2*x[0]*x[0]*x[0]*x[0]*x[0],0,0,2*x[0]*x[0]*x[0]*x[0]*x[1],6*x[0]*x[0]*x[0]*x[0]*x[2],0,0,2*x[0]*x[0]*x[0]*x[1]*x[1],6*x[0]*x[0]*x[0]*x[1]*x[2],12*x[0]*x[0]*x[0]*x[2]*x[2],0,0,2*x[0]*x[0]*x[1]*x[1]*x[1],6*x[0]*x[0]*x[1]*x[1]*x[2],12*x[0]*x[0]*x[1]*x[2]*x[2],20*x[0]*x[0]*x[2]*x[2]*x[2],0,0,2*x[0]*x[1]*x[1]*x[1]*x[1],6*x[0]*x[1]*x[1]*x[1]*x[2],12*x[0]*x[1]*x[1]*x[2]*x[2],20*x[0]*x[1]*x[2]*x[2]*x[2],30*x[0]*x[2]*x[2]*x[2]*x[2],0,0,2*x[1]*x[1]*x[1]*x[1]*x[1],6*x[1]*x[1]*x[1]*x[1]*x[2],12*x[1]*x[1]*x[1]*x[2]*x[2],20*x[1]*x[1]*x[2]*x[2]*x[2],30*x[1]*x[2]*x[2]*x[2]*x[2],42*x[2]*x[2]*x[2]*x[2]*x[2]};
}

//------------------------------------------------------------------------------
// Non-templated frontend helper methods
//------------------------------------------------------------------------------
// Get the corrections size based on the correction order and whether Hessian is included
template<typename Dimension>
inline 
std::pair<RKOrder, bool>
RKCorrectionsSize(const RKOrder correctionOrder,
                  const bool needHessian) {
  switch(correctionOrder) {
  case RKOrder::ZerothOrder:
    return RKUtilities<Dimension, RKOrder::ZerothOrder>::correctionsSize(needHessian);
  case RKOrder::LinearOrder:
    return RKUtilities<Dimension, RKOrder::LinearOrder>::correctionsSize(needHessian);
  case RKOrder::QuadraticOrder:
    return RKUtilities<Dimension, RKOrder::QuadraticOrder>::correctionsSize(needHessian);
  case RKOrder::CubicOrder:
    return RKUtilities<Dimension, RKOrder::CubicOrder>::correctionsSize(needHessian);
  case RKOrder::QuarticOrder:
    return RKUtilities<Dimension, RKOrder::QuarticOrder>::correctionsSize(needHessian);
  case RKOrder::QuinticOrder:
    return RKUtilities<Dimension, RKOrder::QuinticOrder>::correctionsSize(needHessian);
  case RKOrder::SexticOrder:
    return RKUtilities<Dimension, RKOrder::SexticOrder>::correctionsSize(needHessian);
  case RKOrder::SepticOrder:
    return RKUtilities<Dimension, RKOrder::SepticOrder>::correctionsSize(needHessian);
  default:
    VERIFY2("Unknown order passed to computeRKCorrections", false);
    return -1;
  }
}

// Get the correction order and whether Hessian is included based on the size of the corrections
template<typename Dimension>
inline 
std::pair<RKOrder, bool>
RKOrderAndHessian(const int size) {
  switch(size) {
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::ZerothOrder, false);
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::ZerothOrder, true);
  case RKUtilities<Dimension, RKOrder::LinearOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::LinearOrder, false);
  case RKUtilities<Dimension, RKOrder::LinearOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::LinearOrder, true);
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::QuadraticOrder, false);
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::QuadraticOrder, true);
  case RKUtilities<Dimension, RKOrder::CubicOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::CubicOrder, false);
  case RKUtilities<Dimension, RKOrder::CubicOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::CubicOrder, true);
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::QuarticOrder, false);
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::QuarticOrder, true);
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::QuinticOrder, false);
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::QuinticOrder, true);
  case RKUtilities<Dimension, RKOrder::SexticOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::SexticOrder, false);
  case RKUtilities<Dimension, RKOrder::SexticOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::SexticOrder, true);
  case RKUtilities<Dimension, RKOrder::SepticOrder>::gradCorrectionsSize:
    return std::make_pair(RKOrder::SepticOrder, false);
  case RKUtilities<Dimension, RKOrder::SepticOrder>::hessCorrectionsSize:
    return std::make_pair(RKOrder::SepticOrder, true);
  default:
    VERIFY2("Unknown order in RKOrderFromSize", false);
    return std::make_pair(RKOrder::ZerothOrder, false);
  }
}

// RK corrected kernel
template<typename Dimension>
inline
typename Dimension::Scalar
RKKernel(const TableKernel<Dimension>& W,
         const typename Dimension::Vector& x,
         const typename Dimension::SymTensor& H,
         const std::vector<double>& corrections) {
  const int size = corrections.size();
  switch(size) {
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernel(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::LinearOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::LinearOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateKernel(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateKernel(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::CubicOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::CubicOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateKernel(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateKernel(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateKernel(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::SexticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::SexticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateKernel(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::SepticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::SepticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateKernel(W, x, H, corrections);
  default:
    VERIFY2("Unknown order passed to RKKernel", false);
    return 0.0;
  }
}

// RK corrected kernel
template<typename Dimension>
inline
typename Dimension::Vector
RKGradient(const TableKernel<Dimension>& W,
           const typename Dimension::Vector& x,
           const typename Dimension::SymTensor& H,
           const std::vector<double>& corrections) {
  const int size = corrections.size();
  switch(size) {
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateGradient(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::LinearOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::LinearOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateGradient(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateGradient(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::CubicOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::CubicOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateGradient(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateGradient(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateGradient(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::SexticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::SexticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateGradient(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::SepticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::SepticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateGradient(W, x, H, corrections);
  default:
    VERIFY2("Unknown order passed to RKGradient", false);
    return Dimension::Vector::zero;
  }
}

// RK corrected Hessian
template<typename Dimension>
inline
typename Dimension::SymTensor
RKHessian(const TableKernel<Dimension>& W,
          const typename Dimension::Vector& x,
          const typename Dimension::SymTensor& H,
          const std::vector<double>& corrections) {
  const int size = corrections.size();
  switch(size) {
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateHessian(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::LinearOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateHessian(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateHessian(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::CubicOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateHessian(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateHessian(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateHessian(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::SexticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateHessian(W, x, H, corrections);
  case RKUtilities<Dimension, RKOrder::SepticOrder>::hessCorrectionsSize:
    return RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateHessian(W, x, H, corrections);
  default:
    VERIFY2("Unknown order passed to RKHessian or grad corrections not included", false);
    return Dimension::SymTensor::zero;
  }
}

// RK corrected kernel + gradient
template<typename Dimension>
inline
void
RKKernelAndGradient(typename Dimension::Scalar& WRK,
                    typename Dimension::Scalar& gradWSPH,
                    typename Dimension::Vector& gradWRK,
                    const TableKernel<Dimension>& W,
                    const typename Dimension::Vector& x,
                    const typename Dimension::SymTensor& H,
                    const std::vector<double>& corrections) {
  gradWSPH = W.gradValue((H*x).magnitude(), H.Determinant());
  const int size = corrections.size();
  switch(size) {
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::ZerothOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  case RKUtilities<Dimension, RKOrder::LinearOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::LinearOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuadraticOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  case RKUtilities<Dimension, RKOrder::CubicOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::CubicOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuarticOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::QuinticOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  case RKUtilities<Dimension, RKOrder::SexticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::SexticOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  case RKUtilities<Dimension, RKOrder::SepticOrder>::gradCorrectionsSize:
  case RKUtilities<Dimension, RKOrder::SepticOrder>::hessCorrectionsSize:
    std::tie(WRK, gradWRK) = RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateKernelAndGradient(W, x, H, corrections);
    break;
  default:
    VERIFY2("Unknown order passed to RKKernelAndGradient", false);
    break;
  }
}

// RK base kernel
template<typename Dimension>
inline
typename Dimension::Scalar
RKBaseKernel(const TableKernel<Dimension>& W,
             const typename Dimension::Vector& x,
             const typename Dimension::SymTensor& H) {
  return RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseKernel(W, x, H);
}

// RK base gradient
template<typename Dimension>
inline
typename Dimension::Vector
RKBaseGradient(const TableKernel<Dimension>& W,
               const typename Dimension::Vector& x,
               const typename Dimension::SymTensor& H) {
  return RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseGradient(W, x, H);
}

// RK base Hessian
template<typename Dimension>
inline
typename Dimension::SymTensor
RKBaseHessian(const TableKernel<Dimension>& W,
              const typename Dimension::Vector& x,
              const typename Dimension::SymTensor& H) {
  return RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseHessian(W, x, H);
}

// RK corrections
template<typename Dimension>
inline
void
computeRKCorrections(const RKOrder order,
                     const ConnectivityMap<Dimension>& connectivityMap,
                     const TableKernel<Dimension>& kernel,
                     const FieldList<Dimension, typename Dimension::Scalar>& volume,
                     const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     const bool needHessian,
                     FieldList<Dimension, std::vector<double>>& zerothCorrections,
                     FieldList<Dimension, std::vector<double>>& corrections) {
  switch(order) {
  case RKOrder::ZerothOrder:
    RKUtilities<Dimension, RKOrder::ZerothOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  case RKOrder::LinearOrder:
    RKUtilities<Dimension, RKOrder::LinearOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  case RKOrder::QuadraticOrder:
    RKUtilities<Dimension, RKOrder::QuadraticOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  case RKOrder::CubicOrder:
    RKUtilities<Dimension, RKOrder::CubicOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  case RKOrder::QuarticOrder:
    RKUtilities<Dimension, RKOrder::QuarticOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  case RKOrder::QuinticOrder:
    RKUtilities<Dimension, RKOrder::QuinticOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  case RKOrder::SexticOrder:
    RKUtilities<Dimension, RKOrder::SexticOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  case RKOrder::SepticOrder:
    RKUtilities<Dimension, RKOrder::SepticOrder>::computeCorrections(connectivityMap, kernel, volume, position, H, needHessian, zerothCorrections, corrections);
    break;
  default:
    VERIFY2("Unknown order passed to computeRKCorrections", false);
  }
}

// Surface normals
template<typename Dimension>
inline
void
computeRKNormal(const RKOrder order,
                const ConnectivityMap<Dimension>& connectivityMap,
                const TableKernel<Dimension>& kernel,
                const FieldList<Dimension, typename Dimension::Scalar>& volume,
                const FieldList<Dimension, typename Dimension::Vector>& position,
                const FieldList<Dimension, typename Dimension::SymTensor>& H,
                FieldList<Dimension, std::vector<double>>& corrections,
                FieldList<Dimension, typename Dimension::Scalar>& surfaceArea,
                FieldList<Dimension, typename Dimension::Vector>& normal) {
  switch(order) {
  case RKOrder::ZerothOrder:
    RKUtilities<Dimension, RKOrder::ZerothOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  case RKOrder::LinearOrder:
    RKUtilities<Dimension, RKOrder::LinearOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  case RKOrder::QuadraticOrder:
    RKUtilities<Dimension, RKOrder::QuadraticOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  case RKOrder::CubicOrder:
    RKUtilities<Dimension, RKOrder::CubicOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  case RKOrder::QuarticOrder:
    RKUtilities<Dimension, RKOrder::QuarticOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  case RKOrder::QuinticOrder:
    RKUtilities<Dimension, RKOrder::QuinticOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  case RKOrder::SexticOrder:
    RKUtilities<Dimension, RKOrder::SexticOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  case RKOrder::SepticOrder:
    RKUtilities<Dimension, RKOrder::SepticOrder>::computeNormal(connectivityMap, kernel, volume, position, H, corrections, surfaceArea, normal);
    break;
  default:
    VERIFY2("Unknown order passed to computeRKNormal", false);
  }
}

} // end namespace Spheral
