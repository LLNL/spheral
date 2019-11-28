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
SuperiorRKUtilities<Dimension, correctionOrder>::
innerProductRK(const std::vector<double>& x,
               const std::vector<double>& y,
               const int offsetx,
               const int offsety) {
  return std::inner_product(std::begin(x) + offsetx, std::end(x) + offsetx + polynomialSize, std::begin(y) + offsety, 0.0);
}

//------------------------------------------------------------------------------
// Get offsets for coefficients and polynomials to allow inner products
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
inline
int
SuperiorRKUtilities<Dimension, correctionOrder>::
offsetGradC(const int d) {
  return polynomialSize * (1 + d);
}

template<typename Dimension, CRKOrder correctionOrder>
inline
int
SuperiorRKUtilities<Dimension, correctionOrder>::
offsetHessC(const int d1, const int d2) {
  const auto dim = Dimension::nDim;
  const auto k1 = std::min(d1, d2);
  const auto k2 = std::max(d1, d2);
  const auto k = dim * (dim - 1) / 2 - (dim - k1) * (dim - k1 - 1) / 2 + k2;
  return polynomialSize * (1 + dim + k);
}

template<typename Dimension, CRKOrder correctionOrder>
inline
int
SuperiorRKUtilities<Dimension, correctionOrder>::
offsetGradP(const int d) {
  return polynomialSize * d;
}

template<typename Dimension, CRKOrder correctionOrder>
inline
int
SuperiorRKUtilities<Dimension, correctionOrder>::
offsetHessP(const int d1, const int d2) {
  const auto dim = Dimension::nDim;
  const auto k1 = std::min(d1, d2);
  const auto k2 = std::max(d1, d2);
  const auto k = dim * (dim - 1) / 2 - (dim - k1) * (dim - k1 - 1) / 2 + k2;
  return polynomialSize * k;
}


//------------------------------------------------------------------------------
// Get the coefficient sizes for each case
//------------------------------------------------------------------------------
int SuperiorRKUtilities<Dim<1>, CRKOrder::ZerothOrder>::polynomialSize = 1;
int SuperiorRKUtilities<Dim<2>, CRKOrder::ZerothOrder>::polynomialSize = 1;
int SuperiorRKUtilities<Dim<3>, CRKOrder::ZerothOrder>::polynomialSize = 1;
int SuperiorRKUtilities<Dim<1>, CRKOrder::LinearOrder>::polynomialSize = 2;
int SuperiorRKUtilities<Dim<2>, CRKOrder::LinearOrder>::polynomialSize = 3;
int SuperiorRKUtilities<Dim<3>, CRKOrder::LinearOrder>::polynomialSize = 4;
int SuperiorRKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::polynomialSize = 3;
int SuperiorRKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::polynomialSize = 6;
int SuperiorRKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::polynomialSize = 10;
int SuperiorRKUtilities<Dim<1>, CRKOrder::CubicOrder>::polynomialSize = 4;
int SuperiorRKUtilities<Dim<2>, CRKOrder::CubicOrder>::polynomialSize = 10;
int SuperiorRKUtilities<Dim<3>, CRKOrder::CubicOrder>::polynomialSize = 20;
int SuperiorRKUtilities<Dim<1>, CRKOrder::QuarticOrder>::polynomialSize = 5;
int SuperiorRKUtilities<Dim<2>, CRKOrder::QuarticOrder>::polynomialSize = 15;
int SuperiorRKUtilities<Dim<3>, CRKOrder::QuarticOrder>::polynomialSize = 35;
int SuperiorRKUtilities<Dim<1>, CRKOrder::QuinticOrder>::polynomialSize = 6;
int SuperiorRKUtilities<Dim<2>, CRKOrder::QuinticOrder>::polynomialSize = 21;
int SuperiorRKUtilities<Dim<3>, CRKOrder::QuinticOrder>::polynomialSize = 56;

//------------------------------------------------------------------------------
// Get the polynomials
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::ZerothOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::ZerothOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::ZerothOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1};
}

// Linear order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::LinearOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::LinearOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::LinearOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[2]};
}

// Quadratic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2]};
}

// Cubic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::CubicOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::CubicOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::CubicOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2]};
}

// Quartic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuarticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuarticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuarticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]};
}

// Quintic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuinticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[0]*x[0],x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuinticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[0]*x[0],x[0]*x[1],x[1]*x[1],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[1]*x[1],x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuinticOrder>::
getPolynomials(typename Dimension::Vector& x) {
  return {1,x[0],x[1],x[2],x[0]*x[0],x[0]*x[1],x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],x[0]*x[0]*x[0],x[0]*x[0]*x[1],x[0]*x[0]*x[2],x[0]*x[1]*x[1],x[0]*x[1]*x[2],x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[0],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[1],x[0]*x[0]*x[0]*x[2],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[1],x[0]*x[0]*x[1]*x[2],x[0]*x[0]*x[2]*x[2],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[1],x[0]*x[1]*x[1]*x[2],x[0]*x[1]*x[2]*x[2],x[0]*x[2]*x[2]*x[2],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[1],x[1]*x[1]*x[1]*x[2],x[1]*x[1]*x[2]*x[2],x[1]*x[2]*x[2]*x[2],x[2]*x[2]*x[2]*x[2]};
}

//------------------------------------------------------------------------------
// Get the gradients
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::ZerothOrder>::
getDPolynomials(typename Dimension::Vector& x) {
  return {0};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::ZerothOrder>::
getDPolynomials(typename Dimension::Vector& x) {
  return {0,0};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::ZerothOrder>::
getDPolynomials(typename Dimension::Vector& x) {
  return {0,0,0};
}

// Linear order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::LinearOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::LinearOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,0,0,1};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::LinearOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,0,0,1};
}

// Quadratic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,2*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,2*x[0],x[1],0,0,0,1,0,x[0],2*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2]};
}

// Cubic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::CubicOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::CubicOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::CubicOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2]};
}

// Quartic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuarticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuarticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuarticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2]};
}

// Quintic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuinticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,2*x[0],3*x[0]*x[0],4*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuinticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,2*x[0],x[1],0,3*x[0]*x[0],2*x[0]*x[1],x[1]*x[1],0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,4*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],2*x[0]*x[1]*x[1],x[1]*x[1]*x[1],0,0,0,1,0,x[0],2*x[1],0,x[0]*x[0],2*x[0]*x[1],3*x[1]*x[1],0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1],0,0,x[0]*x[0]*x[0],2*x[0]*x[0]*x[1],3*x[0]*x[1]*x[1],4*x[1]*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuinticOrder>::
getDPolynomials(typename Dimension::Vector& x) {
    return {0,1,0,0,2*x[0],x[1],x[2],0,0,0,3*x[0]*x[0],2*x[0]*x[1],2*x[0]*x[2],x[1]*x[1],x[1]*x[2],x[2]*x[2],0,0,0,0,4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,4*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0],4*x[0]*x[0]*x[0],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[1],3*x[0]*x[0]*x[2],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],2*x[0]*x[2]*x[2],x[1]*x[1]*x[1],x[1]*x[1]*x[1],x[1]*x[1]*x[2],x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,0,0,0,0,1,0,0,x[0],0,2*x[1],x[2],0,0,x[0]*x[0],0,2*x[0]*x[1],x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,x[0]*x[0]*x[0],x[0]*x[0]*x[0],0,2*x[0]*x[0]*x[1],2*x[0]*x[0]*x[1],x[0]*x[0]*x[2],0,3*x[0]*x[1]*x[1],3*x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],x[0]*x[2]*x[2],0,4*x[1]*x[1]*x[1],4*x[1]*x[1]*x[1],3*x[1]*x[1]*x[2],2*x[1]*x[2]*x[2],x[2]*x[2]*x[2],0,0,0,0,1,0,0,x[0],0,x[1],2*x[2],0,0,x[0]*x[0],0,x[0]*x[1],2*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,x[0]*x[0]*x[0],0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2],0,0,0,0,0,x[0]*x[0]*x[0],0,0,x[0]*x[0]*x[1],2*x[0]*x[0]*x[2],0,0,x[0]*x[1]*x[1],2*x[0]*x[1]*x[2],3*x[0]*x[2]*x[2],0,0,x[1]*x[1]*x[1],2*x[1]*x[1]*x[2],3*x[1]*x[2]*x[2],4*x[2]*x[2]*x[2]};
}

//------------------------------------------------------------------------------
// Get the hessians
//------------------------------------------------------------------------------

// Zeroth order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::ZerothOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
  return {0};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::ZerothOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::ZerothOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,0,0,0};
}

// Linear order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::LinearOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::LinearOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,0,0,0,0,0,0};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::LinearOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
}

// Quadratic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuadraticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,2};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuadraticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,2};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuadraticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2};
}

// Cubic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::CubicOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,2,6*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::CubicOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::CubicOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2]};
}

// Quartic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuarticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,2,6*x[0],12*x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuarticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuarticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2]};
}

// Quintic order
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<1>, CRKOrder::QuinticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,2,6*x[0],12*x[0]*x[0],12*x[0]*x[0]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<2>, CRKOrder::QuinticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,2,0,0,6*x[0],2*x[1],0,0,12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,12*x[0]*x[0],12*x[0]*x[0],6*x[0]*x[1],2*x[1]*x[1],0,0,0,0,0,0,1,0,0,2*x[0],2*x[1],0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,0,3*x[0]*x[0],4*x[0]*x[1],3*x[1]*x[1],0,0,0,0,0,0,2,0,0,2*x[0],6*x[1],0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1],0,0,0,2*x[0]*x[0],6*x[0]*x[1],12*x[1]*x[1]};
}
template<>
inline
std::vector<double>
SuperiorRKUtilities<Dim<3>, CRKOrder::QuinticOrder>::
getDDPolynomials(typename Dimension::Vector& x) {
    return {0,0,0,0,2,0,0,0,0,0,6*x[0],2*x[1],2*x[2],0,0,0,0,0,0,0,12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,12*x[0]*x[0],12*x[0]*x[0],12*x[0]*x[0],6*x[0]*x[1],6*x[0]*x[1],6*x[0]*x[2],2*x[1]*x[1],2*x[1]*x[1],2*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,2*x[1],x[2],0,0,0,0,0,0,3*x[0]*x[0],0,4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,0,0,3*x[0]*x[0],3*x[0]*x[0],0,4*x[0]*x[1],4*x[0]*x[1],2*x[0]*x[2],0,3*x[1]*x[1],3*x[1]*x[1],2*x[1]*x[2],x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,2*x[0],0,x[1],2*x[2],0,0,0,0,0,0,3*x[0]*x[0],0,2*x[0]*x[1],4*x[0]*x[2],0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,3*x[0]*x[0],0,0,2*x[0]*x[1],4*x[0]*x[2],0,0,x[1]*x[1],2*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,6*x[1],2*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,2*x[0]*x[0],2*x[0]*x[0],0,0,6*x[0]*x[1],6*x[0]*x[1],2*x[0]*x[2],0,0,12*x[1]*x[1],12*x[1]*x[1],6*x[1]*x[2],2*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,x[0],0,0,2*x[1],2*x[2],0,0,0,0,0,x[0]*x[0],0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,x[0]*x[0],0,0,0,2*x[0]*x[1],2*x[0]*x[2],0,0,0,3*x[1]*x[1],4*x[1]*x[2],3*x[2]*x[2],0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,2*x[0],0,0,2*x[1],6*x[2],0,0,0,0,0,2*x[0]*x[0],0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2],0,0,0,0,0,0,0,0,0,2*x[0]*x[0],0,0,0,2*x[0]*x[1],6*x[0]*x[2],0,0,0,2*x[1]*x[1],6*x[1]*x[2],12*x[2]*x[2]};
}
