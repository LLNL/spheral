//---------------------------------Spheral++----------------------------------//
// A collection of useful helper methods to explicitly unroll Dimensional 
// loops for efficiency.
//----------------------------------------------------------------------------//
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct a Tensor with the given diagonal elements, enforcing a minimum.
//------------------------------------------------------------------------------
inline
Dim<1>::Tensor
constructTensorWithMinDiagonal(const Dim<1>::Vector& vec,
                               const double maxval) {
  return Dim<1>::Tensor(std::min(maxval, vec.x()));
}

inline
Dim<2>::Tensor
constructTensorWithMinDiagonal(const Dim<2>::Vector& vec,
                               const double maxval) {
  return Dim<2>::Tensor(std::min(maxval, vec.x()), 0.0,
                        0.0, std::min(maxval, vec.y()));
}

inline
Dim<3>::Tensor
constructTensorWithMinDiagonal(const Dim<3>::Vector& vec,
                               const double maxval) {
  return Dim<3>::Tensor(std::min(maxval, vec.x()), 0.0, 0.0,
                        0.0, std::min(maxval, vec.y()), 0.0,
                        0.0, 0.0, std::min(maxval, vec.z()));
}

inline
Dim<1>::SymTensor
constructSymTensorWithMinDiagonal(const Dim<1>::Vector& vec,
                                  const double maxval) {
  return Dim<1>::SymTensor(std::min(maxval, vec.x()));
}

inline
Dim<2>::SymTensor
constructSymTensorWithMinDiagonal(const Dim<2>::Vector& vec,
                                  const double maxval) {
  return Dim<2>::SymTensor(std::min(maxval, vec.x()), 0.0,
                           0.0, std::min(maxval, vec.y()));
}

inline
Dim<3>::SymTensor
constructSymTensorWithMinDiagonal(const Dim<3>::Vector& vec,
                                  const double maxval) {
  return Dim<3>::SymTensor(std::min(maxval, vec.x()), 0.0, 0.0,
                           0.0, std::min(maxval, vec.y()), 0.0,
                           0.0, 0.0, std::min(maxval, vec.z()));
}

//------------------------------------------------------------------------------
// Construct a Tensor with the given diagonal elements, enforcing a maximum.
//------------------------------------------------------------------------------
inline
Dim<1>::Tensor
constructTensorWithMaxDiagonal(const Dim<1>::Vector& vec,
                               const double minval) {
  return Dim<1>::Tensor(std::max(minval, vec.x()));
}

inline
Dim<2>::Tensor
constructTensorWithMaxDiagonal(const Dim<2>::Vector& vec,
                               const double minval) {
  return Dim<2>::Tensor(std::max(minval, vec.x()), 0.0,
                        0.0, std::max(minval, vec.y()));
}

inline
Dim<3>::Tensor
constructTensorWithMaxDiagonal(const Dim<3>::Vector& vec,
                               const double minval) {
  return Dim<3>::Tensor(std::max(minval, vec.x()), 0.0, 0.0,
                        0.0, std::max(minval, vec.y()), 0.0,
                        0.0, 0.0, std::max(minval, vec.z()));
}

inline
Dim<1>::SymTensor
constructSymTensorWithMaxDiagonal(const Dim<1>::Vector& vec,
                                  const double minval) {
  return Dim<1>::SymTensor(std::max(minval, vec.x()));
}

inline
Dim<2>::SymTensor
constructSymTensorWithMaxDiagonal(const Dim<2>::Vector& vec,
                                  const double minval) {
  return Dim<2>::SymTensor(std::max(minval, vec.x()), 0.0,
                           0.0, std::max(minval, vec.y()));
}

inline
Dim<3>::SymTensor
constructSymTensorWithMaxDiagonal(const Dim<3>::Vector& vec,
                                  const double minval) {
  return Dim<3>::SymTensor(std::max(minval, vec.x()), 0.0, 0.0,
                           0.0, std::max(minval, vec.y()), 0.0,
                           0.0, 0.0, std::max(minval, vec.z()));
}

//------------------------------------------------------------------------------
// Construct a Tensor with the given diagonal elements, constrained to the 
// given range.
//------------------------------------------------------------------------------
inline
Dim<1>::Tensor
constructTensorWithBoundedDiagonal(const Dim<1>::Vector& vec,
                                   const double minval,
                                   const double maxval) {
  return Dim<1>::Tensor(std::max(minval, std::min(maxval, vec.x())));
}

inline
Dim<2>::Tensor
constructTensorWithBoundedDiagonal(const Dim<2>::Vector& vec,
                                   const double minval,
                                   const double maxval) {
  return Dim<2>::Tensor(std::max(minval, std::min(maxval, vec.x())), 0.0,
                        0.0, std::max(minval, std::min(maxval, vec.y())));
}

inline
Dim<3>::Tensor
constructTensorWithBoundedDiagonal(const Dim<3>::Vector& vec,
                                   const double minval,
                                   const double maxval) {
  return Dim<3>::Tensor(std::max(minval, std::min(maxval, vec.x())), 0.0, 0.0,
                        0.0, std::max(minval, std::min(maxval, vec.y())), 0.0,
                        0.0, 0.0, std::max(minval, std::min(maxval, vec.z())));
}

inline
Dim<1>::SymTensor
constructSymTensorWithBoundedDiagonal(const Dim<1>::Vector& vec,
                                      const double minval,
                                      const double maxval) {
  return Dim<1>::SymTensor(std::max(minval, std::min(maxval, vec.x())));
}

inline
Dim<2>::SymTensor
constructSymTensorWithBoundedDiagonal(const Dim<2>::Vector& vec,
                                      const double minval,
                                      const double maxval) {
  return Dim<2>::SymTensor(std::max(minval, std::min(maxval, vec.x())), 0.0,
                           0.0, std::max(minval, std::min(maxval, vec.y())));
}

inline
Dim<3>::SymTensor
constructSymTensorWithBoundedDiagonal(const Dim<3>::Vector& vec,
                                      const double minval,
                                      const double maxval) {
  return Dim<3>::SymTensor(std::max(minval, std::min(maxval, vec.x())), 0.0, 0.0,
                           0.0, std::max(minval, std::min(maxval, vec.y())), 0.0,
                           0.0, 0.0, std::max(minval, std::min(maxval, vec.z())));
}

//------------------------------------------------------------------------------
// Construct a SymTensor with the given diagonal elements, raised to a power.
//------------------------------------------------------------------------------
inline
Dim<1>::SymTensor
constructSymTensorWithPowDiagonal(const Dim<1>::Vector& vec,
                                  const double exponent) {
  return Dim<1>::SymTensor(pow(vec.x(), exponent));
}

inline
Dim<2>::SymTensor
constructSymTensorWithPowDiagonal(const Dim<2>::Vector& vec,
                                  const double exponent) {
  return Dim<2>::SymTensor(pow(vec.x(), exponent), 0.0,
                           0.0, pow(vec.y(), exponent));
}

inline
Dim<3>::SymTensor
constructSymTensorWithPowDiagonal(const Dim<3>::Vector& vec,
                                  const double exponent) {
  return Dim<3>::SymTensor(pow(vec.x(), exponent), 0.0, 0.0,
                           0.0, pow(vec.y(), exponent), 0.0,
                           0.0, 0.0, pow(vec.z(), exponent));
}

//------------------------------------------------------------------------------
// Construct a tensor with the given value in the given row.
//------------------------------------------------------------------------------
template<typename TensorType, size_t col>
inline
TensorType
constructTensorWithColumnValue(const double val) {
  VERIFY(false);
}

template<>
inline
Dim<1>::Tensor
constructTensorWithColumnValue<Dim<1>::Tensor, (size_t) 0>(const double val) {
  return Dim<1>::Tensor(val);
}

template<>
inline
Dim<2>::Tensor
constructTensorWithColumnValue<Dim<2>::Tensor, (size_t) 0>(const double val) {
  return Dim<2>::Tensor(val, 0.0,
                        val, 0.0);
}

template<>
inline
Dim<2>::Tensor
constructTensorWithColumnValue<Dim<2>::Tensor, (size_t) 1>(const double val) {
  return Dim<2>::Tensor(0.0, val,
                        0.0, val);
}

template<>
inline
Dim<3>::Tensor
constructTensorWithColumnValue<Dim<3>::Tensor, (size_t) 0>(const double val) {
  return Dim<3>::Tensor(val, 0.0, 0.0,
                        val, 0.0, 0.0,
                        val, 0.0, 0.0);
}

template<>
inline
Dim<3>::Tensor
constructTensorWithColumnValue<Dim<3>::Tensor, (size_t) 1>(const double val) {
  return Dim<3>::Tensor(0.0, val, 0.0,
                        0.0, val, 0.0,
                        0.0, val, 0.0);
}

template<>
inline
Dim<3>::Tensor
constructTensorWithColumnValue<Dim<3>::Tensor, (size_t) 2>(const double val) {
  return Dim<3>::Tensor(0.0, 0.0, val,
                        0.0, 0.0, val,
                        0.0, 0.0, val);
}

//------------------------------------------------------------------------------
// In place addition of the absolute value of the elements of a tensor.
//------------------------------------------------------------------------------
template<typename TensorType>
inline
void
inPlaceAbsAdd(TensorType& lhs, const TensorType& rhs) {
  VERIFY(false);
}

template<>
inline
void
inPlaceAbsAdd<Dim<1>::Tensor>(Dim<1>::Tensor& lhs, const Dim<1>::Tensor& rhs) {
  lhs.xx(lhs.xx() + std::abs(rhs.xx()));
}

template<>
inline
void
inPlaceAbsAdd<Dim<2>::Tensor>(Dim<2>::Tensor& lhs, const Dim<2>::Tensor& rhs) {
  lhs.xx(lhs.xx() + std::abs(rhs.xx()));
  lhs.xy(lhs.xy() + std::abs(rhs.xy()));
  lhs.yx(lhs.yx() + std::abs(rhs.yx()));
  lhs.yy(lhs.yy() + std::abs(rhs.yy()));
}

template<>
inline
void
inPlaceAbsAdd<Dim<3>::Tensor>(Dim<3>::Tensor& lhs, const Dim<3>::Tensor& rhs) {
  lhs.xx(lhs.xx() + std::abs(rhs.xx()));
  lhs.xy(lhs.xy() + std::abs(rhs.xy()));
  lhs.xz(lhs.xz() + std::abs(rhs.xz()));
  lhs.yx(lhs.yx() + std::abs(rhs.yx()));
  lhs.yy(lhs.yy() + std::abs(rhs.yy()));
  lhs.yz(lhs.yz() + std::abs(rhs.yz()));
  lhs.zx(lhs.zx() + std::abs(rhs.zx()));
  lhs.zy(lhs.zy() + std::abs(rhs.zy()));
  lhs.zz(lhs.zz() + std::abs(rhs.zz()));
}

//------------------------------------------------------------------------------
// In place element wise division of the elements of a tensor.
//------------------------------------------------------------------------------
template<typename TensorType>
inline
void
tensorElementWiseDivide(TensorType& lhs, const TensorType& rhs) {
  VERIFY(false);
}

template<>
inline
void
tensorElementWiseDivide<Dim<1>::Tensor>(Dim<1>::Tensor& lhs, const Dim<1>::Tensor& rhs) {
  REQUIRE(rhs.xx() != 0.0);
  lhs.xx(lhs.xx() / rhs.xx());
}

template<>
inline
void
tensorElementWiseDivide<Dim<2>::Tensor>(Dim<2>::Tensor& lhs, const Dim<2>::Tensor& rhs) {
  REQUIRE(std::find_if(rhs.begin(), rhs.end(), std::bind2nd(std::equal_to<double>(), 0.0)) == rhs.end());
  lhs.xx(lhs.xx() / rhs.xx());
  lhs.xy(lhs.xy() / rhs.xy());
  lhs.yx(lhs.yx() / rhs.yx());
  lhs.yy(lhs.yy() / rhs.yy());
}

template<>
inline
void
tensorElementWiseDivide<Dim<3>::Tensor>(Dim<3>::Tensor& lhs, const Dim<3>::Tensor& rhs) {
  REQUIRE(std::find_if(rhs.begin(), rhs.end(), std::bind2nd(std::equal_to<double>(), 0.0)) == rhs.end());
  lhs.xx(lhs.xx() / rhs.xx());
  lhs.xy(lhs.xy() / rhs.xy());
  lhs.xz(lhs.xz() / rhs.xz());
  lhs.yx(lhs.yx() / rhs.yx());
  lhs.yy(lhs.yy() / rhs.yy());
  lhs.yz(lhs.yz() / rhs.yz());
  lhs.zx(lhs.zx() / rhs.zx());
  lhs.zy(lhs.zy() / rhs.zy());
  lhs.zz(lhs.zz() / rhs.zz());
}

}
