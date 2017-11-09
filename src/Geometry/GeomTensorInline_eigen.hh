#include <algorithm>
#include <limits.h>
#include <string>

#include "GeomVector.hh"
#include "GeomSymmetricTensor.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>::
GeomTensor():
  mTensorData(TensorStorage::Zero()) {
}

//------------------------------------------------------------------------------
// Construct with the given values for the elements.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>::
GeomTensor(const double a11):
  mTensorData(TensorStorage::Constant(a11)) {
}

template<>
inline
GeomTensor<2>::
GeomTensor(const double a11, const double a12, 
           const double a21, const double a22):
  mTensorData() {
  mTensorData << a11, a12,
                 a21, a22;
}

template<>
inline
GeomTensor<3>::
GeomTensor(const double a11, const double a12, const double a13,
           const double a21, const double a22, const double a23,
           const double a31, const double a32, const double a33):
  mTensorData() {
  mTensorData << a11, a12, a13,
                 a21, a22, a23,
                 a31, a32, a33;
}

//------------------------------------------------------------------------------
// Override the generic constructors to throw if they're called in the wrong 
// dimensions.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>::
GeomTensor(const double a11, const double a12,
           const double a21, const double a22):
  mTensorData() {
  VERIFY2(false, "GeomTensor(a11, a12, a21, a22): wrong number of dimensions.");
}

template<int nDim>
inline
GeomTensor<nDim>::
GeomTensor(const double a11, const double a12, const double a13,
           const double a21, const double a22, const double a23,
           const double a31, const double a32, const double a33):
  mTensorData() {
  VERIFY2(false, "GeomTensor(a11, a12, a13, a21, a22, a23, a31, a32, a33): wrong number of dimensions.");
}

//------------------------------------------------------------------------------
// Copy constructors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>::
GeomTensor(const GeomTensor<nDim>& ten):
  mTensorData(ten.mTensorData) {
}

template<int nDim>
inline
GeomTensor<nDim>::
GeomTensor(const GeomSymmetricTensor<nDim>& ten):
  mTensorData(ten.native()) {
}

template<int nDim>
inline
GeomTensor<nDim>::
GeomTensor(const TensorStorage& ten):
  mTensorData(ten) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>::~GeomTensor() {
}

//------------------------------------------------------------------------------
// Assignment operators.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::
operator=(const GeomTensor<nDim>& ten) {
  mTensorData = ten.mTensorData;
  return *this;
}

template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::
operator=(const GeomSymmetricTensor<nDim>& ten) {
  mTensorData = ten.native();
  return *this;
}

template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::
operator=(const TensorStorage& ten) {
  mTensorData = ten;
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::operator()(const typename GeomTensor<nDim>::size_type row,
                             const typename GeomTensor<nDim>::size_type column) const {
  REQUIRE(row < nDim);
  REQUIRE(column < nDim);
  return mTensorData(row,column);
}

template<int nDim>
inline
double&
GeomTensor<nDim>::operator()(const typename GeomTensor<nDim>::size_type row,
                             const typename GeomTensor<nDim>::size_type column) {
  REQUIRE(row < nDim);
  REQUIRE(column < nDim);
  return mTensorData(row,column);
}

//------------------------------------------------------------------------------
// Return the (index) element using the bracket operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::operator[](typename GeomTensor<nDim>::size_type index) const {
  REQUIRE(index < numElements);
  return *(begin() + index);
}

template<int nDim>
inline
double&
GeomTensor<nDim>::operator[](typename GeomTensor<nDim>::size_type index) {
  REQUIRE(index < numElements);
  return *(begin() + index);
}

//------------------------------------------------------------------------------
// Return the individual elements, mapped as:
//    xx, xy, xz     11, 12, 13
//    yx, yy, yz  =  21, 22, 23
//    zx, zy, zz     31, 32, 33
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::xx() const {
  return mTensorData(0,0);
}

template<int nDim>
inline
double
GeomTensor<nDim>::xy() const {
  return mTensorData(0,1);
}

template<int nDim>
inline
double
GeomTensor<nDim>::xz() const {
  return mTensorData(0,2);
}

template<int nDim>
inline
double
GeomTensor<nDim>::yx() const {
  return mTensorData(1,0);
}

template<int nDim>
inline
double
GeomTensor<nDim>::yy() const {
  return mTensorData(1,1);
}

template<int nDim>
inline
double
GeomTensor<nDim>::yz() const {
  return mTensorData(1,2);
}

template<int nDim>
inline
double
GeomTensor<nDim>::zx() const {
  return mTensorData(2,0);
}

template<int nDim>
inline
double
GeomTensor<nDim>::zy() const {
  return mTensorData(2,1);
}

template<int nDim>
inline
double
GeomTensor<nDim>::zz() const {
  return mTensorData(2,2);
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> inline double GeomTensor<1>::xy() const { return 0.0; }
template<> inline double GeomTensor<1>::xz() const { return 0.0; }
template<> inline double GeomTensor<1>::yx() const { return 0.0; }
template<> inline double GeomTensor<1>::yy() const { return 0.0; }
template<> inline double GeomTensor<1>::yz() const { return 0.0; }
template<> inline double GeomTensor<1>::zx() const { return 0.0; }
template<> inline double GeomTensor<1>::zy() const { return 0.0; }
template<> inline double GeomTensor<1>::zz() const { return 0.0; }

//------------------------------------------------------------------------------
// 2D dummy elements
template<> inline double GeomTensor<2>::xz() const { return 0.0; }
template<> inline double GeomTensor<2>::yz() const { return 0.0; }
template<> inline double GeomTensor<2>::zx() const { return 0.0; }
template<> inline double GeomTensor<2>::zy() const { return 0.0; }
template<> inline double GeomTensor<2>::zz() const { return 0.0; }

//------------------------------------------------------------------------------
// Set the individual elements, as above.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomTensor<nDim>::xx(const double val) {
  mTensorData(0,0) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::xy(const double val) {
  mTensorData(0,1) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::xz(const double val) {
  mTensorData(0,2) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::yx(const double val) {
  mTensorData(1,0) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::yy(const double val) {
  mTensorData(1,1) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::yz(const double val) {
  mTensorData(1,2) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::zx(const double val) {
  mTensorData(2,0) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::zy(const double val) {
  mTensorData(2,1) = val;
}

template<int nDim>
inline
void
GeomTensor<nDim>::zz(double val) {
  mTensorData(2,2) = val;
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> inline void GeomTensor<1>::xy(const double val) {}
template<> inline void GeomTensor<1>::xz(const double val) {}
template<> inline void GeomTensor<1>::yx(const double val) {}
template<> inline void GeomTensor<1>::yy(const double val) {}
template<> inline void GeomTensor<1>::yz(const double val) {}
template<> inline void GeomTensor<1>::zx(const double val) {}
template<> inline void GeomTensor<1>::zy(const double val) {}
template<> inline void GeomTensor<1>::zz(const double val) {}

//------------------------------------------------------------------------------
// 2D dummy elements
template<> inline void GeomTensor<2>::xz(const double val) {}
template<> inline void GeomTensor<2>::yz(const double val) {}
template<> inline void GeomTensor<2>::zx(const double val) {}
template<> inline void GeomTensor<2>::zy(const double val) {}
template<> inline void GeomTensor<2>::zz(const double val) {}

//------------------------------------------------------------------------------
// Access the individual rows of the GeomTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomTensor<nDim>::getRow(const GeomTensor<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return GeomVector<nDim>((mTensorData.row(index).transpose()).eval());
}

//------------------------------------------------------------------------------
// Access the individual columns of the GeomTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomTensor<nDim>::getColumn(const GeomTensor<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return GeomVector<nDim>((mTensorData.col(index)).eval());
}

//------------------------------------------------------------------------------
// Set a row of the GeomTensor to a GeomVector.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomTensor<nDim>::setRow(const GeomTensor<nDim>::size_type index,
                         const GeomVector<nDim>& vec) {
  REQUIRE(index < nDim);
  mTensorData.row(index) = vec.native().tranpose();
}

//------------------------------------------------------------------------------
// Set a column of the GeomTensor to a GeomVector.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomTensor<nDim>::setColumn(const GeomTensor<nDim>::size_type index,
                            const GeomVector<nDim>& vec) {
  REQUIRE(index < nDim);
  mTensorData.col(index) = vec.native();
}

//------------------------------------------------------------------------------
// Iterators to the raw data.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomTensor<nDim>::iterator
GeomTensor<nDim>::begin() {
  return mTensorData.data();
}

template<int nDim>
inline
typename GeomTensor<nDim>::iterator
GeomTensor<nDim>::end() {
  return mTensorData.data() + nDim*nDim;
}

template<int nDim>
inline
typename GeomTensor<nDim>::const_iterator
GeomTensor<nDim>::begin() const{
  return mTensorData.data();
}

template<int nDim>
inline
typename GeomTensor<nDim>::const_iterator
GeomTensor<nDim>::end() const {
  return mTensorData.data() + nDim*nDim;
}

//------------------------------------------------------------------------------
// Zero out the GeomTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomTensor<nDim>::Zero() {
  mTensorData = TensorStorage::Zero();
}

//------------------------------------------------------------------------------
// Force the tensor to be the identity tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomTensor<nDim>::Identity() {
  mTensorData = TensorStorage::Identity();
}

//------------------------------------------------------------------------------
// Return the negative of a tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::operator-() const {
  return GeomTensor<nDim>(-mTensorData);
}

//------------------------------------------------------------------------------
// Add two tensors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator+(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData + rhs.mTensorData).eval());
}

template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator+(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData + rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Subtract a tensor from another.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator-(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData - rhs.mTensorData).eval());
}

template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator-(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData - rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply two tensors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator*(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.mTensorData).eval());
}

template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
operator*(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomTensor<nDim>::operator*(const GeomVector<nDim>& rhs) const {
  return GeomVector<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply a tensor by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::operator*(const double rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs).eval());
}

//------------------------------------------------------------------------------
// Divide a tensor by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  return GeomTensor<nDim>((mTensorData / rhs).eval());
}

//------------------------------------------------------------------------------
// Add two tensors in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator+=(const GeomTensor<nDim>& rhs) {
  mTensorData += rhs.mTensorData;
  return *this;
}

template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator+=(const GeomSymmetricTensor<nDim>& rhs) {
  mTensorData += rhs.native();
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a tensor from this one in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator-=(const GeomTensor<nDim>& rhs) {
  mTensorData -= rhs.mTensorData;
  return *this;
}

template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator-=(const GeomSymmetricTensor<nDim>& rhs) {
  mTensorData -= rhs.native();
  return *this;
}

//------------------------------------------------------------------------------
// Multiply by a tensor in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator*=(const GeomTensor<nDim>& rhs) {
  mTensorData *= rhs.mTensorData;
  return *this;
}

template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator*=(const GeomSymmetricTensor<nDim>& rhs) {
  mTensorData *= rhs.native();
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this tensor by a scalar in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator*=(const double rhs) {
  mTensorData *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this tensor by a scalar in place
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>&
GeomTensor<nDim>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  mTensorData /= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Define the equivalence operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomTensor<nDim>::
operator==(const GeomTensor<nDim>& rhs) const {
  return mTensorData == rhs.mTensorData;
}

template<int nDim>
inline
bool
GeomTensor<nDim>::
operator==(const GeomSymmetricTensor<nDim>& rhs) const {
  return mTensorData == rhs.native();
}

//------------------------------------------------------------------------------
// Define the not equivalence than comparitor.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomTensor<nDim>::
operator!=(const GeomTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

template<int nDim>
inline
bool
GeomTensor<nDim>::
operator!=(const GeomSymmetricTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// Define the less than operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomTensor<nDim>::
operator<(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

template<int nDim>
inline
bool
GeomTensor<nDim>::
operator<(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

// template<>
// inline
// bool
// GeomTensor<1>::
// operator<(const double rhs) const {
//   return this->mxx < rhs;
// }

// template<>
// inline
// bool
// GeomTensor<2>::
// operator<(const double rhs) const {
//   return (this->mxx < rhs and
//           this->mxy < rhs and
//           this->myx < rhs and
//           this->myy < rhs);
// }

// template<>
// inline
// bool
// GeomTensor<3>::
// operator<(const double rhs) const {
//   return (this->mxx < rhs and
//           this->mxy < rhs and
//           this->mxz < rhs and
//           this->myx < rhs and
//           this->myy < rhs and
//           this->myz < rhs and
//           this->mzx < rhs and
//           this->mzy < rhs and
//           this->mzz < rhs);
// }

//------------------------------------------------------------------------------
// Define the greater than operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomTensor<nDim>::
operator>(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

template<int nDim>
inline
bool
GeomTensor<nDim>::
operator>(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

// template<int nDim>
// inline
// bool
// GeomTensor<nDim>::
// operator>(const double rhs) const {
//   return this->mxx > rhs;
// }

// template<>
// inline
// bool
// GeomTensor<2>::
// operator>(const double rhs) const {
//   return (this->mxx > rhs and
//           this->mxy > rhs and
//           this->myx > rhs and
//           this->myy > rhs);
// }

// template<>
// inline
// bool
// GeomTensor<3>::
// operator>(const double rhs) const {
//   return (this->mxx > rhs and
//           this->mxy > rhs and
//           this->mxz > rhs and
//           this->myx > rhs and
//           this->myy > rhs and
//           this->myz > rhs and
//           this->mzx > rhs and
//           this->mzy > rhs and
//           this->mzz > rhs);
// }

//------------------------------------------------------------------------------
// Define the less than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomTensor<nDim>::
operator<=(const GeomTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

template<int nDim>
inline
bool
GeomTensor<nDim>::
operator<=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

// template<>
// inline
// bool
// GeomTensor<1>::
// operator<=(const double rhs) const {
//   return this->mxx <= rhs;
// }

// template<>
// inline
// bool
// GeomTensor<2>::
// operator<=(const double rhs) const {
//   return (this->mxx <= rhs and
//           this->mxy <= rhs and
//           this->myx <= rhs and
//           this->myy <= rhs);
// }

// template<>
// inline
// bool
// GeomTensor<3>::
// operator<=(const double rhs) const {
//   return (this->mxx <= rhs and
//           this->mxy <= rhs and
//           this->mxz <= rhs and
//           this->myx <= rhs and
//           this->myy <= rhs and
//           this->myz <= rhs and
//           this->mzx <= rhs and
//           this->mzy <= rhs and
//           this->mzz <= rhs);
// }

//------------------------------------------------------------------------------
// Define the greater than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomTensor<nDim>::
operator>=(const GeomTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

template<int nDim>
inline
bool
GeomTensor<nDim>::
operator>=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

// template<int nDim>
// inline
// bool
// GeomTensor<nDim>::
// operator>=(const double rhs) const {
//   return this->mxx >= rhs;
// }

// template<>
// inline
// bool
// GeomTensor<2>::
// operator>=(const double rhs) const {
//   return (this->mxx >= rhs and
//           this->mxy >= rhs and
//           this->myx >= rhs and
//           this->myy >= rhs);
// }

// template<>
// inline
// bool
// GeomTensor<3>::
// operator>=(const double rhs) const {
//   return (this->mxx >= rhs and
//           this->mxy >= rhs and
//           this->mxz >= rhs and
//           this->myx >= rhs and
//           this->myy >= rhs and
//           this->myz >= rhs and
//           this->mzx >= rhs and
//           this->mzy >= rhs and
//           this->mzz >= rhs);
// }

//------------------------------------------------------------------------------
// Return the symmetric part of a GeomTensor.
//   Bij = 0.5*(Aij + Aji)
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomTensor<nDim>::Symmetric() const {
  return GeomSymmetricTensor<nDim>((0.5*(mTensorData + mTensorData.transpose())).eval());
}

//------------------------------------------------------------------------------
// Return the skew-symmetric part of a GeomTensor.
//   Bij = 0.5*(Aij - Aji)
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::SkewSymmetric() const {
  return GeomTensor<nDim>((0.5*(mTensorData - mTensorData.transpose())).eval());
}

//------------------------------------------------------------------------------
// Return the transpose of the GeomTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
Transpose() const {
  return GeomTensor<nDim>(mTensorData.transpose().eval());
}

//------------------------------------------------------------------------------
// Return the inverse of the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::Inverse() const {
  return GeomTensor<nDim>(mTensorData.inverse().eval());
}

// template<>
// inline
// GeomTensor<1>
// GeomTensor<1>::Inverse() const {
//   REQUIRE(mTensorData(0,0) != 0.0);
//   return GeomTensor<1>(1.0/mTensorData(0,0));
// }

// template<>
// inline
// GeomTensor<2>
// GeomTensor<2>::Inverse() const {
//   return GeomTensor<2>((mTensorData(1,1)), -(mTensorData(0,1)),
//                       -(mTensor(1,0)),  (mTensorData(0,0)))/this->Determinant();
// }

// template<>
// inline
// GeomTensor<3>
// GeomTensor<3>::Inverse() const {
//   REQUIRE(Determinant() != 0.0);
//   return GeomTensor<3>((mTensorData(1,1))*(this->mzz) - (mTensorData(1,2))*(this->mzy), (mTensorData(0,2))*(this->mzy) - (mTensorData(0,1))*(this->mzz), (mTensorData(0,1))*(mTensorData(1,2)) - (mTensorData(0,2))*(mTensorData(1,1)),
//                        (mTensorData(1,2))*(this->mzx) - (mTensor(1,0))*(this->mzz), (mTensorData(0,0))*(this->mzz) - (mTensorData(0,2))*(this->mzx), (mTensorData(0,2))*(mTensor(1,0)) - (mTensorData(0,0))*(mTensorData(1,2)),
//                        (mTensor(1,0))*(this->mzy) - (mTensorData(1,1))*(this->mzx), (mTensorData(0,1))*(this->mzx) - (mTensorData(0,0))*(this->mzy), (mTensorData(0,0))*(mTensorData(1,1)) - (mTensorData(0,1))*(mTensor(1,0)))/this->Determinant();
// }

//------------------------------------------------------------------------------
// Return the diagonal elements of the GeomTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomTensor<nDim>::diagonalElements() const {
  return GeomVector<nDim>(mTensorData.diagonal().eval());
}

//------------------------------------------------------------------------------
// Return the trace of the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::Trace() const {
  return mTensorData.trace();
}

//------------------------------------------------------------------------------
// Return the determinant of the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::Determinant() const {
  return mTensorData.determinant();
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomTensor<nDim>::dot(const GeomVector<nDim>& rhs) const {
  return GeomVector<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply two tensors.  This is just the linear algebra definition for matrix
// multiplication.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::dot(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.mTensorData).eval());
}

template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::dot(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Return the doubledot product.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::
doubledot(const GeomTensor<nDim>& rhs) const {
  return (mTensorData*rhs.mTensorData).trace();
}

template<int nDim>
inline
double
GeomTensor<nDim>::
doubledot(const GeomSymmetricTensor<nDim>& rhs) const {
  return (mTensorData*rhs.native()).trace();
}

//------------------------------------------------------------------------------
// Return the doubledot product with ourself.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::
selfDoubledot() const {
  return (mTensorData*mTensorData).trace();
}

//------------------------------------------------------------------------------
// Return the square of this tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
square() const {
  return GeomTensor<nDim>((mTensorData * mTensorData).eval());
}

//------------------------------------------------------------------------------
// Return a new tensor with the elements of this tensor squared.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomTensor<nDim>::
squareElements() const {
  return GeomTensor<nDim>(mTensorData.array().square().eval());
}

//------------------------------------------------------------------------------
// eigen values
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomTensor<nDim>::
eigenValues() const {
  Eigen::EigenSolver<TensorStorage> eigensolver(mTensorData);
  CHECK(eigensolver.info() == Eigen::Success);
  GeomVector<nDim> result;
  for (unsigned i = 0; i < nDim; ++i) result(i) = std::real(eigensolver.eigenvalues()(i));
  return result;
}

//------------------------------------------------------------------------------
// Apply a rotational transform to this tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomTensor<nDim>::
rotationalTransform(const GeomTensor<nDim>& R) {
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-8));
  mTensorData = R.mTensorData * mTensorData * R.mTensorData.transpose();
}

//------------------------------------------------------------------------------
// Return the maximum absolute value of the elements.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomTensor<nDim>::
maxAbsElement() const {
  return mTensorData.cwiseAbs().maxCoeff();
}

//------------------------------------------------------------------------------
// Access the native Eigen type.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomTensor<nDim>::TensorStorage&
GeomTensor<nDim>::native() {
  return mTensorData;
}

template<int nDim>
inline
const typename GeomTensor<nDim>::TensorStorage&
GeomTensor<nDim>::native() const {
  return mTensorData;
}

//********************************************************************************
// Global Functions.
//********************************************************************************

//------------------------------------------------------------------------------
// Multiplication by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
operator*(double lhs, const GeomTensor<nDim>& rhs) {
  return rhs*lhs;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, GeomTensor<nDim>& ten) {
  std::string parenthesis;
  is >> parenthesis;
  for (typename GeomTensor<nDim>::iterator elementItr = ten.begin();
       elementItr != ten.end();
       ++elementItr) {
    is >> *elementItr;
  }
  is >> parenthesis;
  return is;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::ostream&
operator<<(std::ostream& os, const GeomTensor<nDim>& ten) {
  os << "( ";
  for (typename GeomTensor<nDim>::const_iterator itr = ten.begin();
       itr != ten.end(); ++itr) {
    os << *itr << " ";
  }
  os << ")";
  return os;
}

}

