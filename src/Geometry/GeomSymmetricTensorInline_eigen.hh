#include <algorithm>
#include <string>
#include <limits>
#include "GeomVector.hh"
#include "GeomTensor.hh"
#include "EigenStruct.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor():
  mTensorData(TensorStorage::Zero()) {
}

//------------------------------------------------------------------------------
// Construct with the given values for the elements.
// Unfortunately it does not appear Eigen allows us to specify try symmetric
// matrices, so we have to burn the memory for the symmetric components.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const double a11):
  mTensorData(TensorStorage::Constant(a11)) {
}

template<>
inline
GeomSymmetricTensor<2>::
GeomSymmetricTensor(const double a11, const double a12, 
                    const double a21, const double a22):
  mTensorData() {
  REQUIRE(a12 == a21);
  mTensorData << a11, a12, 
                 a21, a22;
}

template<>
inline
GeomSymmetricTensor<3>::
GeomSymmetricTensor(const double a11, const double a12, const double a13,
                    const double a21, const double a22, const double a23,
                    const double a31, const double a32, const double a33):
  mTensorData() {
  REQUIRE(a12 == a21);
  REQUIRE(a13 == a31);
  REQUIRE(a23 == a32);
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
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const double a11, const double a12,
                    const double a21, const double a22):
  mTensorData() {
  VERIFY2(false, "GeomSymmetricTensor(a11, a12, a21, a22): wrong number of dimensions.");
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const double a11, const double a12, const double a13,
                    const double a21, const double a22, const double a23,
                    const double a31, const double a32, const double a33):
  mTensorData() {
  VERIFY2(false, "GeomSymmetricTensor(a11, a12, a13, a21, a22, a23, a31, a32, a33): wrong number of dimensions.");
}

//------------------------------------------------------------------------------
// Copy constructors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const GeomSymmetricTensor<nDim>& ten):
  mTensorData(ten.mTensorData) {
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const GeomTensor<nDim>& ten):
  mTensorData(0.5*(ten.native() + ten.native().transpose())) {
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const TensorStorage& ten):
  mTensorData(0.5*(ten + ten.transpose())) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>::~GeomSymmetricTensor() {
}

//------------------------------------------------------------------------------
// Assignment operators.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::
operator=(const GeomTensor<nDim>& rhs) {
  mTensorData = 0.5*(rhs.native() + rhs.native().transpose());
  return *this;
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::
operator=(const GeomSymmetricTensor<nDim>& rhs) {
  mTensorData = rhs.mTensorData;
  return *this;
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::
operator=(const TensorStorage& rhs) {
  mTensorData = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::operator()(const typename GeomSymmetricTensor<nDim>::size_type row,
                                      const typename GeomSymmetricTensor<nDim>::size_type column) const {
  REQUIRE(row < nDim);
  REQUIRE(column < nDim);
  return mTensorData(row,column);
}

template<int nDim>
inline
double&
GeomSymmetricTensor<nDim>::operator()(const typename GeomSymmetricTensor<nDim>::size_type row,
                             const typename GeomSymmetricTensor<nDim>::size_type column) {
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
GeomSymmetricTensor<nDim>::operator[](typename GeomSymmetricTensor<nDim>::size_type index) const {
  REQUIRE(index < numElements);
  return *(begin() + index);
}

template<int nDim>
inline
double&
GeomSymmetricTensor<nDim>::operator[](typename GeomSymmetricTensor<nDim>::size_type index) {
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
GeomSymmetricTensor<nDim>::xx() const {
  return mTensorData(0,0);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::xy() const {
  return mTensorData(0,1);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::xz() const {
  return mTensorData(0,2);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::yx() const {
  return mTensorData(1,0);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::yy() const {
  return mTensorData(1,1);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::yz() const {
  return mTensorData(1,2);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::zx() const {
  return mTensorData(2,0);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::zy() const {
  return mTensorData(2,1);
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::zz() const {
  return mTensorData(2,2);
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> inline double GeomSymmetricTensor<1>::xy() const { return 0.0; }
template<> inline double GeomSymmetricTensor<1>::xz() const { return 0.0; }
template<> inline double GeomSymmetricTensor<1>::yx() const { return 0.0; }
template<> inline double GeomSymmetricTensor<1>::yy() const { return 0.0; }
template<> inline double GeomSymmetricTensor<1>::yz() const { return 0.0; }
template<> inline double GeomSymmetricTensor<1>::zx() const { return 0.0; }
template<> inline double GeomSymmetricTensor<1>::zy() const { return 0.0; }
template<> inline double GeomSymmetricTensor<1>::zz() const { return 0.0; }

//------------------------------------------------------------------------------
// 2D dummy elements
template<> inline double GeomSymmetricTensor<2>::xz() const { return 0.0; }
template<> inline double GeomSymmetricTensor<2>::yz() const { return 0.0; }
template<> inline double GeomSymmetricTensor<2>::zx() const { return 0.0; }
template<> inline double GeomSymmetricTensor<2>::zy() const { return 0.0; }
template<> inline double GeomSymmetricTensor<2>::zz() const { return 0.0; }

//------------------------------------------------------------------------------
// Set the individual elements, as above.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::xx(const double val) {
  mTensorData(0,0) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::xy(const double val) {
  mTensorData(0,1) = val;
  mTensorData(1,0) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::xz(const double val) {
  mTensorData(0,2) = val;
  mTensorData(2,0) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::yx(const double val) {
  mTensorData(1,0) = val;
  mTensorData(0,1) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::yy(const double val) {
  mTensorData(1,1) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::yz(const double val) {
  mTensorData(1,2) = val;
  mTensorData(2,1) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::zx(const double val) {
  mTensorData(2,0) = val;
  mTensorData(0,2) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::zy(const double val) {
  mTensorData(2,1) = val;
  mTensorData(1,2) = val;
}

template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::zz(double val) {
  mTensorData(2,2) = val;
}

//------------------------------------------------------------------------------
// 1D dummy elements
template<> inline void GeomSymmetricTensor<1>::xy(const double val) {}
template<> inline void GeomSymmetricTensor<1>::xz(const double val) {}
template<> inline void GeomSymmetricTensor<1>::yx(const double val) {}
template<> inline void GeomSymmetricTensor<1>::yy(const double val) {}
template<> inline void GeomSymmetricTensor<1>::yz(const double val) {}
template<> inline void GeomSymmetricTensor<1>::zx(const double val) {}
template<> inline void GeomSymmetricTensor<1>::zy(const double val) {}
template<> inline void GeomSymmetricTensor<1>::zz(const double val) {}

//------------------------------------------------------------------------------
// 2D dummy elements
template<> inline void GeomSymmetricTensor<2>::xz(const double val) {}
template<> inline void GeomSymmetricTensor<2>::yz(const double val) {}
template<> inline void GeomSymmetricTensor<2>::zx(const double val) {}
template<> inline void GeomSymmetricTensor<2>::zy(const double val) {}
template<> inline void GeomSymmetricTensor<2>::zz(const double val) {}

//------------------------------------------------------------------------------
// Access the individual rows of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::getRow(const GeomSymmetricTensor<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return GeomVector<nDim>((mTensorData.row(index).transpose()).eval());
}

//------------------------------------------------------------------------------
// Access the individual columns of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::getColumn(const GeomSymmetricTensor<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return GeomVector<nDim>((mTensorData.col(index)).eval());
}

//------------------------------------------------------------------------------
// Set a row of the GeomSymmetricTensor to a GeomVector.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::setRow(const GeomSymmetricTensor<nDim>::size_type index,
                         const GeomVector<nDim>& vec) {
  REQUIRE(index < nDim);
  mTensorData.row(index) = vec.native().transpose();
  mTensorData.col(index) = vec.native();
}

//------------------------------------------------------------------------------
// Set a column of the GeomSymmetricTensor to a GeomVector.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::setColumn(const GeomSymmetricTensor<nDim>::size_type index,
                            const GeomVector<nDim>& vec) {
  REQUIRE(index < nDim);
  mTensorData.col(index) = vec.native();
  mTensorData.row(index) = vec.native().transpose();
}

//------------------------------------------------------------------------------
// Iterators to the raw data.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::iterator
GeomSymmetricTensor<nDim>::begin() {
  return mTensorData.data();
}

template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::iterator
GeomSymmetricTensor<nDim>::end() {
  return mTensorData.data() + nDim*nDim;
}

template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::const_iterator
GeomSymmetricTensor<nDim>::begin() const{
  return mTensorData.data();
}

template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::const_iterator
GeomSymmetricTensor<nDim>::end() const {
  return mTensorData.data() + nDim*nDim;
}

//------------------------------------------------------------------------------
// Zero out the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::Zero() {
  mTensorData = TensorStorage::Zero();
}

//------------------------------------------------------------------------------
// Force the tensor to be the identity tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::Identity() {
  mTensorData = TensorStorage::Identity();
}

//------------------------------------------------------------------------------
// Return the negative of a tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::operator-() const {
  return GeomSymmetricTensor<nDim>((-mTensorData).eval());
}

//------------------------------------------------------------------------------
// Add two tensors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator+(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData + rhs.native()).eval());
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
operator+(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomSymmetricTensor<nDim>((mTensorData + rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Subtract a tensor from another.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator-(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData - rhs.native()).eval());
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
operator-(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomSymmetricTensor<nDim>((mTensorData - rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply two tensors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator*(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.native()).eval());
}

template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator*(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::operator*(const GeomVector<nDim>& rhs) const {
  return GeomVector<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply a tensor by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::operator*(const double rhs) const {
  return GeomSymmetricTensor<nDim>((mTensorData * rhs).eval());
}

//------------------------------------------------------------------------------
// Divide a tensor by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  return GeomSymmetricTensor<nDim>((mTensorData / rhs).eval());
}

//------------------------------------------------------------------------------
// Add two tensors in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::operator+=(const GeomSymmetricTensor<nDim>& rhs) {
  mTensorData += rhs.native();
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a tensor from this one in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::operator-=(const GeomSymmetricTensor<nDim>& rhs) {
  mTensorData -= rhs.native();
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this tensor by a scalar in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::operator*=(const double rhs) {
  mTensorData *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this tensor by a scalar in place
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::operator/=(const double rhs) {
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
GeomSymmetricTensor<nDim>::
operator==(const GeomTensor<nDim>& rhs) const {
  return mTensorData == rhs.native();
}

template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator==(const GeomSymmetricTensor<nDim>& rhs) const {
  return mTensorData == rhs.native();
}

//------------------------------------------------------------------------------
// Define the not equivalence than comparitor.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator!=(const GeomTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator!=(const GeomSymmetricTensor<nDim>& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// Define the less than operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator<(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator<(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() < rhs.Determinant();
}

//------------------------------------------------------------------------------
// Define the greater than operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator>(const GeomTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator>(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->Determinant() > rhs.Determinant();
}

//------------------------------------------------------------------------------
// Define the less than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator<=(const GeomTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator<=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this < rhs) or (*this == rhs);
}

//------------------------------------------------------------------------------
// Define the greater than or equal operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator>=(const GeomTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator>=(const GeomSymmetricTensor<nDim>& rhs) const {
  return (*this > rhs) or (*this == rhs);
}

//------------------------------------------------------------------------------
// Return the symmetric part.  A no-op for this class.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::Symmetric() const {
  return *this;
}

//------------------------------------------------------------------------------
// Return the skew-symmetric part of a GeomSymmetricTensor.
//   Bij = 0.5*(Aij - Aji)
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::SkewSymmetric() const {
  return GeomTensor<nDim>::zero;
}

//------------------------------------------------------------------------------
// Return the transpose of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
Transpose() const {
  return *this;
}

//------------------------------------------------------------------------------
// Return the inverse of the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::Inverse() const {
  return GeomSymmetricTensor<nDim>(mTensorData.inverse().eval());
}

//------------------------------------------------------------------------------
// Return the diagonal elements of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::diagonalElements() const {
  return GeomVector<nDim>(mTensorData.diagonal().eval());
}

//------------------------------------------------------------------------------
// Return the trace of the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::Trace() const {
  return mTensorData.trace();
}

//------------------------------------------------------------------------------
// Return the determinant of the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::Determinant() const {
  return mTensorData.determinant();
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::dot(const GeomVector<nDim>& rhs) const {
  return GeomVector<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Multiply two tensors.  This is just the linear algebra definition for matrix
// multiplication.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::dot(const GeomTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.native()).eval());
}

template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::dot(const GeomSymmetricTensor<nDim>& rhs) const {
  return GeomTensor<nDim>((mTensorData * rhs.native()).eval());
}

//------------------------------------------------------------------------------
// Return the doubledot product.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::
doubledot(const GeomTensor<nDim>& rhs) const {
  return (mTensorData*rhs.native()).trace();
}

template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::
doubledot(const GeomSymmetricTensor<nDim>& rhs) const {
  return (mTensorData*rhs.native()).trace();
}

//------------------------------------------------------------------------------
// Return the doubledot product with ourself.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::
selfDoubledot() const {
  return (mTensorData*mTensorData).trace();
}

//------------------------------------------------------------------------------
// Return the square of this tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
square() const {
  return GeomSymmetricTensor<nDim>((mTensorData * mTensorData).eval());
}

//------------------------------------------------------------------------------
// Return the cube of this tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
cube() const {
  return GeomSymmetricTensor<nDim>((mTensorData * mTensorData * mTensorData).eval());
}

//------------------------------------------------------------------------------
// Return the sqrt of this tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
sqrt() const {
  Eigen::SelfAdjointEigenSolver<TensorStorage> eigensolver(mTensorData);
  CHECK(eigensolver.info() == Eigen::Success);
  const TensorStorage vecs = eigensolver.eigenvectors();
  const Eigen::Matrix<double, nDim, 1> vals = eigensolver.eigenvalues().cwiseSqrt();
  return GeomSymmetricTensor<nDim>((vecs*vals.asDiagonal()*vecs.transpose()).eval());
}

//------------------------------------------------------------------------------
// Compute the cube root of the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
cuberoot() const {
  Eigen::SelfAdjointEigenSolver<TensorStorage> eigensolver(mTensorData);
  CHECK(eigensolver.info() == Eigen::Success);
  const TensorStorage vecs = eigensolver.eigenvectors();
  Eigen::Matrix<double, nDim, 1> vals = eigensolver.eigenvalues();
  for (size_t i = 0; i != nDim; ++i) vals(i) = FastMath::CubeRootHalley2(vals(i));
  // const Eigen::Matrix<double, nDim, 1> vals = eigensolver.eigenvalues().array().pow(1.0/3.0);
  return GeomSymmetricTensor<nDim>((vecs*vals.asDiagonal()*vecs.transpose()).eval());
}

//------------------------------------------------------------------------------
// The general version, raise to an arbitrary power.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
pow(const double p) const {
  Eigen::SelfAdjointEigenSolver<TensorStorage> eigensolver(mTensorData);
  CHECK(eigensolver.info() == Eigen::Success);
  const TensorStorage vecs = eigensolver.eigenvectors();
  const Eigen::Matrix<double, nDim, 1> vals = eigensolver.eigenvalues().array().pow(p);
  return GeomSymmetricTensor<nDim>((vecs*vals.asDiagonal()*vecs.transpose()).eval());
}

//------------------------------------------------------------------------------
// Return a new tensor with the elements of this tensor squared.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
squareElements() const {
  return GeomSymmetricTensor<nDim>(mTensorData.array().square().eval());
}

//------------------------------------------------------------------------------
// Apply a rotational transform to this tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::
rotationalTransform(const GeomTensor<nDim>& R) {
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-8));
  mTensorData = R.native() * mTensorData * R.native().transpose();
}

//------------------------------------------------------------------------------
// Return the maximum absolute value of the elements.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::
maxAbsElement() const {
  return mTensorData.cwiseAbs().maxCoeff();
}

//------------------------------------------------------------------------------
// eigen values
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::
eigenValues() const {
  Eigen::SelfAdjointEigenSolver<TensorStorage> eigensolver(mTensorData);
  CHECK(eigensolver.info() == Eigen::Success);
  return GeomVector<nDim>(eigensolver.eigenvalues().eval());
}

//------------------------------------------------------------------------------
// eigen vectors
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::EigenStructType
GeomSymmetricTensor<nDim>::
eigenVectors() const {
  Eigen::SelfAdjointEigenSolver<TensorStorage> eigensolver(mTensorData);
  CHECK(eigensolver.info() == Eigen::Success);
  EigenStructType result;
  result.eigenValues = eigensolver.eigenvalues();
  result.eigenVectors = eigensolver.eigenvectors();
  return result;
}

//------------------------------------------------------------------------------
// Access the native Eigen type.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::TensorStorage&
GeomSymmetricTensor<nDim>::native() {
  return mTensorData;
}

template<int nDim>
inline
const typename GeomSymmetricTensor<nDim>::TensorStorage&
GeomSymmetricTensor<nDim>::native() const {
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
Spheral::GeomSymmetricTensor<nDim>
operator*(double lhs, const Spheral::GeomSymmetricTensor<nDim>& rhs) {
  return rhs*lhs;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, Spheral::GeomSymmetricTensor<nDim>& ten) {
  std::string parenthesis;
  is >> parenthesis;
  for (typename Spheral::GeomSymmetricTensor<nDim>::iterator elementItr = ten.begin();
       elementItr < ten.end();
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
operator<<(std::ostream& os, const Spheral::GeomSymmetricTensor<nDim>& ten) {
  os << "( ";
  for (typename Spheral::GeomSymmetricTensor<nDim>::const_iterator itr = ten.begin();
       itr < ten.end(); ++itr) {
    os << *itr << " ";
  }
  os << ")";
  return os;
}

//------------------------------------------------------------------------------
// Min a symmetric tensor with a scalar -- limit the eigenvalues.
//------------------------------------------------------------------------------
template<int nDim>
inline
Spheral::GeomSymmetricTensor<nDim>
min(const double minValue, const Spheral::GeomSymmetricTensor<nDim>& tensor) {
  typedef typename Spheral::GeomSymmetricTensor<nDim>::TensorStorage TensorStorage;
  Eigen::SelfAdjointEigenSolver<TensorStorage> eigensolver(tensor.native());
  CHECK(eigensolver.info() == Eigen::Success);
  const TensorStorage vecs = eigensolver.eigenvectors();
  const Eigen::Matrix<double, nDim, 1> vals = eigensolver.eigenvalues().array().min(minValue);
  return GeomSymmetricTensor<nDim>(vecs*vals.asDiagonal()*vecs.transpose());
}

template<int nDim>
inline
Spheral::GeomSymmetricTensor<nDim>
min(const Spheral::GeomSymmetricTensor<nDim>& tensor, const double minValue) {
  return min(minValue, tensor);
}

//------------------------------------------------------------------------------
// Max a symmetric tensor with a scalar -- limit the eigenvalues.
//------------------------------------------------------------------------------
template<int nDim>
inline
Spheral::GeomSymmetricTensor<nDim>
max(const double maxValue, const Spheral::GeomSymmetricTensor<nDim>& tensor) {
  typedef typename Spheral::GeomSymmetricTensor<nDim>::TensorStorage TensorStorage;
  Eigen::SelfAdjointEigenSolver<TensorStorage> eigensolver(tensor.native());
  CHECK(eigensolver.info() == Eigen::Success);
  const TensorStorage vecs = eigensolver.eigenvectors();
  const Eigen::Matrix<double, nDim, 1> vals = eigensolver.eigenvalues().array().max(maxValue);
  return GeomSymmetricTensor<nDim>(vecs*vals.asDiagonal()*vecs.transpose());
}

template<int nDim>
inline
Spheral::GeomSymmetricTensor<nDim>
max(const Spheral::GeomSymmetricTensor<nDim>& tensor, double const maxValue) {
  return max(maxValue, tensor);
}

}
