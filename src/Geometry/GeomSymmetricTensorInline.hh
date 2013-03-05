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
// Return the element index corresponding to the given (row,column)
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>::size_type
GeomSymmetricTensor<1>::elementIndex(const GeomSymmetricTensor<1>::size_type row,
                                     const GeomSymmetricTensor<1>::size_type column) const {
  REQUIRE(row < 1);
  REQUIRE(column < 1);
  return 0;
}

template<>
inline
GeomSymmetricTensor<2>::size_type
GeomSymmetricTensor<2>::elementIndex(const GeomSymmetricTensor<2>::size_type row,
                                     const GeomSymmetricTensor<2>::size_type column) const {
  REQUIRE(row < 2);
  REQUIRE(column < 2);
  int i = std::min(row, column);
  int j = std::max(row, column);
  int result = (5 - i)*i/2 + j - i;
  ENSURE(result >= 0 and result < mNumElements);
  return result;
}

template<>
inline
GeomSymmetricTensor<3>::size_type
GeomSymmetricTensor<3>::elementIndex(const GeomSymmetricTensor<3>::size_type row,
                            const GeomSymmetricTensor<3>::size_type column) const {
  REQUIRE(row < 3);
  REQUIRE(column < 3);
  int i = std::min(row, column);
  int j = std::max(row, column);
  int result = (7 - i)*i/2 + j - i;
  ENSURE(result >= 0 and result < mNumElements);
  return result;
}

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor():
  GeomSymmetricTensorBase<nDim>(0.0) {
}

//------------------------------------------------------------------------------
// Construct with the given values for the elements.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const double a11):
  GeomSymmetricTensorBase<nDim>(a11) {
}

template<>
inline
GeomSymmetricTensor<2>::
GeomSymmetricTensor(const double a11, const double a12, 
                    const double a21, const double a22):
  GeomSymmetricTensorBase<2>(a11, a12,
                                  a22) {
  REQUIRE(a12 == a21);
}

template<>
inline
GeomSymmetricTensor<3>::
GeomSymmetricTensor(const double a11, const double a12, const double a13,
                    const double a21, const double a22, const double a23,
                    const double a31, const double a32, const double a33):
  GeomSymmetricTensorBase<3>(a11, a12, a13,
                                  a22, a23,
                                       a33) {
  REQUIRE(a12 == a21);
  REQUIRE(a13 == a31);
  REQUIRE(a23 == a32);
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
  GeomSymmetricTensorBase<nDim>(0.0) {
  VERIFY2(false, "GeomSymmetricTensor(a11, a12, a21, a22): wrong number of dimensions.");
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const double a11, const double a12, const double a13,
                    const double a21, const double a22, const double a23,
                    const double a31, const double a32, const double a33):
  GeomSymmetricTensorBase<nDim>(0.0) {
  VERIFY2(false, "GeomSymmetricTensor(a11, a12, a13, a21, a22, a23, a31, a32, a33): wrong number of dimensions.");
}

//------------------------------------------------------------------------------
// Copy constructors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>::
GeomSymmetricTensor(const GeomSymmetricTensor<nDim>& ten):
  GeomSymmetricTensorBase<nDim>(ten) {
}

template<>
inline
GeomSymmetricTensor<1>::
GeomSymmetricTensor(const GeomTensor<1>& ten):
  GeomSymmetricTensorBase<1>(ten.xx()) {
}

template<>
inline
GeomSymmetricTensor<2>::
GeomSymmetricTensor(const GeomTensor<2>& ten):
  GeomSymmetricTensorBase<2>(ten.xx(), 0.5*(ten.xy() + ten.yx()),
                                       ten.yy()) {
                             
}

template<>
inline
GeomSymmetricTensor<3>::
GeomSymmetricTensor(const GeomTensor<3>& ten):
  GeomSymmetricTensorBase<3>(ten.xx(), 0.5*(ten.xy() + ten.yx()), 0.5*(ten.xz() + ten.zx()),
                                       ten.yy(),                  0.5*(ten.yz() + ten.zy()),
                                                                  ten.zz()) {
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
template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::
operator=(const GeomTensor<1>& rhs) {
  this->mxx = rhs.xx();
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::
operator=(const GeomTensor<2>& rhs) {
  REQUIRE(rhs.xy() == rhs.yx());
  this->mxx = rhs.xx();
  this->mxy = rhs.xy();
  this->myy = rhs.yy();
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::
operator=(const GeomTensor<3>& rhs) {
  this->mxx = rhs.xx();
  this->mxy = rhs.xy();
  this->mxz = rhs.xz();
  this->myy = rhs.yy();
  this->myz = rhs.yz();
  this->mzz = rhs.zz();
  return *this;
}

template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::
operator=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx = rhs.mxx;
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::
operator=(const GeomSymmetricTensor<2>& rhs) {
  this->mxx = rhs.mxx;
  this->mxy = rhs.mxy;
  this->myy = rhs.myy;
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::
operator=(const GeomSymmetricTensor<3>& rhs) {
  this->mxx = rhs.mxx;
  this->mxy = rhs.mxy;
  this->mxz = rhs.mxz;
  this->myy = rhs.myy;
  this->myz = rhs.myz;
  this->mzz = rhs.mzz;
  return *this;
}

//------------------------------------------------------------------------------
// Assignment with a scalar.  This sets the diagonal elements to the given
// value, and the off-diagonals to zero.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>&
GeomSymmetricTensor<nDim>::operator=(const double rhs) {
  *this = one*rhs;
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
  return *(begin() + elementIndex(row, column));
}

template<int nDim>
inline
double&
GeomSymmetricTensor<nDim>::operator()(const typename GeomSymmetricTensor<nDim>::size_type row,
                                      const typename GeomSymmetricTensor<nDim>::size_type column) {
  REQUIRE(row < nDim);
  REQUIRE(column < nDim);
  return *(begin() + elementIndex(row, column));
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
  return this->mxx;
}

//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomSymmetricTensor<nDim>::xy() const {
  return this->mxy;
}

template<>
inline
double
GeomSymmetricTensor<1>::xy() const {
  return 0.0;
}

//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::xz() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<2>::xz() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<3>::xz() const {
  return this->mxz;
}

//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::yx() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<2>::yx() const {
  return this->mxy;
}

template<>
inline
double
GeomSymmetricTensor<3>::yx() const {
  return this->mxy;
}

//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::yy() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<2>::yy() const {
  return this->myy;
}

template<>
inline
double
GeomSymmetricTensor<3>::yy() const {
  return this->myy;
}

//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::yz() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<2>::yz() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<3>::yz() const {
  return this->myz;
}

//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::zx() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<2>::zx() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<3>::zx() const {
  return this->mxz;
}

//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::zy() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<2>::zy() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<3>::zy() const {
  return this->myz;
}

//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::zz() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<2>::zz() const {
  return 0.0;
}

template<>
inline
double
GeomSymmetricTensor<3>::zz() const {
  return this->mzz;
}

//------------------------------------------------------------------------------
// Set the individual elements, as above.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::xx(const double val) {
  this->mxx = val;
}

//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::xy(const double val) {
  this->mxy = val;
}

template<>
inline
void
GeomSymmetricTensor<1>::xy(const double val) {
}

//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::xz(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<2>::xz(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<3>::xz(const double val) {
  this->mxz = val;
}

//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::yx(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<2>::yx(const double val) {
  this->mxy = val;
}

template<>
inline
void
GeomSymmetricTensor<3>::yx(const double val) {
  this->mxy = val;
}

//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::yy(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<2>::yy(const double val) {
  this->myy = val;
}

template<>
inline
void
GeomSymmetricTensor<3>::yy(const double val) {
  this->myy = val;
}

//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::yz(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<2>::yz(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<3>::yz(const double val) {
  this->myz = val;
}

//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::zx(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<2>::zx(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<3>::zx(const double val) {
  this->mxz = val;
}

//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::zy(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<2>::zy(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<3>::zy(const double val) {
  this->myz = val;
}

//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::zz(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<2>::zz(const double val) {
}

template<>
inline
void
GeomSymmetricTensor<3>::zz(const double val) {
  this->mzz = val;
}

//------------------------------------------------------------------------------
// Access the individual rows of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
GeomSymmetricTensor<1>::
getRow(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 1);
  return GeomVector<1>((*this)(index, 0));
}

template<>
inline
GeomVector<2>
GeomSymmetricTensor<2>::
getRow(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 2);
  return GeomVector<2>((*this)(index, 0), (*this)(index, 1));
}

template<>
inline
GeomVector<3>
GeomSymmetricTensor<3>::
getRow(const GeomSymmetricTensor<3>::size_type index) const {
  REQUIRE(index < 3);
  return GeomVector<3>((*this)(index, 0), (*this)(index, 1), (*this)(index, 2));
}

//------------------------------------------------------------------------------
// Access the individual columns of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
GeomSymmetricTensor<1>::
getColumn(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 1);
  return GeomVector<1>((*this)(0, index));
}

template<>
inline
GeomVector<2>
GeomSymmetricTensor<2>::
getColumn(const GeomSymmetricTensor<2>::size_type index) const {
  REQUIRE(index < 2);
  return GeomVector<2>((*this)(0, index), (*this)(1, index));
}

template<>
inline
GeomVector<3>
GeomSymmetricTensor<3>::
getColumn(const GeomSymmetricTensor<3>::size_type index) const {
  REQUIRE(index < 3);
  return GeomVector<3>((*this)(0, index), (*this)(1, index), (*this)(2, index));
}

//------------------------------------------------------------------------------
// Set a row of the GeomSymmetricTensor to a GeomVector.
//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::
setRow(const GeomSymmetricTensor<1>::size_type index,
       const GeomVector<1>& vec) {
  REQUIRE(index < 1);
  (*this)(index, 0) = vec(0);
}

template<>
inline
void
GeomSymmetricTensor<2>::
setRow(const GeomSymmetricTensor<2>::size_type index,
       const GeomVector<2>& vec) {
  REQUIRE(index < 2);
  (*this)(index, 0) = vec(0);
  (*this)(index, 1) = vec(1);
}

template<>
inline
void
GeomSymmetricTensor<3>::
setRow(const GeomSymmetricTensor<3>::size_type index,
       const GeomVector<3>& vec) {
  REQUIRE(index < 3);
  (*this)(index, 0) = vec(0);
  (*this)(index, 1) = vec(1);
  (*this)(index, 2) = vec(2);
}

//------------------------------------------------------------------------------
// Set a column of the GeomSymmetricTensor to a GeomVector.
//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::
setColumn(const GeomSymmetricTensor<1>::size_type index,
          const GeomVector<1>& vec) {
  REQUIRE(index < 1);
  (*this)(0, index) = vec(0);
}

template<>
inline
void
GeomSymmetricTensor<2>::
setColumn(const GeomSymmetricTensor<2>::size_type index,
          const GeomVector<2>& vec) {
  REQUIRE(index < 2);
  (*this)(0, index) = vec(0);
  (*this)(1, index) = vec(1);
}

template<>
inline
void
GeomSymmetricTensor<3>::
setColumn(const GeomSymmetricTensor<3>::size_type index,
          const GeomVector<3>& vec) {
  REQUIRE(index < 3);
  (*this)(0, index) = vec(0);
  (*this)(1, index) = vec(1);
  (*this)(2, index) = vec(2);
}

//------------------------------------------------------------------------------
// Iterators to the raw data.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::iterator
GeomSymmetricTensor<nDim>::begin() {
  return &(this->mxx);
}

template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::iterator
GeomSymmetricTensor<nDim>::end() {
  return &(this->mxx) + mNumElements;
}

template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::const_iterator
GeomSymmetricTensor<nDim>::begin() const{
  return &(this->mxx);
}

template<int nDim>
inline
typename GeomSymmetricTensor<nDim>::const_iterator
GeomSymmetricTensor<nDim>::end() const {
  return &(this->mxx) + mNumElements;
}

//------------------------------------------------------------------------------
// Zero out the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::Zero() {
  *this = zero;
}

//------------------------------------------------------------------------------
// Force the tensor to be the identity tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomSymmetricTensor<nDim>::Identity() {
  *this = one;
}

//------------------------------------------------------------------------------
// Return the negative of a tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator-() const {
  GeomSymmetricTensor<1> result;
  result.mxx = -(this->mxx);
  return result;
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator-() const {
  GeomSymmetricTensor<2> result;
  result.mxx = -(this->mxx);
  result.mxy = -(this->mxy);
  result.myy = -(this->myy);
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator-() const {
  GeomSymmetricTensor<3> result;
  result.mxx = -(this->mxx);
  result.mxy = -(this->mxy);
  result.mxz = -(this->mxz);
  result.myy = -(this->myy);
  result.myz = -(this->myz);
  result.mzz = -(this->mzz);
  return result;
}

//------------------------------------------------------------------------------
// Add two tensors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator+(const GeomTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
  result += rhs;
  return result;
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
operator+(const GeomSymmetricTensor<nDim>& rhs) const {
  GeomSymmetricTensor<nDim> result(*this);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a tensor from another.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator-(const GeomTensor<nDim>& rhs) const {
  GeomTensor<nDim> result(*this);
  result -= rhs;
  return result;
}

template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
operator-(const GeomSymmetricTensor<nDim>& rhs) const {
  GeomSymmetricTensor<nDim> result(*this);
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiply two tensors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator*(const GeomTensor<nDim>& rhs) const {
  return this->dot(rhs);
}

template<int nDim>
inline
GeomTensor<nDim>
GeomSymmetricTensor<nDim>::
operator*(const GeomSymmetricTensor<nDim>& rhs) const {
  return this->dot(rhs);
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomSymmetricTensor<nDim>::operator*(const GeomVector<nDim>& rhs) const {
  return this->dot(rhs);
}

//------------------------------------------------------------------------------
// Add a scalar to a tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator+(const double rhs) const {
  return GeomSymmetricTensor<1>(this->mxx + rhs);
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator+(const double rhs) const {
  GeomSymmetricTensor<2> result(*this);
  result.mxx += rhs;
  result.mxy += rhs;
  result.myy += rhs;
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator+(const double rhs) const {
  GeomSymmetricTensor<3> result(*this);
  result.mxx += rhs;
  result.mxy += rhs;
  result.mxz += rhs;
  result.myy += rhs;
  result.myz += rhs;
  result.mzz += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a scalar from a tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator-(const double rhs) const {
  return GeomSymmetricTensor<1>(this->mxx - rhs);
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator-(const double rhs) const {
  GeomSymmetricTensor<2> result(*this);
  result.mxx -= rhs;
  result.mxy -= rhs;
  result.myy -= rhs;
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator-(const double rhs) const {
  GeomSymmetricTensor<3> result(*this);
  result.mxx -= rhs;
  result.mxy -= rhs;
  result.mxz -= rhs;
  result.myy -= rhs;
  result.myz -= rhs;
  result.mzz -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiply a tensor by a scalar
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator*(const double rhs) const {
  return GeomSymmetricTensor<1>(this->mxx * rhs);
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator*(const double rhs) const {
  GeomSymmetricTensor<2> result(*this);
  result.mxx *= rhs;
  result.mxy *= rhs;
  result.myy *= rhs;
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator*(const double rhs) const {
  GeomSymmetricTensor<3> result(*this);
  result.mxx *= rhs;
  result.mxy *= rhs;
  result.mxz *= rhs;
  result.myy *= rhs;
  result.myz *= rhs;
  result.mzz *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Divide a tensor by a scalar
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  return GeomSymmetricTensor<1>(this->mxx / rhs);
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  GeomSymmetricTensor<2> result(*this);
  result.mxx *= rhsInv;
  result.mxy *= rhsInv;
  result.myy *= rhsInv;
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  GeomSymmetricTensor<3> result(*this);
  result.mxx *= rhsInv;
  result.mxy *= rhsInv;
  result.mxz *= rhsInv;
  result.myy *= rhsInv;
  result.myz *= rhsInv;
  result.mzz *= rhsInv;
  return result;
}

//------------------------------------------------------------------------------
// Add two symmetric tensors in place.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator+=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx += rhs.mxx;
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator+=(const GeomSymmetricTensor<2>& rhs) {
  this->mxx += rhs.mxx;
  this->mxy += rhs.mxy;
  this->myy += rhs.myy;
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator+=(const GeomSymmetricTensor<3>& rhs) {
  this->mxx += rhs.mxx;
  this->mxy += rhs.mxy;
  this->mxz += rhs.mxz;
  this->myy += rhs.myy;
  this->myz += rhs.myz;
  this->mzz += rhs.mzz;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a tensor from this one in place.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator-=(const GeomSymmetricTensor<1>& rhs) {
  this->mxx -= rhs.mxx;
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator-=(const GeomSymmetricTensor<2>& rhs) {
  this->mxx -= rhs.mxx;
  this->mxy -= rhs.mxy;
  this->myy -= rhs.myy;
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator-=(const GeomSymmetricTensor<3>& rhs) {
  this->mxx -= rhs.mxx;
  this->mxy -= rhs.mxy;
  this->mxz -= rhs.mxz;
  this->myy -= rhs.myy;
  this->myz -= rhs.myz;
  this->mzz -= rhs.mzz;
  return *this;
}

//------------------------------------------------------------------------------
// Add a scalar to this tensor in place.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator+=(const double rhs) {
  this->mxx += rhs;
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator+=(const double rhs) {
  this->mxx += rhs;
  this->mxy += rhs;
  this->myy += rhs;
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator+=(const double rhs) {
  this->mxx += rhs;
  this->mxy += rhs;
  this->mxz += rhs;
  this->myy += rhs;
  this->myz += rhs;
  this->mzz += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a scalar from this tensor in place.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator-=(const double rhs) {
  this->mxx -= rhs;
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator-=(const double rhs) {
  this->mxx -= rhs;
  this->mxy -= rhs;
  this->myy -= rhs;
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator-=(const double rhs) {
  this->mxx -= rhs;
  this->mxy -= rhs;
  this->mxz -= rhs;
  this->myy -= rhs;
  this->myz -= rhs;
  this->mzz -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this tensor by a scalar in place.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator*=(const double rhs) {
  this->mxx *= rhs;
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator*=(const double rhs) {
  this->mxx *= rhs;
  this->mxy *= rhs;
  this->myy *= rhs;
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator*=(const double rhs) {
  this->mxx *= rhs;
  this->mxy *= rhs;
  this->mxz *= rhs;
  this->myy *= rhs;
  this->myz *= rhs;
  this->mzz *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this tensor by a scalar in place
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>&
GeomSymmetricTensor<1>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  this->mxx /= rhs;
  return *this;
}

template<>
inline
GeomSymmetricTensor<2>&
GeomSymmetricTensor<2>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  this->mxx *= rhsInv;
  this->mxy *= rhsInv;
  this->myy *= rhsInv;
  return *this;
}

template<>
inline
GeomSymmetricTensor<3>&
GeomSymmetricTensor<3>::operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  const double rhsInv = 1.0/rhs;
  this->mxx *= rhsInv;
  this->mxy *= rhsInv;
  this->mxz *= rhsInv;
  this->myy *= rhsInv;
  this->myz *= rhsInv;
  this->mzz *= rhsInv;
  return *this;
}

//------------------------------------------------------------------------------
// Define the equivalence operator.
//------------------------------------------------------------------------------
template<>
inline
bool
GeomSymmetricTensor<1>::
operator==(const GeomTensor<1>& rhs) const {
  return this->mxx == rhs.xx();
}

template<>
inline
bool
GeomSymmetricTensor<2>::
operator==(const GeomTensor<2>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->mxy == rhs.yx() and
          this->myy == rhs.yy());
}

template<>
inline
bool
GeomSymmetricTensor<3>::
operator==(const GeomTensor<3>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->mxz == rhs.xz() and
          this->myy == rhs.yy() and
          this->myz == rhs.yz() and
          this->mzz == rhs.zz() and 
          this->mxy == rhs.yx() and
          this->mxz == rhs.zx() and
          this->myz == rhs.zy());
}

template<>
inline
bool
GeomSymmetricTensor<1>::
operator==(const GeomSymmetricTensor<1>& rhs) const {
  return this->mxx == rhs.xx();
}

template<>
inline
bool
GeomSymmetricTensor<2>::
operator==(const GeomSymmetricTensor<2>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->myy == rhs.yy());
}

template<>
inline
bool
GeomSymmetricTensor<3>::
operator==(const GeomSymmetricTensor<3>& rhs) const {
  return (this->mxx == rhs.xx() and
          this->mxy == rhs.xy() and
          this->mxz == rhs.xz() and
          this->myy == rhs.yy() and
          this->myz == rhs.yz() and
          this->mzz == rhs.zz() and 
          this->mxy == rhs.yx());
}

template<>
inline
bool
GeomSymmetricTensor<1>::
operator==(const double rhs) const {
  return this->mxx == rhs;
}

template<>
inline
bool
GeomSymmetricTensor<2>::
operator==(const double rhs) const {
  return (this->mxx == rhs and
          this->mxy == rhs and
          this->myy == rhs);
}

template<>
inline
bool
GeomSymmetricTensor<3>::
operator==(const double rhs) const {
  return (this->mxx == rhs and
          this->mxy == rhs and
          this->mxz == rhs and
          this->myy == rhs and
          this->myz == rhs and
          this->mzz == rhs and 
          this->mxy == rhs);
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

template<int nDim>
inline
bool
GeomSymmetricTensor<nDim>::
operator!=(const double rhs) const {
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

template<>
inline
bool
GeomSymmetricTensor<1>::
operator<(const double rhs) const {
  return this->mxx < rhs;
}

template<>
inline
bool
GeomSymmetricTensor<2>::
operator<(const double rhs) const {
  return (this->mxx < rhs and
          this->mxy < rhs and
          this->myy < rhs);
}

template<>
inline
bool
GeomSymmetricTensor<3>::
operator<(const double rhs) const {
  return (this->mxx < rhs and
          this->mxy < rhs and
          this->mxz < rhs and
          this->myy < rhs and
          this->myz < rhs and
          this->mzz < rhs and 
          this->mxy < rhs);
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

template<>
inline
bool
GeomSymmetricTensor<1>::
operator>(const double rhs) const {
  return this->mxx > rhs;
}

template<>
inline
bool
GeomSymmetricTensor<2>::
operator>(const double rhs) const {
  return (this->mxx > rhs and
          this->mxy > rhs and
          this->myy > rhs);
}

template<>
inline
bool
GeomSymmetricTensor<3>::
operator>(const double rhs) const {
  return (this->mxx > rhs and
          this->mxy > rhs and
          this->mxz > rhs and
          this->myy > rhs and
          this->myz > rhs and
          this->mzz > rhs and 
          this->mxy > rhs);
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

template<>
inline
bool
GeomSymmetricTensor<1>::
operator<=(const double rhs) const {
  return this->mxx <= rhs;
}

template<>
inline
bool
GeomSymmetricTensor<2>::
operator<=(const double rhs) const {
  return (this->mxx <= rhs and
          this->mxy <= rhs and
          this->myy <= rhs);
}

template<>
inline
bool
GeomSymmetricTensor<3>::
operator<=(const double rhs) const {
  return (this->mxx <= rhs and
          this->mxy <= rhs and
          this->mxz <= rhs and
          this->myy <= rhs and
          this->myz <= rhs and
          this->mzz <= rhs and 
          this->mxy <= rhs);
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

template<>
inline
bool
GeomSymmetricTensor<1>::
operator>=(const double rhs) const {
  return this->mxx >= rhs;
}

template<>
inline
bool
GeomSymmetricTensor<2>::
operator>=(const double rhs) const {
  return (this->mxx >= rhs and
          this->mxy >= rhs and
          this->myy >= rhs);
}

template<>
inline
bool
GeomSymmetricTensor<3>::
operator>=(const double rhs) const {
  return (this->mxx >= rhs and
          this->mxy >= rhs and
          this->mxz >= rhs and
          this->myy >= rhs and
          this->myz >= rhs and
          this->mzz >= rhs and 
          this->mxy >= rhs);
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
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::Inverse() const {
  CHECK(this->mxx != 0.0);
  return GeomSymmetricTensor<1>(1.0/(this->mxx));
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::Inverse() const {
  REQUIRE(Determinant() != 0.0);
  return GeomSymmetricTensor<2>( (this->myy), -(this->mxy),
                                -(this->mxy),  (this->mxx))/this->Determinant();
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::Inverse() const {
  REQUIRE(Determinant() != 0.0);
  return GeomSymmetricTensor<3>((this->myy)*(this->mzz) - (this->myz)*(this->myz),
                                (this->myz)*(this->mxz) - (this->mxy)*(this->mzz),
                                (this->mxy)*(this->myz) - (this->myy)*(this->mxz),
                                (this->mxz)*(this->myz) - (this->mxy)*(this->mzz),
                                (this->mxx)*(this->mzz) - (this->mxz)*(this->mxz),
                                (this->mxy)*(this->mxz) - (this->mxx)*(this->myz),
                                (this->mxy)*(this->myz) - (this->mxz)*(this->myy),
                                (this->mxz)*(this->mxy) - (this->mxx)*(this->myz),
                                (this->mxx)*(this->myy) - (this->mxy)*(this->mxy))/this->Determinant();
}

//------------------------------------------------------------------------------
// Return the diagonal elements of the GeomSymmetricTensor.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
GeomSymmetricTensor<1>::diagonalElements() const {
  return GeomVector<1>(this->mxx);
}

template<>
inline
GeomVector<2>
GeomSymmetricTensor<2>::diagonalElements() const {
  return GeomVector<2>(this->mxx, this->myy);
}

template<>
inline
GeomVector<3>
GeomSymmetricTensor<3>::diagonalElements() const {
  return GeomVector<3>(this->mxx, this->myy, this->mzz);
}

//------------------------------------------------------------------------------
// Return the trace of the tensor.
//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::Trace() const {
  return this->mxx;
}

template<>
inline
double
GeomSymmetricTensor<2>::Trace() const {
  return this->mxx + this->myy;
}

template<>
inline
double
GeomSymmetricTensor<3>::Trace() const {
  return this->mxx + this->myy + this->mzz;
}

//------------------------------------------------------------------------------
// Return the determinant of the tensor.
//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::Determinant() const {
  return this->mxx;
}

template<>
inline
double
GeomSymmetricTensor<2>::Determinant() const {
  return (this->mxx)*(this->myy) - (this->mxy)*(this->mxy);
}

template<>
inline
double
GeomSymmetricTensor<3>::Determinant() const {
  return ((this->mxx)*(this->myy)*(this->mzz) +
	  (this->mxy)*(this->myz)*(this->mxz) +
	  (this->mxz)*(this->mxy)*(this->myz) -
	  (this->mxx)*(this->myz)*(this->myz) -
	  (this->mxy)*(this->mxy)*(this->mzz) -
	  (this->mxz)*(this->myy)*(this->mxz));
}

//------------------------------------------------------------------------------
// Multiply a tensor with a vector.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
GeomSymmetricTensor<1>::dot(const GeomVector<1>& rhs) const {
  return GeomVector<1>((this->mxx)*(rhs.x()));
}

template<>
inline
GeomVector<2>
GeomSymmetricTensor<2>::dot(const GeomVector<2>& rhs) const {
  return GeomVector<2>((this->mxx)*(rhs.x()) + (this->mxy)*(rhs.y()),
                       (this->mxy)*(rhs.x()) + (this->myy)*(rhs.y()));
}

template<>
inline
GeomVector<3>
GeomSymmetricTensor<3>::dot(const GeomVector<3>& rhs) const {
  return GeomVector<3>((this->mxx)*(rhs.x()) + (this->mxy)*(rhs.y()) + (this->mxz)*(rhs.z()),
                       (this->mxy)*(rhs.x()) + (this->myy)*(rhs.y()) + (this->myz)*(rhs.z()),
                       (this->mxz)*(rhs.x()) + (this->myz)*(rhs.y()) + (this->mzz)*(rhs.z()));
}

//------------------------------------------------------------------------------
// Multiply two tensors.  This is just the linear algebra definition for matrix
// multiplication.
//------------------------------------------------------------------------------
template<>
inline
GeomTensor<1>
GeomSymmetricTensor<1>::dot(const GeomTensor<1>& rhs) const {
  return GeomTensor<1>((this->mxx) * (rhs.xx()));
}

template<>
inline
GeomTensor<2>
GeomSymmetricTensor<2>::dot(const GeomTensor<2>& rhs) const {
  return GeomTensor<2>((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()),
                       (this->mxx)*(rhs.xy()) + (this->mxy)*(rhs.yy()),
                       (this->mxy)*(rhs.xx()) + (this->myy)*(rhs.yx()),
                       (this->mxy)*(rhs.xy()) + (this->myy)*(rhs.yy()));
}

template<>
inline
GeomTensor<3>
GeomSymmetricTensor<3>::dot(const GeomTensor<3>& rhs) const {
  return GeomTensor<3>((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()) + (this->mxz)*(rhs.zx()),
                       (this->mxx)*(rhs.xy()) + (this->mxy)*(rhs.yy()) + (this->mxz)*(rhs.zy()),
                       (this->mxx)*(rhs.xz()) + (this->mxy)*(rhs.yz()) + (this->mxz)*(rhs.zz()),
                       (this->mxy)*(rhs.xx()) + (this->myy)*(rhs.yx()) + (this->myz)*(rhs.zx()),
                       (this->mxy)*(rhs.xy()) + (this->myy)*(rhs.yy()) + (this->myz)*(rhs.zy()),
                       (this->mxy)*(rhs.xz()) + (this->myy)*(rhs.yz()) + (this->myz)*(rhs.zz()),
                       (this->mxz)*(rhs.xx()) + (this->myz)*(rhs.yx()) + (this->mzz)*(rhs.zx()),
                       (this->mxz)*(rhs.xy()) + (this->myz)*(rhs.yy()) + (this->mzz)*(rhs.zy()),
                       (this->mxz)*(rhs.xz()) + (this->myz)*(rhs.yz()) + (this->mzz)*(rhs.zz()));
}

template<>
inline
GeomTensor<1>
GeomSymmetricTensor<1>::dot(const GeomSymmetricTensor<1>& rhs) const {
  return GeomTensor<1>((this->mxx) * (rhs.xx()));
}

template<>
inline
GeomTensor<2>
GeomSymmetricTensor<2>::dot(const GeomSymmetricTensor<2>& rhs) const {
  return GeomTensor<2>((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()),
                       (this->mxx)*(rhs.xy()) + (this->mxy)*(rhs.yy()),
                       (this->mxy)*(rhs.xx()) + (this->myy)*(rhs.yx()),
                       (this->mxy)*(rhs.xy()) + (this->myy)*(rhs.yy()));
}

template<>
inline
GeomTensor<3>
GeomSymmetricTensor<3>::dot(const GeomSymmetricTensor<3>& rhs) const {
  return GeomTensor<3>((this->mxx)*(rhs.xx()) + (this->mxy)*(rhs.yx()) + (this->mxz)*(rhs.zx()),
                       (this->mxx)*(rhs.xy()) + (this->mxy)*(rhs.yy()) + (this->mxz)*(rhs.zy()),
                       (this->mxx)*(rhs.xz()) + (this->mxy)*(rhs.yz()) + (this->mxz)*(rhs.zz()),
                       (this->mxy)*(rhs.xx()) + (this->myy)*(rhs.yx()) + (this->myz)*(rhs.zx()),
                       (this->mxy)*(rhs.xy()) + (this->myy)*(rhs.yy()) + (this->myz)*(rhs.zy()),
                       (this->mxy)*(rhs.xz()) + (this->myy)*(rhs.yz()) + (this->myz)*(rhs.zz()),
                       (this->mxz)*(rhs.xx()) + (this->myz)*(rhs.yx()) + (this->mzz)*(rhs.zx()),
                       (this->mxz)*(rhs.xy()) + (this->myz)*(rhs.yy()) + (this->mzz)*(rhs.zy()),
                       (this->mxz)*(rhs.xz()) + (this->myz)*(rhs.yz()) + (this->mzz)*(rhs.zz()));
}


//------------------------------------------------------------------------------
// Return the doubledot product.
//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::
doubledot(const GeomTensor<1>& rhs) const {
  return (this->mxx)*(rhs.xx());
}

template<>
inline
double
GeomSymmetricTensor<2>::
doubledot(const GeomTensor<2>& rhs) const {
  return ((this->mxx)*(rhs.xx()) +
          (this->mxy)*(rhs.xy()) +
          (this->mxy)*(rhs.yx()) +
          (this->myy)*(rhs.yy()));
}

template<>
inline
double
GeomSymmetricTensor<3>::
doubledot(const GeomTensor<3>& rhs) const {
  return ((this->mxx)*(rhs.xx()) +
          (this->mxy)*(rhs.xy()) +
          (this->mxz)*(rhs.xz()) +
          (this->mxy)*(rhs.yx()) +
          (this->myy)*(rhs.yy()) +
          (this->myz)*(rhs.yz()) +
          (this->mxz)*(rhs.zx()) +
          (this->myz)*(rhs.zy()) +
          (this->mzz)*(rhs.zz()));
}

//------------------------------------------------------------------------------
// Return the doubledot product with a symmetric tensor.
//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::
doubledot(const GeomSymmetricTensor<1>& rhs) const {
  return (this->mxx)*(rhs.xx());
}

template<>
inline
double
GeomSymmetricTensor<2>::
doubledot(const GeomSymmetricTensor<2>& rhs) const {
  return ((this->mxx)*(rhs.xx()) +
          (this->mxy)*(rhs.xy()) +
          (this->mxy)*(rhs.xy()) +
          (this->myy)*(rhs.yy()));
}

template<>
inline
double
GeomSymmetricTensor<3>::
doubledot(const GeomSymmetricTensor<3>& rhs) const {
  return ((this->mxx)*(rhs.xx()) +
          (this->mxy)*(rhs.xy()) +
          (this->mxz)*(rhs.xz()) +
          (this->mxy)*(rhs.yx()) +
          (this->myy)*(rhs.yy()) +
          (this->myz)*(rhs.yz()) +
          (this->mxz)*(rhs.zx()) +
          (this->myz)*(rhs.zy()) +
          (this->mzz)*(rhs.zz()));
}

//------------------------------------------------------------------------------
// Return the doubledot product with ourself.
//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::
selfDoubledot() const {
  return FastMath::square(this->mxx);
}

template<>
inline
double
GeomSymmetricTensor<2>::
selfDoubledot() const {
  return (FastMath::square(this->mxx) +
          2.0*FastMath::square(this->mxy) +
          FastMath::square(this->myy));
}

template<>
inline
double
GeomSymmetricTensor<3>::
selfDoubledot() const {
  return (FastMath::square(this->mxx) +
          2.0*FastMath::square(this->mxy) +
          2.0*FastMath::square(this->mxz) +
          FastMath::square(this->myy) +
          2.0*FastMath::square(this->myz) +
          FastMath::square(this->mzz));
}

//------------------------------------------------------------------------------
// Return the square of this tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
square() const {
  return GeomSymmetricTensor<1>((this->mxx)*(this->mxx));
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
square() const {
  GeomSymmetricTensor<2> result;
  result.mxx = (this->mxx)*(this->mxx) + (this->mxy)*(this->mxy);
  result.mxy = (this->mxy)*(this->mxx + this->myy);
  result.myy = (this->mxy)*(this->mxy) + (this->myy)*(this->myy);
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
square() const {
  GeomSymmetricTensor<3> result;
  result.mxx = (this->mxx)*(this->mxx) + (this->mxy)*(this->mxy) + (this->mxz)*(this->mxz);
  result.mxy = (this->mxy)*(this->mxx + this->myy) + (this->mxz)*(this->myz);
  result.mxz = (this->mxz)*(this->mxx + this->mzz) + (this->mxy)*(this->myz);
  result.myy = (this->mxy)*(this->mxy) + (this->myy)*(this->myy) + (this->myz)*(this->myz);
  result.myz = (this->myz)*(this->myy + this->mzz) + (this->mxy)*(this->mxz);
  result.mzz = (this->mxz)*(this->mxz) + (this->myz)*(this->myz) + (this->mzz)*(this->mzz);
  return result;
}

//------------------------------------------------------------------------------
// Return the cube of this tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
cube() const {
  return GeomSymmetricTensor<1>((this->mxx)*(this->mxx)*(this->mxx));
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
cube() const {
  GeomSymmetricTensor<2> result;
  result.mxx = (this->mxx)*((this->mxx)*(this->mxx) + (this->mxy)*(this->mxy)) + (this->mxy)*((this->mxx)*(this->mxy) + (this->mxy)*(this->myy));
  result.mxy = (this->mxy)*((this->mxx)*(this->mxx) + (this->mxy)*(this->mxy)) + (this->myy)*((this->mxx)*(this->mxy) + (this->mxy)*(this->myy));
  result.myy = (this->mxy)*((this->mxx)*(this->mxy) + (this->mxy)*(this->myy)) + (this->myy)*((this->mxy)*(this->mxy) + (this->myy)*(this->myy));
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
cube() const {
  GeomSymmetricTensor<3> result;
  result.mxx = (this->mxx)*(this->mxx)*(this->mxx) + 2.0*(this->mxx)*((this->mxy)*(this->mxy) + (this->mxz)*(this->mxz)) +
    (this->mxy)*(this->mxy)*(this->myy) + 2.0*(this->mxy)*(this->mxz)*(this->myz) + (this->mxz)*(this->mxz)*(this->mzz);
  result.mxy = (this->mxx)*(this->mxx)*(this->mxy) + (this->mxy)*(this->mxy)*(this->mxy) + (this->mxx)*((this->mxy)*(this->myy) + (this->mxz)*(this->myz)) + 
    (this->mxy)*((this->mxz)*(this->mxz) + (this->myy)*(this->myy) + (this->myz)*(this->myz)) + (this->mxz)*(this->myz)*((this->myy) + (this->mzz));
  result.mxz = (this->mxz)*((this->mxx)*(this->mxx) + (this->mxy)*(this->mxy) + (this->mxz)*(this->mxz)) + (this->myz)*((this->mxx)*(this->mxy) + (this->mxy)*(this->myy) + (this->mxz)*(this->myz)) +
    (this->mzz)*((this->mxx)*(this->mxz) + (this->mxy)*(this->myz) + (this->mxz)*(this->mzz));
  result.myy = (this->mxx)*(this->mxy)*(this->mxy) + 2.0*(this->mxy)*(this->mxy)*(this->myy) + (this->myy)*(this->myy)*(this->myy) + 
    2.0*(this->mxy)*(this->mxz)*(this->myz) + 2.0*(this->myy)*(this->myz)*(this->myz) + (this->myz)*(this->myz)*(this->mzz);
  result.myz = (this->mxx)*(this->mxy)*(this->mxz) + (this->mxy)*(this->mxy)*(this->myz) + (this->mxy)*(this->mxz)*((this->myy) + (this->mzz)) +
    (this->myz)*((this->mxz)*(this->mxz) + (this->myy)*(this->myy) + (this->myz)*(this->myz) + (this->myy)*(this->mzz) + (this->mzz)*(this->mzz));
  result.mzz = (this->mxx)*(this->mxz)*(this->mxz) + 2.0*(this->mxy)*(this->mxz)*(this->myz) + (this->myy)*(this->myz)*(this->myz) + 
    2.0*(this->mxz)*(this->mxz)*(this->mzz) + 2.0*(this->myz)*(this->myz)*(this->mzz) + (this->mzz)*(this->mzz)*(this->mzz);
  return result;
}

//------------------------------------------------------------------------------
// Return a new tensor with the elements of this tensor squared.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomSymmetricTensor<1>::
squareElements() const {
  return GeomSymmetricTensor<1>((this->mxx)*(this->mxx));
}

template<>
inline
GeomSymmetricTensor<2>
GeomSymmetricTensor<2>::
squareElements() const {
  GeomSymmetricTensor<2> result;
  result.mxx = (this->mxx)*(this->mxx);
  result.mxy = (this->mxy)*(this->mxy);
  result.myy = (this->myy)*(this->myy);
  return result;
}

template<>
inline
GeomSymmetricTensor<3>
GeomSymmetricTensor<3>::
squareElements() const {
  GeomSymmetricTensor<3> result;
  result.mxx = (this->mxx)*(this->mxx);
  result.mxy = (this->mxy)*(this->mxy);
  result.mxz = (this->mxz)*(this->mxz);
  result.myy = (this->myy)*(this->myy);
  result.myz = (this->myz)*(this->myz);
  result.mzz = (this->mzz)*(this->mzz);
  return result;
}

//------------------------------------------------------------------------------
// Apply a rotational transform to this tensor.
//------------------------------------------------------------------------------
template<>
inline
void
GeomSymmetricTensor<1>::
rotationalTransform(const GeomTensor<1>& R) {
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-5));
}

template<>
inline
void
GeomSymmetricTensor<2>::
rotationalTransform(const GeomTensor<2>& R) {
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-5));

  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->myy;

  const double R0 = R.xx();
  const double R1 = R.xy();
  const double R2 = R.yx();
  const double R3 = R.yy();

  const double T0 = A0*R0 + A1*R1;
  const double T1 = A1*R0 + A2*R1;

  this->mxx = R0*T0 + R1*T1;
  this->mxy = R2*T0 + R3*T1;
  this->myy = R2*(A0*R2 + A1*R3) + R3*(A1*R2 + A2*R3);
}

template<>
inline
void
GeomSymmetricTensor<3>::
rotationalTransform(const GeomTensor<3>& R) {
  REQUIRE(fuzzyEqual(std::abs(R.Determinant()), 1.0, 1.0e-5));

  const double A0 = this->mxx;
  const double A1 = this->mxy;
  const double A2 = this->mxz;
  const double A3 = this->myy;
  const double A4 = this->myz;
  const double A5 = this->mzz;

  const double R0 = R.xx();
  const double R1 = R.xy();
  const double R2 = R.xz();
  const double R3 = R.yx();
  const double R4 = R.yy();
  const double R5 = R.yz();
  const double R6 = R.zx();
  const double R7 = R.zy();
  const double R8 = R.zz();

  const double T0 = A0*R0 + A1*R1 + A2*R2;
  const double T1 = A0*R3 + A1*R4 + A2*R5;
  const double T2 = A1*R0 + A3*R1 + A4*R2;
  const double T3 = A1*R3 + A3*R4 + A4*R5;
  const double T4 = A2*R0 + A4*R1 + A5*R2;
  const double T5 = A2*R3 + A4*R4 + A5*R5;

  this->mxx = R0*T0 + R1*T2 + R2*T4;
  this->mxy = R3*T0 + R4*T2 + R5*T4;
  this->mxz = R6*T0 + R7*T2 + R8*T4;
  this->myy = R3*T1 + R4*T3 + R5*T5;
  this->myz = R6*T1 + R7*T3 + R8*T5;
  this->mzz = R6*(A0*R6 + A1*R7 + A2*R8) + R7*(A1*R6 + A3*R7 + A4*R8) + R8*(A2*R6 + A4*R7 + A5*R8);
}

//------------------------------------------------------------------------------
// Return the maximum absolute value of the elements.
//------------------------------------------------------------------------------
template<>
inline
double
GeomSymmetricTensor<1>::
maxAbsElement() const {
  return std::abs(this->mxx);
}

template<>
inline
double
GeomSymmetricTensor<2>::
maxAbsElement() const {
  return std::max(std::abs(this->mxx),
                  std::max(std::abs(this->mxy), 
                           std::abs(this->myy)));
}

template<>
inline
double
GeomSymmetricTensor<3>::
maxAbsElement() const {
  return std::max(std::abs(this->mxx), 
                  std::max(std::abs(this->mxy), 
                           std::max(std::abs(this->mxz), 
                                    std::max(std::abs(this->myy), 
                                             std::max(std::abs(this->myz), 
                                                      std::abs(this->mzz))))));
}

//------------------------------------------------------------------------------
// Find the eigenvalues of a tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
GeomSymmetricTensor<1>::eigenValues() const {
  return GeomVector<1>(mxx);
}

//----------------------------------------------------------------------
template<>
inline
GeomVector<2>
GeomSymmetricTensor<2>::eigenValues() const {
  if (std::abs(xy()) < 1.0e-50) {
    return diagonalElements();
  } else {
    const double b = Trace();
    const double c = Determinant();
    const double q = 0.5*(b + sgn(b)*std::sqrt(std::max(0.0, b*b - 4.0*c)));
    CHECK(q != 0.0);
    return GeomVector<2>(q, c/q); // (q + 1.0e-50*sgn(q)));
  }
}

//----------------------------------------------------------------------
// This 3-D version is based on the ideas from David Eberly at
// www.geometrictools.com
//----------------------------------------------------------------------
template<>
inline
GeomVector<3>
GeomSymmetricTensor<3>::eigenValues() const {
  const double fscale = std::max(10.0*std::numeric_limits<double>::epsilon(), this->maxAbsElement()); 
  CHECK(fscale > 0.0);
  const double fscalei = 1.0/fscale;
  const double a00 = this->xx()*fscalei;
  const double a01 = this->xy()*fscalei;
  const double a02 = this->xz()*fscalei;
  const double a11 = this->yy()*fscalei;
  const double a12 = this->yz()*fscalei;
  const double a22 = this->zz()*fscalei;
  const double c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
  const double c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12;
  const double c2 = a00 + a11 + a22;
  const double c2Div3 = c2*onethird;
  const double aDiv3 = std::min(0.0, onethird*(c1 - c2*c2Div3));
  const double mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1));
  const double q = std::min(0.0, mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3);
  CHECK(-aDiv3 >= 0.0);
  CHECK(-q >= 0.0);
  const double mag = std::sqrt(-aDiv3);
  const double angle = atan2(std::sqrt(-q), mbDiv2)*onethird;
  const double cs = cos(angle);
  const double sn = sin(angle);
  return GeomVector<3>(fscale*(c2Div3 + 2.0*mag*cs),
                       fscale*(c2Div3 - mag*(cs + sqrt3*sn)),
                       fscale*(c2Div3 - mag*(cs - sqrt3*sn)));
}

//------------------------------------------------------------------------------
// Return the eigen values and eigen vectors of a symmetric tensor
//------------------------------------------------------------------------------
// 1-D.
template<>
inline
EigenStruct<1>
GeomSymmetricTensor<1>::eigenVectors() const {
  EigenStruct<1> result;
  result.eigenValues.x(mxx);
  result.eigenVectors.xx(1.0);
  return result;
}

//------------------------------------------------------------------------------
// 2-D.
template<>
inline
EigenStruct<2>
GeomSymmetricTensor<2>::eigenVectors() const {
  const double fscale = std::max(10.0*std::numeric_limits<double>::epsilon(), this->maxAbsElement()); 
  CHECK(fscale > 0.0);
  const double fscalei = 1.0/fscale;
  const double axx = mxx*fscalei;
  const double axy = mxy*fscalei;
  const double ayy = myy*fscalei;
  EigenStruct<2> result;
  if (std::abs(axy) < 1.0e-50) {
    result.eigenValues = diagonalElements();
    result.eigenVectors = one;
  } else {
    const double theta = 0.5*atan2(2.0*axy, ayy - axx);
    const double xhat = cos(theta);
    const double yhat = sin(theta);
    result.eigenValues.x(xhat*(axx*xhat - axy*yhat) - yhat*(axy*xhat - ayy*yhat));
    result.eigenValues.y(yhat*(axx*yhat + axy*xhat) + xhat*(axy*yhat + ayy*xhat));
    result.eigenValues *= fscale;
    result.eigenVectors = GeomTensor<2>( xhat, yhat,
                                        -yhat, xhat);
  }
  return result;
}

//------------------------------------------------------------------------------
// Compute the square root of the tensor.
//------------------------------------------------------------------------------
template <int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
sqrt() const {
  const typename GeomSymmetricTensor<nDim>::EigenStructType eigen = this->eigenVectors();
  GeomSymmetricTensor<nDim> result;
  for (int i = 0; i != nDim; ++i) {
    REQUIRE(eigen.eigenValues(i) >= 0.0);
    result(i,i) = std::sqrt(eigen.eigenValues(i));
  }
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// Compute the cube root of the tensor.
//------------------------------------------------------------------------------
template <int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
cuberoot() const {
  const typename GeomSymmetricTensor<nDim>::EigenStructType eigen = this->eigenVectors();
  GeomSymmetricTensor<nDim> result;
  const double onethird = 1.0/3.0;
  for (int i = 0; i != nDim; ++i) result(i,i) = FastMath::CubeRootHalley2(eigen.eigenValues(i));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// The general version, raise to an arbitrary power.
//------------------------------------------------------------------------------
template <int nDim>
inline
GeomSymmetricTensor<nDim>
GeomSymmetricTensor<nDim>::
pow(const double p) const {
  const typename GeomSymmetricTensor<nDim>::EigenStructType eigen = this->eigenVectors();
  GeomSymmetricTensor<nDim> result;
  for (int i = 0; i != nDim; ++i)
    result(i,i) = std::pow(std::abs(eigen.eigenValues(i)), p) * sgn(eigen.eigenValues(i));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
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

  typedef Spheral::GeomVector<nDim> Vector;
  typedef Spheral::GeomTensor<nDim> Tensor;
  typedef Spheral::GeomSymmetricTensor<nDim> SymTensor;


  // Get the eigen values and eigen vectors.
  EigenStruct<nDim> eigen = tensor.eigenVectors();

  // Limit the eigen values if necessary.
  if (eigen.eigenValues.maxElement() < minValue) {
    return tensor;
  } else {
    SymTensor result;
    for (int i = 0; i != nDim; ++i) {
      result(i,i) = std::min(minValue, eigen.eigenValues(i));
    }
    result.rotationalTransform(eigen.eigenVectors);
    return result;
  }
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

  typedef Spheral::GeomVector<nDim> Vector;
  typedef Spheral::GeomTensor<nDim> Tensor;
  typedef Spheral::GeomSymmetricTensor<nDim> SymTensor;

  // Get the eigen values and eigen vectors.
  EigenStruct<nDim> eigen = tensor.eigenVectors();

  // Limit the eigen values if necessary.
  if (eigen.eigenValues.minElement() > maxValue) {
    return tensor;
  } else {
    SymTensor result;
    for (int i = 0; i != nDim; ++i) {
      result(i,i) = std::max(maxValue, eigen.eigenValues(i));
    }
    result.rotationalTransform(eigen.eigenVectors);
    return result;
  }
}

template<int nDim>
inline
Spheral::GeomSymmetricTensor<nDim>
max(const Spheral::GeomSymmetricTensor<nDim>& tensor, double const maxValue) {
  return max(maxValue, tensor);
}

}
