#include <algorithm>
#include <limits.h>
#include <math.h>
#include <cfloat>
#include <string>

#include "GeomTensor.hh"
#include "GeomSymmetricTensor.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given values for the elements.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>::GeomVector(const double x,
                          const double y,
                          const double z):
  mVecData(x) {
}

template<>
inline
GeomVector<2>::GeomVector(const double x,
                          const double y,
                          const double z):
  mVecData(x, y) {
}

template<>
inline
GeomVector<3>::GeomVector(const double x,
                          const double y,
                          const double z):
  mVecData(x, y, z) {
}

//------------------------------------------------------------------------------
// Copy constructors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>::GeomVector(const GeomVector<nDim>& vec):
  mVecData(vec.mVecData) {
}

template<int nDim>
inline
GeomVector<nDim>::GeomVector(const typename GeomVector<nDim>::VectorStorage& vec):
  mVecData(vec) {
}

//------------------------------------------------------------------------------
// The assignment operators.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator=(const GeomVector<nDim>& vec) {
  this->mVecData = vec.mVecData;
  return *this;
}

template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator=(const VectorStorage& vec) {
  this->mVecData = vec;
  return *this;
}

template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator=(const double val) {
  mVecData = VectorStorage::Constant(val);
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>::~GeomVector() {}

//------------------------------------------------------------------------------
// Return the (index) element using the parenthesis operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::operator()(typename GeomVector<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return mVecData(index);
}

template<int nDim>
inline
double&
GeomVector<nDim>::operator()(typename GeomVector<nDim>::size_type index) {
  REQUIRE(index < nDim);
  return mVecData(index);
}

//------------------------------------------------------------------------------
// Return the (index) element using the bracket operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::operator[](typename GeomVector<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return mVecData[nDim];
}

template<int nDim>
inline
double&
GeomVector<nDim>::operator[](typename GeomVector<nDim>::size_type index) {
  REQUIRE(index < nDim);
  return mVecData[nDim];
}

//------------------------------------------------------------------------------
// Return the x (first) element.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::x() const {
  return mVecData.x();
}

//------------------------------------------------------------------------------
// Return the y (second) element
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::y() const {
  return mVecData.y();
}

template<>
inline
double
GeomVector<1>::y() const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Return the z (third) element
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::z() const {
  return mVecData.z();
}

template<>
inline
double
GeomVector<1>::z() const {
  return 0.0;
}

template<>
inline
double
GeomVector<2>::z() const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Set the x (first) element.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomVector<nDim>::x(const double val) {
  mVecData(0) = val;
}

//------------------------------------------------------------------------------
// Set the y (second) element
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomVector<nDim>::y(const double val) {
  mVecData(1) = val;
}

template<>
inline
void
GeomVector<1>::y(const double val) {
}

//------------------------------------------------------------------------------
// Set the z (third) element
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomVector<nDim>::z(const double val) {
  mVecData(2) = val;
}

template<>
inline
void
GeomVector<1>::z(const double val) {
}

template<>
inline
void
GeomVector<2>::z(const double val) {
}

//------------------------------------------------------------------------------
// Provide begin/end iterators over the elements of the Vector.
//------------------------------------------------------------------------------
// Non-const versions.
template<int nDim>
inline
typename GeomVector<nDim>::iterator
GeomVector<nDim>::begin() {
  return mVecData.data();
}

template<int nDim>
inline
typename GeomVector<nDim>::iterator
GeomVector<nDim>::end() {
  return mVecData.data() + nDim;
}

// Const versions.
template<int nDim>
inline
typename GeomVector<nDim>::const_iterator
GeomVector<nDim>::begin() const {
  return mVecData.data();
}

template<int nDim>
inline
typename GeomVector<nDim>::const_iterator
GeomVector<nDim>::end() const {
  return mVecData.data() + nDim;
}

//------------------------------------------------------------------------------
// Zero out the Vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomVector<nDim>::Zero() {
  mVecData = VectorStorage::Zero();
}

//------------------------------------------------------------------------------
// Return the negative of a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator-() const {
  return GeomVector<nDim>((-mVecData).eval());
}

//------------------------------------------------------------------------------
// Add two vectors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator+(const GeomVector<nDim>& rhs) const {
  return GeomVector<nDim>((mVecData + rhs.mVecData).eval());
}

//------------------------------------------------------------------------------
// Subtract a vector from another.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator-(const GeomVector<nDim>& rhs) const {
  return GeomVector<nDim>((mVecData - rhs.mVecData).eval());
}

//------------------------------------------------------------------------------
// Mutiply two vectors.  For our purposes this returns the dyad.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomVector<nDim>::operator*(const GeomVector<nDim>& rhs) const {
  return GeomTensor<nDim>((mVecData * rhs.mVecData.transpose()).eval());
}

//------------------------------------------------------------------------------
// Multiply a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator*(const double val) const {
  return GeomVector<nDim>((mVecData * val).eval());
}

//------------------------------------------------------------------------------
// Divide a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator/(const double val) const {
  CHECK(val != 0.0);
  return GeomVector<nDim>((mVecData / val).eval());
}

//------------------------------------------------------------------------------
// Add two vectors in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator+=(const GeomVector<nDim>& vec) {
  mVecData += vec.mVecData;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a vector in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator-=(const GeomVector<nDim>& vec) {
  mVecData -= vec.mVecData;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this vector by a scalar in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator*=(const double val) {
  mVecData *= val;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this vector by a scalar in place.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator/=(const double val) {
  REQUIRE(val != 0.0);
  mVecData /= val;
  return *this;
}

//------------------------------------------------------------------------------
// Return (-1, 0, 1) if this vector is (less than, equal to, greater than)
// the given vector.
//------------------------------------------------------------------------------
template<>
inline
int
GeomVector<1>::compare(const GeomVector<1>& vec) const {
  return (mVecData(0) < vec.mVecData(0) ? -1 :
          mVecData(0) > vec.mVecData(0) ?  1 :
          0);
}

template<>
inline
int
GeomVector<2>::compare(const GeomVector<2>& vec) const {
  return (mVecData(1) < vec.mVecData(1) ? -1 :
          mVecData(1) > vec.mVecData(1) ?  1 :
          mVecData(0) < vec.mVecData(0) ? -1 :
          mVecData(0) > vec.mVecData(0) ?  1 :
          0);
}

template<>
inline
int
GeomVector<3>::compare(const GeomVector<3>& vec) const {
  return (mVecData(2) < vec.mVecData(2) ? -1 :
          mVecData(2) > vec.mVecData(2) ?  1 :
          mVecData(1) < vec.mVecData(1) ? -1 :
          mVecData(1) > vec.mVecData(1) ?  1 :
          mVecData(0) < vec.mVecData(0) ? -1 :
          mVecData(0) > vec.mVecData(0) ?  1 :
          0);
}

//------------------------------------------------------------------------------
// Return (-1, 0, 1) if this vector is (less than, equal to, greater than)
// the given double.
//------------------------------------------------------------------------------
template<>
inline
int
GeomVector<1>::compare(const double val) const {
  return (mVecData(0) < val ? -1 :
          mVecData(0) > val ?  1 :
          0);
}

template<>
inline
int
GeomVector<2>::compare(const double val) const {
  return (mVecData(1) < val ? -1 :
          mVecData(1) > val ?  1 :
          mVecData(0) < val ? -1 :
          mVecData(0) > val ?  1 :
          0);
}

template<>
inline
int
GeomVector<3>::compare(const double val) const {
  return (mVecData(2) < val ? -1 :
          mVecData(2) > val ?  1 :
          mVecData(1) < val ? -1 :
          mVecData(1) > val ?  1 :
          mVecData(0) < val ? -1 :
          mVecData(0) > val ?  1 :
          0);
}

//------------------------------------------------------------------------------
// The equivalence comparator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator==(const GeomVector<nDim>& rhs) const {
  return this->mVecData == rhs.mVecData;
}

//------------------------------------------------------------------------------
// The non-equivalence comparator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator!=(const GeomVector<nDim>& rhs) const {
  return this->mVecData != rhs.mVecData;
}

//------------------------------------------------------------------------------
// The less than comparator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator<(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) == -1;
}

//------------------------------------------------------------------------------
// The greater than comparator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator>(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) == 1;
}

//------------------------------------------------------------------------------
// The less than or equal comparator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator<=(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) <= 0;
}

//------------------------------------------------------------------------------
// The greater than or equal comparator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator>=(const GeomVector<nDim>& rhs) const {
  return this->compare(rhs) >= 0;
}

//------------------------------------------------------------------------------
// The equivalence comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator==(const double val) const {
  return (mVecData.array() == val).all();
}

//------------------------------------------------------------------------------
// The non-equivalence comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator!=(const double val) const {
  return not (*this == val);
}

//------------------------------------------------------------------------------
// The less than comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator<(const double val) const {
  return this->compare(val) == -1;
}

//------------------------------------------------------------------------------
// The greater than comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator>(const double val) const {
  return this->compare(val) == 1;
}

//------------------------------------------------------------------------------
// The less than or equal comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator<=(const double val) const {
  return this->compare(val) <= 0;
}

//------------------------------------------------------------------------------
// The greater than or equal comparator (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator>=(const double val) const {
  return this->compare(val) >= 0;
}

//------------------------------------------------------------------------------
// Dot two vectors.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::dot(const GeomVector<nDim>& vec) const {
  return this->mVecData.dot(vec.mVecData);
}

//------------------------------------------------------------------------------
// Cross two vectors.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<3>
GeomVector<3>::cross(const GeomVector<3>& vec) const {
  return GeomVector<3>((this->mVecData.cross(vec.mVecData)).eval());
}

template<>
inline
GeomVector<3>
GeomVector<2>::cross(const GeomVector<2>& vec) const {
  return GeomVector<3>(0.0, 0.0,
                       this->mVecData(0)*vec.mVecData(1) - this->mVecData(1)*vec.mVecData(0));
}

template<>
inline
GeomVector<3>
GeomVector<1>::cross(const GeomVector<1>& vec) const {
  return GeomVector<3>(0.0, 0.0, 0.0);
}

//------------------------------------------------------------------------------
// Perform the dyad operation with another vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomVector<nDim>::dyad(const GeomVector<nDim>& rhs) const {
  return GeomTensor<nDim>((mVecData * rhs.mVecData.transpose()).eval());
}

//------------------------------------------------------------------------------
// Perform the dyad operation with ourself, resulting in a symmetric tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomSymmetricTensor<nDim>
GeomVector<nDim>::selfdyad() const {
  return GeomSymmetricTensor<nDim>((mVecData*(mVecData.transpose())).eval());
}

//------------------------------------------------------------------------------
// Return a unit vector with the direction of this one.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::unitVector() const {
  const double mag = this->magnitude();
  return mag > 1.0e-50 ? (*this)/mag : GeomVector<nDim>(1.0, 0.0, 0.0);
}

//------------------------------------------------------------------------------
// Return the magnitude of the Vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::magnitude() const {
  return this->mVecData.norm();
}

//------------------------------------------------------------------------------
// Return the square of the magnitude.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::magnitude2() const {
  return this->mVecData.squaredNorm();
}

//------------------------------------------------------------------------------
// Return the minimum element.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::minElement() const {
  return this->mVecData.minCoeff();
}

//------------------------------------------------------------------------------
// Return the maximum element.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::maxElement() const {
  return this->mVecData.maxCoeff();
}

//------------------------------------------------------------------------------
// Return the maximum element by absolute value.
//------------------------------------------------------------------------------
template<>
inline
double
GeomVector<1>::maxAbsElement() const {
  return std::abs(this->mVecData(0));
}

template<>
inline
double
GeomVector<2>::maxAbsElement() const {
  return std::max(std::abs(this->mVecData(0)), 
                  std::abs(this->mVecData(1)));
}

template<>
inline
double
GeomVector<3>::maxAbsElement() const {
  return std::max(std::abs(this->mVecData(0)),
                  std::max(std::abs(this->mVecData(1)),
                           std::abs(this->mVecData(2))));
}

//------------------------------------------------------------------------------
// Return the sum of the elements.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::sumElements() const {
  return this->mVecData.sum();
}

//******************************************************************************
// Global functions
//******************************************************************************

//------------------------------------------------------------------------------
// Multiply a scalar by a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
operator*(const double val, const GeomVector<nDim>& vec) {
  return vec*val;
}

//------------------------------------------------------------------------------
// Element wise minimum comparison.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
elementWiseMin(const GeomVector<1>& lhs, const GeomVector<1>& rhs) {
  return GeomVector<1>(std::min(lhs.x(), rhs.x()));
}

template<>
inline
GeomVector<2>
elementWiseMin(const GeomVector<2>& lhs, const GeomVector<2>& rhs) {
  return GeomVector<2>(std::min(lhs.x(), rhs.x()),
                       std::min(lhs.y(), rhs.y()));
}

template<>
inline
GeomVector<3>
elementWiseMin(const GeomVector<3>& lhs, const GeomVector<3>& rhs) {
  return GeomVector<3>(std::min(lhs.x(), rhs.x()),
                       std::min(lhs.y(), rhs.y()),
                       std::min(lhs.z(), rhs.z()));
}

//------------------------------------------------------------------------------
// Element wise maximum comparison.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
elementWiseMax(const GeomVector<1>& lhs, const GeomVector<1>& rhs) {
  return GeomVector<1>(std::max(lhs.x(), rhs.x()));
}

template<>
inline
GeomVector<2>
elementWiseMax(const GeomVector<2>& lhs, const GeomVector<2>& rhs) {
  return GeomVector<2>(std::max(lhs.x(), rhs.x()),
                       std::max(lhs.y(), rhs.y()));
}

template<>
inline
GeomVector<3>
elementWiseMax(const GeomVector<3>& lhs, const GeomVector<3>& rhs) {
  return GeomVector<3>(std::max(lhs.x(), rhs.x()),
                       std::max(lhs.y(), rhs.y()),
                       std::max(lhs.z(), rhs.z()));
}

//------------------------------------------------------------------------------
// Access the native Eigen type.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomVector<nDim>::VectorStorage&
GeomVector<nDim>::native() {
  return mVecData;
}

template<int nDim>
inline
const typename GeomVector<nDim>::VectorStorage&
GeomVector<nDim>::native() const {
  return mVecData;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, GeomVector<nDim>& vec) {
  std::string parenthesis;
  is >> parenthesis;
  for (typename GeomVector<nDim>::iterator elementItr = vec.begin();
       elementItr < vec.end();
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
operator<<(std::ostream& os, const GeomVector<nDim>& vec) {
  os << "( ";
  for (typename GeomVector<nDim>::const_iterator elementItr = vec.begin();
       elementItr < vec.end();
       ++elementItr) {
    os << *elementItr << " ";
  }
  os << ")";
  return os;
}

}
