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
  mData(new double[1]) {
  mData[0] = x;
}

template<>
inline
GeomVector<2>::GeomVector(const double x,
                          const double y,
                          const double z):
  mData(new double[2]) {
  mData[0] = x;
  mData[1] = y;
}

template<>
inline
GeomVector<3>::GeomVector(const double x,
                          const double y,
                          const double z):
  mData(new double[3]) {
  mData[0] = x;
  mData[1] = y;
  mData[2] = z;
}

//------------------------------------------------------------------------------
// Copy constructors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>::GeomVector(const GeomVector<nDim>& vec):
  mData(new double[nDim]) {
  std::copy(vec.mData, vec.mData + nDim, mData);
}

//------------------------------------------------------------------------------
// The assignment operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator=(const GeomVector<nDim>& vec) {
  std::copy(vec.mData, vec.mData + nDim, mData);
  return *this;
}

//------------------------------------------------------------------------------
// Set the vector elements to a constant scalar value.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>&
GeomVector<nDim>::operator=(const double val) {
  std::fill(mData, mData + nDim, val);
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>::~GeomVector() {
  delete [] mData;
}

//------------------------------------------------------------------------------
// Return the (index) element using the parenthesis operator.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::operator()(typename GeomVector<nDim>::size_type index) const {
  REQUIRE(index < nDim);
  return mData[index];
}

template<int nDim>
inline
double&
GeomVector<nDim>::operator()(typename GeomVector<nDim>::size_type index) {
  REQUIRE(index < nDim);
  return mData[index];
}

//------------------------------------------------------------------------------
// Return the x (first) element.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::x() const {
  return mData[0];
}

//------------------------------------------------------------------------------
// Return the y (second) element
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomVector<nDim>::y() const {
  REQUIRE(nDim > 1);
  return mData[1];
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
  REQUIRE(nDim > 2);
  return mData[2];
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
  mData[0] = val;
}

//------------------------------------------------------------------------------
// Set the y (second) element
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomVector<nDim>::y(const double val) {
  REQUIRE(nDim > 1);
  mData[1] = val;
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
  REQUIRE(nDim > 2);
  mData[2] = val;
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
  return &mData[0];
}

template<int nDim>
inline
typename GeomVector<nDim>::iterator
GeomVector<nDim>::end() {
  return &mData[nDim];
}

// Const versions.
template<int nDim>
inline
typename GeomVector<nDim>::const_iterator
GeomVector<nDim>::begin() const {
  return &mData[0];
}

template<int nDim>
inline
typename GeomVector<nDim>::const_iterator
GeomVector<nDim>::end() const {
  return &mData[nDim];
}

//------------------------------------------------------------------------------
// Zero out the Vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomVector<nDim>::Zero() {
  std::fill(mData, mData + nDim, 0.0);
}

//------------------------------------------------------------------------------
// Return the negative of a vector.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>
GeomVector<1>::operator-() const {
  return GeomVector<1>(-(mData[0]));
}

template<>
inline
GeomVector<2>
GeomVector<2>::operator-() const {
  return GeomVector<2>(-(mData[0]), -(mData[1]));
}

template<>
inline
GeomVector<3>
GeomVector<3>::operator-() const {
  return GeomVector<3>(-(mData[0]), -(mData[1]), -(mData[2]));
}

//------------------------------------------------------------------------------
// Add two vectors.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator+(const GeomVector<nDim>& vec) const {
  GeomVector<nDim> result(*this);
  result += vec;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a vector from another.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator-(const GeomVector<nDim>& vec) const {
  GeomVector<nDim> result(*this);
  result -= vec;
  return result;
}

//------------------------------------------------------------------------------
// Mutiply two vectors.  For our purposes this returns the dyad.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomTensor<nDim>
GeomVector<nDim>::operator*(const GeomVector<nDim>& vec) const {
  return this->dyad(vec);
}

//------------------------------------------------------------------------------
// Add a scalar to a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator+(const double val) const {
  GeomVector<nDim> result(*this);
  result += val;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a scalar from a vector.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator-(const double val) const {
  GeomVector<nDim> result(*this);
  result -= val;
  return result;
}

//------------------------------------------------------------------------------
// Multiply a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator*(const double val) const {
  GeomVector<nDim> result(*this);
  result *= val;
  return result;
}

//------------------------------------------------------------------------------
// Divide a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
GeomVector<nDim>::operator/(const double val) const {
  CHECK(val != 0.0);
  GeomVector<nDim> result(*this);
  result /= val;
  return result;
}

//------------------------------------------------------------------------------
// Add two vectors in place.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>&
GeomVector<1>::operator+=(const GeomVector<1>& vec) {
  mData[0] += vec.mData[0];
  return *this;
}

template<>
inline
GeomVector<2>&
GeomVector<2>::operator+=(const GeomVector<2>& vec) {
  mData[0] += vec.mData[0];
  mData[1] += vec.mData[1];
  return *this;
}

template<>
inline
GeomVector<3>&
GeomVector<3>::operator+=(const GeomVector<3>& vec) {
  mData[0] += vec.mData[0];
  mData[1] += vec.mData[1];
  mData[2] += vec.mData[2];
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a vector in place.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>& 
GeomVector<1>::operator-=(const GeomVector<1>& vec) {
  mData[0] -= vec.mData[0];
  return *this;
}

template<>
inline
GeomVector<2>&
GeomVector<2>::operator-=(const GeomVector<2>& vec) {
  mData[0] -= vec.mData[0];
  mData[1] -= vec.mData[1];
  return *this;
}

template<>
inline
GeomVector<3>&
GeomVector<3>::operator-=(const GeomVector<3>& vec) {
  mData[0] -= vec.mData[0];
  mData[1] -= vec.mData[1];
  mData[2] -= vec.mData[2];
  return *this;
}

//------------------------------------------------------------------------------
// Add a scalar to this vector in place.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>&
GeomVector<1>::operator+=(const double val) {
  mData[0] += val;
  return *this;
}

template<>
inline
GeomVector<2>&
GeomVector<2>::operator+=(const double val) {
  mData[0] += val;
  mData[1] += val;
  return *this;
}

template<>
inline
GeomVector<3>&
GeomVector<3>::operator+=(const double val) {
  mData[0] += val;
  mData[1] += val;
  mData[2] += val;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a scalar from this vector in place.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>&
GeomVector<1>::operator-=(const double val) {
  mData[0] -= val;
  return *this;
}

template<>
inline
GeomVector<2>&
GeomVector<2>::operator-=(const double val) {
  mData[0] -= val;
  mData[1] -= val;
  return *this;
}

template<>
inline
GeomVector<3>&
GeomVector<3>::operator-=(const double val) {
  mData[0] -= val;
  mData[1] -= val;
  mData[2] -= val;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this vector by a scalar in place.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>&
GeomVector<1>::operator*=(const double val) {
  mData[0] *= val;
  return *this;
}

template<>
inline
GeomVector<2>&
GeomVector<2>::operator*=(const double val) {
  mData[0] *= val;
  mData[1] *= val;
  return *this;
}

template<>
inline
GeomVector<3>&
GeomVector<3>::operator*=(const double val) {
  mData[0] *= val;
  mData[1] *= val;
  mData[2] *= val;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this vector by a scalar in place.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<1>&
GeomVector<1>::operator/=(const double val) {
  REQUIRE(val != 0.0);
  mData[0] /= val;
  return *this;
}

template<>
inline
GeomVector<2>&
GeomVector<2>::operator/=(const double val) {
  REQUIRE(val != 0.0);
  const double valInv = 1.0/val;
  mData[0] *= valInv;
  mData[1] *= valInv;
  return *this;
}

template<>
inline
GeomVector<3>&
GeomVector<3>::operator/=(const double val) {
  REQUIRE(val != 0.0);
  const double valInv = 1.0/val;
  mData[0] *= valInv;
  mData[1] *= valInv;
  mData[2] *= valInv;
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
  return (mData[0] < vec.mData[0] ? -1 :
          mData[0] > vec.mData[0] ?  1 :
          0);
}

template<>
inline
int
GeomVector<2>::compare(const GeomVector<2>& vec) const {
  return (mData[1] < vec.mData[1] ? -1 :
          mData[1] > vec.mData[1] ?  1 :
          mData[0] < vec.mData[0] ? -1 :
          mData[0] > vec.mData[0] ?  1 :
          0);
}

template<>
inline
int
GeomVector<3>::compare(const GeomVector<3>& vec) const {
  return (mData[2] < vec.mData[2] ? -1 :
          mData[2] > vec.mData[2] ?  1 :
          mData[1] < vec.mData[1] ? -1 :
          mData[1] > vec.mData[1] ?  1 :
          mData[0] < vec.mData[0] ? -1 :
          mData[0] > vec.mData[0] ?  1 :
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
  return (mData[0] < val ? -1 :
          mData[0] > val ?  1 :
          0);
}

template<>
inline
int
GeomVector<2>::compare(const double val) const {
  return (mData[1] < val ? -1 :
          mData[1] > val ?  1 :
          mData[0] < val ? -1 :
          mData[0] > val ?  1 :
          0);
}

template<>
inline
int
GeomVector<3>::compare(const double val) const {
  return (mData[2] < val ? -1 :
          mData[2] > val ?  1 :
          mData[1] < val ? -1 :
          mData[1] > val ?  1 :
          mData[0] < val ? -1 :
          mData[0] > val ?  1 :
          0);
}

//------------------------------------------------------------------------------
// The equivalence comparator.
//------------------------------------------------------------------------------
template<>
inline
bool
GeomVector<1>::operator==(const GeomVector<1>& vec) const {
  return mData[0] == vec.mData[0];
}

template<>
inline
bool
GeomVector<2>::operator==(const GeomVector<2>& vec) const {
  return (mData[0] == vec.mData[0]) and (mData[1] == vec.mData[1]);
}

template<>
inline
bool
GeomVector<3>::operator==(const GeomVector<3>& vec) const {
  return (mData[0] == vec.mData[0]) and (mData[1] == vec.mData[1]) and (mData[2] == vec.mData[2]);
}

//------------------------------------------------------------------------------
// The non-equivalence comparator.
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomVector<nDim>::operator!=(const GeomVector<nDim>& vec) const {
  return not (*this == vec);
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
template<>
inline
bool
GeomVector<1>::operator==(const double val) const {
  return mData[0] == val;
}

template<>
inline
bool
GeomVector<2>::operator==(const double val) const {
  return (mData[0] == val) and (mData[1] == val);
}

template<>
inline
bool
GeomVector<3>::operator==(const double val) const {
  return (mData[0] == val) and (mData[1] == val) and (mData[2] == val);
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
template<>
inline
double
GeomVector<1>::dot(const GeomVector<1>& vec) const {
  return mData[0]*vec.mData[0];
}

template<>
inline
double
GeomVector<2>::dot(const GeomVector<2>& vec) const {
  return mData[0]*vec.mData[0] + mData[1]*vec.mData[1];
}

template<>
inline
double
GeomVector<3>::dot(const GeomVector<3>& vec) const {
  return mData[0]*vec.mData[0] + mData[1]*vec.mData[1] + mData[2]*vec.mData[2];
}

//------------------------------------------------------------------------------
// Cross two vectors.
//------------------------------------------------------------------------------
template<>
inline
GeomVector<3>
GeomVector<3>::cross(const GeomVector<3>& vec) const {
  return GeomVector<3>(mData[1]*vec.mData[2] - mData[2]*vec.mData[1],
		       mData[2]*vec.mData[0] - mData[0]*vec.mData[2],
		       mData[0]*vec.mData[1] - mData[1]*vec.mData[0]);
}

template<>
inline
GeomVector<3>
GeomVector<2>::cross(const GeomVector<2>& vec) const {
  return GeomVector<3>(0.0, 0.0,
		       mData[0]*vec.mData[1] - mData[1]*vec.mData[0]);
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
template<>
inline
GeomTensor<1>
GeomVector<1>::dyad(const GeomVector<1>& rhs) const {
  return GeomTensor<1>(mData[0]*rhs(0));
}

template<>
inline
GeomTensor<2>
GeomVector<2>::dyad(const GeomVector<2>& rhs) const {
  return GeomTensor<2>(mData[0]*rhs(0), mData[0]*rhs(1),
                       mData[1]*rhs(0), mData[1]*rhs(1));
}

template<>
inline
GeomTensor<3>
GeomVector<3>::dyad(const GeomVector<3>& rhs) const {
  return GeomTensor<3>(mData[0]*rhs(0), mData[0]*rhs(1), mData[0]*rhs(2),
                       mData[1]*rhs(0), mData[1]*rhs(1), mData[1]*rhs(2),
                       mData[2]*rhs(0), mData[2]*rhs(1), mData[2]*rhs(2));
}

//------------------------------------------------------------------------------
// Perform the dyad operation with ourself, resulting in a symmetric tensor.
//------------------------------------------------------------------------------
template<>
inline
GeomSymmetricTensor<1>
GeomVector<1>::selfdyad() const {
  return GeomSymmetricTensor<1>((mData[0])*(mData[0]));
}

template<>
inline
GeomSymmetricTensor<2>
GeomVector<2>::selfdyad() const {
  const double a = (mData[0])*(mData[0]);
  const double b = (mData[0])*(mData[1]);
  const double c = (mData[1])*(mData[1]);
  return GeomSymmetricTensor<2>(a, b,
                                b, c);
}

template<>
inline
GeomSymmetricTensor<3>
GeomVector<3>::selfdyad() const {
  const double a = (mData[0])*(mData[0]);
  const double b = (mData[0])*(mData[1]);
  const double c = (mData[0])*(mData[2]);
  const double d = (mData[1])*(mData[1]);
  const double e = (mData[1])*(mData[2]);
  const double f = (mData[2])*(mData[2]);
  return GeomSymmetricTensor<3>(a, b, c,
                                b, d, e,
                                c, e, f);
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
template<>
inline
double
GeomVector<1>::magnitude() const {
  return std::abs(mData[0]);
}

template<>
inline
double
GeomVector<2>::magnitude() const {
  return sqrt((mData[0])*(mData[0]) + (mData[1])*(mData[1]));
}

template<>
inline
double
GeomVector<3>::magnitude() const {
  return sqrt((mData[0])*(mData[0]) + (mData[1])*(mData[1]) + (mData[2])*(mData[2]));
}

//------------------------------------------------------------------------------
// Return the square of the magnitude.
//------------------------------------------------------------------------------
template<>
inline
double
GeomVector<1>::magnitude2() const {
  return (mData[0])*(mData[0]);
}

template<>
inline
double
GeomVector<2>::magnitude2() const {
  return (mData[0])*(mData[0]) + (mData[1])*(mData[1]);
}

template<>
inline
double
GeomVector<3>::magnitude2() const {
  return (mData[0])*(mData[0]) + (mData[1])*(mData[1]) + (mData[2])*(mData[2]);
}

//------------------------------------------------------------------------------
// Return the minimum element.
//------------------------------------------------------------------------------
template<>
inline
double
GeomVector<1>::minElement() const {
  return mData[0];
}

template<>
inline
double
GeomVector<2>::minElement() const {
  return std::min(mData[0], mData[1]);
}

template<>
inline
double
GeomVector<3>::minElement() const {
  return std::min(mData[0], std::min(mData[1], mData[2]));
}

//------------------------------------------------------------------------------
// Return the maximum element.
//------------------------------------------------------------------------------
template<>
inline
double
GeomVector<1>::maxElement() const {
  return mData[0];
}

template<>
inline
double
GeomVector<2>::maxElement() const {
  return std::max(mData[0], mData[1]);
}

template<>
inline
double
GeomVector<3>::maxElement() const {
  return std::max(mData[0], std::max(mData[1], mData[2]));
}

//------------------------------------------------------------------------------
// Return the maximum element by absolute value.
//------------------------------------------------------------------------------
template<>
inline
double
GeomVector<1>::maxAbsElement() const {
  return std::abs(mData[0]);
}

template<>
inline
double
GeomVector<2>::maxAbsElement() const {
  return std::max(std::abs(mData[0]), 
                  std::abs(mData[1]));
}

template<>
inline
double
GeomVector<3>::maxAbsElement() const {
  return std::max(std::abs(mData[0]),
                  std::max(std::abs(mData[1]),
                           std::abs(mData[2])));
}

//------------------------------------------------------------------------------
// Return the sum of the elements.
//------------------------------------------------------------------------------
template<>
inline
double
GeomVector<1>::sumElements() const {
  return mData[0];
}

template<>
inline
double
GeomVector<2>::sumElements() const {
  return mData[0] + mData[1];
}

template<>
inline
double
GeomVector<3>::sumElements() const {
  return mData[0] + mData[1] + mData[2];
}

//******************************************************************************
// Global functions
//******************************************************************************

//------------------------------------------------------------------------------
// Add a vector to a scalar.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
operator+(const double val, const GeomVector<nDim>& vec) {
  return vec + val;
}

//------------------------------------------------------------------------------
// Subtract a vector from a scalar.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim>
operator-(const double val, const GeomVector<nDim>& vec) {
  return -(vec - val);
}

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
