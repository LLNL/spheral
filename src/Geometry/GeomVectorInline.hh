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
GeomVector<1, true>::
GeomVector(const double x,
           const double y,
           const double z):
  mData(new double[1]) {
  mData[0] = x;
}

template<>
inline
GeomVector<2, true>::
GeomVector(const double x,
           const double y,
           const double z):
  mData(new double[2]) {
  mData[0] = x;
  mData[1] = y;
}

template<>
inline
GeomVector<3, true>::
GeomVector(const double x,
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
template<bool otherMemory> 
inline
GeomVector<nDim, true>::
GeomVector(const GeomVector<nDim, otherMemory>& vec):
  mData(new double[nDim]) {
  std::copy(vec.mData, vec.mData + nDim, mData);
}

//------------------------------------------------------------------------------
// The assignment operator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
GeomVector<nDim, ownMemory>&
GeomVector<nDim, ownMemory>::operator=(const GeomVector<nDim, otherMemory>& vec) {
  std::copy(vec.mData, vec.mData + nDim, mData);
  return *this;
}

//------------------------------------------------------------------------------
// Set the vector elements to a constant scalar value.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim, ownMemory>&
GeomVector<nDim, ownMemory>::operator=(const double val) {
  std::fill(mData, mData + nDim, val);
  return *this;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomVector<nDim, true>::
~GeomVector() {
  delete [] mData;
}

template<int nDim>
inline
GeomVector<nDim, false>::
~GeomVector() {
}

//------------------------------------------------------------------------------
// Return the (index) element using the parenthesis operator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
double
GeomVector<nDim, ownMemory>::
operator()(typename GeomVector<nDim, ownMemory>::size_type index) const {
  REQUIRE(index < nDim);
  return mData[index];
}

template<int nDim, bool ownMemory>
inline
double&
GeomVector<nDim, ownMemory>::
operator()(typename GeomVector<nDim, ownMemory>::size_type index) {
  REQUIRE(index < nDim);
  return mData[index];
}

//------------------------------------------------------------------------------
// Return the x (first) element.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
double
GeomVector<nDim, ownMemory>::
x() const {
  return mData[0];
}

//------------------------------------------------------------------------------
// Return the y (second) element
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
double
GeomVector<nDim, ownMemory>::
y() const {
  REQUIRE(nDim > 1);
  return mData[1];
}

template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
y() const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Return the z (third) element
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
double
GeomVector<nDim, ownMemory>::
z() const {
  REQUIRE(nDim > 2);
  return mData[2];
}

template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
z() const {
  return 0.0;
}

template<bool ownMemory>
inline
double
GeomVector<2, ownMemory>::
z() const {
  return 0.0;
}

//------------------------------------------------------------------------------
// Set the x (first) element.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
void
GeomVector<nDim, ownMemory>::
x(const double val) {
  mData[0] = val;
}

//------------------------------------------------------------------------------
// Set the y (second) element
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
void
GeomVector<nDim, ownMemory>::
y(const double val) {
  REQUIRE(nDim > 1);
  mData[1] = val;
}

template<bool ownMemory>
inline
void
GeomVector<1, ownMemory>::
y(const double val) {
}

//------------------------------------------------------------------------------
// Set the z (third) element
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
void
GeomVector<nDim, ownMemory>::
z(const double val) {
  REQUIRE(nDim > 2);
  mData[2] = val;
}

template<bool ownMemory>
inline
void
GeomVector<1, ownMemory>::
z(const double val) {
}

template<bool ownMemory>
inline
void
GeomVector<2, ownMemory>::
z(const double val) {
}

//------------------------------------------------------------------------------
// Provide begin/end iterators over the elements of the Vector.
//------------------------------------------------------------------------------
// Non-const versions.
template<int nDim, bool ownMemory>
inline
typename GeomVector<nDim, ownMemory>::iterator
GeomVector<nDim>::
begin() {
  return &mData[0];
}

template<int nDim, bool ownMemory>
inline
typename GeomVector<nDim, ownMemory>::iterator
GeomVector<nDim>::
end() {
  return &mData[nDim];
}

// Const versions.
template<int nDim, bool ownMemory>
inline
typename GeomVector<nDim, ownMemory>::const_iterator
GeomVector<nDim>::
begin() const {
  return &mData[0];
}

template<int nDim, bool ownMemory>
inline
typename GeomVector<nDim, ownMemory>::const_iterator
GeomVector<nDim>::
end() const {
  return &mData[nDim];
}

//------------------------------------------------------------------------------
// Zero out the Vector.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
void
GeomVector<nDim, ownMemory>::
Zero() {
  std::fill(mData, mData + nDim, 0.0);
}

//------------------------------------------------------------------------------
// Return the negative of a vector.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
GeomVector<1>
GeomVector<1, ownMemory>::
operator-() const {
  return GeomVector<1>(-(mData[0]));
}

template<bool ownMemory>
inline
GeomVector<2>
GeomVector<2, ownMemory>::
operator-() const {
  return GeomVector<2>(-(mData[0]), -(mData[1]));
}

template<bool ownMemory>
inline
GeomVector<3>
GeomVector<3, ownMemory>::
operator-() const {
  return GeomVector<3>(-(mData[0]), -(mData[1]), -(mData[2]));
}

//------------------------------------------------------------------------------
// Add two vectors.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
GeomVector<nDim>
GeomVector<nDim, ownMemory>::
operator+(const GeomVector<nDim, otherMemory>& vec) const {
  GeomVector<nDim> result(*this);
  result += vec;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a vector from another.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
GeomVector<nDim>
GeomVector<nDim, ownMemory>::
operator-(const GeomVector<nDim, otherMemory>& vec) const {
  GeomVector<nDim> result(*this);
  result -= vec;
  return result;
}

//------------------------------------------------------------------------------
// Mutiply two vectors.  For our purposes this returns the dyad.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
GeomTensor<nDim>
GeomVector<nDim, ownMemory>::
operator*(const GeomVector<nDim, otherMemory>& vec) const {
  return this->dyad(vec);
}

//------------------------------------------------------------------------------
// Add a scalar to a vector.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
GeomVector<nDim, ownMemory>::
operator+(const double val) const {
  GeomVector<nDim> result(*this);
  result += val;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a scalar from a vector.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
GeomVector<nDim, ownMemory>::
operator-(const double val) const {
  GeomVector<nDim> result(*this);
  result -= val;
  return result;
}

//------------------------------------------------------------------------------
// Multiply a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
GeomVector<nDim, ownMemory>::
operator*(const double val) const {
  GeomVector<nDim> result(*this);
  result *= val;
  return result;
}

//------------------------------------------------------------------------------
// Divide a vector by a scalar
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
GeomVector<nDim, ownMemory>::
operator/(const double val) const {
  CHECK(val != 0.0);
  GeomVector<nDim> result(*this);
  result /= val;
  return result;
}

//------------------------------------------------------------------------------
// Add two vectors in place.
//------------------------------------------------------------------------------
template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<1, ownMemory>&
GeomVector<1, ownMemory>::
operator+=(const GeomVector<1, otherMemory>& vec) {
  mData[0] += vec.mData[0];
  return *this;
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<2, ownMemory>&
GeomVector<2, ownMemory>::
operator+=(const GeomVector<2, otherMemory>& vec) {
  mData[0] += vec.mData[0];
  mData[1] += vec.mData[1];
  return *this;
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<3, ownMemory>&
GeomVector<3, ownMemory>::
operator+=(const GeomVector<3, otherMemory>& vec) {
  mData[0] += vec.mData[0];
  mData[1] += vec.mData[1];
  mData[2] += vec.mData[2];
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a vector in place.
//------------------------------------------------------------------------------
template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<1, bool ownMemory>& 
GeomVector<1, bool ownMemory>::
operator-=(const GeomVector<1, otherMemory>& vec) {
  mData[0] -= vec.mData[0];
  return *this;
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<2, ownMemory>&
GeomVector<2, ownMemory>::
operator-=(const GeomVector<2, otherMemory>& vec) {
  mData[0] -= vec.mData[0];
  mData[1] -= vec.mData[1];
  return *this;
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<3, ownMemory>&
GeomVector<3, ownMemory>::
operator-=(const GeomVector<3, otherMemory>& vec) {
  mData[0] -= vec.mData[0];
  mData[1] -= vec.mData[1];
  mData[2] -= vec.mData[2];
  return *this;
}

//------------------------------------------------------------------------------
// Add a scalar to this vector in place.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
GeomVector<1, bool ownMemory>&
GeomVector<1, bool ownMemory>::
operator+=(const double val) {
  mData[0] += val;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<2, ownMemory>&
GeomVector<2, ownMemory>::
operator+=(const double val) {
  mData[0] += val;
  mData[1] += val;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<3, ownMemory>&
GeomVector<3, ownMemory>::
operator+=(const double val) {
  mData[0] += val;
  mData[1] += val;
  mData[2] += val;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a scalar from this vector in place.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
GeomVector<1, ownMemory>&
GeomVector<1, ownMemory>::
operator-=(const double val) {
  mData[0] -= val;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<2, ownMemory>&
GeomVector<2, ownMemory>::
operator-=(const double val) {
  mData[0] -= val;
  mData[1] -= val;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<3, ownMemory>&
GeomVector<3, ownMemory>::
operator-=(const double val) {
  mData[0] -= val;
  mData[1] -= val;
  mData[2] -= val;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this vector by a scalar in place.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
GeomVector<1, ownMemory>&
GeomVector<1, ownMemory>::
operator*=(const double val) {
  mData[0] *= val;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<2, ownMemory>&
GeomVector<2, ownMemory>::
operator*=(const double val) {
  mData[0] *= val;
  mData[1] *= val;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<3, ownMemory>&
GeomVector<3, ownMemory>::
operator*=(const double val) {
  mData[0] *= val;
  mData[1] *= val;
  mData[2] *= val;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this vector by a scalar in place.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
GeomVector<1, ownMemory>&
GeomVector<1, ownMemory>::
operator/=(const double val) {
  REQUIRE(val != 0.0);
  mData[0] /= val;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<2, ownMemory>&
GeomVector<2, ownMemory>::
operator/=(const double val) {
  REQUIRE(val != 0.0);
  const double valInv = 1.0/val;
  mData[0] *= valInv;
  mData[1] *= valInv;
  return *this;
}

template<bool ownMemory>
inline
GeomVector<3, ownMemory>&
GeomVector<3, ownMemory>::
operator/=(const double val) {
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
template<bool ownMemory>
template<bool otherMemory>
inline
int
GeomVector<1, ownMemory>::
compare(const GeomVector<1, otherMemory>& vec) const {
  return (mData[0] < vec.mData[0] ? -1 :
          mData[0] > vec.mData[0] ?  1 :
          0);
}

template<bool ownMemory>
template<bool otherMemory>
inline
int
GeomVector<2, ownMemory>::
compare(const GeomVector<2, otherMemory>& vec) const {
  return (mData[1] < vec.mData[1] ? -1 :
          mData[1] > vec.mData[1] ?  1 :
          mData[0] < vec.mData[0] ? -1 :
          mData[0] > vec.mData[0] ?  1 :
          0);
}

template<bool ownMemory>
template<bool otherMemory>
inline
int
GeomVector<3, ownMemory>::
compare(const GeomVector<3, otherMemory>& vec) const {
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
template<bool ownMemory>
inline
int
GeomVector<1, ownMemory>::
compare(const double val) const {
  return (mData[0] < val ? -1 :
          mData[0] > val ?  1 :
          0);
}

template<bool ownMemory>
inline
int
GeomVector<2, ownMemory>::
compare(const double val) const {
  return (mData[1] < val ? -1 :
          mData[1] > val ?  1 :
          mData[0] < val ? -1 :
          mData[0] > val ?  1 :
          0);
}

template<bool ownMemory>
inline
int
GeomVector<3, ownMemory>::
compare(const double val) const {
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
template<bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<1, ownMemory>::
operator==(const GeomVector<1, otherMemory>& vec) const {
  return mData[0] == vec.mData[0];
}

template<bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<2, ownMemory>::
operator==(const GeomVector<2, otherMemory>& vec) const {
  return (mData[0] == vec.mData[0]) and (mData[1] == vec.mData[1]);
}

template<bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<3, ownMemory>::
operator==(const GeomVector<3, otherMemory>& vec) const {
  return (mData[0] == vec.mData[0]) and (mData[1] == vec.mData[1]) and (mData[2] == vec.mData[2]);
}

//------------------------------------------------------------------------------
// The non-equivalence comparator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator!=(const GeomVector<nDim, otherMemory>& vec) const {
  return not (*this == vec);
}

//------------------------------------------------------------------------------
// The less than comparator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator<(const GeomVector<nDim, otherMemory>& rhs) const {
  return this->compare(rhs) == -1;
}

//------------------------------------------------------------------------------
// The greater than comparator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator>(const GeomVector<nDim, otherMemory>& rhs) const {
  return this->compare(rhs) == 1;
}

//------------------------------------------------------------------------------
// The less than or equal comparator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator<=(const GeomVector<nDim, otherMemory>& rhs) const {
  return this->compare(rhs) <= 0;
}

//------------------------------------------------------------------------------
// The greater than or equal comparator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
template<bool otherMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator>=(const GeomVector<nDim, otherMemory>& rhs) const {
  return this->compare(rhs) >= 0;
}

//------------------------------------------------------------------------------
// The equivalence comparator (double).
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
bool
GeomVector<1, ownMemory>::
operator==(const double val) const {
  return mData[0] == val;
}

template<bool ownMemory>
inline
bool
GeomVector<2, ownMemory>::
operator==(const double val) const {
  return (mData[0] == val) and (mData[1] == val);
}

template<bool ownMemory>
inline
bool
GeomVector<3, ownMemory>::
operator==(const double val) const {
  return (mData[0] == val) and (mData[1] == val) and (mData[2] == val);
}

//------------------------------------------------------------------------------
// The non-equivalence comparator (double).
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator!=(const double val) const {
  return not (*this == val);
}

//------------------------------------------------------------------------------
// The less than comparator (double).
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator<(const double val) const {
  return this->compare(val) == -1;
}

//------------------------------------------------------------------------------
// The greater than comparator (double).
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator>(const double val) const {
  return this->compare(val) == 1;
}

//------------------------------------------------------------------------------
// The less than or equal comparator (double).
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator<=(const double val) const {
  return this->compare(val) <= 0;
}

//------------------------------------------------------------------------------
// The greater than or equal comparator (double).
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
bool
GeomVector<nDim, ownMemory>::
operator>=(const double val) const {
  return this->compare(val) >= 0;
}

//------------------------------------------------------------------------------
// Dot two vectors.
//------------------------------------------------------------------------------
template<bool ownMemory>
template<bool otherMemory>
inline
double
GeomVector<1, ownMemory>::
dot(const GeomVector<1, otherMemory>& vec) const {
  return mData[0]*vec.mData[0];
}

template<bool ownMemory>
template<bool otherMemory>
inline
double
GeomVector<2, ownMemory>::
dot(const GeomVector<2, otherMemory>& vec) const {
  return mData[0]*vec.mData[0] + mData[1]*vec.mData[1];
}

template<bool ownMemory>
template<bool otherMemory>
inline
double
GeomVector<3, ownMemory>::
dot(const GeomVector<3, otherMemory>& vec) const {
  return mData[0]*vec.mData[0] + mData[1]*vec.mData[1] + mData[2]*vec.mData[2];
}

//------------------------------------------------------------------------------
// Cross two vectors.
//------------------------------------------------------------------------------
template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<3>
GeomVector<3, ownMemory>::
cross(const GeomVector<3, otherMemory>& vec) const {
  return GeomVector<3>(mData[1]*vec.mData[2] - mData[2]*vec.mData[1],
		       mData[2]*vec.mData[0] - mData[0]*vec.mData[2],
		       mData[0]*vec.mData[1] - mData[1]*vec.mData[0]);
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<3>
GeomVector<2, ownMemory>::
cross(const GeomVector<2, otherMemory>& vec) const {
  return GeomVector<3>(0.0,
                       0.0,
		       mData[0]*vec.mData[1] - mData[1]*vec.mData[0]);
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomVector<3>
GeomVector<1, ownMemory>::
cross(const GeomVector<1, otherMemory>& vec) const {
  return GeomVector<3>(0.0, 0.0, 0.0);
}

//------------------------------------------------------------------------------
// Perform the dyad operation with another vector.
//------------------------------------------------------------------------------
template<bool ownMemory>
template<bool otherMemory>
inline
GeomTensor<1>
GeomVector<1, ownMemory>::
dyad(const GeomVector<1, otherMemory>& vec) const {
  return GeomTensor<1>(mData[0]*rhs(0));
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomTensor<2>
GeomVector<2, ownMemory>::
dyad(const GeomVector<2, otherMemory>& rhs) const {
  return GeomTensor<2>(mData[0]*rhs(0), mData[0]*rhs(1),
                       mData[1]*rhs(0), mData[1]*rhs(1));
}

template<bool ownMemory>
template<bool otherMemory>
inline
GeomTensor<3>
GeomVector<3, ownMemory>::
dyad(const GeomVector<3, otherMemory>& rhs) const {
  return GeomTensor<3>(mData[0]*rhs(0), mData[0]*rhs(1), mData[0]*rhs(2),
                       mData[1]*rhs(0), mData[1]*rhs(1), mData[1]*rhs(2),
                       mData[2]*rhs(0), mData[2]*rhs(1), mData[2]*rhs(2));
}

//------------------------------------------------------------------------------
// Perform the dyad operation with ourself, resulting in a symmetric tensor.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
GeomSymmetricTensor<1>
GeomVector<1, ownMemory>::
selfdyad() const {
  return GeomSymmetricTensor<1>((mData[0])*(mData[0]));
}

template<bool ownMemory>
inline
GeomSymmetricTensor<2>
GeomVector<2, ownMemory>::
selfdyad() const {
  const double a = (mData[0])*(mData[0]);
  const double b = (mData[0])*(mData[1]);
  const double c = (mData[1])*(mData[1]);
  return GeomSymmetricTensor<2>(a, b,
                                b, c);
}

template<bool ownMemory>
inline
GeomSymmetricTensor<3>
GeomVector<3, ownMemory>::
selfdyad() const {
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
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
GeomVector<nDim, ownMemory>::
unitVector() const {
  const double mag = this->magnitude();
  return mag > 1.0e-50 ? (*this)/mag : GeomVector<nDim>(1.0, 0.0, 0.0);
}

//------------------------------------------------------------------------------
// Return the magnitude of the Vector.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
magnitude() const {
  return std::abs(mData[0]);
}

template<bool ownMemory>
inline
double
GeomVector<2, ownMemory>::
magnitude() const {
  return sqrt((mData[0])*(mData[0]) + (mData[1])*(mData[1]));
}

template<bool ownMemory>
inline
double
GeomVector<3, ownMemory>::
magnitude() const {
  return sqrt((mData[0])*(mData[0]) + (mData[1])*(mData[1]) + (mData[2])*(mData[2]));
}

//------------------------------------------------------------------------------
// Return the square of the magnitude.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
magnitude2() const {
  return (mData[0])*(mData[0]);
}

template<bool ownMemory>
inline
double
GeomVector<2, ownMemory>::
magnitude2() const {
  return (mData[0])*(mData[0]) + (mData[1])*(mData[1]);
}

template<bool ownMemory>
inline
double
GeomVector<3, ownMemory>::
magnitude2() const {
  return (mData[0])*(mData[0]) + (mData[1])*(mData[1]) + (mData[2])*(mData[2]);
}

//------------------------------------------------------------------------------
// Return the minimum element.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
minElement() const {
  return mData[0];
}

template<bool ownMemory>
inline
double
GeomVector<2, ownMemory>::
minElement() const {
  return std::min(mData[0], mData[1]);
}

template<bool ownMemory>
inline
double
GeomVector<3, ownMemory>::
minElement() const {
  return std::min(mData[0], std::min(mData[1], mData[2]));
}

//------------------------------------------------------------------------------
// Return the maximum element.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
maxElement() const {
  return mData[0];
}

template<bool ownMemory>
inline
double
GeomVector<2, ownMemory>::
maxElement() const {
  return std::max(mData[0], mData[1]);
}

template<bool ownMemory>
inline
double
GeomVector<3, ownMemory>::
maxElement() const {
  return std::max(mData[0], std::max(mData[1], mData[2]));
}

//------------------------------------------------------------------------------
// Return the maximum element by absolute value.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
maxAbsElement() const {
  return std::abs(mData[0]);
}

template<bool ownMemory>
inline
double
GeomVector<2, ownMemory>::
maxAbsElement() const {
  return std::max(std::abs(mData[0]), 
                  std::abs(mData[1]));
}

template<bool ownMemory>
inline
double
GeomVector<3, ownMemory>::
maxAbsElement() const {
  return std::max(std::abs(mData[0]),
                  std::max(std::abs(mData[1]),
                           std::abs(mData[2])));
}

//------------------------------------------------------------------------------
// Return the sum of the elements.
//------------------------------------------------------------------------------
template<bool ownMemory>
inline
double
GeomVector<1, ownMemory>::
sumElements() const {
  return mData[0];
}

template<bool ownMemory>
inline
double
GeomVector<2, ownMemory>::
sumElements() const {
  return mData[0] + mData[1];
}

template<bool ownMemory>
inline
double
GeomVector<3, ownMemory>::
sumElements() const {
  return mData[0] + mData[1] + mData[2];
}

//******************************************************************************
// Global functions
//******************************************************************************

//------------------------------------------------------------------------------
// Add a vector to a scalar.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
operator+(const double val, const GeomVector<nDim, ownMemory>& vec) {
  return vec + val;
}

//------------------------------------------------------------------------------
// Subtract a vector from a scalar.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
operator-(const double val, const GeomVector<nDim, ownMemory>& vec) {
  return -(vec - val);
}

//------------------------------------------------------------------------------
// Multiply a scalar by a vector.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
GeomVector<nDim>
operator*(const double val, const GeomVector<nDim, ownMemory>& vec) {
  return vec*val;
}

//------------------------------------------------------------------------------
// Element wise minimum comparison.
//------------------------------------------------------------------------------
template<bool ownMemory, bool otherMemory>
inline
GeomVector<1>
elementWiseMin(const GeomVector<1, ownMemory>& lhs, 
               const GeomVector<1, otherMemory>& rhs) {
  return GeomVector<1>(std::min(lhs.x(), rhs.x()));
}

template<bool ownMemory, bool otherMemory>
inline
GeomVector<2>
elementWiseMin(const GeomVector<2, ownMemory>& lhs, 
               const GeomVector<2, otherMemory>& rhs) {
  return GeomVector<2>(std::min(lhs.x(), rhs.x()),
                       std::min(lhs.y(), rhs.y()));
}

template<bool ownMemory, bool otherMemory>
inline
GeomVector<3>
elementWiseMin(const GeomVector<3, ownMemory>& lhs, 
               const GeomVector<3, otherMemory>& rhs) {
  return GeomVector<3>(std::min(lhs.x(), rhs.x()),
                       std::min(lhs.y(), rhs.y()),
                       std::min(lhs.z(), rhs.z()));
}

//------------------------------------------------------------------------------
// Element wise maximum comparison.
//------------------------------------------------------------------------------
template<bool ownMemory, bool otherMemory>
inline
GeomVector<1>
elementWiseMax(const GeomVector<1, ownMemory>& lhs, 
               const GeomVector<1, otherMemory>& rhs) {
  return GeomVector<1>(std::max(lhs.x(), rhs.x()));
}

template<bool ownMemory, bool otherMemory>
inline
GeomVector<2>
elementWiseMax(const GeomVector<2, ownMemory>& lhs, 
               const GeomVector<2, otherMemory>& rhs) {
  return GeomVector<2>(std::max(lhs.x(), rhs.x()),
                       std::max(lhs.y(), rhs.y()));
}

template<bool ownMemory, bool otherMemory>
inline
GeomVector<3>
elementWiseMax(const GeomVector<3, ownMemory>& lhs, 
               const GeomVector<3, otherMemory>& rhs) {
  return GeomVector<3>(std::max(lhs.x(), rhs.x()),
                       std::max(lhs.y(), rhs.y()),
                       std::max(lhs.z(), rhs.z()));
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<int nDim, bool ownMemory>
inline
std::istream&
operator>>(std::istream& is, GeomVector<nDim, ownMemory>& vec) {
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
template<int nDim, bool ownMemory>
inline
std::ostream&
operator<<(std::ostream& os, const GeomVector<nDim, ownMemory>& vec) {
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
