#include <algorithm>
#include <limits.h>
#include <string>

#include "GeomVector.hh"
#include "GeomTensor.hh"
#include "GeomSymmetricTensor.hh"
#include "GeomThirdRankTensor.hh"
#include "Infrastructure/SpheralError.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
GeomThirdRankTensor():
  mElements(new double[numElements]) {
  std::fill(mElements, mElements + numElements, 0.0);
}

//------------------------------------------------------------------------------
// Construct with the given value filling the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
GeomThirdRankTensor(const double val):
  mElements(new double[numElements]) {
  std::fill(mElements, mElements + numElements, val);
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
GeomThirdRankTensor(const GeomThirdRankTensor& ten):
  mElements(new double[numElements]) {
  std::copy(ten.begin(), ten.begin() + numElements, this->begin());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>::
~GeomThirdRankTensor() {
  delete [] mElements;
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator=(const GeomThirdRankTensor& rhs) {
  if (this != &rhs) std::copy(rhs.begin(), rhs.begin() + numElements, this->begin());
  return *this;
}

//------------------------------------------------------------------------------
// Assignment (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator=(const double rhs) {
  std::fill(mElements, mElements + numElements, rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Access the elements by indicies.
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomThirdRankTensor<nDim>::
operator()(const GeomThirdRankTensor::size_type i,
           const GeomThirdRankTensor::size_type j,
           const GeomThirdRankTensor::size_type k) const {
  return mElements[elementIndex(i, j, k)];
}

template<int nDim>
inline
double&
GeomThirdRankTensor<nDim>::
operator()(const GeomThirdRankTensor::size_type i,
           const GeomThirdRankTensor::size_type j,
           const GeomThirdRankTensor::size_type k) {
  return mElements[elementIndex(i, j, k)];
}

//------------------------------------------------------------------------------
// Iterators.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomThirdRankTensor<nDim>::iterator
GeomThirdRankTensor<nDim>::
begin() {
  return mElements;
}

template<int nDim>
inline
typename GeomThirdRankTensor<nDim>::iterator
GeomThirdRankTensor<nDim>::
end() {
  return mElements + numElements;
}

template<int nDim>
inline
typename GeomThirdRankTensor<nDim>::const_iterator
GeomThirdRankTensor<nDim>::
begin() const {
  return mElements;
}

template<int nDim>
inline
typename GeomThirdRankTensor<nDim>::const_iterator
GeomThirdRankTensor<nDim>::
end() const {
  return mElements + numElements;
}

//------------------------------------------------------------------------------
// Zero out the tensor.
//------------------------------------------------------------------------------
template<int nDim>
inline
void
GeomThirdRankTensor<nDim>::
Zero() {
  std::fill(mElements, mElements + numElements, 0.0);
}

//------------------------------------------------------------------------------
// Negative.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
GeomThirdRankTensor<nDim>::
operator-() const {
  GeomThirdRankTensor<nDim> result;
  for (size_type i = 0; i != numElements; ++i) result.mElements[i] = -mElements[i];
  return result;
}

//------------------------------------------------------------------------------
// In place addition.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator+=(const GeomThirdRankTensor& rhs) {
  for (size_type i = 0; i != numElements; ++i) mElements[i] += rhs.mElements[i];
  return *this;
}

//------------------------------------------------------------------------------
// In place subtraction.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator-=(const GeomThirdRankTensor& rhs) {
  for (size_type i = 0; i != numElements; ++i) mElements[i] -= rhs.mElements[i];
  return *this;
}

//------------------------------------------------------------------------------
// Addition.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
GeomThirdRankTensor<nDim>::
operator+(const GeomThirdRankTensor& rhs) const {
  GeomThirdRankTensor<nDim> result(*this);
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtraction.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
GeomThirdRankTensor<nDim>::
operator-(const GeomThirdRankTensor& rhs) const {
  GeomThirdRankTensor<nDim> result(*this);
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// In place multiplication (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator*=(const double rhs) {
  for (size_type i = 0; i != numElements; ++i) mElements[i] *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// In place division (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>&
GeomThirdRankTensor<nDim>::
operator/=(const double rhs) {
  REQUIRE(rhs != 0.0);
  for (size_type i = 0; i != numElements; ++i) mElements[i] /= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
GeomThirdRankTensor<nDim>::
operator*(const double rhs) const {
  GeomThirdRankTensor<nDim> result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division (double).
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
GeomThirdRankTensor<nDim>::
operator/(const double rhs) const {
  REQUIRE(rhs != 0.0);
  GeomThirdRankTensor<nDim> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// ==
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator==(const GeomThirdRankTensor& rhs) const {
  bool result = mElements[0] == rhs.mElements[0];
  size_type i = 1;
  while (i != numElements and result) {
    result = result and (mElements[i] == rhs.mElements[i]);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// !=
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator!=(const GeomThirdRankTensor& rhs) const {
  return not this->operator==(rhs);
}

//------------------------------------------------------------------------------
// <
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator<(const GeomThirdRankTensor& rhs) const {
  return this->doubledot(*this) < rhs.doubledot(rhs);
}

//------------------------------------------------------------------------------
// >
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator>(const GeomThirdRankTensor& rhs) const {
  return this->doubledot(*this) > rhs.doubledot(rhs);
}

//------------------------------------------------------------------------------
// <=
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator<=(const GeomThirdRankTensor& rhs) const {
  return this->operator<(rhs) or this->operator==(rhs);
}

//------------------------------------------------------------------------------
// >=
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator>=(const GeomThirdRankTensor& rhs) const {
  return this->operator>(rhs) or this->operator==(rhs);
}

//------------------------------------------------------------------------------
// == (double)
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator==(const double rhs) const {
  bool result = mElements[0] == rhs;
  size_type i = 1;
  while (i != numElements and result) {
    result = result and (mElements[i] == rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// != (double)
//------------------------------------------------------------------------------
template<int nDim>
inline
bool
GeomThirdRankTensor<nDim>::
operator!=(const double rhs) const {
  return not this->operator==(rhs);
}

//------------------------------------------------------------------------------
// doubledot
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomThirdRankTensor<nDim>::
doubledot(const GeomThirdRankTensor& rhs) const {
  double result = mElements[0] * rhs.mElements[0];
  for (size_type i = 1; i != rhs.numElements; ++i) 
    result += mElements[i] * rhs.mElements[i];
  return result;
}

//------------------------------------------------------------------------------
// squareElements
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
GeomThirdRankTensor<nDim>::
squareElements() const {
  GeomThirdRankTensor result(*this);
  for (size_type i = 1; i != this->numElements; ++i) 
    result.mElements[i] *= mElements[i];
  return result;
}

//------------------------------------------------------------------------------
// maxAbsElement
//------------------------------------------------------------------------------
template<int nDim>
inline
double
GeomThirdRankTensor<nDim>::
maxAbsElement() const {
  double result = mElements[0];
  for (size_type i = 1; i != numElements; ++i) result = std::max(result, mElements[i]);
  return result;
}

//------------------------------------------------------------------------------
// Return the element index corresponding to the given (i,j,k) triple.
//------------------------------------------------------------------------------
template<int nDim>
inline
typename GeomThirdRankTensor<nDim>::size_type
GeomThirdRankTensor<nDim>::
elementIndex(const typename GeomThirdRankTensor<nDim>::size_type i,
             const typename GeomThirdRankTensor<nDim>::size_type j,
             const typename GeomThirdRankTensor<nDim>::size_type k) const {
  REQUIRE(i < nDim);
  REQUIRE(j < nDim);
  REQUIRE(k < nDim);
  const size_type result = i*nDim2 + j*nDim + k;
  ENSURE(result < numElements);
  return result;
}

//******************************************************************************
// Global methods.  (unbound to class)
//******************************************************************************

//------------------------------------------------------------------------------
// Multiplication with a double.
//------------------------------------------------------------------------------
template<int nDim>
inline
GeomThirdRankTensor<nDim>
operator*(const double lhs,
          const GeomThirdRankTensor<nDim>& rhs) {
  GeomThirdRankTensor<nDim> result(rhs);
  result *= lhs;
  return result;
}

//------------------------------------------------------------------------------
// operator>> (input)
//------------------------------------------------------------------------------
template<int nDim>
inline
std::istream&
operator>>(std::istream& is, GeomThirdRankTensor<nDim>& ten) {
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
// operator<< (output)
//------------------------------------------------------------------------------
template<int nDim>
inline
std::ostream&
operator<<(std::ostream& os, const GeomThirdRankTensor<nDim>& ten) {
  os << "( ";
  for (typename GeomTensor<nDim>::const_iterator itr = ten.begin();
       itr != ten.end(); ++itr) {
    os << *itr << " ";
  }
  os << ")";
  return os;
}

}
