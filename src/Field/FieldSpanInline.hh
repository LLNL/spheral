#include "Field/Field.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/safeInv.hh"
#include "Distributed/allReduce.hh"
#include "Distributed/Communicator.hh"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <limits>

#ifdef USE_MPI
extern "C" {
#include <mpi.h>
}
#endif

// Inlined methods.
namespace Spheral {

//------------------------------------------------------------------------------
// Construct from a Field
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>::
FieldSpan(Field<Dimension, DataType>& field):
  mDataSpan(field.mDataArray),
  mNumInternalElements(field.numInternalElements()),
  mNumGhostElements(field.numGhostElements()) {
}

//------------------------------------------------------------------------------
// Assignment operator with a constant value of DataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto& x: mDataSpan) x = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Element access by integer index.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldSpan<Dimension, DataType>::operator()(size_t index) {
  CHECK2(index < mDataSpan.size(), "FieldSpan index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldSpan<Dimension, DataType>::operator()(size_t index) const {
  CHECK2(index < mDataSpan.size(), "FieldSpan index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with STL interface.
// Since std::span doesn't support at() until C++26, we emulate it with VERIFY
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldSpan<Dimension, DataType>::at(size_t index) {
  VERIFY2(index < mDataSpan.size(), "FieldSpan index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldSpan<Dimension, DataType>::at(size_t index) const {
  VERIFY2(index < mDataSpan.size(), "FieldSpan index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

//------------------------------------------------------------------------------
// Index operators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldSpan<Dimension, DataType>::
operator[](const size_t index) {
  CHECK2(index < mDataSpan.size(), "FieldSpan index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldSpan<Dimension, DataType>::
operator[](const size_t index) const {
  CHECK2(index < mDataSpan.size(), "FieldSpan index out of range: " << index << " " << mDataSpan.size());
  return mDataSpan[index];
}

//------------------------------------------------------------------------------
// Apply a minimum value to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpan<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (auto& x: mDataSpan) x = std::max(x, dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum value to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpan<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (auto& x: mDataSpan) x = std::min(x, dataMax);
}

//------------------------------------------------------------------------------
// Apply a scalar minimum value  to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpan<Dimension, DataType>::applyScalarMin(const Scalar& dataMin) {
  for (auto& x: mDataSpan) x = std::max(x, dataMin);
}

//------------------------------------------------------------------------------
// Apply a scalar maximum value to the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpan<Dimension, DataType>::applyScalarMax(const Scalar& dataMax) {
  for (auto& x: mDataSpan) x = std::min(x, dataMax);
}

//------------------------------------------------------------------------------
// Addition with another FieldSpan in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator+=(const FieldSpan<Dimension, DataType>& rhs) {
  const auto n = this->numElements();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) += rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Subtract another FieldSpan from this one in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator-=(const FieldSpan<Dimension, DataType>& rhs) {
  const auto n = this->numElements();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) -= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Addition with a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator+=(const DataType& rhs) {
  for (auto& x: mDataSpan) x += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator-=(const DataType& rhs) {
  for (auto& x: mDataSpan) x -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar FieldSpan in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator*=(const FieldSpan<Dimension, Scalar>& rhs) {
  const auto n = this->numElements();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) *= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar FieldSpan in place.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator/=(const FieldSpan<Dimension, typename Dimension::Scalar>& rhs) {
  const auto n = this->numElements();
  REQUIRE(rhs.numElements() == n);
  for (auto i = 0u; i < n; ++i) (*this)(i) *= safeInvVar(rhs(i), 1.0e-60);
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator*=(const Scalar& rhs) {
  for (auto& x: mDataSpan) x *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar value in place
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpan<Dimension, DataType>&
FieldSpan<Dimension, DataType>::
operator/=(const Scalar& rhs) {
  const auto rhsInv = safeInvVar(rhs, 1.0e-60);
  for (auto& x: mDataSpan) x *= rhsInv;
  return *this;
}

//------------------------------------------------------------------------------
// Sum the elements of the FieldSpan (assumes the DataType::operator+= is 
// available).
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpan<Dimension, DataType>::
sumElements() const {
  return allReduce(this->localSumElements(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Minimum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpan<Dimension, DataType>::
min() const {
  return allReduce(this->localMin(), SPHERAL_OP_MIN);
}

//------------------------------------------------------------------------------
// Maximum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpan<Dimension, DataType>::
max() const {
  return allReduce(this->localMax(), SPHERAL_OP_MAX);
}

//------------------------------------------------------------------------------
// Sum the elements of the FieldSpan (assumes the DataType::operator+= is 
// available).  
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpan<Dimension, DataType>::
localSumElements() const {
  const auto* start = &mDataSpan.front();
  return std::accumulate(start, start + numInternalElements(), DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Minimum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpan<Dimension, DataType>::
localMin() const {
  const auto* start = &mDataSpan.front();
  return (mDataSpan.empty() ?
          std::numeric_limits<DataType>::max() :
          *std::min_element(start, start + numInternalElements()));
}

//------------------------------------------------------------------------------
// Maximum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpan<Dimension, DataType>::
localMax() const {
  const auto* start = &mDataSpan.front();
  return (mDataSpan.empty() ?
          std::numeric_limits<DataType>::lowest() :
          *std::max_element(start, start + numInternalElements()));
}

//------------------------------------------------------------------------------
// operator==(FieldSpan)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator==(const FieldSpan<Dimension, DataType>& rhs) const {
  const auto n = this->numElements();
  if (rhs.numElements() != n) return false;
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] == rhs[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(FieldSpan)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator!=(const FieldSpan<Dimension, DataType>& rhs) const {
  return !((*this) == rhs);
}

// //------------------------------------------------------------------------------
// // operator>(FieldSpan)
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// bool
// FieldSpan<Dimension, DataType>::
// operator>(const FieldSpan<Dimension, DataType>& rhs) const {
//   const auto n = this->numElements();
//   if (rhs.numElements() != n) return false;
//   auto result = true;
//   size_t i = 0u;
//   while (i < n and result) {
//     result = (*this)[i] > rhs[i];
//     ++i;
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // operator<(FieldSpan)
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// bool
// FieldSpan<Dimension, DataType>::
// operator<(const FieldSpan<Dimension, DataType>& rhs) const {
//   const auto n = this->numElements();
//   if (rhs.numElements() != n) return false;
//   auto result = true;
//   size_t i = 0u;
//   while (i < n and result) {
//     result = (*this)[i] < rhs[i];
//     ++i;
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // operator>=(FieldSpan)
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// bool
// FieldSpan<Dimension, DataType>::
// operator>=(const FieldSpan<Dimension, DataType>& rhs) const {
//   const auto n = this->numElements();
//   if (rhs.numElements() != n) return false;
//   auto result = true;
//   size_t i = 0u;
//   while (i < n and result) {
//     result = (*this)[i] >= rhs[i];
//     ++i;
//   }
// }

// //------------------------------------------------------------------------------
// // operator<=(Field)
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// bool
// FieldSpan<Dimension, DataType>::
// operator<=(const FieldSpan<Dimension, DataType>& rhs) const {
//   const auto n = this->numElements();
//   if (rhs.numElements() != n) return false;
//   auto result = true;
//   size_t i = 0u;
//   while (i < n and result) {
//     result = (*this)[i] <= rhs[i];
//     ++i;
//   }
//   return result;
// }

//------------------------------------------------------------------------------
// operator==(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator==(const DataType& rhs) const {
  const auto n = this->numElements();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] == rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator>(const DataType& rhs) const {
  const auto n = this->numElements();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] > rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator<(const DataType& rhs) const {
  const auto n = this->numElements();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] < rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  const auto n = this->numElements();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] >= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  const auto n = this->numElements();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] <= rhs;
    ++i;
  }
  return result;
}

//****************************** Global Functions ******************************
//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::ostream&
operator<<(std::ostream& os, const FieldSpan<Dimension, DataType>& fieldSpan) {

  // Write the number of internal elements.
  os << fieldSpan.numInternalElements() << " ";

  // Write the internal elements.
  for (auto itr = fieldSpan.internalBegin(); itr < fieldSpan.internalEnd(); ++itr) {
    os << *itr << " ";
  }
//   os << endl;
  return os;
}

} // namespace Spheral
