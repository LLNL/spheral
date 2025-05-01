#include "Field/FieldSpanBase.hh"
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
  FieldSpanBase<Dimension>(),
  mDataSpan(field.mDataArray),
  mNumInternalElements(field.numInternalElements()),
  mNumGhostElements(field.numGhostElements()) {
}

//------------------------------------------------------------------------------
// Assignment operator with FieldSpanBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanBase<Dimension>&
FieldSpan<Dimension, DataType>::
operator=(FieldSpanBase<Dimension>& rhs) {
  if (this != &rhs) {
    try {
      const FieldSpan<Dimension, DataType>* rhsPtr = dynamic_cast<const FieldSpan<Dimension, DataType>*>(&rhs);
      CHECK2(rhsPtr != nullptr, "Passed incorrect Field to operator=!");
      mDataSpan = rhsPtr->mDataSpan;
    } catch (const std::bad_cast &) {
      VERIFY2(false, "Attempt to assign a FieldSpan to an incompatible field type.");
    }
  }
  return *this;
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
// Test equivalence with a FieldSpanBase.
//------------------------------------------------------------------------------
// template<typename Value>
// struct CrappyFieldCompareMethod {
//   static bool compare(const std::vector<Value,DataAllocator<Value>>& lhs, 
//                       const std::vector<Value,DataAllocator<Value>>& rhs) {
//     return lhs == rhs;
//   }
// };

template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::operator==(const FieldSpanBase<Dimension>& rhs) const {
  try {
    const FieldSpan<Dimension, DataType>* rhsPtr = dynamic_cast<const FieldSpan<Dimension, DataType>*>(&rhs);
    if (rhsPtr == nullptr) return false;
    const auto n = this->size();
    if (rhsPtr->size() != n) return false;
    size_t i = 0u;
    auto result = true;
    while (i < n and result) {
      result = (*this)[i] == (*rhsPtr)[i];
      i++;
    }
    return result;
  } catch (const std::bad_cast &) {
    return false;
  }
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
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldSpan<Dimension, DataType>::at(size_t index) {
  return mDataSpan.at(index);
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldSpan<Dimension, DataType>::at(size_t index) const {
  return mDataSpan.at(index);
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
// Number of elements in the FieldSpan
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
FieldSpan<Dimension, DataType>::numElements() const {
  return mDataSpan.size();
}

template<typename Dimension, typename DataType>
inline
size_t
FieldSpan<Dimension, DataType>::numInternalElements() const {
  return mNumInternalElements;
}

template<typename Dimension, typename DataType>
inline
size_t
FieldSpan<Dimension, DataType>::numGhostElements() const {
  return mNumGhostElements;
}

template<typename Dimension, typename DataType>
inline
size_t
FieldSpan<Dimension, DataType>::size() const {
  return mDataSpan.size();
}

//------------------------------------------------------------------------------
// Zero out the elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpan<Dimension, DataType>::Zero() {
  std::fill(this->begin(), this->end(), DataTypeTraits<DataType>::zero());
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
  const auto n = this->size();
  if (rhs.size() != n) return false;
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

//------------------------------------------------------------------------------
// operator>(FieldSpan)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator>(const FieldSpan<Dimension, DataType>& rhs) const {
  const auto n = this->size();
  if (rhs.size() != n) return false;
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] > rhs[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<(FieldSpan)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator<(const FieldSpan<Dimension, DataType>& rhs) const {
  const auto n = this->size();
  if (rhs.size() != n) return false;
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] < rhs[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=(FieldSpan)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator>=(const FieldSpan<Dimension, DataType>& rhs) const {
  const auto n = this->size();
  if (rhs.size() != n) return false;
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] >= rhs[i];
    ++i;
  }
}

//------------------------------------------------------------------------------
// operator<=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator<=(const FieldSpan<Dimension, DataType>& rhs) const {
  const auto n = this->size();
  if (rhs.size() != n) return false;
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] <= rhs[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpan<Dimension, DataType>::
operator==(const DataType& rhs) const {
  const auto n = this->size();
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
  const auto n = this->size();
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
  const auto n = this->size();
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
  const auto n = this->size();
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
  const auto n = this->size();
  auto result = true;
  size_t i = 0u;
  while (i < n and result) {
    result = (*this)[i] <= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the span
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldSpan<Dimension, DataType>::iterator
FieldSpan<Dimension, DataType>::begin() {
  return mDataSpan.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the span
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldSpan<Dimension, DataType>::iterator
FieldSpan<Dimension, DataType>::end() {
  return mDataSpan.end();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the internal span values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldSpan<Dimension, DataType>::iterator
FieldSpan<Dimension, DataType>::internalBegin() {
  return mDataSpan.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the internal span values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldSpan<Dimension, DataType>::iterator
FieldSpan<Dimension, DataType>::internalEnd() {
  REQUIRE(mDataSpan.size() >= mNumInternalElements);
  return mDataSpan.begin() + mNumInternalElements;
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the ghost values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldSpan<Dimension, DataType>::iterator
FieldSpan<Dimension, DataType>::ghostBegin() {
  return this->internalEnd();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the ghost values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldSpan<Dimension, DataType>::iterator
FieldSpan<Dimension, DataType>::ghostEnd() {
  return mDataSpan.end();
}

// //------------------------------------------------------------------------------
// // Const_iterator pointing to the beginning of the values.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// typename FieldSpan<Dimension, DataType>::const_iterator
// FieldSpan<Dimension, DataType>::begin() const {
//   return mDataSpan.begin();
// }

// //------------------------------------------------------------------------------
// // Const_iterator pointing to the end of the values.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// typename FieldSpan<Dimension, DataType>::const_iterator
// FieldSpan<Dimension, DataType>::end() const {
//   return mDataSpan.end();
// }

// //------------------------------------------------------------------------------
// // Const_iterator pointing to the beginning of the internal values.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// typename FieldSpan<Dimension, DataType>::const_iterator
// FieldSpan<Dimension, DataType>::internalBegin() const {
//   return mDataSpan.begin();
// }

// //------------------------------------------------------------------------------
// // Const_iterator pointing to the end of the internal values.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// typename FieldSpan<Dimension, DataType>::const_iterator
// FieldSpan<Dimension, DataType>::internalEnd() const {
//   REQUIRE(mDataSpan.size() >= mNumInternalElements);
//   return mDataSpan.begin() + mNumInternalElements;
// }

// //------------------------------------------------------------------------------
// // Const_iterator pointing to the beginning of the ghost values.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// typename FieldSpan<Dimension, DataType>::const_iterator
// FieldSpan<Dimension, DataType>::ghostBegin() const {
//   return this->internalEnd();
// }

// //------------------------------------------------------------------------------
// // Const_iterator pointing to the end of the ghost values.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// typename FieldSpan<Dimension, DataType>::const_iterator
// FieldSpan<Dimension, DataType>::ghostEnd() const {
//   return mDataSpan.end();
// }

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
