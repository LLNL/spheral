// Includes.
#include "Geometry/MathTraits.hh"

#include "Field/FieldSpan.hh"
#include "Distributed/allReduce.hh"

#include <algorithm>
#include <limits>

namespace Spheral {

//------------------------------------------------------------------------------
// Construct from a span of FieldSpan
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>::
FieldSpanList(std::span<FieldSpan<Dimension, DataType>>& rhs):
  mFieldSpans(rhs) {
}

//------------------------------------------------------------------------------
// Assignment with a constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto& fspan: mFieldSpans) fspan = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Set the FieldSpans of this FieldSpanList equal to those of another.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpanList<Dimension, DataType>::
assignFields(FieldSpanList<Dimension, DataType>& rhs) {
#pragma omp critical (FieldSpanList_assignFields)
  {
    const auto n = this->size();
    CHECK(rhs.size() == n);
    for (size_t k = 0u; k < n; ++k) {
      CHECK((*this)[k].size() == rhs[k].size());
      (*this)[k] = rhs[k];
    }
  } // OMP critical
}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
const typename FieldSpanList<Dimension, DataType>::value_type&
FieldSpanList<Dimension, DataType>::
operator[](const size_t index) const {
  REQUIRE2(index < this->size(), "FieldSpanList index ERROR: out of bounds " << index << " !< " << this->size());
  return mFieldSpans[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldSpanList<Dimension, DataType>::value_type&
FieldSpanList<Dimension, DataType>::
operator[](const size_t index) {
  REQUIRE2(index < this->size(), "FieldSpanList index ERROR: out of bounds " << index << " !< " << this->size());
  return mFieldSpans[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
const typename FieldSpanList<Dimension, DataType>::value_type&
FieldSpanList<Dimension, DataType>::
at(const size_t index) const {
  return (*this)[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldSpanList<Dimension, DataType>::value_type&
FieldSpanList<Dimension, DataType>::
at(const size_t index) {
  return (*this)[index];
}

//------------------------------------------------------------------------------
// Provide direct access to Field elements
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldSpanList<Dimension, DataType>::
operator()(const size_t fieldIndex,
           const size_t nodeIndex) {
  REQUIRE2(fieldIndex < mFieldSpans.size(), "FieldSpanList index ERROR: out of bounds " << fieldIndex << " !< " << mFieldSpans.size());
  REQUIRE2(nodeIndex < mFieldSpans[fieldIndex].size(), "FieldSpanList node index ERROR: out of bounds " << nodeIndex << " !< " << mFieldSpans[fieldIndex].size());
  return mFieldSpans[fieldIndex][nodeIndex];
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldSpanList<Dimension, DataType>::
operator()(const size_t fieldIndex,
           const size_t nodeIndex) const {
  REQUIRE2(fieldIndex < mFieldSpans.size(), "FieldSpanList index ERROR: out of bounds " << fieldIndex << " !< " << mFieldSpans.size());
  REQUIRE2(nodeIndex < mFieldSpans[fieldIndex].size(), "FieldSpanList node index ERROR: out of bounds " << nodeIndex << " !< " << mFieldSpans[fieldIndex].size());
  return mFieldSpans[fieldIndex][nodeIndex];
}

//------------------------------------------------------------------------------
// Zero out the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpanList<Dimension, DataType>::Zero() {
  for (auto& x: mFieldSpans) x.Zero();
}

//------------------------------------------------------------------------------
// Apply a minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpanList<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (auto& x: mFieldSpans) x.applyMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpanList<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (auto& x: mFieldSpans) x.applyMax(dataMax);
}

//------------------------------------------------------------------------------
// Apply a (scalar) minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpanList<Dimension, DataType>::applyScalarMin(const Scalar dataMin) {
  for (auto& x: mFieldSpans) x.applyScalarMin(dataMin);
}

//------------------------------------------------------------------------------
// Apply a (scalar) maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldSpanList<Dimension, DataType>::applyScalarMax(const Scalar dataMax) {
  for (auto& x: mFieldSpans) x.applyScalarMax(dataMax);
}

//------------------------------------------------------------------------------
// Add two FieldSpanLists in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator+=(const FieldSpanList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] += rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a FieldList from another in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator-=(const FieldSpanList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] -= rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to the Field in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator+=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value from the FieldSpanList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator-=(const DataType& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldSpanList by a Scalar FieldSpanList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::
operator*=(const FieldSpanList<Dimension, typename Dimension::Scalar>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] *= rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldSpanList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator*=(const Scalar& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldSpanList by a Scalar FieldSpanList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::
operator/=(const FieldSpanList<Dimension, typename Dimension::Scalar>& rhs) {

  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  for (size_t i = 0u; i < n; ++i) (*this)[i] /= rhs[i];
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldSpanList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldSpanList<Dimension, DataType>&
FieldSpanList<Dimension, DataType>::operator/=(const typename Dimension::Scalar& rhs) {
  const auto n = this->size();
  for (size_t i = 0u; i < n; ++i) (*this)[i] /= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Sum the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpanList<Dimension, DataType>::
sumElements() const {
  return allReduce(this->localSumElements(), SPHERAL_OP_SUM);
}

//------------------------------------------------------------------------------
// Find the minimum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpanList<Dimension, DataType>::
min() const {
  return allReduce(this->localMin(), SPHERAL_OP_MIN);
}

//------------------------------------------------------------------------------
// Find the maximum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpanList<Dimension, DataType>::
max() const {
  return allReduce(this->localMax(), SPHERAL_OP_MAX);
}

//------------------------------------------------------------------------------
// Sum the field elements.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpanList<Dimension, DataType>::
localSumElements() const {
  auto result = DataTypeTraits<DataType>::zero();
  for (auto& x: mFieldSpans) result += x.localSumElements();
  return result;
}

//------------------------------------------------------------------------------
// Find the minimum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpanList<Dimension, DataType>::
localMin() const {
  auto result = std::numeric_limits<DataType>::max();
  for (auto& x: mFieldSpans) result = std::min(result, x.localMin());
  return result;
}

//------------------------------------------------------------------------------
// Find the maximum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldSpanList<Dimension, DataType>::
localMax() const {
  auto result = std::numeric_limits<DataType>::lowest();
  for (auto& x: mFieldSpans) result = std::max(result, x.localMax());
  return result;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator==(const FieldSpanList<Dimension, DataType>& rhs) const {
  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] == rhs.mFieldSpans[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator!=(const FieldSpanList<Dimension, DataType>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator>(const FieldSpanList<Dimension, DataType>& rhs) const {
  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] > rhs.mFieldSpans[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator<(const FieldSpanList<Dimension, DataType>& rhs) const {
  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] < rhs.mFieldSpans[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator>=(const FieldSpanList<Dimension, DataType>& rhs) const {
  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] >= rhs.mFieldSpans[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator<=(const FieldSpanList<Dimension, DataType>& rhs) const {
  // Pre-conditions.
  const auto n = this->size();
  BEGIN_CONTRACT_SCOPE
  {
    REQUIRE(rhs.size() == n);
    for (size_t i = 0u; i < n; ++i) REQUIRE(mFieldSpans[i].size() == rhs.mFieldSpans[i].size());
  }
  END_CONTRACT_SCOPE

  auto result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] <= rhs.mFieldSpans[i];
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator==(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] == rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator>(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] > rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator<(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] < rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] >= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldSpanList<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  const auto n = this->size();
  bool result = true;
  size_t i = 0u;
  while (result and i < n) {
    result = mFieldSpans[i] <= rhs;
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// numNodes
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
FieldSpanList<Dimension, DataType>::
numElements() const {
  size_t result = 0u;
  for (const auto& x: mFieldSpans) result += x.size();
  return result;
}

//------------------------------------------------------------------------------
// numInternalNodes
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
FieldSpanList<Dimension, DataType>::
numInternalElements() const {
  size_t result = 0u;
  for (const auto& x: mFieldSpans) result += x.numInternalElements();
  return result;
}

//------------------------------------------------------------------------------
// numGhostNodes
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
size_t
FieldSpanList<Dimension, DataType>::
numGhostElements() const {
  size_t result = 0u;
  for (const auto& x: mFieldSpans) result += x.numGhostElements();
  return result;
}

}
