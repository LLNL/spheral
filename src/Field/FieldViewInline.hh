#include "Field/FieldBase.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/allReduce.hh"
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
// Default Constructor
// node list.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::FieldView():
  mDataArray() {}

//------------------------------------------------------------------------------
// Field Constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::FieldView(FieldType const& field):
  mDataArray(field.mDataArray),
  mFieldPtr(const_cast<FieldType*>(&field)) {}

//------------------------------------------------------------------------------
// Copy Constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::FieldView(const FieldView& field):
  mDataArray(field.mDataArray),
  mFieldPtr(field.mFieldPtr) {}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>::~FieldView() {}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FieldType&
FieldView<Dimension, DataType>::
operator*() const{
  return *mFieldPtr;
}

template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::FieldType*
FieldView<Dimension, DataType>::
operator->() const {
  return mFieldPtr;
}

//------------------------------------------------------------------------------
// Assignment operator to another FieldView.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::operator=(const FieldView<Dimension, DataType>& rhs) {
  if (this != &rhs) {
    mDataArray = rhs.mDataArray;
    mFieldPtr = rhs.mFieldPtr;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assigment operator with a ManagedVector.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::operator=(const ContainerType& rhs) {
  mDataArray = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Assignment operator with a constant value of DataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::operator=(const DataType& rhs) {
  std::fill(begin(), end(), rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Element access by integer index.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldView<Dimension, DataType>::operator()(int index) {
  CHECK(index >= 0 && index < (int)numElements());
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldView<Dimension, DataType>::operator()(int index) const {
  CHECK(index >= 0 && index < (int)numElements());
  return mDataArray[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with STL interface.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldView<Dimension, DataType>::at(int index) {
  return this->operator()(index);
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldView<Dimension, DataType>::at(int index) const {
  return this->operator()(index);
}

//------------------------------------------------------------------------------
// Number of elements in the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
unsigned 
FieldView<Dimension, DataType>::numElements() const {
  return mDataArray.size();
}

template<typename Dimension, typename DataType>
inline
unsigned 
FieldView<Dimension, DataType>::size() const {
  return numElements();
}

//------------------------------------------------------------------------------
// Zero out the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldView<Dimension, DataType>::Zero() {
  std::fill(begin(), end(), DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Apply a minimum value to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldView<Dimension, DataType>::applyMin(const DataType& dataMin) {
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) mDataArray[i] = std::max(mDataArray[i], dataMin);
}

//------------------------------------------------------------------------------
// Apply a maximum value to the FieldView elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldView<Dimension, DataType>::applyMax(const DataType& dataMax) {
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) mDataArray[i] = std::min(mDataArray[i], dataMax);
}

//------------------------------------------------------------------------------
// Apply a (double) minimum value  to the FieldView elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldView<Dimension, DataType>::applyScalarMin(const double dataMin) {
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) mDataArray[i] = std::max(mDataArray[i], dataMin);
}

//------------------------------------------------------------------------------
// Apply a (double) maximum value to the FieldView elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldView<Dimension, DataType>::applyScalarMax(const double dataMax) {
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) mDataArray[i] = std::min(mDataArray[i], dataMax);
}

//------------------------------------------------------------------------------
// Addition with another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::operator+(const FieldView<Dimension, DataType>& rhs) const {
  REQUIRE(this->size() == rhs.size());
  FieldView<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) += rhs(i);
  return result;
}

//------------------------------------------------------------------------------
// Subtract another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::operator-(const FieldView<Dimension, DataType>& rhs) const {
  REQUIRE(this->size() == rhs.size());
  FieldView<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) -= rhs(i);
  return result;
}

//------------------------------------------------------------------------------
// Addition with another field in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::operator+=(const FieldView<Dimension, DataType>& rhs) {
  REQUIRE(this->size() == rhs.size());
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) (*this)(i) += rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Subtract another field from this one in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::operator-=(const FieldView<Dimension, DataType>& rhs) {
  REQUIRE(this->size() == rhs.size());
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) (*this)(i) -= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::operator+(const DataType& rhs) const {
  FieldView<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a single value from every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::operator-(const DataType& rhs) const {
  FieldView<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Addition with a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::operator+=(const DataType& rhs) {
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) (*this)(i) += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::operator-=(const DataType& rhs) {
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) (*this)(i) -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar FieldView
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::
operator*(const FieldView<Dimension, Scalar>& rhs) const {
  REQUIRE(this->size() == rhs.size());
  FieldView<Dimension, DataType> result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division by a Scalar FieldView
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::
operator/(const FieldView<Dimension, Scalar>& rhs) const {
  REQUIRE(this->size() == rhs.size());
  FieldView<Dimension, DataType> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar FieldView in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator*=(const FieldView<Dimension, Scalar>& rhs) {
  REQUIRE(this->size() == rhs.size());
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) (*this)(i) *= rhs(i);
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar FieldView in place.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator/=(const FieldView<Dimension, typename Dimension::Scalar>& rhs) {
  REQUIRE(this->size() == rhs.size());
  const unsigned n = this->numElements();
  for (auto i = 0u; i < n; ++i) {
    (*this)(i) *= safeInvVar(rhs(i), 1.0e-60);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::
operator*(const Scalar& rhs) const {
  FieldView<Dimension, DataType> result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division by a Scalar
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>
FieldView<Dimension, DataType>::
operator/(const Scalar& rhs) const {
  REQUIRE(rhs != 0.0);
  FieldView<Dimension, DataType> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator*=(const Scalar& rhs) {
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) (*this)(i) *= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar value in place
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldView<Dimension, DataType>&
FieldView<Dimension, DataType>::
operator/=(const Scalar& rhs) {
  REQUIRE(rhs != 0.0);
  const unsigned n = this->numElements();
  for (int i = 0; i < (int)n; ++i) (*this)(i) /= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Sum the elements of the field (assumes the DataType::operator+= is 
// available).
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldView<Dimension, DataType>::
sumElements() const {
  return allReduce(this->localSumElements(), MPI_SUM, Communicator::communicator());
}

//------------------------------------------------------------------------------
// Minimum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldView<Dimension, DataType>::
min() const {
  return allReduce(this->localMin(), MPI_MIN, Communicator::communicator());
}

//------------------------------------------------------------------------------
// Maximum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldView<Dimension, DataType>::
max() const {
  return allReduce(this->localMax(), MPI_MAX, Communicator::communicator());
}

//------------------------------------------------------------------------------
// Sum the elements of the field (assumes the DataType::operator+= is 
// available).  
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldView<Dimension, DataType>::
localSumElements() const {
  return std::accumulate(begin(), end(), DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Minimum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldView<Dimension, DataType>::
localMin() const {
  DataType result;
  if (size() == 0) {
    result = std::numeric_limits<DataType>::max();
  } else {
    result = *std::min_element(begin(), end());
  }
  return result;
}

//------------------------------------------------------------------------------
// Maximum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldView<Dimension, DataType>::
localMax() const {
  DataType result;
  if (size() == 0) {
    result = std::numeric_limits<DataType>::lowest();
  } else {
    result = *std::max_element(begin(), end());
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator==(const FieldView<Dimension, DataType>& rhs) const {
  if (&rhs == this) return true;
  //return mDataArray == rhs.mDataArray;
  return mFieldPtr == rhs.mFieldPtr;
}

//------------------------------------------------------------------------------
// operator!=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator!=(const FieldView<Dimension, DataType>& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(FieldView)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator>(const FieldView<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr > *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<(FieldView)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator<(const FieldView<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr < *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=(FieldView)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator>=(const FieldView<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr >= *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=(FieldView)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator<=(const FieldView<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr <= *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator==(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr == rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator>(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr > rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator<(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr < rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr >= rhs;
    ++lhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<=(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldView<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  bool result = true;
  const_iterator lhsItr = this->begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr <= rhs;
    ++lhsItr;
  }
  return result;
}


//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::iterator
FieldView<Dimension, DataType>::begin() {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::iterator
FieldView<Dimension, DataType>::end() {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::const_iterator
FieldView<Dimension, DataType>::begin() const {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldView<Dimension, DataType>::const_iterator
FieldView<Dimension, DataType>::end() const {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Index operators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldView<Dimension, DataType>::
operator[](const unsigned int index) {
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
DataType&
FieldView<Dimension, DataType>::
operator[](const unsigned int index) const {
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
void
FieldView<Dimension, DataType>::
move (chai::ExecutionSpace space, bool touch) const {
  mDataArray.move(space, touch);
}

//****************************** Global Functions ******************************
//------------------------------------------------------------------------------
// Multiplication by another FieldView.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
FieldView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const FieldView<Dimension, DataType>& lhs,
          const FieldView<Dimension, OtherDataType>& rhs) {
  CHECK(lhs.size() == rhs.size());
  FieldView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<FieldView<Dimension, DataType>&>(lhs).nodeList());
  for (auto i = 0u; i < result.numElements(); ++i) {
    result(i) = lhs(i) * rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a single value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
FieldView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const FieldView<Dimension, DataType>& lhs,
          const OtherDataType& rhs) {
  FieldView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<FieldView<Dimension, DataType>&>(lhs).nodeList());
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = lhs(i) * rhs;
  }
  return result;
}

template<typename Dimension, typename DataType, typename OtherDataType>
FieldView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const DataType& lhs,
          const FieldView<Dimension, OtherDataType>& rhs) {
  FieldView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<FieldView<Dimension, OtherDataType>&>(rhs).nodeList());
  for (auto i = 0u; i < result.numElements(); ++i) {
    result(i) = lhs * rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::ostream&
operator<<(std::ostream& os, const FieldView<Dimension, DataType>& field) {

  // Write the number of internal elements.
  os << field.size() << " ";

  // Write the internal elements.
  for (typename FieldView<Dimension, DataType>::const_iterator itr = field.begin();
       itr < field.end();
       ++itr) {
    os << *itr << " ";
  }
//   os << endl;
  return os;
}




} // namespace Spheral
