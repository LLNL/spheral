#include <cmath>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <limits>

#include "Field/FieldBase.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/packElement.hh"
#include "Utilities/removeElements.hh"
#include "Distributed/Communicator.hh"

#ifdef USE_MPI
extern "C" {
#include "mpi.h"
}
#endif

// Inlined methods.
namespace Spheral {
namespace FieldSpace {

//------------------------------------------------------------------------------
// Construct with name.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name):
  FieldBase<Dimension>(name),
  mDataArray(),
  mValid(false) {}

//------------------------------------------------------------------------------
// Construct with name and field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const Field<Dimension, DataType>& field):
  FieldBase<Dimension>(name, *field.nodeListPtr()),
  mDataArray(field.mDataArray),
  mValid(field.mValid) {}

//------------------------------------------------------------------------------
// Construct with the given name and NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const NodeSpace::NodeList<Dimension>& nodeList):
  FieldBase<Dimension>(name, nodeList),
  mDataArray((size_t) nodeList.numNodes(), DataType()),
  mValid(true) {
  REQUIRE(numElements() == nodeList.numNodes());
}

template<>
inline
Field<Dim<1>, Dim<1>::Scalar>::
Field(FieldBase<Dim<1> >::FieldName name,
      const NodeSpace::NodeList<Dim<1> >& nodeList):
  FieldBase<Dim<1> >(name, nodeList),
  mDataArray((size_t) nodeList.numNodes(), 0.0),
  mValid(true) {
  REQUIRE(numElements() == nodeList.numNodes());
}

template<>
inline
Field<Dim<2>, Dim<2>::Scalar>::
Field(FieldBase<Dim<2> >::FieldName name,
      const NodeSpace::NodeList<Dim<2> >& nodeList):
  FieldBase<Dim<2> >(name, nodeList),
  mDataArray((size_t) nodeList.numNodes(), 0.0),
  mValid(true) {
  REQUIRE(numElements() == nodeList.numNodes());
}

template<>
inline
Field<Dim<3>, Dim<3>::Scalar>::
Field(FieldBase<Dim<3> >::FieldName name,
      const NodeSpace::NodeList<Dim<3> >& nodeList):
  FieldBase<Dim<3> >(name, nodeList),
  mDataArray((size_t) nodeList.numNodes(), 0.0),
  mValid(true) {
  REQUIRE(numElements() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Construct with given name, NodeList, and value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const NodeSpace::NodeList<Dimension>& nodeList,
      DataType value):
  FieldBase<Dimension>(name, nodeList),
  mDataArray((size_t) nodeList.numNodes(), value),
  mValid(true) {
  REQUIRE(numElements() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Construct for a given name, NodeList, and vector of values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name, 
      const NodeSpace::NodeList<Dimension>& nodeList,
      const std::vector<DataType>& array):
  FieldBase<Dimension>(name, nodeList),
  mDataArray((size_t) nodeList.numNodes()),
  mValid(true) {
  REQUIRE(numElements() == nodeList.numNodes());
  REQUIRE(numElements() == array.size());
  mDataArray = array;
}

//------------------------------------------------------------------------------
// Construct by copying the values of another Field, but using a different
// node list.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::Field(const NodeSpace::NodeList<Dimension>& nodeList,
                                  const Field<Dimension, DataType>& field):
  FieldBase<Dimension>(field.name(), nodeList),
  mDataArray(field.mDataArray),
  mValid(true) {
  ENSURE(numElements() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Copy Constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::Field(const Field& field):
  FieldBase<Dimension>(const_cast<Field<Dimension, DataType>&>(field)),
  mDataArray(field.mDataArray),
  mValid(field.valid()) {
}

//------------------------------------------------------------------------------
// The virtual clone method, allowing us to duplicate fields with just 
// FieldBase*.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
boost::shared_ptr<FieldBase<Dimension> >
Field<Dimension, DataType>::clone() const {
  return boost::shared_ptr<FieldBase<Dimension> >(new Field<Dimension, DataType>(*this));
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::~Field() {
}

//------------------------------------------------------------------------------
// Assignment operator with FieldBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldBase<Dimension>&
Field<Dimension, DataType>::operator=(const FieldBase<Dimension>& rhs) {
  if (this != &rhs) {
    try {
      const Field<Dimension, DataType>* rhsPtr = dynamic_cast<const Field<Dimension, DataType>*>(&rhs);
      CHECK2(rhsPtr != 0, "Passed incorrect Field to operator=!");
      FieldBase<Dimension>::operator=(rhs);
      mDataArray = rhsPtr->mDataArray;
      mValid = rhsPtr->mValid;
    } catch (std::bad_cast) {
      VERIFY2(false, "Attempt to assign a field to an incompatible field type.");
    }
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assignment operator to another Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const Field<Dimension, DataType>& rhs) {
  REQUIRE(rhs.valid());
  if (this != &rhs) {
    FieldBase<Dimension>::operator=(rhs);
    mDataArray = rhs.mDataArray;
    mValid = rhs.mValid;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assigment operator with a vector.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const std::vector<DataType>& rhs) {
  REQUIRE(mValid);
  REQUIRE(this->nodeList().numNodes() == rhs.size());
  mDataArray = rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Assignment operator with a constant value of DataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const DataType& rhs) {
  REQUIRE(mValid);
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    *fieldItr = rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Test equivalence with a FieldBase.
//------------------------------------------------------------------------------
template<typename Value>
struct CrappyFieldCompareMethod {
  static bool compare(const std::vector<Value>& lhs, 
                      const std::vector<Value>& rhs) {
    return lhs == rhs;
  }
};

template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::operator==(const FieldBase<Dimension>& rhs) const {
  if (this->name() != rhs.name()) return false;
  if (this->nodeListPtr() != rhs.nodeListPtr()) return false;
  try {
    const Field<Dimension, DataType>* rhsPtr = dynamic_cast<const Field<Dimension, DataType>*>(&rhs);
    if (rhsPtr == 0) return false;
    return CrappyFieldCompareMethod<DataType>::compare(mDataArray, rhsPtr->mDataArray);
  } catch (std::bad_cast) {
    return false;
  }
}

//------------------------------------------------------------------------------
// Element access by integer index.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::operator()(int index) {
  CHECK(index >= 0 && index < numElements());
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::operator()(int index) const {
  CHECK(index >= 0 && index < numElements());
  return mDataArray[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with STL interface.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::at(int index) {
  return this->operator()(index);
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::at(int index) const {
  return this->operator()(index);
}

//------------------------------------------------------------------------------
// Element access by Node ID iterator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) {
  CHECK(itr.nodeListPtr() == this->nodeListPtr());
  CHECK(itr.nodeID() >= 0 && itr.nodeID() < numElements());
  return mDataArray[itr.nodeID()];
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) const {
  CHECK(itr.nodeListPtr() == this->nodeListPtr());
  CHECK(itr.nodeID() >= 0 && itr.nodeID() < numElements());
  return mDataArray[itr.nodeID()];
}

// //------------------------------------------------------------------------------
// // Element access by CoarseNodeIterator.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// DataType&
// Field<Dimension, DataType>::
// operator()(const CoarseNodeIterator<Dimension>& itr) {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewCoarseNodes) {
//     cacheCoarseValues();
//     mNewCoarseNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mCoarseCache.size());
//   return mCoarseCache[itr.cacheID()];
// }

// template<typename Dimension, typename DataType>
// inline
// const DataType&
// Field<Dimension, DataType>::
// operator()(const CoarseNodeIterator<Dimension>& itr) const {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewCoarseNodes) {
//     cacheCoarseValues();
//     mNewCoarseNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mCoarseCache.size());
//   return mCoarseCache[itr.cacheID()];
// }

// //------------------------------------------------------------------------------
// // Element access by RefineNodeIterator.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// DataType&
// Field<Dimension, DataType>::
// operator()(const RefineNodeIterator<Dimension>& itr) {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewRefineNodes) {
//     cacheRefineValues();
//     mNewRefineNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mRefineCache.size());
//   return mRefineCache[itr.cacheID()];
// }

// template<typename Dimension, typename DataType>
// inline
// const DataType&
// Field<Dimension, DataType>::
// operator()(const RefineNodeIterator<Dimension>& itr) const {
//   CHECK(itr.nodeListPtr() == nodeListPtr());
//   if (mNewRefineNodes) {
//     cacheRefineValues();
//     mNewRefineNodes = false;
//   }
//   CHECK(itr.cacheID() >= 0 && itr.cacheID() < mRefineCache.size());
//   return mRefineCache[itr.cacheID()];
// }

//------------------------------------------------------------------------------
// Number of elements in the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
int
Field<Dimension, DataType>::numElements() const {
  return mDataArray.size();
}

template<typename Dimension, typename DataType>
inline
int
Field<Dimension, DataType>::numInternalElements() const {
  return this->nodeList().numInternalNodes();
}

template<typename Dimension, typename DataType>
inline
int
Field<Dimension, DataType>::size() const {
  return numElements();
}

//------------------------------------------------------------------------------
// Zero out the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::Zero() {
  REQUIRE(mValid);
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    *fieldItr = DataTypeTraits<DataType>::zero();
  }
}

//------------------------------------------------------------------------------
// Apply a minimum value to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyMin(const DataType& dataMin) {
  REQUIRE(mValid);
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    *fieldItr = std::max(dataMin, *fieldItr);
  }
}

//------------------------------------------------------------------------------
// Apply a maximum value to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyMax(const DataType& dataMax) {
  REQUIRE(mValid);
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    *fieldItr = std::min(dataMax, *fieldItr);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) minimum value  to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyScalarMin(const double dataMin) {
  REQUIRE(mValid);
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    *fieldItr = Spheral::max(dataMin, *fieldItr);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) maximum value to the Field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::applyScalarMax(const double dataMax) {
  REQUIRE(mValid);
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    *fieldItr = Spheral::min(dataMax, *fieldItr);
  }
}

//------------------------------------------------------------------------------
// Addition with another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator+(const Field<Dimension, DataType>& rhs) const {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  for (int i = 0; i < result.size(); ++i) {
    result(i) += rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Subtract another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator-(const Field<Dimension, DataType>& rhs) const {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  for (int i = 0; i < result.size(); ++i) {
    result(i) -= rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Addition with another field in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator+=(const Field<Dimension, DataType>& rhs) {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  for (int i = 0; i < size(); ++i) {
    (*this)(i) += rhs(i);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Subtract another field from this one in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator-=(const Field<Dimension, DataType>& rhs) {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  for (int i = 0; i < size(); ++i) {
    (*this)(i) -= rhs(i);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator+(const DataType& rhs) const {
  REQUIRE(valid());
  Field<Dimension, DataType> result(*this);
  for (int i = 0; i < result.size(); ++i) {
    result(i) += rhs;
  }
  return result;
}

//------------------------------------------------------------------------------
// Subtract a single value from every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator-(const DataType& rhs) const {
  REQUIRE(valid());
  Field<Dimension, DataType> result(*this);
  for (int i = 0; i < result.size(); ++i) {
    result(i) -= rhs;
  }
  return result;
}

//------------------------------------------------------------------------------
// Addition with a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator+=(const DataType& rhs) {
  REQUIRE(valid());
  for (int i = 0; i < size(); ++i) {
    (*this)(i) += rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value in place
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator-=(const DataType& rhs) {
  REQUIRE(valid());
  for (int i = 0; i < size(); ++i) {
    (*this)(i) -= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by another Field in place.  Only meaningful when multiplying
// by a scalar field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator*=(const Field<Dimension, Scalar>& rhs) {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  for (int i = 0; i < size(); ++i) {
    (*this)(i) *= rhs(i);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator*=(const Scalar& rhs) {
  REQUIRE(valid());
  for (int i = 0; i < size(); ++i) {
    (*this)(i) *= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Division by another Field.
// Only meaningful for Scalar Fields!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator/(const Field<Dimension, typename Dimension::Scalar>& rhs) const {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  for (int i = 0; i < result.size(); ++i) {
    CHECK(rhs(i) != 0.0);
    result(i) /= rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Division by another Field in place.
// Only meaningful for Scalar Fields!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator/=(const Field<Dimension, typename Dimension::Scalar>& rhs) {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  for (int i = 0; i < size(); ++i) {
    REQUIRE(rhs(i) != 0.0);
    (*this)(i) /= rhs(i);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Division by a Scalar value.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator/(const Scalar& rhs) const {
  REQUIRE(valid());
  REQUIRE(rhs != 0.0);
  Field<Dimension, DataType> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division by a Scalar value in place.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::
operator/=(const Scalar& rhs) {
  REQUIRE(valid());
  REQUIRE(rhs != 0.0);
  for (int i = 0; i < size(); ++i) {
    CHECK(rhs != 0.0);
    (*this)(i) /= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Sum the elements of the field (assumes the DataType::operator+= is 
// available).
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
sumElements() const {
  DataType result = localSumElements();
#ifdef USE_MPI
  {
    // In parallel, do the reduction.
    // Note this will fail to compile for types that do not have corresponding MPI type!
    DataType tmp = result;
    MPI_Allreduce(&tmp, &result, 1, DataTypeTraits<DataType>::MpiDataType(), MPI_SUM, Communicator::communicator());
  }
#endif
  return result;
}

//------------------------------------------------------------------------------
// Minimum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
min() const {
  DataType result = localMin();
#ifdef USE_MPI
  {
    // In parallel, do the reduction.
    // Note this will fail to compile for types that do not have corresponding MPI type!
    DataType tmp = result;
    MPI_Allreduce(&tmp, &result, 1, DataTypeTraits<DataType>::MpiDataType(), MPI_MIN, Communicator::communicator());
  }
#endif
  return result;
}

//------------------------------------------------------------------------------
// Maximum.
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
max() const {
  DataType result = localMax();
#ifdef USE_MPI
  {
    // In parallel, do the reduction.
    // Note this will fail to compile for types that do not have corresponding MPI type!
    DataType tmp = result;
    MPI_Allreduce(&tmp, &result, 1, DataTypeTraits<DataType>::MpiDataType(), MPI_MAX, Communicator::communicator());
  }
#endif
  return result;
}

//------------------------------------------------------------------------------
// Sum the elements of the field (assumes the DataType::operator+= is 
// available).  
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
localSumElements() const {
  DataType result = DataTypeTraits<DataType>::zero();
  for (const_iterator elementItr = begin(); elementItr < begin() + numInternalElements(); ++elementItr) {
    result += *elementItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// Minimum.
// LOCAL to processor!
//-------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
Field<Dimension, DataType>::
localMin() const {
  DataType result;
  if (size() == 0) {
    result = std::numeric_limits<DataType>::max();
  } else {
    result = *std::min_element(begin(), begin() + numInternalElements());
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
Field<Dimension, DataType>::
localMax() const {
  DataType result;
  if (size() == 0) {
    result = -std::numeric_limits<DataType>::max() < std::numeric_limits<DataType>::min() ? -std::numeric_limits<DataType>::max() : std::numeric_limits<DataType>::min();
  } else {
    result = *std::max_element(begin(), begin() + numInternalElements());
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator==(const Field<Dimension, DataType>& rhs) const {
  if (this->size() != rhs.size()) return false;
  bool result = true;
  const_iterator lhsItr = this->begin();
  const_iterator rhsItr = rhs.begin();
  while (lhsItr < this->end() && result) {
    result = *lhsItr == *rhsItr;
    ++lhsItr;
    ++rhsItr;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator!=(const Field<Dimension, DataType>& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator>(const Field<Dimension, DataType>& rhs) const {
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
// operator<(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator<(const Field<Dimension, DataType>& rhs) const {
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
// operator>=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator>=(const Field<Dimension, DataType>& rhs) const {
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
// operator<=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator<=(const Field<Dimension, DataType>& rhs) const {
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
Field<Dimension, DataType>::
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
Field<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// operator>(value)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
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
Field<Dimension, DataType>::
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
Field<Dimension, DataType>::
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
Field<Dimension, DataType>::
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
// Test if the field is valid.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::valid() const {
  return mValid && this->nodeListPtr();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::begin() {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::end() {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::internalBegin() {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::internalEnd() {
  CHECK(this->nodeList().firstGhostNode() >= 0 &&
         this->nodeList().firstGhostNode() <= this->nodeList().numNodes());
  return mDataArray.begin() + this->nodeList().firstGhostNode();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::ghostBegin() {
  return this->internalEnd();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::ghostEnd() {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::begin() const {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::end() const {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::internalBegin() const {
  return mDataArray.begin();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::internalEnd() const {
  REQUIRE(this->nodeList().firstGhostNode() >= 0 &&
          this->nodeList().firstGhostNode() <= this->nodeList().numNodes());
  return mDataArray.begin() + this->nodeList().firstGhostNode();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::ghostBegin() const {
  return this->internalEnd();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::ghostEnd() const {
  return mDataArray.end();
}

//------------------------------------------------------------------------------
// Index operators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::
operator[](const unsigned int index) {
  return mDataArray[index];
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::
operator[](const unsigned int index) const {
  return mDataArray[index];
}

//------------------------------------------------------------------------------
// Set the NodeList with which this field is defined.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::setNodeList(const NodeSpace::NodeList<Dimension>& nodeList) {
  int oldSize = this->size();
  this->setFieldBaseNodeList(nodeList);
  mDataArray.resize(nodeList.numNodes());
  if (this->size() > oldSize) {
    for (int i = oldSize; i < this->size(); ++i) {
      (*this)(i) = DataTypeTraits<DataType>::zero();
    }
  }
  mValid = true;
}

//------------------------------------------------------------------------------
// Resize the field to the given number of nodes.  This operation ignores
// the distinction between internal and ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeField(int size) {
  REQUIRE(size == this->nodeList().numNodes());
  int oldSize = this->size();
  mDataArray.resize(size);
  if (oldSize < size) {
    std::fill(mDataArray.begin() + oldSize,
              mDataArray.end(),
              DataTypeTraits<DataType>::zero());
  }
  mValid = true;
}

//------------------------------------------------------------------------------
// Delete the given element id.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::deleteElement(int nodeID) {
  const int originalSize = this->size();
  REQUIRE(nodeID >= 0 && nodeID < originalSize);
  mDataArray.erase(mDataArray.begin() + nodeID);
  ENSURE(mDataArray.size() == originalSize - 1);
}

//------------------------------------------------------------------------------
// Delete the given set of elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::deleteElements(const std::vector<int>& nodeIDs) {
  // The standalone method does the actual work.
  removeElements(mDataArray, nodeIDs);
}

//------------------------------------------------------------------------------
// Pack the given Field values by appending the elements decomposed into simple
// Scalars onto the given array.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<char>
Field<Dimension, DataType>::
packValues(const std::vector<int>& nodeIDs) const {
  return packFieldValues(*this, nodeIDs);
}

//------------------------------------------------------------------------------
// Unpack the given buffer into the requested field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
unpackValues(const int numElements,
             const int beginInsertionIndex,
             const std::vector<char>& buffer) {

  REQUIRE(numElements >= 0);
  REQUIRE(beginInsertionIndex >= 0 &&
          beginInsertionIndex < this->size());
  const int endIndex = beginInsertionIndex + numElements;
  REQUIRE(endIndex <= this->size());

  // The unpackFieldValues method requires the insertion indicies explicitly.
  std::vector<int> indicies;
  indicies.reserve(numElements);
  for (int i = 0; i != numElements; ++i) indicies.push_back(beginInsertionIndex + i);

  // Now we're ready to do the deed...
  unpackFieldValues(*this, indicies, buffer);
}

//------------------------------------------------------------------------------
// Resize the field to the given number of internal nodes, preserving any ghost
// values at the end.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeFieldInternal(const int size,
                                                const int oldFirstGhostNode) {
  const int currentSize = this->size();
  const int currentInternalSize = oldFirstGhostNode;
  const int numGhostNodes = this->nodeList().numGhostNodes();
  const int newSize = size + numGhostNodes;
  REQUIRE(numGhostNodes == currentSize - oldFirstGhostNode);
  REQUIRE(newSize == this->nodeList().numNodes());

  // If there is ghost data, we must preserve it.
  std::vector<DataType> oldGhostValues(numGhostNodes);
  if (numGhostNodes > 0) {
    for (int i = 0; i != numGhostNodes; ++i) {
      const int j = oldFirstGhostNode + i;
      CHECK(i >= 0 && i < numGhostNodes);
      CHECK(j >= 0 && j < this->size());
      oldGhostValues[i] = (*this)(j);
    }
  }

  // Resize the field data.
  mDataArray.resize(newSize);

  // Fill in any new internal values.
  if (newSize > currentSize) {
    CHECK(currentInternalSize < this->nodeList().firstGhostNode());
    std::fill(mDataArray.begin() + currentInternalSize,
              mDataArray.begin() + this->nodeList().firstGhostNode(),
              DataTypeTraits<DataType>::zero());
  }

  // Fill the ghost data back in.
  if (numGhostNodes > 0) {
    for (int i = 0; i != numGhostNodes; ++i) {
      const int j = this->nodeList().firstGhostNode() + i;
      CHECK(i >= 0 && i < oldGhostValues.size());
      CHECK(j >= 0 && j < this->size());
      (*this)(j) = oldGhostValues[i];
    }
  }

  mValid = true;
}

//------------------------------------------------------------------------------
// Resize the field to the given number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeFieldGhost(const int size) {
  const int currentSize = this->size();
  const int numInternalNodes = this->nodeList().numInternalNodes();
  const int currentNumGhostNodes = currentSize - numInternalNodes;
  REQUIRE(currentNumGhostNodes >= 0);
  const int newSize = numInternalNodes + size;
  REQUIRE(newSize == this->nodeList().numNodes());

  // Resize the field data.
  mDataArray.resize(newSize);
  CHECK(this->size() == newSize);

  // Fill in any new ghost values.
  if (newSize > currentSize) {
    std::fill(mDataArray.begin() + numInternalNodes + currentNumGhostNodes,
              mDataArray.end(),
              DataTypeTraits<DataType>::zero());
  }

  mValid = true;
}

//------------------------------------------------------------------------------
// Pack the Field into a string.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::string
Field<Dimension, DataType>::
string(const int precision) const {
  const int n = numInternalElements();
  std::vector<int> indicies;
  indicies.reserve(n);
  for (int i = 0; i != n; ++i) indicies.push_back(i);
  CHECK(indicies.size() == n);
  const std::vector<char> packedValues = packFieldValues(*this, indicies);
  return std::string(this->name()) + "|" + std::string(packedValues.begin(), packedValues.end());
}

//------------------------------------------------------------------------------
// Unpack the values from a string into this Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
string(const std::string& s) {
  const int n = numInternalElements();
  std::vector<int> indicies;
  indicies.reserve(n);
  for (int i = 0; i != n; ++i) indicies.push_back(i);
  CHECK(indicies.size() == n);
  const size_t j = s.find("|");
  CHECK(j != std::string::npos and
        j < s.size());
  this->name(s.substr(0, j));
  const std::vector<char> packedValues(s.begin() + j + 1, s.end());
  unpackFieldValues(*this, indicies, packedValues);
}

//------------------------------------------------------------------------------
// Construct std::vectors of pointers to the values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
internalValues() const {
  std::vector<DataType> result;
  result.reserve(this->nodeList().numInternalNodes());
  for (const_iterator itr = internalBegin();
       itr != internalEnd();
       ++itr) result.push_back(*itr);
  ENSURE(result.size() == this->nodeList().numInternalNodes());
  return result;
}

template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
ghostValues() const {
  std::vector<DataType> result;
  result.reserve(this->nodeList().numGhostNodes());
  for (const_iterator itr = ghostBegin();
       itr != ghostEnd();
       ++itr) result.push_back(*itr);
  ENSURE(result.size() == this->nodeList().numGhostNodes());
  return result;
}

template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
allValues() const {
  std::vector<DataType> result;
  result.reserve(this->nodeList().numNodes());
  for (const_iterator itr = begin();
       itr != end();
       ++itr) result.push_back(*itr);
  ENSURE(result.size() == this->nodeList().numNodes());
  return result;
}

//****************************** Global Functions ******************************
//------------------------------------------------------------------------------
// Multiplication by another Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const Field<Dimension, DataType>& lhs,
          const Field<Dimension, OtherDataType>& rhs) {
  CHECK(lhs.valid() && rhs.valid());
  CHECK(lhs.nodeList().numNodes() == rhs.nodeList().numNodes());
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, DataType>&>(lhs).nodeList());
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = lhs(i) * rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a single value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const Field<Dimension, DataType>& lhs,
          const OtherDataType& rhs) {
  CHECK(lhs.valid());
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, DataType>&>(lhs).nodeList());
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = lhs(i) * rhs;
  }
  return result;
}

template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const DataType& lhs,
          const Field<Dimension, OtherDataType>& rhs) {
  CHECK(rhs.valid());
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, OtherDataType>&>(rhs).nodeList());
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = lhs * rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Absolute value.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
abs(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::abs(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Inverse cosine
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
acos(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::acos(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Inverse sine
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
asin(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::asin(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Inverse tangent
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
atan(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::atan(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Inverse tangent2.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
atan2(const Field<Dimension, typename Dimension::Scalar>& field1,
      const Field<Dimension, typename Dimension::Scalar>& field2) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field1.valid() && field2.valid());
  CHECK(field1.nodeListPtr() == field2.nodeListPtr());
  Field<Dimension, Scalar> result(field1);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::atan2(field1(i), field2(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Ceiling -- smallest floating-point integer value not less than the argument.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
ceil(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::ceil(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Cosine
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
cos(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::cos(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Hyperbolic cosine
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
cosh(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::cosh(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Exponential e^x
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
exp(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::exp(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// fabs -- same as abs, absolute value
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
fabs(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::abs(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Floor -- largest floating-point integer value not greater than the argument
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
floor(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::floor(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Natural logarithm
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
log(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::log(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Log base 10
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
log10(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::log10(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// pow -- raise each element to an arbitrary power
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow(const Field<Dimension, typename Dimension::Scalar>& field,
    const double exponent) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::pow(result(i), exponent);
  }
  return result;
}

//------------------------------------------------------------------------------
// powN -- raise each element to the power N
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow2(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = FastMath::pow2(result(i));
  }
  return result;
}

template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow3(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = FastMath::pow3(result(i));
  }
  return result;
}

template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow4(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = FastMath::pow4(result(i));
  }
  return result;
}

template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow5(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = FastMath::pow5(result(i));
  }
  return result;
}

template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow6(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = FastMath::pow6(result(i));
  }
  return result;
}

template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow7(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = FastMath::pow7(result(i));
  }
  return result;
}

template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
pow8(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = FastMath::pow8(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Sine
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
sin(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::sin(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Hyperbolic sine
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
sinh(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = std::sinh(result(i));
  }
  return result;
}

//------------------------------------------------------------------------------
// Sqr -- square, same as pow2()
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
sqr(const Field<Dimension, typename Dimension::Scalar>& field) {
  return field.pow2();
}

//------------------------------------------------------------------------------
// Sqrt
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
sqrt(const Field<Dimension, typename Dimension::Scalar>& field) {
  typedef typename Dimension::Scalar Scalar;
  CHECK(field.valid());
  Field<Dimension, Scalar> result(field);
  for (int i = 0; i < result.numElements(); ++i) {
    CHECK(result(i) >= 0.0);
    result(i) = std::sqrt(result(i));
  }
  return result;
}

}

//------------------------------------------------------------------------------
// Minimum
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldSpace::Field<Dimension, DataType>
min(const FieldSpace::Field<Dimension, DataType>& field1,
    const FieldSpace::Field<Dimension, DataType>& field2) {
  CHECK(field1.valid() && field2.valid());
  CHECK(field1.numElements() == field2.numElements());
  CHECK(field1.nodeListPtr() == field2.nodeListPtr());
  FieldSpace::Field<Dimension, DataType> result("min", const_cast<NodeSpace::NodeList<Dimension>&>(field1.nodeList()));
  for (int i = 0; i < field1.numElements(); ++i) {
    result(i) = std::min(field1(i), field2(i));
  }
  return result;
}

template<typename Dimension, typename DataType>
FieldSpace::Field<Dimension, DataType>
min(const DataType& value,
    const FieldSpace::Field<Dimension, DataType>& field) {
  CHECK(field.valid());
  FieldSpace::Field<Dimension, DataType> result("min", const_cast<NodeSpace::NodeList<Dimension>&>(field.nodeList()));
  for (int i = 0; i < field.numElements(); ++i) {
    result(i) = std::min(value, field(i));
  }
  return result;
}

template<typename Dimension, typename DataType>
FieldSpace::Field<Dimension, DataType>
min(const FieldSpace::Field<Dimension, DataType>& field, 
    const DataType& value) {
  CHECK(field.valid());
  FieldSpace::Field<Dimension, DataType> result("min", const_cast<NodeSpace::NodeList<Dimension>&>(field.nodeList()));
  for (int i = 0; i < field.numElements(); ++i) {
    result(i) = std::min(field(i), value);
  }
  return result;
}

//------------------------------------------------------------------------------
// Maximum
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldSpace::Field<Dimension, DataType>
max(const FieldSpace::Field<Dimension, DataType>& field1,
    const FieldSpace::Field<Dimension, DataType>& field2) {
  CHECK(field1.valid() && field2.valid());
  CHECK(field1.numElements() == field2.numElements());
  CHECK(field1.nodeListPtr() == field2.nodeListPtr());
  FieldSpace::Field<Dimension, DataType> result("max", const_cast<NodeSpace::NodeList<Dimension>&>(field1.nodeList()));
  for (int i = 0; i < field1.numElements(); ++i) {
    result(i) = std::max(field1(i), field2(i));
  }
  return result;
}

template<typename Dimension, typename DataType>
FieldSpace::Field<Dimension, DataType>
max(const DataType& value,
    const FieldSpace::Field<Dimension, DataType>& field) {
  CHECK(field.valid());
  FieldSpace::Field<Dimension, DataType> result("max", const_cast<NodeSpace::NodeList<Dimension>&>(field.nodeList()));
  for (int i = 0; i < field.numElements(); ++i) {
    result(i) = std::max(value, field(i));
  }
  return result;
}

template<typename Dimension, typename DataType>
FieldSpace::Field<Dimension, DataType>
max(const FieldSpace::Field<Dimension, DataType>& field, 
    const DataType& value) {
  CHECK(field.valid());
  FieldSpace::Field<Dimension, DataType> result("max", const_cast<NodeSpace::NodeList<Dimension>&>(field.nodeList()));
  for (int i = 0; i < field.numElements(); ++i) {
    result(i) = std::max(field(i), value);
  }
  return result;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::istream&
operator>>(std::istream& is, FieldSpace::Field<Dimension, DataType>& field) {

  // Start by reading the number of elements.
  int numElementsInStream;
  is >> numElementsInStream;
  CHECK(numElementsInStream == field.nodeList().numInternalNodes());

  // Read in the elements.
  for (typename FieldSpace::Field<Dimension, DataType>::iterator itr = field.internalBegin();
       itr < field.internalEnd();
       ++itr) {
    is >> *itr;
  }
  return is;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::ostream&
operator<<(std::ostream& os, const FieldSpace::Field<Dimension, DataType>& field) {

  // Write the number of internal elements.
  os << field.nodeList().numInternalNodes() << " ";

  // Write the internal elements.
  for (typename FieldSpace::Field<Dimension, DataType>::const_iterator itr = field.internalBegin();
       itr < field.internalEnd();
       ++itr) {
    os << *itr << " ";
  }
//   os << endl;
  return os;
}

}

