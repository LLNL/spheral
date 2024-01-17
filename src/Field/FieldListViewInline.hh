// Includes.
#include "Geometry/MathTraits.hh"
#include "NodeIterators.hh"

#include "NodeList/FluidNodeList.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Neighbor/Neighbor.hh"
#include "Field/Field.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/allReduce.hh"

#ifdef USE_MPI
#include <mpi.h>
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/packElement.hh"
#include "Distributed/Communicator.hh"
#endif

#include <algorithm>
#include <limits.h>
#include <float.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>::FieldListView():
  mFieldPtrs(0) {}

//------------------------------------------------------------------------------
// FieldList constructor.
// Note that the copy constructor copies the data by reference by default.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>::
FieldListView(const FieldList<Dimension, DataType>& rhs):
  mFieldPtrs(rhs.mFieldPtrs) {}

//------------------------------------------------------------------------------
// Copy constructor.
// Note that the copy constructor copies the data by reference by default.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>::
FieldListView(const FieldListView<Dimension, DataType>& rhs):
  mFieldPtrs(rhs.mFieldPtrs) {}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>::~FieldListView() {}

//------------------------------------------------------------------------------
// Assignment with another FieldListView.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator=(const FieldListView<Dimension, DataType>& rhs) {
  if (this != &rhs) {
    mFieldPtrs = rhs.mFieldPtrs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assignment with a constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldListView<Dimension, DataType>&
FieldListView<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (auto fieldPtrItr = begin(); fieldPtrItr < end(); ++fieldPtrItr) {
    (*fieldPtrItr)->operator=(rhs);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Standard iterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::iterator
FieldListView<Dimension, DataType>::
begin() {
  return mFieldPtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::iterator
FieldListView<Dimension, DataType>::
end() {
  return mFieldPtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::const_iterator
FieldListView<Dimension, DataType>::
begin() const {
  return mFieldPtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::const_iterator
FieldListView<Dimension, DataType>::
end() const {
  return mFieldPtrs.end();
}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::ElementType
FieldListView<Dimension, DataType>::
operator[](const unsigned int index) const {
  REQUIRE2(index < mFieldPtrs.size(), "FieldListView index ERROR: out of bounds " << index << " !< " << mFieldPtrs.size());
  return mFieldPtrs[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::ElementType
FieldListView<Dimension, DataType>::
operator[](const unsigned int index) {
  REQUIRE2(index < mFieldPtrs.size(), "FieldListView index ERROR: out of bounds " << index << " !< " << mFieldPtrs.size());
  return mFieldPtrs[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::ElementType
FieldListView<Dimension, DataType>::
at(const unsigned int index) const {
  return (*this)[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldListView<Dimension, DataType>::ElementType
FieldListView<Dimension, DataType>::
at(const unsigned int index) {
  return (*this)[index];
}

//------------------------------------------------------------------------------
// Provide a more primitive access to Field elements, based on the index of the Field
// and the node index within that field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldListView<Dimension, DataType>::
operator()(const unsigned int fieldIndex,
           const unsigned int nodeIndex) {
  REQUIRE2(fieldIndex < mFieldPtrs.size(), "FieldListView index ERROR: out of bounds " << fieldIndex << " !< " << mFieldPtrs.size());
  REQUIRE2(nodeIndex < mFieldPtrs[fieldIndex]->size(), "FieldListView node index ERROR: out of bounds " << nodeIndex << " !< " << mFieldPtrs[fieldIndex]->size());
  return mFieldPtrs[fieldIndex]->operator()(nodeIndex);
}

template<typename Dimension, typename DataType>
inline
DataType&
FieldListView<Dimension, DataType>::
operator()(const unsigned int fieldIndex,
           const unsigned int nodeIndex) const {
  REQUIRE2(fieldIndex < mFieldPtrs.size(), "FieldListView index ERROR: out of bounds " << fieldIndex << " !< " << mFieldPtrs.size());
  REQUIRE2(nodeIndex < mFieldPtrs[fieldIndex]->size(), "FieldListView node index ERROR: out of bounds " << nodeIndex << " !< " << mFieldPtrs[fieldIndex]->size());
  return mFieldPtrs[fieldIndex]->operator()(nodeIndex);
}

//------------------------------------------------------------------------------
// Zero out the FieldListView.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListView<Dimension, DataType>::Zero() {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->Zero();
  }
}

//------------------------------------------------------------------------------
// Apply a minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListView<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyMin(dataMin);
  }
}

//------------------------------------------------------------------------------
// Apply a maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListView<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyMax(dataMax);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListView<Dimension, DataType>::applyScalarMin(const double dataMin) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyScalarMin(dataMin);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldListView<Dimension, DataType>::applyScalarMax(const double dataMax) {
  for (auto fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyScalarMax(dataMax);
  }
}

////------------------------------------------------------------------------------
//// Add two FieldListViews.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::operator+(const FieldListView<Dimension, DataType>& rhs) const {
//
//  // Pre-conditions.
//  BEGIN_CONTRACT_SCOPE
//  REQUIRE(numFields() == rhs.numFields());
//  for (int i = 0; i != (int)numFields(); ++i) 
//    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
//  END_CONTRACT_SCOPE
//
//  FieldListView<Dimension, DataType> result(*this);
//  result.copyFields();
//  result += rhs;
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Subtract a FieldListView from another.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::operator-(const FieldListView<Dimension, DataType>& rhs) const {
//
//  // Pre-conditions.
//  BEGIN_CONTRACT_SCOPE
//  REQUIRE(numFields() == rhs.numFields());
//  for (int i = 0; i != (int)numFields(); ++i) 
//    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
//  END_CONTRACT_SCOPE
//
//  FieldListView<Dimension, DataType> result(*this);
//  result.copyFields();
//  result -= rhs;
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Add two FieldListViews in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::operator+=(const FieldListView<Dimension, DataType>& rhs) {
//
//  // Pre-conditions.
//  BEGIN_CONTRACT_SCOPE
//  REQUIRE(numFields() == rhs.numFields());
//  for (auto i = 0u; i != numFields(); ++i) 
//    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
//  END_CONTRACT_SCOPE
//
//  for (auto i = 0u; i < numFields(); ++i) {
//    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
//    *((*this)[i]) += *(rhs[i]);
//  }
//  return *this;
//}
//
////------------------------------------------------------------------------------
//// Subtract a FieldListView from another in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::operator-=(const FieldListView<Dimension, DataType>& rhs) {
//
//  // Pre-conditions.
//  BEGIN_CONTRACT_SCOPE
//  REQUIRE(numFields() == rhs.numFields());
//  for (int i = 0; i != (int)numFields(); ++i) 
//    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
//  END_CONTRACT_SCOPE
//
//  for (auto i = 0u; i < numFields(); ++i) {
//    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
//    *((*this)[i]) -= *(rhs[i]);
//  }
//  return *this;
//}
//
////------------------------------------------------------------------------------
//// Add a single value to a FieldListView.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::operator+(const DataType& rhs) const {
//  FieldListView<Dimension, DataType> result(*this);
//  result.copyFields();
//  result += rhs;
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Subtract a single value from a Field.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::operator-(const DataType& rhs) const {
//  FieldListView<Dimension, DataType> result(*this);
//  result.copyFields();
//  result -= rhs;
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Add a single value to the Field in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::operator+=(const DataType& rhs) {
//  for (auto i = 0u; i < numFields(); ++i) {
//    *((*this)[i]) += rhs;
//  }
//  return *this;
//}
//
////------------------------------------------------------------------------------
//// Subtract a single value from the FieldListView in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::operator-=(const DataType& rhs) {
//  for (auto i = 0u; i < numFields(); ++i) {
//    *((*this)[i]) -= rhs;
//  }
//  return *this;
//}
//
////------------------------------------------------------------------------------
//// Multiply this FieldListView by a Scalar FieldListView in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::
//operator*=(const FieldListView<Dimension, typename Dimension::Scalar>& rhs) {
//  REQUIRE(this->numFields() == rhs.numFields());
//  for (auto i = 0u; i < numFields(); ++i) {
//    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
//    *((*this)[i]) *= *(rhs[i]);
//  }
//  return *this;
//}
//
////------------------------------------------------------------------------------
//// Multiply this FieldListView by a Scalar in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::operator*=(const Scalar& rhs) {
//  for (auto i = 0u; i < numFields(); ++i) {
//    *((*this)[i]) *= rhs;
//  }
//  return *this;
//}
//
////------------------------------------------------------------------------------
//// Divide this FieldListView by a Scalar FieldListView.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::
//operator/(const FieldListView<Dimension, typename Dimension::Scalar>& rhs) const {
//  REQUIRE(this->numFields() == rhs.numFields());
//  FieldListView<Dimension, DataType> result(*this);
//  result.copyFields();
//  result /= rhs;
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Divide this FieldListView by a Scalar value.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::operator/(const Scalar& rhs) const {
//  REQUIRE(rhs != 0.0);
//  FieldListView<Dimension, DataType> result(*this);
//  result.copyFields();
//  result /= rhs;
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Divide this FieldListView by a Scalar FieldListView in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::
//operator/=(const FieldListView<Dimension, typename Dimension::Scalar>& rhs) {
//  REQUIRE(this->numFields() == rhs.numFields());
//  for (auto i = 0u; i < numFields(); ++i) {
//    CHECK2((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr(), (*this)[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
//    *((*this)[i]) /= *(rhs[i]);
//  }
//  return *this;
//}
//
////------------------------------------------------------------------------------
//// Divide this FieldListView by a Scalar in place.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>&
//FieldListView<Dimension, DataType>::operator/=(const typename Dimension::Scalar& rhs) {
//  REQUIRE(rhs != 0.0);
//  for (auto i = 0u; i < numFields(); ++i) {
//    *((*this)[i]) /= rhs;
//  }
//  return *this;
//}

//------------------------------------------------------------------------------
// Sum the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldListView<Dimension, DataType>::
sumElements() const {
  return allReduce(this->localSumElements(), MPI_SUM, Communicator::communicator());
}

//------------------------------------------------------------------------------
// Find the minimum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldListView<Dimension, DataType>::
min() const {
  return allReduce(this->localMin(), MPI_MIN, Communicator::communicator());
}

//------------------------------------------------------------------------------
// Find the maximum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldListView<Dimension, DataType>::
max() const {
  return allReduce(this->localMax(), MPI_MAX, Communicator::communicator());
}

//------------------------------------------------------------------------------
// Sum the field elements.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldListView<Dimension, DataType>::
localSumElements() const {
  auto result = DataTypeTraits<DataType>::zero();
  for (auto itr = begin(); itr < end(); ++itr) result += (*itr)->localSumElements();
  return result;
}

//------------------------------------------------------------------------------
// Find the minimum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldListView<Dimension, DataType>::
localMin() const {
  auto result = std::numeric_limits<DataType>::max();
  for (auto itr = begin(); itr != end(); ++itr) result = std::min(result, (*itr)->localMin());
  return result;
}

//------------------------------------------------------------------------------
// Find the maximum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldListView<Dimension, DataType>::
localMax() const {
  auto result = std::numeric_limits<DataType>::lowest();
  for (auto itr = begin(); itr < end(); ++itr) result = std::max(result, (*itr)->localMax());
  return result;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListView<Dimension, DataType>::
operator==(const FieldListView<Dimension, DataType>& rhs) const {
  if (&rhs == this) return true;
  return mFieldPtrs == rhs.mFieldPtrs;

  //// Pre-conditions.
  //BEGIN_CONTRACT_SCOPE
  //REQUIRE(numFields() == rhs.numFields());
  //for (int i = 0; i != (int)numFields(); ++i) 
  //  REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  //END_CONTRACT_SCOPE

  //bool result = true;
  //int i = 0;
  //while (result && i != (int)numFields()) {
  //  result = result && (*(mFieldPtrs[i]) == *(rhs[i]));
  //  ++i;
  //}
  //return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListView<Dimension, DataType>::
operator!=(const FieldListView<Dimension, DataType>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListView<Dimension, DataType>::
operator>(const FieldListView<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) > *(rhs[i]));
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
FieldListView<Dimension, DataType>::
operator<(const FieldListView<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) < *(rhs[i]));
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
FieldListView<Dimension, DataType>::
operator>=(const FieldListView<Dimension, DataType>& rhs) const {
  return operator>(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListView<Dimension, DataType>::
operator<=(const FieldListView<Dimension, DataType>& rhs) const {
  return operator<(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListView<Dimension, DataType>::
operator==(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != (int)numFields()) {
    result = result && (*(mFieldPtrs[i]) == rhs);
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
FieldListView<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListView<Dimension, DataType>::
operator>(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != (int)numFields()) {
    result = result && (*(mFieldPtrs[i]) > rhs);
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
FieldListView<Dimension, DataType>::
operator<(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != (int)numFields()) {
    result = result && (*(mFieldPtrs[i]) < rhs);
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
FieldListView<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  return operator>(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldListView<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  return operator<(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// Return the number of Fields stored in this Field List.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
unsigned 
FieldListView<Dimension, DataType>::numFields() const {
  return mFieldPtrs.size();
}

template<typename Dimension, typename DataType>
inline
unsigned 
FieldListView<Dimension, DataType>::size() const {
  return numFields();
}

////------------------------------------------------------------------------------
//// Make a thread local copy of the FieldListView -- assumes we're already in a
//// threaded region.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::
//threadCopy(const ThreadReduction reductionType,
//           const bool copy) {
//  FieldListView<Dimension, DataType> result;
//#pragma omp critical (FieldListView_threadCopy)
//  {
//    if (omp_get_num_threads() == 1) {
//
//      // In serial we can skip all the work copying
//      result.referenceFields(*this);
//
//    } else if (copy or
//               reductionType == ThreadReduction::MIN or
//               reductionType == ThreadReduction::MAX) {
//
//      // For min/max operations, we need to copy the original data
//      result.copyFields(*this);
//
//    } else {    
//
//      // Otherwise make standalone Fields of zeros
//      result = FieldListView<Dimension, DataType>(FieldStorageType::CopyFields);
//      for (auto fitr = this->begin(); fitr < this->end(); ++fitr) result.appendNewField((*fitr)->name(),
//                                                                                        (*fitr)->nodeList(),
//                                                                                        DataTypeTraits<DataType>::zero());
//    }
//    result.reductionType = reductionType;
//    result.threadMasterPtr = this;
//  }
//  return result;
//}
//
//template<typename Dimension, typename DataType>
//inline
//FieldListView<Dimension, DataType>
//FieldListView<Dimension, DataType>::
//threadCopy(typename SpheralThreads<Dimension>::FieldListViewStack& stack,
//           const ThreadReduction reductionType,
//           const bool copy) {
//  FieldListView<Dimension, DataType> result = this->threadCopy(reductionType, copy);
//  stack.push_back(&result);
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Reduce the values in the FieldListView with the passed thread-local values.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType>
//inline
//void
//FieldListView<Dimension, DataType>::
//threadReduce() const {
//  REQUIRE(threadMasterPtr != NULL);
//  REQUIRE(threadMasterPtr->size() == this->size());
//  if (omp_get_num_threads() > 1) {
//    const auto numNL = this->size();
//    for (auto k = 0u; k < numNL; ++k) {
//      const auto n = mFieldPtrs[k]->numInternalElements();
//      for (auto i = 0u; i < n; ++i) {
//
//        switch (reductionType) {
//
//        case ThreadReduction::SUM:
//          (*threadMasterPtr)(k,i) += (*this)(k,i);
//          break;
//
//        case ThreadReduction::MIN:
//          (*threadMasterPtr)(k,i) = std::min((*threadMasterPtr)(k,i), (*this)(k,i));
//          break;
//
//        case ThreadReduction::MAX:
//          (*threadMasterPtr)(k,i) = std::max((*threadMasterPtr)(k,i), (*this)(k,i));
//          break;
//
//        }
//      }
//    }
//  }
//}

//****************************** Global Functions ******************************

////------------------------------------------------------------------------------
//// Multiply a FieldListView by another FieldListView.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType, typename OtherDataType>
//inline
//FieldListView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
//operator*(const FieldListView<Dimension, DataType>& lhs,
//          const FieldListView<Dimension, OtherDataType>& rhs) {
//  REQUIRE(lhs.numFields() == rhs.numFields());
//  FieldListView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
//  result.copyFields();
//  for (auto i = 0u; i < lhs.numFields(); ++i) {
//    CHECK2(lhs[i]->nodeListPtr() == rhs[i]->nodeListPtr(), lhs[i]->nodeListPtr()->name() << " != " << rhs[i]->nodeListPtr()->name());
//    result.appendField((*(lhs[i])) * (*(rhs[i])));
//  }
//  return result;
//}
//
////------------------------------------------------------------------------------
//// Multiply a FieldListView by a single value.
////------------------------------------------------------------------------------
//template<typename Dimension, typename DataType, typename OtherDataType>
//inline
//FieldListView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
//operator*(const FieldListView<Dimension, DataType>& lhs,
//          const OtherDataType& rhs) {
//  FieldListView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
//  result.copyFields();
//  for (int i = 0; i < lhs.numFields(); ++i) {
//    result.appendField((*(lhs[i])) * rhs);
//  }
//  return result;
//}
//
//template<typename Dimension, typename DataType, typename OtherDataType>
//inline
//FieldListView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
//operator*(const DataType& lhs,
//          const FieldListView<Dimension, OtherDataType>& rhs) {
//  FieldListView<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
//  result.copyFields();
//  for (auto i = 0u; i < rhs.numFields(); ++i) {
//    result.appendField(lhs * (*(rhs[i])));
//  }
//  return result;
//}

} // namespace Spheral
