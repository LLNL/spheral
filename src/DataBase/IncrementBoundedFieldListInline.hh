//---------------------------------Spheral++----------------------------------//
// IncrementBoundedFieldList -- An implementation of FieldListUpdatePolicyBase 
// appropriate for when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt.
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "UpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>() {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementBoundedState<Dimension, Value>>(minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementBoundedState<Dimension, Value>>(depend0,
                                                                                                        minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0,
                          const std::string& depend1,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementBoundedState<Dimension, Value>>(depend0,
                                                                                                        depend1,
                                                                                                        minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementBoundedState<Dimension, Value>>(depend0,
                                                                                                        depend1,
                                                                                                        depend2,
                                                                                                        minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const std::string& depend3,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementBoundedState<Dimension, Value>>(depend0,
                                                                                                        depend1,
                                                                                                        depend2,
                                                                                                        depend3,
                                                                                                        minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const std::string& depend3,
                          const std::string& depend4,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementBoundedState<Dimension, Value>>(depend0,
                                                                                                        depend1,
                                                                                                        depend2,
                                                                                                        depend3,
                                                                                                        depend4,
                                                                                                        minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const FieldList<Dimension, ValueType>& fieldList,
                          const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const std::string& depend3,
                          const std::string& depend4,
                          const std::string& depend5,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4, depend5) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementBoundedState<Dimension, Value>>(depend0,
                                                                                                        depend1,
                                                                                                        depend2,
                                                                                                        depend3,
                                                                                                        depend4,
                                                                                                        depend5,
                                                                                                        minValue, maxValue));
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
~IncrementBoundedFieldList() {
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
bool
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>* rhsPtr = dynamic_cast<const IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>*>(&rhs);
  if (rhsPtr == 0) return false;

  // Ok, now do we agree on min & max?
  return (minValue() == rhsPtr->minValue()) && (maxValue() == rhsPtr->maxValue());
}

}
