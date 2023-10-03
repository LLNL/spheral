//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//
#include "ReplaceBoundedFieldList.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const FieldList<Dimension, Value>& fieldList,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>() {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceBoundedState<Dimension, Value>>(minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const FieldList<Dimension, Value>& fieldList,
                        const std::string& depend0,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceBoundedState<Dimension, Value>>(depend0,
                                                                                                      minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const FieldList<Dimension, Value>& fieldList,
                        const std::string& depend0,
                        const std::string& depend1,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceBoundedState<Dimension, Value>>(depend0, depend1,
                                                                                                      minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const FieldList<Dimension, Value>& fieldList,
                        const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceBoundedState<Dimension, Value>>(depend0, depend1, depend2,
                                                                                                      minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const FieldList<Dimension, Value>& fieldList,
                        const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2,
                        const std::string& depend3,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceBoundedState<Dimension, Value>>(depend0, depend1, depend2, depend3,
                                                                                                      minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const FieldList<Dimension, Value>& fieldList,
                        const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2,
                        const std::string& depend3,
                        const std::string& depend4,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceBoundedState<Dimension, Value>>(depend0, depend1, depend2, depend3, depend4,
                                                                                                      minValue, maxValue));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const FieldList<Dimension, Value>& fieldList,
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
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceBoundedState<Dimension, Value>>(depend0, depend1, depend2, depend3, depend4, depend5,
                                                                                                      minValue, maxValue));
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
~ReplaceBoundedFieldList() {
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
bool
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>* rhsPtr = dynamic_cast<const ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>*>(&rhs);
  return rhsPtr != 0;
}

}
