//---------------------------------Spheral++----------------------------------//
// ReplaceFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "ReplaceFieldList.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const FieldList<Dimension, Value>& fieldList):
  FieldListUpdatePolicyBase<Dimension, Value>() {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceState<Dimension, Value>>());
  }
}

template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const FieldList<Dimension, Value>& fieldList,
                 const std::string& depend0):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceState<Dimension, Value>>(depend0));
  }
}

template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const FieldList<Dimension, Value>& fieldList,
                 const std::string& depend0,
                 const std::string& depend1):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceState<Dimension, Value>>(depend0, depend1));
  }
}

template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const FieldList<Dimension, Value>& fieldList,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceState<Dimension, Value>>(depend0, depend1, depend2));
  }
}

template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const FieldList<Dimension, Value>& fieldList,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceState<Dimension, Value>>(depend0, depend1, depend2, depend3));
  }
}

template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const FieldList<Dimension, Value>& fieldList,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceState<Dimension, Value>>(depend0, depend1, depend2, depend3, depend4));
  }
}

template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const FieldList<Dimension, Value>& fieldList,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4,
                 const std::string& depend5):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<ReplaceState<Dimension, Value>>(depend0, depend1, depend2, depend3, depend4, depend5));
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
ReplaceFieldList<Dimension, Value>::
~ReplaceFieldList() {
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
bool
ReplaceFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const ReplaceFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const ReplaceFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

