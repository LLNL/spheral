//---------------------------------Spheral++----------------------------------//
// IncrementFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <regex>
#include <vector>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>() {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(wildCardDerivs));
  }
}

template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const std::string& depend0,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(depend0,
                                                                                                 wildCardDerivs));
  }
}

template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const std::string& depend0,
                   const std::string& depend1,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(depend0,
                                                                                                 depend1,
                                                                                                 wildCardDerivs));
  }
}

template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(depend0,
                                                                                                 depend1,
                                                                                                 depend2,
                                                                                                 wildCardDerivs));
  }
}

template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(depend0,
                                                                                                 depend1,
                                                                                                 depend2,
                                                                                                 depend3,
                                                                                                 wildCardDerivs));
  }
}

template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(depend0,
                                                                                                 depend1,
                                                                                                 depend2,
                                                                                                 depend3,
                                                                                                 depend4,
                                                                                                 wildCardDerivs));
  }
}

template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4,
                   const std::string& depend5,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(depend0,
                                                                                                 depend1,
                                                                                                 depend2,
                                                                                                 depend3,
                                                                                                 depend4,
                                                                                                 depend5,
                                                                                                 wildCardDerivs));
  }
}

template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const FieldList<Dimension, Value>& fieldList,
                   const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4,
                   const std::string& depend5,
                   const std::string& depend6,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5) {
  for (const auto* fieldPtr: fieldList) {
    this->enroll(fieldPtr->nodeList().name(), std::make_shared<IncrementState<Dimension, Value>>(depend0,
                                                                                                 depend1,
                                                                                                 depend2,
                                                                                                 depend3,
                                                                                                 depend4,
                                                                                                 depend5,
                                                                                                 depend6,
                                                                                                 wildCardDerivs));
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
IncrementFieldList<Dimension, Value>::
~IncrementFieldList() {
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
bool
IncrementFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const IncrementFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0 and mWildCardDerivs == rhsPtr->mWildCardDerivs;
}

}

