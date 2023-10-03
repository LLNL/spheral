//---------------------------------Spheral++----------------------------------//
// CopyFieldList -- An implementation of UpdatePolicyBase appropriate for
// copying one state FieldList to another.
// Assumes that the copied state is dependent upon the master state, *and* that
// that is the only dependency.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "FieldListUpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
CopyFieldList<Dimension, ValueType>::
CopyFieldList(const FieldList<Dimension, Value>& fieldList,
              const std::string& masterState,
              const std::string& copyState):
  FieldListUpdatePolicyBase<Dimension, ValueType>(masterState),
  mMasterStateName(masterState),
  mCopyStateName(copyState) {
  for (const auto& field: fieldList) {
    const KeyType key = StateBase::key(field);
    KeyType fieldKey, nodeListKey;
    StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
    REQUIRE(fieldKey == copyState);
    this->enroll(field.nodeList().name(),
                 std::make_shared<CopyState<Dimension, Field<Dimension, Value>>>(StateBase::buildFieldKey(masterState, nodeListKey), key));
  }
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
CopyFieldList<Dimension, ValueType>::
~CopyFieldList() {
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
bool
CopyFieldList<Dimension, ValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const CopyFieldList<Dimension, ValueType>* rhsPtr = dynamic_cast<const CopyFieldList<Dimension, ValueType>*>(&rhs);
  if (rhsPtr == 0) return false;
  return (mMasterStateName == rhsPtr->mMasterStateName && 
          mCopyStateName == rhsPtr->mCopyStateName);
}

}

