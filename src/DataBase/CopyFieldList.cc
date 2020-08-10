//---------------------------------Spheral++----------------------------------//
// CopyFieldList -- An implementation of UpdatePolicyBase appropriate for
// copying one state field to another.
// Assumes that the copied state is dependent upon the master state, *and* that
// that is the only dependency.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "CopyFieldList.hh"
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
CopyFieldList<Dimension, ValueType>::
CopyFieldList(const std::string& masterState,
              const std::string& copyState):
  FieldListUpdatePolicyBase<Dimension, ValueType>(),
  mMasterStateName(masterState),
  mCopyStateName(copyState) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
CopyFieldList<Dimension, ValueType>::
~CopyFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
void
CopyFieldList<Dimension, ValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == mCopyStateName);
  REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the field for this key.
  FieldList<Dimension, ValueType> f = state.fields(key, ValueType());

  // Find the master state field from the State.
  const FieldList<Dimension, ValueType> masterField = state.fields(mMasterStateName, ValueType());

  // Copy the master state.
  f.assignFields(masterField);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
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

