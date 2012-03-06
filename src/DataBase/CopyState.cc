//---------------------------------Spheral++----------------------------------//
// CopyState -- An implementation of UpdatePolicyBase appropriate for
// copying one state field to another.
// Assumes that the copied state is dependent upon the master state, *and* that
// that is the only dependency.
//
// Created by JMO, Tue Oct 5 11:08:48 2004
//----------------------------------------------------------------------------//

#include "CopyState.hh"
#include "FieldUpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
CopyState<Dimension, ValueType>::
CopyState(const std::string& masterState,
          const std::string& copyState):
  FieldUpdatePolicyBase<Dimension, ValueType>(),
  mMasterStateName(masterState),
  mCopyStateName(copyState) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
CopyState<Dimension, ValueType>::
~CopyState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
void
CopyState<Dimension, ValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  
  // Find the field for this key.
  Field<Dimension, ValueType>& f = state.field(key, ValueType());

  // Find the master state field from the State.
  const Field<Dimension, ValueType>& masterField = state.field(mMasterStateName, ValueType());

  // Copy the master state.
  f = masterField;
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
bool
CopyState<Dimension, ValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const CopyState<Dimension, ValueType>* rhsPtr = dynamic_cast<const CopyState<Dimension, ValueType>*>(&rhs);
  if (rhsPtr == 0) return false;
  return (mMasterStateName == rhsPtr->mMasterStateName && 
          mCopyStateName == rhsPtr->mCopyStateName);
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class CopyState<Dim<1>, Dim<1>::Scalar>;
  template class CopyState<Dim<1>, Dim<1>::Vector>;
  template class CopyState<Dim<1>, Dim<1>::Vector3d>;
  template class CopyState<Dim<1>, Dim<1>::Tensor>;
  template class CopyState<Dim<1>, Dim<1>::SymTensor>;
                 
  template class CopyState<Dim<2>, Dim<2>::Scalar>;
  template class CopyState<Dim<2>, Dim<2>::Vector>;
  template class CopyState<Dim<2>, Dim<2>::Vector3d>;
  template class CopyState<Dim<2>, Dim<2>::Tensor>;
  template class CopyState<Dim<2>, Dim<2>::SymTensor>;
                 
  template class CopyState<Dim<3>, Dim<3>::Scalar>;
  template class CopyState<Dim<3>, Dim<3>::Vector>;
  template class CopyState<Dim<3>, Dim<3>::Vector3d>;
  template class CopyState<Dim<3>, Dim<3>::Tensor>;
  template class CopyState<Dim<3>, Dim<3>::SymTensor>;
}
