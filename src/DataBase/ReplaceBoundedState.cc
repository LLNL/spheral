//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//

#include "ReplaceBoundedState.hh"
#include "IncrementBoundedState.hh"
#include "FieldUpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Infrastructure/SpheralFunctions.hh"

namespace Spheral {

using Spheral::FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const std::string& depend3,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const std::string& depend3,
                    const std::string& depend4,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const std::string& depend3,
                    const std::string& depend4,
                    const std::string& depend5,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4, depend5),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
~ReplaceBoundedState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
void
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Find the matching derivative field from the StateDerivatives.
  KeyType replaceKey = prefix() + key;
  Field<Dimension, ValueType>& f = state.field(key, ValueType());
  const Field<Dimension, ValueType>& df = derivs.field(replaceKey, ValueType());

  // Loop over the internal values of the field.
  for (int i = 0; i != f.nodeList().numInternalNodes(); ++i) {
    f(i) = min(mMaxValue, max(mMinValue, df(i)));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
void
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {
  IncrementBoundedState<Dimension, ValueType, BoundValueType> mIncrementStatePolicy(mMinValue, mMaxValue);
  mIncrementStatePolicy.update(key, state, derivs, multiplier, t, dt);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
bool
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ReplaceBoundedState<Dimension, ValueType, BoundValueType>* rhsPtr = dynamic_cast<const ReplaceBoundedState<Dimension, ValueType, BoundValueType>*>(&rhs);
  return rhsPtr != 0;
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Scalar>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Vector>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Vector3d>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::Tensor>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::SymTensor>;
  template class ReplaceBoundedState<Dim<1>, Dim<1>::SymTensor, Dim<1>::Scalar>;
                 
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Scalar>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Vector>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Vector3d>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::Tensor>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::SymTensor>;
  template class ReplaceBoundedState<Dim<2>, Dim<2>::SymTensor, Dim<2>::Scalar>;
                 
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Scalar>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Vector>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Vector3d>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::Tensor>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::SymTensor>;
  template class ReplaceBoundedState<Dim<3>, Dim<3>::SymTensor, Dim<3>::Scalar>;
}
