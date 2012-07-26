//---------------------------------Spheral++----------------------------------//
// IncrementBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//

#include "IncrementBoundedState.hh"
#include "UpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using Spheral::FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const std::string& depend1,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const std::string& depend1,
                      const std::string& depend2,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
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
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
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
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
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
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
~IncrementBoundedState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
void
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Find the matching derivative field from the StateDerivatives.
  KeyType incrementKey = prefix() + key;
  Field<Dimension, ValueType>& f = state.field(key, ValueType());
  const Field<Dimension, ValueType>& df = derivs.field(incrementKey, ValueType());

  // Loop over the internal values of the field.
  for (int i = 0; i != f.nodeList().numInternalNodes(); ++i) {
    f(i) = min(mMaxValue, max(mMinValue, f(i) + multiplier*(df(i))));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
bool
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementBoundedState<Dimension, ValueType, BoundValueType>* rhsPtr = dynamic_cast<const IncrementBoundedState<Dimension, ValueType, BoundValueType>*>(&rhs);
  if (rhsPtr == 0) return false;

  // Ok, now do we agree on min & max?
  return (minValue() == rhsPtr->minValue()) && (maxValue() == rhsPtr->maxValue());
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class IncrementBoundedState<Dim<1>, Dim<1>::Scalar>;
  template class IncrementBoundedState<Dim<1>, Dim<1>::Vector>;
  template class IncrementBoundedState<Dim<1>, Dim<1>::Vector3d>;
  template class IncrementBoundedState<Dim<1>, Dim<1>::Tensor>;
  template class IncrementBoundedState<Dim<1>, Dim<1>::SymTensor>;
  template class IncrementBoundedState<Dim<1>, Dim<1>::SymTensor, Dim<1>::Scalar>;

  template class IncrementBoundedState<Dim<2>, Dim<2>::Scalar>;
  template class IncrementBoundedState<Dim<2>, Dim<2>::Vector>;
  template class IncrementBoundedState<Dim<2>, Dim<2>::Vector3d>;
  template class IncrementBoundedState<Dim<2>, Dim<2>::Tensor>;
  template class IncrementBoundedState<Dim<2>, Dim<2>::SymTensor>;
  template class IncrementBoundedState<Dim<2>, Dim<2>::SymTensor, Dim<2>::Scalar>;

  template class IncrementBoundedState<Dim<3>, Dim<3>::Scalar>;
  template class IncrementBoundedState<Dim<3>, Dim<3>::Vector>;
  template class IncrementBoundedState<Dim<3>, Dim<3>::Vector3d>;
  template class IncrementBoundedState<Dim<3>, Dim<3>::Tensor>;
  template class IncrementBoundedState<Dim<3>, Dim<3>::SymTensor>;
  template class IncrementBoundedState<Dim<3>, Dim<3>::SymTensor, Dim<3>::Scalar>;
}
