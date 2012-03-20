//---------------------------------Spheral++----------------------------------//
// ReplaceState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//

#include "ReplaceState.hh"
#include "IncrementState.hh"
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
ReplaceState<Dimension, ValueType>::
ReplaceState():
  FieldUpdatePolicyBase<Dimension, ValueType>() {
}

template<typename Dimension, typename ValueType>
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0) {
}

template<typename Dimension, typename ValueType>
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1) {
}

template<typename Dimension, typename ValueType>
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2) {
}

template<typename Dimension, typename ValueType>
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2,
             const std::string& depend3):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3) {
}

template<typename Dimension, typename ValueType>
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2,
             const std::string& depend3,
             const std::string& depend4):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4) {
}

template<typename Dimension, typename ValueType>
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2,
             const std::string& depend3,
             const std::string& depend4,
             const std::string& depend5):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4, depend5) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
ReplaceState<Dimension, ValueType>::
~ReplaceState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
void
ReplaceState<Dimension, ValueType>::
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
    f(i) = df(i);
  }
}

template<typename Dimension, typename ValueType>
void
ReplaceState<Dimension, ValueType>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {
  IncrementState<Dimension, ValueType> mIncrementStatePolicy;
  mIncrementStatePolicy.update(key, state, derivs, multiplier, t, dt);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
bool
ReplaceState<Dimension, ValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ReplaceState<Dimension, ValueType>* rhsPtr = dynamic_cast<const ReplaceState<Dimension, ValueType>*>(&rhs);
  return rhsPtr != 0;
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class ReplaceState<Dim<1>, Dim<1>::Scalar>;
  template class ReplaceState<Dim<1>, Dim<1>::Vector>;
  template class ReplaceState<Dim<1>, Dim<1>::Vector3d>;
  template class ReplaceState<Dim<1>, Dim<1>::Tensor>;
  template class ReplaceState<Dim<1>, Dim<1>::SymTensor>;
                 
  template class ReplaceState<Dim<2>, Dim<2>::Scalar>;
  template class ReplaceState<Dim<2>, Dim<2>::Vector>;
  template class ReplaceState<Dim<2>, Dim<2>::Vector3d>;
  template class ReplaceState<Dim<2>, Dim<2>::Tensor>;
  template class ReplaceState<Dim<2>, Dim<2>::SymTensor>;
                 
  template class ReplaceState<Dim<3>, Dim<3>::Scalar>;
  template class ReplaceState<Dim<3>, Dim<3>::Vector>;
  template class ReplaceState<Dim<3>, Dim<3>::Vector3d>;
  template class ReplaceState<Dim<3>, Dim<3>::Tensor>;
  template class ReplaceState<Dim<3>, Dim<3>::SymTensor>;
}
