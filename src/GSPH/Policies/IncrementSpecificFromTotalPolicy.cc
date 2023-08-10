//---------------------------------Spheral++----------------------------------//
// IncrementSpecificFromTotalPolicy -- replaces one fieldlist with the ratio 
//                                     of two fieldlists from the state.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/Policies/IncrementSpecificFromTotalPolicy.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"

#include <limits.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivsKey):
  FieldListUpdatePolicyBase<Dimension, Value>(),
  mStateKey(stateKey),
  mDerivativeKey(derivsKey){
}

template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string&  derivsKey, const std::string& depend0):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0),
  mStateKey(stateKey),
  mDerivativeKey(derivsKey) {
}

template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(const std::string& stateKey,
                              const std::string& derivsKey,
                              const std::string& depend0,
                              const std::string& depend1):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1),
  mStateKey(stateKey),
  mDerivativeKey(derivsKey) {
}

template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(const std::string& stateKey,
                              const std::string& derivsKey,
                              const std::string& depend0,
                              const std::string& depend1,
                              const std::string& depend2):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2),
  mStateKey(stateKey),
  mDerivativeKey(derivsKey){
}

template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(const std::string& stateKey,
                              const std::string& derivsKey,
                              const std::string& depend0,
                              const std::string& depend1,
                              const std::string& depend2,
                              const std::string& depend3):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3),
  mStateKey(stateKey),
  mDerivativeKey(derivsKey){
}

template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(const std::string& stateKey,
                              const std::string& derivsKey,
                              const std::string& depend0,
                              const std::string& depend1,
                              const std::string& depend2,
                              const std::string& depend3,
                              const std::string& depend4):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4),
  mStateKey(stateKey),
  mDerivativeKey(derivsKey){
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
~IncrementSpecificFromTotalPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
IncrementSpecificFromTotalPolicy<Dimension, Value>::
update(const KeyType& /*key*/,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();

  const auto  m = state.fields(HydroFieldNames::mass, Scalar());
        auto  q = state.fields(mStateKey, Value());

  const auto  DmDt = derivs.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::mass, Scalar());
  const auto  DQDt = derivs.fields(mDerivativeKey, Value());

  // Loop over the internal values of the field.
  const unsigned numNodeLists = q.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = q[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto m1 = m(k,i)+DmDt(k,i)*multiplier;
      if (m1 > tiny) q(k, i) += (DQDt(k, i) - DmDt(k, i)*q(k, i)) * multiplier * safeInv(m1);
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
IncrementSpecificFromTotalPolicy<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const IncrementSpecificFromTotalPolicy<Dimension, Value>* rhsPtr = dynamic_cast<const IncrementSpecificFromTotalPolicy<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

