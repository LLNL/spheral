//---------------------------------Spheral++----------------------------------//
// IncrementSpecificFromTotalPolicy -- replaces one fieldlist with the ratio 
//                                     of two fieldlists from the state.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/Policies/IncrementSpecificFromTotalPolicy.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Hydro/HydroFieldNames.hh"

#include <limits.h>

namespace Spheral {
//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(std::initializer_list<std::string> depends,const std::string& stateKey, const std::string& derivKey):
  UpdatePolicyBase<Dimension>(depends),
  mStateKey(stateKey),
  mDerivativeKey(derivKey){
}

template<typename Dimension, typename Value>
IncrementSpecificFromTotalPolicy<Dimension, Value>::
IncrementSpecificFromTotalPolicy(const std::string& stateKey, const std::string& derivKey):
  UpdatePolicyBase<Dimension>(),
  mStateKey(stateKey),
  mDerivativeKey(derivKey){
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
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  const auto massKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const auto stateFieldKey = StateBase<Dimension>::buildFieldKey(mStateKey, nodeListKey);
  const auto derivFieldKey = StateBase<Dimension>::buildFieldKey(mDerivativeKey, nodeListKey);

  const auto  m = state.field(massKey,       Scalar());
        auto  q = state.field(stateFieldKey, Value());

  const auto  DmDt = derivs.field(prefix() + massKey, Scalar());
  const auto  DQDt = derivs.field(derivFieldKey,      Value());

// Loop over the internal values of the field.
  const auto n = m.numInternalElements();
#pragma omp parallel for
  for (unsigned i = 0; i != n; ++i) {
    const auto m1 = m(i)+DmDt(i)*multiplier;
    if (m1 > tiny) q(i) += (DQDt(i) - DmDt(i)*q(i)) * multiplier * safeInv(m1);
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
  const auto* rhsPtr = dynamic_cast<const IncrementSpecificFromTotalPolicy<Dimension, Value>*>(&rhs);
  return rhsPtr != nullptr;
}

}

