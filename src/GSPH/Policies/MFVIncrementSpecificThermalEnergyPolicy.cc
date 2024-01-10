//---------------------------------Spheral++----------------------------------//
// MFVIncrementSpecificThermalEnergyPolicy -- replaces one fieldlist with the ratio 
//                                     of two fieldlists from the state.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/Policies/MFVIncrementSpecificThermalEnergyPolicy.hh"
#include "GSPH/GSPHFieldNames.hh"
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
template<typename Dimension>
MFVIncrementSpecificThermalEnergyPolicy<Dimension>::
MFVIncrementSpecificThermalEnergyPolicy(std::initializer_list<std::string> depends):
  FieldUpdatePolicy<Dimension>(depends){
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MFVIncrementSpecificThermalEnergyPolicy<Dimension>::
~MFVIncrementSpecificThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MFVIncrementSpecificThermalEnergyPolicy<Dimension>::
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
  const auto derivFieldKey = StateBase<Dimension>::buildFieldKey(prefix() + GSPHFieldNames::thermalEnergy, nodeListKey);

  const auto&  m   = state.field(massKey, Scalar());
        auto&  eps = state.field(key,     Scalar());

  const auto&  DmDt = derivs.field(prefix() + massKey, Scalar());
  const auto&  DmepsDt = derivs.field(derivFieldKey,   Scalar());

// Loop over the internal values of the field.
  const auto n = m.numInternalElements();
#pragma omp parallel for
  for (unsigned i = 0; i != n; ++i) {
    const auto m1 = m(i)+DmDt(i)*multiplier;
    if (m1 > tiny) eps(i) += (DmepsDt(i) - DmDt(i)*eps(i)) * multiplier * safeInv(m1);
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MFVIncrementSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const auto* rhsPtr = dynamic_cast<const MFVIncrementSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  return rhsPtr != nullptr;
}

}

