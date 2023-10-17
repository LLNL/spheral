//---------------------------------Spheral++----------------------------------//
// DamagedPressurePolicy -- Override the default pressure policy in the 
// presence of damage.
//----------------------------------------------------------------------------//
#include "DamagedPressurePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/SolidEquationOfState.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamagedPressurePolicy<Dimension>::
DamagedPressurePolicy():
  PressurePolicy<Dimension>() {
  this->addDependency(SolidFieldNames::tensorDamage);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
DamagedPressurePolicy<Dimension>::
~DamagedPressurePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DamagedPressurePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey, Dkey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::pressure);

  // Get our state
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  auto&       pressure = state.field(key, 0.0);
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero);

  // Have the base class set the initial pressure.
  PressurePolicy<Dimension>::update(key, state, derivs, multiplier, t, dt);

  // Scale by the damage.
  const auto ni = pressure.numInternalElements();
  auto Pmin = 0.0;
  bool enforceDamagedPmin = false;
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(pressure.nodeListPtr());
  if (solidNodeListPtr != nullptr) {
    const auto* eosPtr = dynamic_cast<const SolidEquationOfState<Dimension>*>(&(solidNodeListPtr->equationOfState()));
    if (eosPtr != nullptr) {
      Pmin = eosPtr->minimumPressureDamage();
      enforceDamagedPmin = true;
    }
  }
  if (enforceDamagedPmin) {
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto Di = std::max(0.0, std::min(1.0, D(i).eigenValues().maxElement()));
      CHECK(Di >= 0.0 and Di <= 1.0);
      pressure(i) = (1.0 - Di)*pressure(i) + Di*std::max(pressure(i), Pmin);
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
DamagedPressurePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const DamagedPressurePolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}
