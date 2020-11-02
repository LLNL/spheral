//---------------------------------Spheral++----------------------------------//
// MeltEnergyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent melt energy state.
//
// Created by JMO, Fri Jul 13 12:57:54 PDT 2018
//----------------------------------------------------------------------------//
#include "MeltEnergyPolicy.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MeltEnergyPolicy<Dimension>::
MeltEnergyPolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
                                                                   HydroFieldNames::specificThermalEnergy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MeltEnergyPolicy<Dimension>::
~MeltEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MeltEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::meltSpecificEnergy);
  auto& melt = state.field(key, 0.0);

  // Get the mass density and specific thermal energy fields from the state.
  const auto rhoKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
  const auto epsKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  CHECK(state.registered(rhoKey));
  CHECK(state.registered(epsKey));
  const auto& rho = state.field(rhoKey, 0.0);
  const auto& eps = state.field(epsKey, 0.0);

  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(melt.nodeListPtr());
  CHECK(solidNodeListPtr != 0);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Set the melt.
  strengthModel.meltSpecificEnergy(melt, rho, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MeltEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a bulk modulus operator.
  const MeltEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const MeltEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

