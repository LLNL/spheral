//---------------------------------Spheral++----------------------------------//
// StrengthSoundSpeedPolicy -- Override the default sound speed policy in the 
// presence of strength.
//----------------------------------------------------------------------------//
#include "StrengthSoundSpeedPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NodeSpace::SolidNodeList;
using FieldSpace::FieldList;
using SolidMaterial::StrengthModel;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrengthSoundSpeedPolicy<Dimension>::
StrengthSoundSpeedPolicy():
  SoundSpeedPolicy<Dimension>() {
  this->addDependency(HydroFieldNames::pressure);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrengthSoundSpeedPolicy<Dimension>::
~StrengthSoundSpeedPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrengthSoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> stateFields = state.fields(fieldKey, Scalar());
  const unsigned numFields = stateFields.numFields();

  // Have the base class set the initial sound speed.
  SoundSpeedPolicy<Dimension>::update(key, state, derivs, multiplier, t, dt);

  // Get the density, energy, and pressure fields from the state.
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);

  // Walk the individual fields.
  for (unsigned k = 0; k != numFields; ++k) {

    // Get the strength model.  This cast is ugly, but is a work-around for now.
    const SolidNodeList<Dimension>* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(stateFields[k]->nodeListPtr());
    CHECK(solidNodeListPtr != 0);
    const SolidMaterial::StrengthModel<Dimension>& strengthModel = solidNodeListPtr->strengthModel();

    // Set the full sound speed.
    strengthModel.soundSpeed(*stateFields[k], *rho[k], *eps[k], *P[k], *stateFields[k]);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StrengthSoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const StrengthSoundSpeedPolicy<Dimension>* rhsPtr = dynamic_cast<const StrengthSoundSpeedPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

