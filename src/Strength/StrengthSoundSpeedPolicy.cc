//---------------------------------Spheral++----------------------------------//
// StrengthSoundSpeedPolicy -- Override the default sound speed policy in the 
// presence of strength.
//----------------------------------------------------------------------------//
#include "StrengthSoundSpeedPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Strength/SolidNodeList.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using SolidMaterial::SolidNodeList;
using FieldSpace::Field;
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
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed);

  // Get the density, energy, and pressure fields from the state.
  const KeyType rhoKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
  const KeyType epsKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  CHECK(state.registered(rhoKey));
  CHECK(state.registered(epsKey));
  CHECK(state.registered(PKey));
  Field<Dimension, Scalar>& soundSpeed = state.field(key, 0.0);
  const Field<Dimension, Scalar>& rho = state.field(rhoKey, 0.0);
  const Field<Dimension, Scalar>& eps = state.field(epsKey, 0.0);
  const Field<Dimension, Scalar>& P = state.field(PKey, 0.0);

  // Get the solid node list and strength model.
  const SolidNodeList<Dimension>* nodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(soundSpeed.nodeListPtr());
  REQUIRE(nodeListPtr != 0);
  const StrengthModel& strengthModel = nodeListPtr->strengthModel();

  // Have the base class set the initial sound speed.
  SoundSpeedPolicy<Dimension>::update(key, state, derivs, multiplier, t, dt);

  // Iterate over the nodes, and augment the sound speed by the value including
  // strength.
  for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
    soundSpeed(i) = strengthModel.soundSpeed(rho(i), eps(i), P(i), soundSpeed(i));
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

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class StrengthSoundSpeedPolicy<Dim<1> >;
  template class StrengthSoundSpeedPolicy<Dim<2> >;
  template class StrengthSoundSpeedPolicy<Dim<3> >;
}

