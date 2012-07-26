//---------------------------------Spheral++----------------------------------//
// SoundSpeedPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the node weight as a dependent quantity.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//

#include "SoundSpeedPolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SoundSpeedPolicy<Dimension>::
SoundSpeedPolicy():
  FieldUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
                                                               HydroFieldNames::specificThermalEnergy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SoundSpeedPolicy<Dimension>::
~SoundSpeedPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed);

  // Get the mass density and specific thermal energy fields from the state.
  const KeyType massDensityKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
  const KeyType energyKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  CHECK(state.registered(massDensityKey));
  CHECK(state.registered(energyKey));
  Field<Dimension, Scalar>& soundSpeed = state.field(key, 0.0);
  const Field<Dimension, Scalar>& massDensity = state.field(massDensityKey, 0.0);
  const Field<Dimension, Scalar>& energy = state.field(energyKey, 0.0);

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const FluidNodeList<Dimension>* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(soundSpeed.nodeListPtr());
  CHECK(fluidNodeListPtr != 0);
  const Material::EquationOfState<Dimension>& eos = fluidNodeListPtr->equationOfState();

  // Now set the sound speed.
  eos.setSoundSpeed(soundSpeed, massDensity, energy);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const SoundSpeedPolicy<Dimension>* rhsPtr = dynamic_cast<const SoundSpeedPolicy<Dimension>*>(&rhs);
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
  template class SoundSpeedPolicy<Dim<1> >;
  template class SoundSpeedPolicy<Dim<2> >;
  template class SoundSpeedPolicy<Dim<3> >;
}

