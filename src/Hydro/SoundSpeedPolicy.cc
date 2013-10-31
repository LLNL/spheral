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
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::Field;
using FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SoundSpeedPolicy<Dimension>::
SoundSpeedPolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
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
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> soundSpeed = state.fields(fieldKey, Scalar());
  const unsigned numFields = soundSpeed.numFields();

  // Get the mass density and specific thermal energy fields from the state.
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, Scalar());
  const FieldList<Dimension, Scalar> energy = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  CHECK(massDensity.numFields() == numFields);
  CHECK(energy.numFields() == numFields);

  // Walk the fields.
  for (unsigned i = 0; i != numFields; ++i) {

    // Get the eos.  This cast is ugly, but is a work-around for now.
    const FluidNodeList<Dimension>* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(soundSpeed[i]->nodeListPtr());
    CHECK(fluidNodeListPtr != 0);
    const Material::EquationOfState<Dimension>& eos = fluidNodeListPtr->equationOfState();

    // Now set the soundSpeed for this field.
    eos.setSoundSpeed(*soundSpeed[i], *massDensity[i], *energy[i]);
  }
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

