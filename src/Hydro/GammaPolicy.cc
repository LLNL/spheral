//---------------------------------Spheral++----------------------------------//
// GammaPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent gamma state.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//

#include "GammaPolicy.hh"
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


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GammaPolicy<Dimension>::
GammaPolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
                                                                   HydroFieldNames::specificThermalEnergy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GammaPolicy<Dimension>::
~GammaPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::gamma and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> gamma = state.fields(fieldKey, Scalar());
  const unsigned numFields = gamma.numFields();

  // Get the mass density and specific thermal energy fields from the state.
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, Scalar());
  const FieldList<Dimension, Scalar> energy = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  CHECK(massDensity.numFields() == numFields);
  CHECK(energy.numFields() == numFields);

  // Walk the fields.
  for (unsigned i = 0; i != numFields; ++i) {

    // Get the eos.  This cast is ugly, but is a work-around for now.
    const FluidNodeList<Dimension>* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(gamma[i]->nodeListPtr());
    CHECK(fluidNodeListPtr != 0);
    const EquationOfState<Dimension>& eos = fluidNodeListPtr->equationOfState();

    // Now set the gamma for this field.
    eos.setGammaField(*gamma[i], *massDensity[i], *energy[i]);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GammaPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const GammaPolicy<Dimension>* rhsPtr = dynamic_cast<const GammaPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

