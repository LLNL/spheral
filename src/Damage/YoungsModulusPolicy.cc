//---------------------------------Spheral++----------------------------------//
// YoungsModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent bulk modulus state.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#include "YoungsModulusPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using std::abs;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YoungsModulusPolicy<Dimension>::
YoungsModulusPolicy():
  UpdatePolicyBase<Dimension>(SolidFieldNames::bulkModulus,
                              SolidFieldNames::shearModulus) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YoungsModulusPolicy<Dimension>::
~YoungsModulusPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
YoungsModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::YoungsModulus);
  Field<Dimension, Scalar>& stateField = state.field(key, 0.0);

  // Get the bulk and shear modulus fields from the state.
  const KeyType KKey = State<Dimension>::buildFieldKey(SolidFieldNames::bulkModulus, nodeListKey);
  const KeyType muKey = State<Dimension>::buildFieldKey(SolidFieldNames::shearModulus, nodeListKey);
  CHECK(state.registered(KKey));
  CHECK(state.registered(muKey));
  const Field<Dimension, Scalar>& K = state.field(KKey, 0.0);
  const Field<Dimension, Scalar>& mu = state.field(muKey, 0.0);

  // Now set Youngs modulus.
  for (auto i = 0u; i != stateField.numInternalElements(); ++i) {
    const double thpt = 9.0*std::abs(K(i)*mu(i));
    const double ack = 3.0*std::abs(K(i)) + std::abs(mu(i)) + 1.0e-50*std::max(1.0, thpt);
    CHECK(ack > 0.0);
    stateField(i) = thpt/ack;
    CHECK2(stateField(i) >= 0.0, stateField(i) << " " << K(i) << " " << mu(i));
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
YoungsModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a Youngs modulus operator.
  const YoungsModulusPolicy<Dimension>* rhsPtr = dynamic_cast<const YoungsModulusPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

