//---------------------------------Spheral++----------------------------------//
// LongitudinalSoundSpeedPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent longitudinal sound speed.
//
// Created by JMO, Sun Jan 23 14:14:24 PST 2005
//----------------------------------------------------------------------------//
#include "LongitudinalSoundSpeedPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
LongitudinalSoundSpeedPolicy<Dimension>::
LongitudinalSoundSpeedPolicy():
  UpdatePolicyBase<Dimension>(SolidFieldNames::YoungsModulus,
                              SolidFieldNames::bulkModulus,
                              SolidFieldNames::shearModulus,
                              HydroFieldNames::massDensity) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
LongitudinalSoundSpeedPolicy<Dimension>::
~LongitudinalSoundSpeedPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LongitudinalSoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::longitudinalSoundSpeed);
  Field<Dimension, Scalar>& stateField = state.field(key, 0.0);

  // Get the bulk, shear, and Youngs modulus, and the mass density fields 
  // from the state.
  const KeyType rhoKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
  const KeyType EKey = State<Dimension>::buildFieldKey(SolidFieldNames::YoungsModulus, nodeListKey);
  const KeyType KKey = State<Dimension>::buildFieldKey(SolidFieldNames::bulkModulus, nodeListKey);
  const KeyType muKey = State<Dimension>::buildFieldKey(SolidFieldNames::shearModulus, nodeListKey);
  CHECK(state.registered(rhoKey));
  CHECK(state.registered(EKey));
  CHECK(state.registered(KKey));
  CHECK(state.registered(muKey));
  const Field<Dimension, Scalar>& rho = state.field(rhoKey, 0.0);
  const Field<Dimension, Scalar>& E = state.field(EKey, 0.0);
  const Field<Dimension, Scalar>& K = state.field(KKey, 0.0);
  const Field<Dimension, Scalar>& mu = state.field(muKey, 0.0);

  // Now set the longitudinal sound speed.
  for (auto i = 0u; i != stateField.numInternalElements(); ++i) {
    const double ack = 3.0*K(i) + mu(i) + 1.0e-30*std::max(1.0, K(i));
    CHECK2(ack > 0.0, nodeListKey << " : " << i << " " << ack << " " << K(i) << " " << mu(i));
    const double nu = min(0.5, max(0.0, 0.5*(3.0*K(i) - 2.0*mu(i))/ack));
    CHECK(nu >= 0.0 && nu <= 0.5);
    const double barf = (1.0 + nu)*(1.0 - 2.0*nu) + 1.0e-10;
    CHECK(distinctlyGreaterThan(barf, 0.0));
    CHECK(distinctlyGreaterThan(rho(i), 0.0));
    stateField(i) = sqrt(abs(E(i)*(1.0 - nu)/(rho(i)*barf)));
    CHECK(stateField(i) >= 0.0);
  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
LongitudinalSoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a Youngs modulus operator.
  const LongitudinalSoundSpeedPolicy<Dimension>* rhsPtr = dynamic_cast<const LongitudinalSoundSpeedPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

