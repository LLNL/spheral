//---------------------------------Spheral++----------------------------------//
// CompatibleFaceSpecificThermalEnergyPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#include <vector>

#include "CompatibleFaceSpecificThermalEnergyPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/FieldUpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
CompatibleFaceSpecificThermalEnergyPolicy():
  IncrementState<Dimension, typename Dimension::Scalar>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
~CompatibleFaceSpecificThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy);

  // Get the state fields.
  // const KeyType eps0Key = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy + "0", nodeListKey);
  const KeyType massKey = State<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const KeyType velfKey = State<Dimension>::buildFieldKey(HydroFieldNames::faceVelocity, nodeListKey);
  const KeyType forcefKey = State<Dimension>::buildFieldKey(HydroFieldNames::faceForce, nodeListKey);
  const KeyType massfKey = State<Dimension>::buildFieldKey(HydroFieldNames::faceMass, nodeListKey);
  // CHECK(state.registered(eps0Key));
  CHECK(state.registered(massKey));
  CHECK(derivs.registered(velfKey));
  CHECK(derivs.registered(forcefKey));
  CHECK(derivs.registered(massfKey));
  Field<Dimension, Scalar>& eps = state.field(key, 0.0);
  // const Field<Dimension, Scalar>& eps0 = state.field(eps0Key, 0.0);
  const Field<Dimension, Scalar>& mass = state.field(massKey, 0.0);
  const Field<Dimension, vector<Vector> >& vel_f = derivs.field(velfKey, vector<Vector>());
  const Field<Dimension, vector<Vector> >& force_f = derivs.field(forcefKey, vector<Vector>());
  const Field<Dimension, vector<Scalar> >& mass_f = derivs.field(massfKey, vector<Scalar>());
  
  const unsigned n = eps.numInternalElements();
  const double hdt = 0.5*multiplier;
  double depsi, mcheck;
  for (unsigned i = 0; i != n; ++i) {
    depsi = 0.0;
    mcheck = 0.0;
    const vector<Vector>& veli_f = vel_f[i];
    const vector<Vector>& fi_f = force_f[i];
    const vector<Scalar>& massi_f = mass_f[i];
    const unsigned nfaces = veli_f.size();
    CHECK(fi_f.size() == nfaces);
    CHECK(massi_f.size() == nfaces);
    for (unsigned k = 0; k != nfaces; ++k) {
      CHECK(massi_f[k] > 0.0);
      depsi += (veli_f[k] + hdt*fi_f[k]/massi_f[k]).dot(fi_f[k]);
      mcheck += massi_f[k];
    }
    CHECK(fuzzyEqual(mcheck, mass[i]));
    CHECK(mass[i] > 0.0);
    depsi *= multiplier/mass[i];
    eps[i] += depsi;
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a CompatibleSpecificThermalEnergy operator.
  const CompatibleFaceSpecificThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const CompatibleFaceSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

