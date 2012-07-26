//---------------------------------Spheral++----------------------------------//
// ShearModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent shear modulus state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#include "ShearModulusPolicy.hh"
#include "SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using FieldSpace::Field;
using NodeSpace::NodeList;
using SolidMaterial::SolidNodeList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ShearModulusPolicy<Dimension>::
ShearModulusPolicy():
  UpdatePolicyBase<Dimension>(HydroFieldNames::massDensity,
                              HydroFieldNames::specificThermalEnergy,
                              HydroFieldNames::pressure) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ShearModulusPolicy<Dimension>::
~ShearModulusPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ShearModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::shearModulus);
  Field<Dimension, Scalar>& stateField = state.field(key, 0.0);

  // We only do this if this is a solid node list.
  const SolidNodeList<Dimension>* nodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(stateField.nodeListPtr());
  if (nodeListPtr != 0) {

    // Get the mass density, specific thermal energy, and pressure fields
    // from the state.
    const KeyType massDensityKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
    const KeyType energyKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
    const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
    CHECK(state.registered(massDensityKey));
    CHECK(state.registered(energyKey));
    CHECK(state.registered(PKey));
    const Field<Dimension, Scalar>& massDensity = state.field(massDensityKey, 0.0);
    const Field<Dimension, Scalar>& energy = state.field(energyKey, 0.0);
    const Field<Dimension, Scalar>& P = state.field(PKey, 0.0);
    
    // Get the strength model.
    const SolidMaterial::StrengthModel& strengthModel = nodeListPtr->strengthModel();

    // Now set the shear modulus.
    for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
      stateField(i) = strengthModel.shearModulus(massDensity(i), energy(i), P(i));
    }

//     // Is there a scalar damage field for this NodeList?
//     {
//       const KeyType DKey(nodeListPtr, SolidFieldNames::scalarDamage);
//       if (state.scalarFieldRegistered(DKey)) {
//         const Field<Dimension, Scalar>& D = state.scalarField(DKey);
//         for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
//           CHECK(D(i) >= 0.0 && D(i) <= 1.0);
//           stateField(i) *= 1.0 - D(i);
//         }
//       }
//     }

//     // Is there a tensor damage field for this NodeList?
//     {
//       typedef typename Dimension::SymTensor SymTensor;
//       const KeyType DKey(nodeListPtr, SolidFieldNames::tensorDamage);
//       if (state.symTensorFieldRegistered(DKey)) {
//         const Field<Dimension, SymTensor>& D = state.symTensorField(DKey);
//         for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
//           const Scalar Di = std::max(0.0, std::min(1.0, D(i).eigenValues().maxElement() + 1.0e-5));
//           CHECK(Di >= 0.0 && Di <= 1.0);
//           stateField(i) *= 1.0 - Di;
//         }
//       }
//     }

  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ShearModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ShearModulusPolicy<Dimension>* rhsPtr = dynamic_cast<const ShearModulusPolicy<Dimension>*>(&rhs);
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
  template class ShearModulusPolicy<Dim<1> >;
  template class ShearModulusPolicy<Dim<2> >;
  template class ShearModulusPolicy<Dim<3> >;
}

