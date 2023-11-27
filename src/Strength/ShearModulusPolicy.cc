//---------------------------------Spheral++----------------------------------//
// ShearModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent shear modulus state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#include "ShearModulusPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ShearModulusPolicy<Dimension>::
ShearModulusPolicy():
  FieldUpdatePolicy<Dimension>({HydroFieldNames::massDensity,
                                HydroFieldNames::specificThermalEnergy,
                                HydroFieldNames::pressure,
                                SolidFieldNames::tensorDamage}) {
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
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::shearModulus);
  auto& mu = state.field(key, 0.0);

  // Get the mass density, specific thermal energy, and pressure fields
  // from the state.
  const auto buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& massDensity = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
  const auto& energy = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero);
    
  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(mu.nodeListPtr());
  CHECK(solidNodeListPtr != 0);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Now set the shear modulus.
  strengthModel.shearModulus(mu, massDensity, energy, P, D);

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

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ShearModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  const auto* rhsPtr = dynamic_cast<const ShearModulusPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

