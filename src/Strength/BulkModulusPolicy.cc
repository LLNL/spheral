//---------------------------------Spheral++----------------------------------//
// BulkModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent bulk modulus state.
//
// Created by JMO, Tue Oct 5 16:40:54 2004
//----------------------------------------------------------------------------//
#include "BulkModulusPolicy.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "NodeList/SolidNodeList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
BulkModulusPolicy<Dimension>::
BulkModulusPolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
                                                                   HydroFieldNames::specificThermalEnergy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
BulkModulusPolicy<Dimension>::
~BulkModulusPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
BulkModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::bulkModulus and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> stateFields = state.fields(fieldKey, Scalar());
  const unsigned numFields = stateFields.numFields();

  // Get the mass density and specific thermal energy fields from the state.
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> energy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  CHECK(massDensity.size() == numFields);
  CHECK(energy.size() == numFields);

  // Walk the individual fields.
  for (unsigned k = 0; k != numFields; ++k) {

    // Check if we have a FluidNodeList or SolidNodeList.  Has to be at least fluid!
    const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(stateFields[k]->nodeListPtr());
    const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(stateFields[k]->nodeListPtr());
    CHECK(fluidNodeListPtr != NULL);

    // Is there a strength model that wants to set the bulk modulus?
    if (solidNodeListPtr != NULL and
        solidNodeListPtr->strengthModel().providesBulkModulus()) {
      solidNodeListPtr->strengthModel().bulkModulus(*stateFields[k], *massDensity[k], *energy[k]);
    } else {
      fluidNodeListPtr->equationOfState().setBulkModulus(*stateFields[k], *massDensity[k], *energy[k]);
    }
  }

//   // Is there a scalar damage field for this NodeList?
//   {
//     const KeyType DKey(nodeListPtr, SolidFieldNames::scalarDamage);
//     if (state.scalarFieldRegistered(DKey)) {
//       const Field<Dimension, Scalar>& D = state.scalarField(DKey);
//       for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
//         CHECK(D(i) >= 0.0 && D(i) <= 1.0);
//         stateField(i) *= 1.0 - D(i);
//       }
//     }
//   }

//   // Is there a tensor damage field for this NodeList?
//   {
//     typedef typename Dimension::SymTensor SymTensor;
//     const KeyType DKey(nodeListPtr, SolidFieldNames::tensorDamage);
//     if (state.symTensorFieldRegistered(DKey)) {
//       const Field<Dimension, SymTensor>& D = state.symTensorField(DKey);
//       for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
//         const Scalar Di = std::max(0.0, std::min(1.0, D(i).eigenValues().maxElement() + 1.0e-5));
//         CHECK(Di >= 0.0 && Di <= 1.0);
//         stateField(i) *= 1.0 - Di;
//       }
//     }
//   }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
BulkModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a bulk modulus operator.
  const BulkModulusPolicy<Dimension>* rhsPtr = dynamic_cast<const BulkModulusPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

