//---------------------------------Spheral++----------------------------------//
// TensorStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the tensor strain.
//
// Created by JMO, Mon Oct 17 10:56:28 PDT 2005
//----------------------------------------------------------------------------//
#include <vector>

#include "TensorStrainPolicy.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

using namespace std;

using FieldSpace::Field;
using NodeSpace::NodeList;
using NodeSpace::SolidNodeList;
using KernelSpace::TableKernel;
using PhysicsSpace::TensorStrainAlgorithm;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorStrainPolicy<Dimension>::
TensorStrainPolicy(const TensorStrainAlgorithm strainType):
  UpdatePolicyBase<Dimension>(HydroFieldNames::position,
                              HydroFieldNames::H,
                              SolidFieldNames::YoungsModulus,
                              SolidFieldNames::bulkModulus,
                              SolidFieldNames::shearModulus,
                              HydroFieldNames::pressure,
                              SolidFieldNames::deviatoricStress),
  mStrainType(strainType) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorStrainPolicy<Dimension>::
~TensorStrainPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorStrainPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::effectiveStrainTensor);
  Field<Dimension, SymTensor>& stateField = state.field(key, SymTensor::zero);

  const double tiny = 1.0e-15;

  // Get the state fields.
  const auto eKey = State<Dimension>::buildFieldKey(SolidFieldNames::strainTensor, nodeListKey);
  const auto EKey = State<Dimension>::buildFieldKey(SolidFieldNames::YoungsModulus, nodeListKey);
  const auto KKey = State<Dimension>::buildFieldKey(SolidFieldNames::bulkModulus, nodeListKey);
  const auto muKey = State<Dimension>::buildFieldKey(SolidFieldNames::shearModulus, nodeListKey);
  const auto PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const auto psKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrain, nodeListKey);
  const auto stressKey = State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, nodeListKey);
  const auto gradvKey = State<Dimension>::buildFieldKey(HydroFieldNames::internalVelocityGradient, nodeListKey);
  const auto dSKey = State<Dimension>::buildFieldKey(IncrementState<Dimension, Field<Dimension, SymTensor> >::prefix() + SolidFieldNames::deviatoricStress, nodeListKey);
  CHECK(state.registered(eKey));
  CHECK(state.registered(EKey));
  CHECK(state.registered(KKey));
  CHECK(state.registered(muKey));
  CHECK(state.registered(PKey));
  CHECK(state.registered(psKey));
  CHECK(state.registered(stressKey));
  CHECK(derivs.registered(gradvKey));
  CHECK(derivs.registered(dSKey));

  auto&       strain = state.field(eKey, SymTensor::zero);
  const auto& E = state.field(EKey, 0.0);
  const auto& K = state.field(KKey, 0.0);
  const auto& mu = state.field(muKey, 0.0);
  const auto& P = state.field(PKey, 0.0);
  const auto& plasticStrain = state.field(psKey, 0.0);
  const auto& S = state.field(stressKey, SymTensor::zero);
  const auto& gradv = derivs.field(gradvKey, Tensor::zero);
  const auto& DSDt = derivs.field(dSKey, SymTensor::zero);

  // Iterate over the internal nodes.
  const auto ni = stateField.numInternalElements();
#pragma omp parallel for
  for (int i = 0; i < ni; ++i) {

    // Begin the big bonanza of options!

    // PseudoPlasticStrain.
    if (mStrainType == PhysicsSpace::TensorStrainAlgorithm::PseudoPlasticStrain) {

      strain(i) += multiplier*safeInv(mu(i), 1.0e-10)*DSDt(i);
      stateField(i) = strain(i);

    } else {

      // First apply the rotational term to the current strain history.
      const auto spin = gradv(i).SkewSymmetric();
      strain(i) += multiplier*(spin*strain(i) + strain(i)*spin).Symmetric();

      // Update the strain history with the current instantaneous deformation.
      const auto eigenv = gradv(i).Symmetric().eigenVectors();
      auto       gradvi = constructSymTensorWithMaxDiagonal(eigenv.eigenValues, 0.0);
      gradvi.rotationalTransform(eigenv.eigenVectors);
      strain(i) += multiplier*gradvi;

      const auto volstrain = strain(i).Trace();

      // Update the effective strain according to the specified algorithm.
      switch(mStrainType) {

      case(PhysicsSpace::TensorStrainAlgorithm::BenzAsphaugStrain):
        CHECK(E(i) >= 0.0);
        stateField(i) = (S(i) - P(i)*SymTensor::one)/(E(i) + tiny); // thpt);
        break;

      case(PhysicsSpace::TensorStrainAlgorithm::StrainHistory):
        stateField(i) = strain(i);
        break;

      case(PhysicsSpace::TensorStrainAlgorithm::MeloshRyanAsphaugStrain):
        stateField(i) = ((K(i) - 2.0*mu(i)/Dimension::nDim)*volstrain*SymTensor::one + 2.0*mu(i)*strain(i))/(E(i) + tiny); // thpt);
        break;

      case(PhysicsSpace::TensorStrainAlgorithm::PlasticStrain):
        stateField(i) = plasticStrain(i)*SymTensor::one;
        break;

      default:
        VERIFY2(false, "TensorStrainPolicy ERROR:  no update for case " << static_cast<int>(mStrainType) << "!");
        break;

      }
    }

    // Apply limiting to the effective strain.
    stateField(i) = max(1.0e-7*max(1.0, std::abs(stateField(i).Trace())/Dimension::nDim), stateField(i));
    ENSURE2(fuzzyGreaterThanOrEqual(stateField(i).eigenValues().minElement(), 0.0, 1.0e-5),
            "Effective strain bad eigenvalues!  " << stateField(i).eigenValues());

  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TensorStrainPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const TensorStrainPolicy<Dimension>* rhsPtr = dynamic_cast<const TensorStrainPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

