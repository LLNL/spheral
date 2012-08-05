//---------------------------------Spheral++----------------------------------//
// TensorStrainPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the tensor strain.
//
// Created by JMO, Mon Oct 17 10:56:28 PDT 2005
//----------------------------------------------------------------------------//
#include <vector>

#include "TensorStrainPolicy.hh"
#include "NodeList/NodeList.hh"
#include "Strength/SolidNodeList.hh"
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
using SolidMaterial::SolidNodeList;
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
  const KeyType eKey = State<Dimension>::buildFieldKey(SolidFieldNames::strainTensor, nodeListKey);
  const KeyType EKey = State<Dimension>::buildFieldKey(SolidFieldNames::YoungsModulus, nodeListKey);
  const KeyType KKey = State<Dimension>::buildFieldKey(SolidFieldNames::bulkModulus, nodeListKey);
  const KeyType muKey = State<Dimension>::buildFieldKey(SolidFieldNames::shearModulus, nodeListKey);
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType psKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrain, nodeListKey);
  const KeyType stressKey = State<Dimension>::buildFieldKey(SolidFieldNames::deviatoricStress, nodeListKey);
  const KeyType gradvKey = State<Dimension>::buildFieldKey(HydroFieldNames::internalVelocityGradient, nodeListKey);
  const KeyType dSKey = State<Dimension>::buildFieldKey(IncrementState<Dimension, Field<Dimension, SymTensor> >::prefix() + SolidFieldNames::deviatoricStress, nodeListKey);
  CHECK(state.registered(eKey));
  CHECK(state.registered(EKey));
  CHECK(state.registered(KKey));
  CHECK(state.registered(muKey));
  CHECK(state.registered(PKey));
  CHECK(state.registered(psKey));
  CHECK(state.registered(stressKey));
  CHECK(derivs.registered(gradvKey));
  CHECK(derivs.registered(dSKey));

  Field<Dimension, SymTensor>& strain = state.field(eKey, SymTensor::zero);
  const Field<Dimension, Scalar>& E = state.field(EKey, 0.0);
  const Field<Dimension, Scalar>& K = state.field(KKey, 0.0);
  const Field<Dimension, Scalar>& mu = state.field(muKey, 0.0);
  const Field<Dimension, Scalar>& P = state.field(PKey, 0.0);
  const Field<Dimension, Scalar>& plasticStrain = state.field(psKey, 0.0);
  const Field<Dimension, SymTensor>& S = state.field(stressKey, SymTensor::zero);
  const Field<Dimension, Tensor>& gradv = derivs.field(gradvKey, Tensor::zero);
  const Field<Dimension, SymTensor>& DSDt = derivs.field(dSKey, SymTensor::zero);

  // Iterate over the internal nodes.
  for (int i = 0; i != stateField.numInternalElements(); ++i) {

    // Begin the big bonanza of options!

    // PseudoPlasticStrain.
    if (mStrainType == PhysicsSpace::PseudoPlasticStrain) {

      strain(i) += multiplier*safeInv(mu(i), 1.0e-10)*DSDt(i);
      stateField(i) = strain(i);

    } else {

      // First apply the rotational term to the current strain history.
      const Tensor spin = gradv(i).SkewSymmetric();
      strain(i) += multiplier*(spin*strain(i) + strain(i)*spin).Symmetric();

      // Update the strain history with the current instantaneous deformation.
      const typename SymTensor::EigenStructType eigenv = gradv(i).Symmetric().eigenVectors();
      SymTensor gradvi = constructSymTensorWithMaxDiagonal(eigenv.eigenValues, 0.0);
      gradvi.rotationalTransform(eigenv.eigenVectors);
      strain(i) += multiplier*gradvi;

      const double volstrain = strain(i).Trace();

      // Update the effective strain according to the specified algorithm.
      switch(mStrainType) {

      case(PhysicsSpace::BenzAsphaug):
        CHECK(E(i) >= 0.0);
        stateField(i) = (S(i) - P(i)*SymTensor::one)/(E(i) + tiny); // thpt);
        break;

      case(PhysicsSpace::StrainHistory):
        stateField(i) = strain(i);
        break;

      case(PhysicsSpace::MeloshRyanAsphaug):
        stateField(i) = ((K(i) - 2.0*mu(i)/Dimension::nDim)*volstrain*SymTensor::one + 2.0*mu(i)*strain(i))/(E(i) + tiny); // thpt);
        break;

      case(PhysicsSpace::PlasticStrain):
        stateField(i) = plasticStrain(i)*SymTensor::one;
        break;

      default:
        VERIFY2(false, "TensorStrainPolicy ERROR:  no update for case " << mStrainType << "!");
        break;

      }
    }

    // Apply limiting to the effective strain.
    stateField(i) = max(1.0e-7*max(1.0, std::abs(stateField(i).Trace())/Dimension::nDim), stateField(i));
    ENSURE2(fuzzyGreaterThanOrEqual(stateField(i).eigenValues().minElement(), 0.0, 1.0e-7),
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

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class Spheral::TensorStrainPolicy<Dim<1> >;
  template class Spheral::TensorStrainPolicy<Dim<2> >;
  template class Spheral::TensorStrainPolicy<Dim<3> >;
}
