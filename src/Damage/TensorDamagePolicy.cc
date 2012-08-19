//---------------------------------Spheral++----------------------------------//
// TensorDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent scalar damage state.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//   Randles, P.W. & Libersky, L.D., 1996, "Comput. Methods Appl. Engrg, 
//     139, 375-408
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#include <algorithm>

#include "TensorDamagePolicy.hh"
#include "TensorDamageModel.hh"
#include "oneMinusDamage.hh"
#include "Strength/SolidNodeList.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/GeometricUtilities.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;

using FieldSpace::Field;
using SolidMaterial::SolidNodeList;
using PhysicsSpace::TensorDamageModel;
using Material::EquationOfState;

//------------------------------------------------------------------------------
// Sort the eigen values (& associated eigen vectors) in decreasing order.
//------------------------------------------------------------------------------
inline
void
sortEigen(Dim<1>::SymTensor::EigenStructType& eigeni) {
}

inline
void
sortEigen(Dim<2>::SymTensor::EigenStructType& eigeni) {
  if (eigeni.eigenValues(0) < eigeni.eigenValues(1)) {
    swap(eigeni.eigenValues(0), eigeni.eigenValues(1));
    eigeni.eigenVectors = Dim<2>::Tensor(eigeni.eigenVectors.xy(), eigeni.eigenVectors.xx(),
                                         eigeni.eigenVectors.yy(), eigeni.eigenVectors.yx());
  }
  ENSURE(eigeni.eigenValues(0) >= eigeni.eigenValues(1));
}

inline
void
sortEigen(Dim<3>::SymTensor::EigenStructType& eigeni) {
  int i = 0;
  int j = 1;
  int k = 2;
  bool flipped = true;
  while (flipped) {
    flipped = false;
    if (eigeni.eigenValues(i) < eigeni.eigenValues(j)) {
      flipped = true;
      swap(i, j);
    }
    if (eigeni.eigenValues(j) < eigeni.eigenValues(k)) {
      flipped = true;
      swap(j, k);
    }
  }
  eigeni.eigenValues = Dim<3>::Vector(eigeni.eigenValues(i),
                                      eigeni.eigenValues(j),
                                      eigeni.eigenValues(k));
  eigeni.eigenVectors = Dim<3>::Tensor(eigeni.eigenVectors(0,i), eigeni.eigenVectors(0,j), eigeni.eigenVectors(0,k),
                                       eigeni.eigenVectors(1,i), eigeni.eigenVectors(1,j), eigeni.eigenVectors(1,k),
                                       eigeni.eigenVectors(2,i), eigeni.eigenVectors(2,j), eigeni.eigenVectors(2,k));
  ENSURE((eigeni.eigenValues(0) >= eigeni.eigenValues(1)) &&
         (eigeni.eigenValues(1) >= eigeni.eigenValues(2)));
}

//------------------------------------------------------------------------------
// Determine the effective rotational transformation from the given velocity
// gradient.
//------------------------------------------------------------------------------
inline
Dim<1>::Tensor
effectiveRotation(const Dim<1>::Tensor& DvDx) {
  return Dim<1>::Tensor::one;
}

inline
Dim<2>::Tensor
effectiveRotation(const Dim<2>::Tensor& DvDx) {
  const double theta = DvDx.xy() - DvDx.yx();
  return Dim<2>::Tensor(cos(theta), sin(theta),
                        -sin(theta), cos(theta));
}

inline
Dim<3>::Tensor
effectiveRotation(const Dim<3>::Tensor& DvDx) {
  const Dim<3>::Vector theta(DvDx.yz() - DvDx.zy(),
                             DvDx.zx() - DvDx.xz(),
                             DvDx.xy() - DvDx.yx());
  return (Dim<3>::Tensor(cos(theta.x()), sin(theta.x()), 0.0,
                         -sin(theta.x()), cos(theta.x()), 0.0,
                         0.0, 0.0, 1.0)*
          Dim<3>::Tensor(cos(theta.y()), sin(theta.y()), 0.0,
                         0.0, 1.0, 0.0,
                         -sin(theta.y()), 0.0, cos(theta.y()))*
          Dim<3>::Tensor(1.0, 0.0, 0.0,
                         0.0, cos(theta.z()), sin(theta.z()),
                         0.0, -sin(theta.z()), cos(theta.z())));
}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorDamagePolicy<Dimension>::
TensorDamagePolicy(const TensorDamageModel<Dimension>& damageModel):
  UpdatePolicyBase<Dimension>(SolidFieldNames::strain,
                              SolidFieldNames::effectiveFlaws),
  mDamageModelPtr(&damageModel) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
TensorDamagePolicy<Dimension>::
~TensorDamagePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TensorDamagePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::tensorDamage);
  Field<Dimension, SymTensor>& stateField = state.field(key, SymTensor::zero);

  const double tiny = 1.0e-30;
  const double tol = 1.0e-5;

  // Get the state fields.
  const KeyType strainKey = State<Dimension>::buildFieldKey(SolidFieldNames::effectiveStrainTensor, nodeListKey);
  const KeyType psKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrain, nodeListKey);
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType clKey = State<Dimension>::buildFieldKey(SolidFieldNames::longitudinalSoundSpeed, nodeListKey);
  const KeyType HKey = State<Dimension>::buildFieldKey(HydroFieldNames::H, nodeListKey);
  const KeyType DdamageDtKey = State<Dimension>::buildFieldKey(this->prefix() + SolidFieldNames::scalarDamage, nodeListKey);
  const KeyType DvDxKey = State<Dimension>::buildFieldKey(HydroFieldNames::internalVelocityGradient, nodeListKey);
  CHECK(state.registered(strainKey));
  CHECK(state.registered(psKey));
  CHECK(state.registered(PKey));
  CHECK(state.registered(clKey));
  CHECK(state.registered(HKey));
  CHECK(derivs.registered(DdamageDtKey));
  CHECK(derivs.registered(DvDxKey));

  const Field<Dimension, SymTensor>& strain = state.field(strainKey, SymTensor::zero);
  const Field<Dimension, Scalar>& plasticStrain = state.field(psKey, 0.0);
  const Field<Dimension, Scalar>& pressure = state.field(PKey, 0.0);
  const Field<Dimension, Scalar>& cl = state.field(clKey, 0.0);
  const Field<Dimension, SymTensor>& H = state.field(HKey, SymTensor::zero);
  const Field<Dimension, Scalar>& DDDt = derivs.field(DdamageDtKey, 0.0);
  const Field<Dimension, Tensor>& localDvDx = derivs.field(DvDxKey, Tensor::zero);

  // Iterate over the internal nodes.
  for (int i = 0; i != stateField.numInternalElements(); ++i) {

    SymTensor& Di = stateField(i);
    CHECK(Di.eigenValues().minElement() >= 0.0 and
          fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    
    // First apply the rotational term to the current damage.
    const Tensor spin = localDvDx(i).SkewSymmetric();
    const SymTensor spinCorrection = (spin*Di - Di*spin).Symmetric();
    Di += multiplier*spinCorrection;
    Di = max(1.0e-5, min(1.0 - 2.0e-5, Di));
    {
      const double maxValue = Di.eigenValues().maxElement();
      if (maxValue > 1.0) Di /= maxValue;
    }
    CHECK(Di.eigenValues().minElement() >= 0.0 and
          fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));

    // The flaw population for this node.
    const vector<double> flaws = mDamageModelPtr->flawsForNode(i);
    const int totalCracks = flaws.size();
    if (totalCracks > 0) {

      // The tensor strain on this node.
      typename SymTensor::EigenStructType eigeni = strain(i).eigenVectors();

      // Boost the effective strain by the damage.
//       const Scalar fDi = max(0.0, 1.0 - Di.eigenValues().maxElement());
//       CHECK(fDi >= 0.0 and fDi <= 1.0);
//       eigeni.eigenValues *= fDi/(fDi*fDi + 1.0e-10);
    
      // We want to go over these from max to min.
      sortEigen(eigeni);

      // Iterate over the eigen values/vectors of the strain.
      for (int j = 0; j != Dimension::nDim; ++j) {

        // Only increment the damage in this direction if there is a positive strain.
        const Scalar straini = (eigeni.eigenValues(j)); //   + plasticStrain(i))/(fDi*fDi + 1.0e-20);
        if (straini > 0.0) {

          //         // The total damage on this node thus far.
          //         CHECK(fuzzyGreaterThanOrEqual(Di.Trace(), 0.0) &&
          //               fuzzyLessThanOrEqual(Di.Trace(), 1.0));
          //         const Scalar D0max = max(0.0, min(1.0, Di.Trace()));
          //         CHECK(D0max >= 0.0 && D0max <= 1.0);

          // The direction of this strain, and projected damage.
          const Vector strainDirection = eigeni.eigenVectors.getColumn(j);
          CHECK(fuzzyEqual(strainDirection.magnitude2(), 1.0, 1.0e-10));
          const double D0 = max(0.0, min(1.0, (Di * strainDirection).magnitude()));
          CHECK(D0 >= 0.0 && D0 <= 1.0);

          // Count how many flaws are remaining, and how many have completely failed.
          const int numFailedCracks = int(D0*double(totalCracks));
          const int numRemainingCracks = totalCracks - numFailedCracks;
          CHECK(numFailedCracks >=0 && numFailedCracks <= totalCracks);
          CHECK(numRemainingCracks >= 0 && numRemainingCracks <= totalCracks);
          CHECK(numRemainingCracks + numFailedCracks == totalCracks);

          // Count how many cracks are currently active.
          int numActiveCracks = 0;
          //         if (weightedNeighborSum(i) < wSumCutoff) {
          //           numActiveCracks = numRemainingCracks;
          //         } else {
          for (int k = numFailedCracks; k < totalCracks; ++k) {
            CHECK(k < flaws.size());
            if (straini >= flaws[k]) ++numActiveCracks;
          }
          //         }
          CHECK(numActiveCracks >= 0 && numActiveCracks <= totalCracks);

          // Choose the allowed range of D.
          double Dmin, Dmax;
          if (multiplier >= 0.0) {
            Dmin = D0;
            Dmax = max(D0, double(numFailedCracks + numActiveCracks)/double(totalCracks));
          } else {
            Dmin = 0.0;
            Dmax = D0;
          }
          CHECK(Dmax >= Dmin);
          CHECK(Dmin >= 0.0 && Dmax <= 1.0);

          // Increment the damage.
          const double D013 = FastMath::CubeRootHalley2(D0);
          const double D113 = D013 + multiplier*double(numActiveCracks + numFailedCracks)/double(totalCracks)*DDDt(i);
          const double D1 = max(Dmin, min(Dmax, D113*D113*D113));
          CHECK((D1 - D0)*multiplier >= 0.0);

          // Increment the damage tensor.
          Di += (D1 - D0)*strainDirection.selfdyad();

        }
      }
    }

    // Enforce bounds on the damage.
    const Vector Dvals = Di.eigenValues();
    if (Dvals.minElement() < 0.0 or
        Dvals.maxElement() > 1.0) {
      const typename SymTensor::EigenStructType Deigen = Di.eigenVectors();
      Di = constructSymTensorWithBoundedDiagonal(Deigen.eigenValues,
                                                 1.0e-5,
                                                 1.0);
      Di.rotationalTransform(Deigen.eigenVectors);
    }
    ENSURE(Di.eigenValues().minElement() >= 0.0 and
           fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    //     ENSURE(fuzzyGreaterThanOrEqual(Di.eigenValues().minElement(), 0.0, 1.0e-5) &&
    //            fuzzyLessThanOrEqual(Di.eigenValues().maxElement(), 1.0, 1.0e-5));
    //     ENSURE(fuzzyGreaterThanOrEqual(Di.eigenValues().minElement(), 0.0) &&
    //            fuzzyLessThanOrEqual(Di.Trace(), 1.0));

  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
TensorDamagePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const TensorDamagePolicy<Dimension>* rhsPtr = dynamic_cast<const TensorDamagePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

