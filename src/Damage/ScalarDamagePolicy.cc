//---------------------------------Spheral++----------------------------------//
// ScalarDamagePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent scalar damage state.
//
// Created by JMO, Sun Sep 26 16:15:17 PDT 2004
//----------------------------------------------------------------------------//
#include "ScalarDamagePolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidFieldNames.hh"
#include "DamageModel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Kernel/TableKernel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ScalarDamagePolicy<Dimension>::
ScalarDamagePolicy(const DamageModel<Dimension>& damageModel):
  UpdatePolicyBase<Dimension, FieldType>(SolidFieldNames::strain),
  mDamageModelPtr(&damageModel) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ScalarDamagePolicy<Dimension>::
~ScalarDamagePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ScalarDamagePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  REQUIRE(key.second == SolidFieldNames::scalarDamage);
  Field<Dimension, Scalar>& stateField = state.scalarField(key);

  const double tiny = 1.0e-30;

  // Get the effective weightedNeighborSum cutoff implied by the nodes per 
  // smoothing scale cutoff.
  const SolidNodeList<Dimension>* nodeListPtr = &(mDamageModelPtr->nodeList());
  const double wSumCutoff = nodeListPtr->kernel().equivalentWsum(mDamageModelPtr->criticalNodesPerSmoothingScale());

  // Get the state fields.
  const KeyType strainKey(nodeListPtr, SolidFieldNames::strain);
  const KeyType DdamageDtKey(nodeListPtr, this->prefix() + SolidFieldNames::scalarDamage);
//   const KeyType wsumKey(nodeListPtr, HydroFieldNames::weightedNeighborSum);

  Field<Dimension, Scalar>& strain = state.scalarField(strainKey);
  const Field<Dimension, Scalar>& DDDt = derivs.scalarField(DdamageDtKey);
//   const Field<Dimension, Scalar>& weightedNeighborSum = derivs.scalarField(wsumKey);

  // Iterate over the internal nodes.
  for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
    
    // Remember the initial damage.
    const double D0 = stateField(i);
    CHECK(D0 >= 0.0 && D0 <= 1.0);

    // The strain on this node.
    const Scalar straini = strain(i)/(1.0 - D0 + tiny*max(1.0, strain(i)));

    // The flaw population for this node.
    const vector<double>& flaws = mDamageModelPtr->flawsForNode(i);

    // Count how many flaws are remaining, and how many have completely failed.
    const int totalCracks = flaws.size();
    int numRemainingCracks = int((1.0 - D0)*double(totalCracks));
    if (D0 < 1.0) numRemainingCracks = max(1, numRemainingCracks);
    const int numFailedCracks = totalCracks - numRemainingCracks;
    CHECK(numRemainingCracks + numFailedCracks == totalCracks);

    // Count how many cracks are currently active.
    int numActiveCracks = 0;
//     if (weightedNeighborSum(i) < wSumCutoff) {
//       numActiveCracks = numRemainingCracks;
//     } else {
      for (int k = numFailedCracks; k < totalCracks; ++k) {
        CHECK(k < flaws.size());
        if (straini >= flaws[k]) ++numActiveCracks;
      }
//     }
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

    // Increment the damage.
    CHECK(DDDt(i) >= 0.0);
    const double D013 = FastMath::CubeRootHalley2(D0);
    const double D113 = D013 + multiplier*double(numActiveCracks)/double(totalCracks)*DDDt(i);
    stateField(i) = max(Dmin, min(Dmax, D113*D113*D113));
    CHECK((multiplier >= 0.0 && stateField(i) >= D0) || (multiplier < 0.0 && stateField(i) <= D0));
    CHECK(stateField(i) >= 0.0 && stateField(i) <= 1.0);

  }

}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ScalarDamagePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension, FieldType>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ScalarDamagePolicy<Dimension>* rhsPtr = dynamic_cast<const ScalarDamagePolicy<Dimension>*>(&rhs);
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
  template class ScalarDamagePolicy<Dim<1> >;
  template class ScalarDamagePolicy<Dim<2> >;
  template class ScalarDamagePolicy<Dim<3> >;
}
