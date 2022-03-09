//---------------------------------Spheral++----------------------------------//
// CompatibleDifferenceSpecificThermalEnergyPolicy -- An implementation of UpdatePolicyBase 
// specialized for the updating the specific thermal energy as a dependent 
// quantity.
// 
// This version is specialized for materials with different properties. A 
// compatible energy discretization in which pairwise work allows for opposite
// sign pair-wise work. DepsDti and  DepsDtj are used as weights and the 
// difference between the conservative  and consistent formulations is added 
// back in.
//----------------------------------------------------------------------------//
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/FieldUpdatePolicyBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/SpheralFunctions.hh"

#include <vector>
#include <limits>
using std::vector;
using std::numeric_limits;
using std::abs;
using std::min;
using std::max;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>::
CompatibleDifferenceSpecificThermalEnergyPolicy(const DataBase<Dimension>& dataBase):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(),
  mDataBasePtr(&dataBase) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>::
~CompatibleDifferenceSpecificThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();

  // constant we'll need for the weighting scheme
  const auto tiny = numeric_limits<Scalar>::epsilon();

  // Get the state fields.
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto& pairAccelerations = derivs.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  const auto& pairDepsDt = derivs.getAny(HydroFieldNames::pairWork, vector<Scalar>());
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  CHECK(pairAccelerations.size() == npairs);
  CHECK(pairDepsDt.size() == 2*npairs);

  auto  DepsDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  DepsDt.Zero();

  const auto hdt = 0.5*multiplier;
  
  // Walk all pairs and figure out the discrete work for each point
#pragma omp parallel
  {
    // Thread private variables
    auto DepsDt_thread = DepsDt.threadCopy();

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto i = pairs[kk].i_node;
      const auto j = pairs[kk].j_node;
      const auto nodeListi = pairs[kk].i_list;
      const auto nodeListj = pairs[kk].j_list;

      const auto& paccij = pairAccelerations[kk];
      const auto& DepsDt0i = pairDepsDt[2*kk];
      const auto& DepsDt0j = pairDepsDt[2*kk+1];

      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);

      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& aj = acceleration(nodeListj, j);

      // half-step velocity
      const auto vi12 = vi + ai*hdt;
      const auto vj12 = vj + aj*hdt;
      const auto vij = vi12 - vj12;

      // weighting scheme
      const auto weighti = abs(DepsDt0i) + tiny;
      const auto weightj = abs(DepsDt0j) + tiny;
      const auto wi = weighti/(weighti+weightj);

      // difference between assessed derivs and conserative ones
      const Scalar delta_duij = vij.dot(paccij)-DepsDt0i-DepsDt0j;

      CHECK(wi >= 0.0 and wi <= 1.0);

      // make conservative
      DepsDt_thread(nodeListi, i) += mj*(wi*delta_duij+DepsDt0i);
      DepsDt_thread(nodeListj, j) += mi*((1.0-wi)*delta_duij+DepsDt0j);

    }

#pragma omp critical
    {
      DepsDt_thread.threadReduce();
    }
  }

  // Now we can update the energy.
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto n = eps[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      eps(nodeListi, i) += DepsDt(nodeListi, i)*multiplier;
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

