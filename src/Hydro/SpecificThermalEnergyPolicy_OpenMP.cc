//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the specific thermal energy as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method, *and* also for OpenMP.
//
//  Created by JMO, Tue Sep 14 22:27:08 2004
// Modified by JMO, Fri Nov  3 15:44:41 PDT 2017
//----------------------------------------------------------------------------//
#include "SpecificThermalEnergyPolicy.hh"
#include "HydroFieldNames.hh"
#include "entropyWeightingFunction.hh"
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
SpecificThermalEnergyPolicy<Dimension>::
SpecificThermalEnergyPolicy(const DataBase<Dimension>& dataBase):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(),
  mDataBasePtr(&dataBase) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SpecificThermalEnergyPolicy<Dimension>::
~SpecificThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  typedef typename Dimension::SymTensor SymTensor;

//   // HACK!
//   std::cerr.setf(std::ios::scientific, std::ios::floatfield);
//   std::cerr.precision(15);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const unsigned numFields = eps.numFields();

  // Get the state fields.
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  // const auto  entropy = state.fields(HydroFieldNames::entropy, Scalar());
  const auto  eps0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const auto  pairAccelerations = derivs.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  const auto  DepsDt0 = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  CHECK(nodeLists.size() == numFields);

  // Prepare a counter to keep track of how we go through the pair-accelerations.
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");

  // Walk all the NodeLists.
  const auto hdt = 0.5*multiplier;
  for (size_t nodeListi = 0; nodeListi != numFields; ++nodeListi) {
    const auto ni = connectivityMap.numNodes(nodeListi);

    // Iterate over the internal nodes of this NodeList.
#pragma omp parallel for
    for (auto k = 0; k < ni; ++k) {
      const auto i = connectivityMap.ithNode(nodeListi, k);

      // State for node i.
      auto& DepsDti = DepsDt(nodeListi, i);
      const auto  weighti = abs(DepsDt0(nodeListi, i)) + std::numeric_limits<Scalar>::epsilon();
      // const auto  si = entropy(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  ui = eps0(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      const auto& pacci = pairAccelerations(nodeListi, i);
      CHECK(pacci.size() == connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i) or
            NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());
      size_t offseti = 0;

      // Get the connectivity (neighbor set) for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      // Iterate over the neighbor NodeLists.
      for (size_t nodeListj = 0; nodeListj != numFields; ++nodeListj) {

        // The set of neighbors from this NodeList.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {

          // Iterate over the neighbors, and accumulate the specific energy
          // change.
          for (auto jitr = connectivity.begin();
               jitr != connectivity.end();
               ++jitr) {
            const int j = *jitr;

            const auto  weightj = abs(DepsDt0(nodeListj, j)) + std::numeric_limits<Scalar>::epsilon();
            // const auto  sj = entropy(nodeListj, j);
            const auto& vj = velocity(nodeListj, j);
            const auto& aj = acceleration(nodeListj, j);
            const auto  vj12 = vj + aj*hdt;
            const auto  vji12 = vj12 - vi12;
            const auto& pai = pacci[offseti++];
            const Scalar duij = vji12.dot(pai);
            const Scalar wi = weighti/(weighti + weightj);      // Du/Dt weighting
            // const Scalar wi = entropyWeighting(si, sj, duij);   // entropy weighting
            CHECK(wi >= 0.0 and wi <= 1.0);
            DepsDti += wi*duij;
          }
        }
      }

      // Now we can update the energy.
      CHECK2(offseti == pacci.size() or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent(),
             "Bad sizing : (" << nodeListi << " " << i << ") " << offseti << " " << pacci.size());
      eps(nodeListi, i) += DepsDti*multiplier;
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const SpecificThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const SpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

