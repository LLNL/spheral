//---------------------------------Spheral++----------------------------------//
// NonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method in the case where we do *not* assume that pairwise forces are
// equal and opposite.
//
// Created by JMO, Sat Aug 10 23:03:39 PDT 2013
//----------------------------------------------------------------------------//
#include "NonSymmetricSpecificThermalEnergyPolicy.hh"
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
using std::numeric_limits;
using std::abs;
using std::min;
using std::max;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NonSymmetricSpecificThermalEnergyPolicy<Dimension>::
NonSymmetricSpecificThermalEnergyPolicy(const DataBase<Dimension>& dataBase):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(),
  mDataBasePtr(&dataBase) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
NonSymmetricSpecificThermalEnergyPolicy<Dimension>::
~NonSymmetricSpecificThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
NonSymmetricSpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

//   // HACK!
//   std::cerr.setf(std::ios::scientific, std::ios::floatfield);
//   std::cerr.precision(15);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto eps = state.fields(fieldKey, Scalar());
  const auto numFields = eps.numFields();

  // Get the state fields.
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto  eps0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const auto& pairAccelerations = derivs.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  const auto  DepsDt0 = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  const auto  nint = mDataBasePtr->numInternalNodes();
  CHECK(pairAccelerations.size() == 2*npairs or
        pairAccelerations.size() == 2*npairs + nint);
  const bool selfInteraction = (pairAccelerations.size() == 2*npairs + nint);

  // // Check if there is a surface point flag field registered.  If so, we use non-compatible energy evolution 
  // // on such points.
  // const bool surface = state.fieldNameRegistered(HydroFieldNames::surfacePoint);
  // FieldList<Dimension, int> surfacePoint;
  // if (surface) {
  //   surfacePoint = state.fields(HydroFieldNames::surfacePoint, int(0));
  // }

  const auto hdt = 0.5*multiplier;
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");
  // auto poisoned = mDataBasePtr->newFluidFieldList(0, "poisoned flag");

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

      // State for node i.
      const auto  weighti = 0.5; // abs(DepsDt0(nodeListi, i)) + numeric_limits<Scalar>::epsilon();
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      const auto& pacci = pairAccelerations[2*kk];

      // State for node j.
      const auto  weightj = 0.5; // abs(DepsDt0(nodeListj, j)) + numeric_limits<Scalar>::epsilon();
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& aj = acceleration(nodeListj, j);
      const auto  vj12 = vj + aj*hdt;
      const auto& paccj = pairAccelerations[2*kk+1];

      const auto dEij = -(mi*vi12.dot(pacci) + mj*vj12.dot(paccj));
      const auto wi = weighti/(weighti + weightj);
      CHECK(wi >= 0.0 and wi <= 1.0);

      DepsDt_thread(nodeListi, i) += wi*dEij/mi;
      DepsDt_thread(nodeListj, j) += (1.0 - wi)*dEij/mj;

      // // Check if either of these points was advanced non-conservatively.
      // if (surface) {
      //   poisoned(nodeListi, i) |= (surfacePoint(nodeListi, i) > 1 or surfacePoint(nodeListj, j) > 1 ? 1 : 0);
      //   poisoned(nodeListj, j) |= (surfacePoint(nodeListi, i) > 1 or surfacePoint(nodeListj, j) > 1 ? 1 : 0);
      // }
    }

#pragma omp critical
    {
      DepsDt_thread.threadReduce();
    }
  }

  // Now we can update the energy.
  auto offset = 2*npairs;
  for (auto nodeListi = 0u; nodeListi < numFields; ++nodeListi) {
    const auto n = eps[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {

      // Add the self-contribution if any (RZ with strength does this for instance).
      if (selfInteraction) {
        const auto& vi = velocity(nodeListi, i);
        const auto& ai = acceleration(nodeListi, i);
        const auto  vi12 = vi + ai*hdt;
        const auto duii = -2.0*vi12.dot(pairAccelerations[offset + i]);
        DepsDt(nodeListi, i) += duii;
      }

      eps(nodeListi, i) += DepsDt(nodeListi, i)*multiplier;
    }
    offset += n;
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
NonSymmetricSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const NonSymmetricSpecificThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const NonSymmetricSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

