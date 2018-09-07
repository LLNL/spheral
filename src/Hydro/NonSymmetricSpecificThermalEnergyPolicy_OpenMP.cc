//---------------------------------Spheral++----------------------------------//
// NonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method in the case where we do *not* assume that pairwise forces are
// equal and opposite.
//
// Also specialized for OpenMP.
//
//  Created by JMO, Sat Aug 10 23:03:39 PDT 2013
// Modified by JMO, Wed Feb 28 15:46:04 PST 2018
//----------------------------------------------------------------------------//
#include "NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/FieldListUpdatePolicyBase.hh"
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

// namespace {

// inline double weighting(const double wi,
//                         const double wj,
//                         const double Eij) {
//   const int si = isgn0(wi);
//   const int sj = isgn0(wj);
//   const int sij = isgn0(wij);
//   if (sij == 0 or (si == 0 and sj == 0))             return 0.5;
//   if (si == sj and sj == sij)                        return wi/(wi + wj);  // All the same sign and non-zero
//   if (si == sj)                                      return 0.5;
//   if (si == sij)                                     return 1.0;
//   return 0.0;
// }

// }

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
  const auto numFields = eps.numFields();

  // Get the state fields.
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto  eps0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const auto  pairAccelerations = derivs.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  const auto  DepsDt0 = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  CHECK(nodeLists.size() == numFields);

  // Check if there is a surface point flag field registered.  If so, we use non-compatible energy evolution 
  // on such points.
  const bool surface = state.fieldNameRegistered(HydroFieldNames::surfacePoint);
  FieldList<Dimension, int> surfacePoint;
  if (surface) {
    surfacePoint = state.fields(HydroFieldNames::surfacePoint, int(0));
  }

  // Prepare a counter to keep track of how we go through the pair-accelerations.
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");
  // auto poisoned = mDataBasePtr->newFluidFieldList(0, "poisoned flag");

  // Walk all the NodeLists and compute the energy change.
  const auto hdt = 0.5*multiplier;
  for (size_t nodeListi = 0; nodeListi != numFields; ++nodeListi) {
    const auto ni = connectivityMap.numNodes(nodeListi);

    // Iterate over the internal nodes of this NodeList.
// #pragma omp parallel for reduction(+:DepsDt)
    for (auto k = 0; k < ni; ++k) {
      const auto i = connectivityMap.ithNode(nodeListi, k);

      // State for node i.
      auto&       DepsDti = DepsDt(nodeListi, i);
      const auto  weighti = abs(DepsDt0(nodeListi, i)) + std::numeric_limits<Scalar>::epsilon();
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  ui = eps0(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      const auto& pacci = pairAccelerations(nodeListi, i);
      CHECK(pacci.size() == 2*connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i) or
            pacci.size() == 2*connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i) + 1);
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
            const auto        j = *jitr;
            auto&       DepsDtj = DepsDt(nodeListj, j);
            const auto  weightj = abs(DepsDt0(nodeListj, j)) + std::numeric_limits<Scalar>::epsilon();
            const auto  mj = mass(nodeListj, j);
            const auto& vj = velocity(nodeListj, j);
            const auto& aj = acceleration(nodeListj, j);
            const auto  vj12 = vj + aj*hdt;

            CHECK(offseti < pacci.size());
            const auto& pai = pacci[offseti++];
            const auto& paj = pacci[offseti++];
            const auto dEij = -(mi*vi12.dot(pai) + mj*vj12.dot(paj));
            const auto wi = weighti/(weighti + weightj);
            // const auto wi = entropyWeighting(si, sj, duij);
            // CHECK2(fuzzyEqual(wi + entropyWeighting(sj, si, dEij/mj), 1.0, 1.0e-10),
            //        wi << " " << entropyWeighting(sj, si, dEij/mj) << " " << (wi + entropyWeighting(sj, si, dEij/mj)));
            CHECK(wi >= 0.0 and wi <= 1.0);
            DepsDti += wi*dEij/mi;

            // // Check if either of these points was advanced non-conservatively.
            // if (surface) {
            //   poisoned(nodeListi, i) |= (surfacePoint(nodeListi, i) > 1 or surfacePoint(nodeListj, j) > 1 ? 1 : 0);
            // }
          }
        }
      }
      CHECK(offseti == pacci.size() or offseti == pacci.size() - 1);

      // Add the self-contribution if any (RZ does this for instance).
      if (offseti == pacci.size() - 1) {
        const auto duii = -2.0*vi12.dot(pacci.back());
        DepsDti += duii;
      }

      // if (poisoned(nodeListi, i) == 0) {
        eps(nodeListi, i) += DepsDti*multiplier;
      // } else {
      //   eps(nodeListi, i) += DepsDt0(nodeListi, i)*multiplier;
      // }
    }
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

