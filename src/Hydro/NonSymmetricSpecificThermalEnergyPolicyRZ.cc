//---------------------------------Spheral++----------------------------------//
// NonSymmetricSpecificThermalEnergyPolicy -- An implementation of 
// UpdatePolicyBase specialized for the updating the specific thermal energy
// as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method in RZ coordinates, which implicitly means we have to use the
// non-symmetric (pairwise) form.
//
// Created by JMO, Wed May  4 16:49:59 PDT 2016
//----------------------------------------------------------------------------//
#include <vector>

#include "NonSymmetricSpecificThermalEnergyPolicyRZ.hh"
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

namespace Spheral {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
NonSymmetricSpecificThermalEnergyPolicyRZ::
NonSymmetricSpecificThermalEnergyPolicyRZ(const DataBase<Dim<2> >& dataBase):
  IncrementFieldList<Dim<2>, Dim<2>::Scalar>(),
  mDataBasePtr(&dataBase) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
NonSymmetricSpecificThermalEnergyPolicyRZ::
~NonSymmetricSpecificThermalEnergyPolicyRZ() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
NonSymmetricSpecificThermalEnergyPolicyRZ::
update(const KeyType& key,
       State<Dim<2> >& state,
       StateDerivatives<Dim<2> >& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  typedef Dimension::SymTensor SymTensor;

//   // HACK!
//   std::cerr.setf(std::ios::scientific, std::ios::floatfield);
//   std::cerr.precision(15);

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> eps = state.fields(fieldKey, Scalar());
  const unsigned numFields = eps.numFields();

  // Get the state fields.
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  mass = state.fields(HydroFieldNames::mass, Scalar());
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  acceleration = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto  eps0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const auto  pairAccelerations = derivs.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  const auto  DepsDt0 = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  const auto& connectivityMap = mDataBasePtr->connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  CHECK(nodeLists.size() == numFields);

  // Prepare a counter to keep track of how we go through the pair-accelerations.
  auto DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");
  auto offset = mDataBasePtr->newFluidFieldList(0, "offset");

  // Walk all the NodeLists.
  const auto hdt = 0.5*multiplier;
  for (size_t nodeListi = 0; nodeListi != numFields; ++nodeListi) {

    // Iterate over the internal nodes of this NodeList.
    for (auto iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // State for node i.
      auto&       DepsDti = DepsDt(nodeListi, i);
      const auto  weighti = abs(DepsDt0(nodeListi, i)) + numeric_limits<Scalar>::epsilon();
      const auto  ri = abs(position(nodeListi, i).y());
      const auto  zi = position(nodeListi, i).x();
      const auto  mi = mass(nodeListi, i);
      const auto  mRZi = mi/(2.0*M_PI*ri);
      const auto& vi = velocity(nodeListi, i);
      const auto  ui = eps0(nodeListi, i);
      const auto& ai = acceleration(nodeListi, i);
      const auto  vi12 = vi + ai*hdt;
      // vi12.x(vi12.x()/(2.0*M_PI*ri));
      const auto& pacci = pairAccelerations(nodeListi, i);
      CHECK(ri > 0.0);
      CHECK(pacci.size() == connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i) + 1);

      // Get the connectivity (neighbor set) for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      // Iterate over the neighbor NodeLists.
      for (size_t nodeListj = 0; nodeListj != numFields; ++nodeListj) {

        // The set of neighbors from this NodeList.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors, and accumulate the specific energy
          // change.
          for (auto jitr = connectivity.begin();
               jitr != connectivity.end();
               ++jitr) {
            const int j = *jitr;

            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              auto&       DepsDtj = DepsDt(nodeListj, j);
              const auto  weightj = abs(DepsDt0(nodeListj, j)) + numeric_limits<Scalar>::epsilon();
              const auto  rj = abs(position(nodeListj, j).y());
              const auto  zj = position(nodeListj, j).x();
              const auto  mj = mass(nodeListj, j);
              const auto  mRZj = mj/(2.0*M_PI*rj);
              const auto& vj = velocity(nodeListj, j);
              const auto  uj = eps0(nodeListj, j);
              const auto& aj = acceleration(nodeListj, j);
              const auto  vj12 = vj + aj*hdt;
              // vj12.x(vj12.x()/(2.0*M_PI*rj));
              const auto& paccj = pairAccelerations(nodeListj, j);
              CHECK(rj > 0.0);
              CHECK(j >= firstGhostNodej or paccj.size() == (connectivityMap.numNeighborsForNode(nodeLists[nodeListj], j) + 1));

              CHECK(offset(nodeListi, i) < pacci.size());
              const auto& pai = pacci[offset(nodeListi, i)];
              ++offset(nodeListi, i);

              CHECK(offset(nodeListj, j) < paccj.size());
              const auto& paj = paccj[offset(nodeListj, j)];
              ++offset(nodeListj, j);

              const auto dEij = -(mi*vi12.dot(pai) + mj*vj12.dot(paj));
              const auto duij = dEij/mi;
              const auto wi = weighti/(weighti + weightj);      // Du/Dt weighting

              CHECK(wi >= 0.0 and wi <= 1.0);
              DepsDti += wi*duij;
              DepsDtj += (1.0 - wi)*dEij/mj;
            }
          }
        }
      }

      // Add the self-contribution.
      const Vector& pai = pacci[offset(nodeListi, i)];
      const Scalar duij = -2.0*vi12.dot(pai);
      DepsDti += duij;

      // Now we can update the energy.
      CHECK(offset(nodeListi, i) == pacci.size() - 1);
      eps(nodeListi, i) += DepsDti*multiplier;
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
NonSymmetricSpecificThermalEnergyPolicyRZ::
operator==(const UpdatePolicyBase<Dim<2> >& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const NonSymmetricSpecificThermalEnergyPolicyRZ* rhsPtr = dynamic_cast<const NonSymmetricSpecificThermalEnergyPolicyRZ*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

