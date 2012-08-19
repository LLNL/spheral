//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the specific thermal energy as a dependent quantity.
// 
// This version is specialized for the compatible energy discretization 
// method.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//

#include <vector>

#include "SpecificThermalEnergyPolicy.hh"
#include "HydroFieldNames.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/FieldUpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/safeInv.hh"
#include "Infrastructure/SpheralFunctions.hh"

namespace Spheral {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Define the weighting function deciding how to divvy up the work between
// nodes.
//------------------------------------------------------------------------------
// The old method.
// inline
// double standardWeighting(const double& ui,
//                          const double& uj,
//                          const double& mi,
//                          const double& mj,
//                          const double& duij) {
//   const double uji = uj - ui;
//   const double result = 0.5*(1.0 + uji/(abs(uji) + 1.0/(1.0 + abs(uji)))*sgn0(duij));
//   ENSURE(result >= 0.0 and result <= 1.0);
//   return result;
// }

inline
double standardWeighting(const double& ui,
                         const double& uj,
                         const double& mi,
                         const double& mj,
                         const double& duij) {
  const double uji = uj - ui;
  const double result = 0.5*(1.0 + uji*sgn0(duij)/(abs(uji) + 0.5));
  ENSURE(result >= 0.0 and result <= 1.0);
  return result;
}

// A strictly monotonic scheme that *will* pull the specific thermal energies
// together.
inline
double monotonicWeighting(const double& ui,
                          const double& uj,
                          const double& mi,
                          const double& mj,
                          const double& DEij) {
  REQUIRE(mi > 0.0);
  REQUIRE(mj > 0.0);
  const double uji = uj - ui;
  const double miInv = 1.0/mi;
  const double mjInv = 1.0/mj;

  double result;
  if (uji*DEij >= 0.0) {
    const double A = min(1.0, max(0.0, uji*mi*safeInv(DEij, 1.0e-50)));
    const double B = min(1.0, max(0.0, mjInv/(miInv + mjInv)));
    CHECK(A >= 0.0 and A <= 1.0);
    CHECK(B >= 0.0 and B <= 1.0);
    result = A + B*(1.0 - A);
  } else {
    const double A = min(1.0, max(0.0, -uji*mj*safeInv(DEij, 1.0e-50)));
    const double B = min(1.0, max(0.0, miInv/(miInv + mjInv)));
    CHECK(A >= 0.0 and A <= 1.0);
    CHECK(B >= 0.0 and B <= 1.0);
    result = 1.0 - A - B*(1.0 - A);
  }
  ENSURE(result >= 0.0 and result <= 1.0);
  return result;
}

inline
double weighting(const double& ui,
                 const double& uj,
                 const double& mi,
                 const double& mj,
                 const double& duij,
                 const double& dt) {

  // First use our standard weighting algorithm.
  const double fi = standardWeighting(ui, uj, mi, mj, duij);
  const double fj = 1.0 - fi;
  CHECK(fi >= 0.0 and fi <= 1.0);
  CHECK(fj >= 0.0 and fj <= 1.0);
  CHECK(fuzzyEqual(fi + fj, 1.0));

  // Now the monotonic weighting.
  const double mfi = monotonicWeighting(ui, uj, mi, mj, mi*duij*dt);

  // Combine them for the final answer.
  const double chi = abs(ui - uj)/(abs(ui) + abs(uj) + 1.0e-50);
  CHECK(chi >= 0.0 and chi <= 1.0);
  return chi*mfi + (1.0 - chi)*fi;
}

// inline
// double weighting(const double& ui,
//                  const double& uj,
//                  const double& mi,
//                  const double& mj,
//                  const double& duij,
//                  const double& dt) {

//   // First use our standard weighting algorithm.
//   const double fi = standardWeighting(ui, uj, mi, mj, duij);
//   const double fj = 1.0 - fi;
//   CHECK(fi >= 0.0 and fi <= 1.0);
//   CHECK(fj >= 0.0 and fj <= 1.0);
//   CHECK(fuzzyEqual(fi + fj, 1.0));

//   // Check if this would result in one (but not both) energies
//   // flipping sign.  If so, we revert to the strictly monotonic, but not
//   // quite as accurate scheme.
//   const double ui1 = ui + fi*duij*dt;
//   const double uj1 = uj + fj*duij*dt;
//   if (sgn(ui*uj) == sgn(ui1*uj1)) {
//     return fi;
//   } else {
//     return monotonicWeighting(ui, uj, mi, mj, mi*duij*dt);
//   }
// }

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SpecificThermalEnergyPolicy<Dimension>::
SpecificThermalEnergyPolicy(const DataBase<Dimension>& dataBase):
  IncrementState<Dimension, typename Dimension::Scalar>(),
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
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy);

  typedef typename Dimension::SymTensor SymTensor;

//   // HACK!
//   std::cerr.setf(std::ios::scientific, std::ios::floatfield);
//   std::cerr.precision(15);

  // Grab the state field.
  Field<Dimension, Scalar>& eps = state.field(key, Scalar());

  // Get the FluidNodeList and check if we're enforcing compatible energy 
  // evolution or not.
  const FluidNodeList<Dimension>* nodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(eps.nodeListPtr());
  CHECK(nodeListPtr != 0);

  // Get the state fields.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, Scalar());
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Vector> acceleration = derivs.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> specificEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const FieldList<Dimension, vector<Vector> > pairAccelerations = derivs.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  const ConnectivityMap<Dimension>& connectivityMap = mDataBasePtr->connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();

  // Find this NodeList.
  const int nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), nodeListPtr));
  const int numNodeLists = nodeLists.size();
  const Field<Dimension, Scalar>& massi = *mass[nodeListi];
  const Field<Dimension, Vector>& velocityi = *velocity[nodeListi];
  const Field<Dimension, Vector>& accelerationi = *acceleration[nodeListi];
  const Field<Dimension, vector<Vector> >& pairAccelerationsi = *pairAccelerations[nodeListi];
  const Field<Dimension, Scalar>& eps0i = *specificEnergy0[nodeListi];

  // Iterate over the internal nodes of this NodeList.
  const double hdt = 0.5*multiplier;
  for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
       iItr != connectivityMap.end(nodeListi);
       ++iItr) {
    const int i = *iItr;
    if (nodeListPtr->nodeType(i) == NodeSpace::InternalNode) {

      // State for node i.
      const Scalar& mi = massi(i);
      const Vector& vi = velocityi(i);
      const Scalar& ui = eps0i(i);
      const Vector& ai = accelerationi(i);
      const Vector vi12 = vi + ai*hdt;
      const vector<Vector>& pacci = pairAccelerationsi(i);
      CHECK(pacci.size() == connectivityMap.numNeighborsForNode(nodeListPtr, i));

      // Get the connectivity (neighbor set) for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListPtr, i);

      // Iterate over the neighbor NodeLists.
      int jj = 0;
      Scalar DepsDti = 0.0;

      // The following bit of two loop insanity is the result of an optimization we employ in 
      // iterating over the neighbors in SPHFluidDerivatives::calculateDerivatives.  The problem is that
      // in that method we only hit each neighbor pair once, updating the the derivative state
      // for each node pair simultaneously.  This means that we build up the vector of pairwise
      // accelerations in a rather non-intuitive manner, and here we must reproduce that ordering
      // exactly.   This really raises it's head when you multiple NodeLists and ghost nodes.
      // Hence the following two loops reproduce the order that each node pair for this node
      // was computed precisely.
    
      // First we do all the nodes that in SPHFluidDerivatives will have
      // been calculated before we did the loop on this node, since they
      // will be first in the list of partial accelerations.
      for (int nodeListj = 0; nodeListj <= nodeListi; ++nodeListj) {

        // The set of neighbors from this NodeList.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const Field<Dimension, Scalar>& massj = *mass[nodeListj];
          const Field<Dimension, Vector>& velocityj = *velocity[nodeListj];
          const Field<Dimension, Vector>& accelerationj = *acceleration[nodeListj];
          const Field<Dimension, Scalar>& eps0j = *specificEnergy0[nodeListj];
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors, and accumulate the specific energy
          // change.
          for (vector<int>::const_iterator jitr = connectivity.begin();
               jitr != connectivity.end();
               ++jitr) {
            const int j = *jitr;

            if (not connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                             nodeListj, j,
                                                             firstGhostNodej)) {
              CHECK(jj < pacci.size());
              const Vector& pai = pacci[jj];
              const Scalar& mj = massj(j);
              const Vector& vj = velocityj(j);
              const Scalar& uj = eps0j(j);
              const Vector& aj = accelerationj(j);
              const Vector vji12 = vj + aj*hdt - vi12;
              const Scalar duij = vji12.dot(pai);
              const Scalar wi = weighting(ui, uj, mi, mj, duij, dt);
              CHECK(wi >= 0.0 and wi <= 1.0);
              CHECK(fuzzyEqual(wi + weighting(uj, ui, mj, mi, duij*mi/mj, dt), 1.0, 1.0e-10));
              DepsDti += wi*duij;
              ++jj;
            }
          }
        }
      }

      // Now hit all the nodes that will have been computed in the loop
      // over this nodes neighbors.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // The set of neighbors from this NodeList.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const Field<Dimension, Scalar>& massj = *mass[nodeListj];
          const Field<Dimension, Vector>& velocityj = *velocity[nodeListj];
          const Field<Dimension, Vector>& accelerationj = *acceleration[nodeListj];
          const Field<Dimension, Scalar>& eps0j = *specificEnergy0[nodeListj];
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors, and accumulate the specific energy
          // change.
          for (vector<int>::const_iterator jitr = connectivity.begin();
               jitr != connectivity.end();
               ++jitr) {
            const int j = *jitr;

            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              CHECK(jj < pacci.size());
              const Vector& pai = pacci[jj];
              const Scalar& mj = massj(j);
              const Vector& vj = velocityj(j);
              const Scalar& uj = eps0j(j);
              const Vector& aj = accelerationj(j);
              const Vector vji12 = vj + aj*hdt - vi12;
              const Scalar duij = vji12.dot(pai);
              //             const Vector vj12 = vj + aj*hdt;
              //             const Scalar duij = (vj12 - vi12).dot(pai);
              const Scalar wi = weighting(ui, uj, mi, mj, duij, dt);
              CHECK(wi >= 0.0 and wi <= 1.0);
              DepsDti += wi*duij;
              ++jj;
            }
          }
        }
      }

      CHECK(jj == pacci.size());

      // Now we can update the energy.
      eps(i) += DepsDti*multiplier;
      //     VERIFY(eps(i) >= 0.0);
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

