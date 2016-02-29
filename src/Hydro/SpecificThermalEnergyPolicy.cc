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

namespace {

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
//   const double result = 0.5*(1.0 + uji/(std::abs(uji) + 1.0/(1.0 + std::abs(uji)))*sgn0(duij));
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
  const double result = 0.5*(1.0 + uji*sgn0(duij)/(std::abs(uji) + 0.5));
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
    const double A = std::min(1.0, std::max(0.0, uji*mi*safeInv(DEij, 1.0e-50)));
    const double B = std::min(1.0, std::max(0.0, mjInv/(miInv + mjInv)));
    CHECK(A >= 0.0 and A <= 1.0);
    CHECK(B >= 0.0 and B <= 1.0);
    result = A + B*(1.0 - A);
  } else {
    const double A = std::min(1.0, std::max(0.0, -uji*mj*safeInv(DEij, 1.0e-50)));
    const double B = std::min(1.0, std::max(0.0, miInv/(miInv + mjInv)));
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
  CHECK(fi >= 0.0 and fi <= 1.0);

  // Now the monotonic weighting.
  const double mfi = monotonicWeighting(ui, uj, mi, mj, mi*duij*dt);

  // Combine them for the final answer.
  const double chi = std::abs(ui - uj)/(std::abs(ui) + std::abs(uj) + 1.0e-50);
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

}

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
  FieldList<Dimension, Scalar> eps = state.fields(fieldKey, Scalar());
  const unsigned numFields = eps.numFields();

  // Get the state fields.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, Scalar());
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Vector> acceleration = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> eps0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, Scalar());
  const FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, Scalar());
  const FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, Scalar());
  const FieldList<Dimension, vector<Vector> > pairAccelerations = derivs.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  const ConnectivityMap<Dimension>& connectivityMap = mDataBasePtr->connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  CHECK(nodeLists.size() == numFields);

  // Prepare a counter to keep track of how we go through the pair-accelerations.
  FieldList<Dimension, Scalar> DepsDt = mDataBasePtr->newFluidFieldList(0.0, "delta E");
  FieldList<Dimension, int> offset = mDataBasePtr->newFluidFieldList(0, "offset");

  // Walk all the NodeLists.
  const double hdt = 0.5*multiplier;
  for (size_t nodeListi = 0; nodeListi != numFields; ++nodeListi) {

    // Iterate over the internal nodes of this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // State for node i.
      Scalar& DepsDti = DepsDt(nodeListi, i);
      const Scalar mi = mass(nodeListi, i);
      const Scalar Ai = std::abs(P(nodeListi, i)/std::pow(rho(nodeListi, i), gamma(nodeListi, i)));
      const Vector& vi = velocity(nodeListi, i);
      const Scalar ui = eps0(nodeListi, i);
      const Vector& ai = acceleration(nodeListi, i);
      const Vector vi12 = vi + ai*hdt;
      const vector<Vector>& pacci = pairAccelerations(nodeListi, i);
      CHECK(pacci.size() == connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i) or
            NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());

      // Get the connectivity (neighbor set) for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      // Iterate over the neighbor NodeLists.
      for (size_t nodeListj = 0; nodeListj != numFields; ++nodeListj) {

        // The set of neighbors from this NodeList.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
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
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              const Scalar mj = mass(nodeListj, j);
              const Scalar Aj = std::abs(P(nodeListj, j)/std::pow(rho(nodeListj, j), gamma(nodeListj, j)));
              const Vector& vj = velocity(nodeListj, j);
              const Scalar uj = eps0(nodeListj, j);
              const Vector& aj = acceleration(nodeListj, j);
              const Vector vj12 = vj + aj*hdt;
              const Vector vji12 = vj12 - vi12;
              const vector<Vector>& paccj = pairAccelerations(nodeListj, j);
              CHECK(j >= firstGhostNodej or 
                    paccj.size() == connectivityMap.numNeighborsForNode(nodeLists[nodeListj], j) or
                    NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent());

              CHECK(offset(nodeListi, i) < pacci.size());
              const Vector& pai =  pacci[offset(nodeListi, i)];
              ++offset(nodeListi, i);

              CHECK(offset(nodeListj, j) < paccj.size());
              const Vector& paj =  paccj[offset(nodeListj, j)];
              ++offset(nodeListj, j);

              CHECK2(fuzzyEqual(mi*mj*pai.dot(paj) + mi*mi*pai.dot(pai), 0.0, 1.0e-10),
                     "Symmetric forces?  (" << nodeListi << " " << i << ") (" << nodeListj << " " << j << ") " << mi << " " << mj << " " << pai << " " << paj << " " << mi*pai << " " << mj*paj);

              const Scalar duij = vji12.dot(pai);
              const Scalar wi = (Ai + Aj == 0 ? 0.5 :
                                 (duij >= 0.0 ? Aj : Ai)/(Ai + Aj));
              // const Scalar wi = (duij >= 0.0 ? 
              //                    safeInvVar(Ai)/(safeInvVar(Ai) + safeInvVar(Aj)) :
              //                    (Ai + Aj == 0.0 ? 0.5 : Ai/(Ai + Aj)));
              // const Scalar wi = (duij >= 0.0 ? 0.5 : weighting(ui, uj, mi, mj, duij, dt));

              CHECK(wi >= 0.0 and wi <= 1.0);
              // CHECK(fuzzyEqual(wi + weighting(uj, ui, mj, mi, duij*mi/mj, dt), 1.0, 1.0e-10));
              DepsDti += wi*duij;
              DepsDtj += (1.0 - wi)*duij*mi/mj;
            }
          }
        }
      }

      // Now we can update the energy.
      CHECK2(offset(nodeListi, i) == pacci.size() or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent(),
             "Bad sizing : (" << nodeListi << " " << i << ") " << offset(nodeListi, i) << " " << pacci.size());
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

