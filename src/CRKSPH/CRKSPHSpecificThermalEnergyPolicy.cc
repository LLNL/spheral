//---------------------------------Spheral++----------------------------------//
// CRKSPHSpecificThermalEnergyPolicy -- An implementation of UpdatePolicyBase
// specialized for the updating the specific thermal energy as a dependent 
// quantity.
// 
// This version is specialized for the compatible energy discretization 
// for use with CRKSPH.
//
// Created by JMO, Sat Aug  2 06:44:50 PDT 2014
//----------------------------------------------------------------------------//
#include <vector>

#include "CRKSPHSpecificThermalEnergyPolicy.hh"
#include "computeCRKSPHCorrections.hh"
#include "Hydro/HydroFieldNames.hh"
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
double PoverRho2Weighting(const double Pi,
                          const double rhoi,
                          const double Pj,
                          const double rhoj) {
  const double wi = max(1.0e-100, Pi/max(1.0e-100, rhoi*rhoi));
  const double wj = max(1.0e-100, Pj/max(1.0e-100, rhoj*rhoj));
  return wi/(wi + wj);
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

}

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHSpecificThermalEnergyPolicy<Dimension>::
CRKSPHSpecificThermalEnergyPolicy(const DataBase<Dimension>& dataBase,
                                const TableKernel<Dimension>& W):
  IncrementFieldList<Dimension, typename Dimension::Scalar>(),
  mDataBasePtr(&dataBase),
  mWT(W) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHSpecificThermalEnergyPolicy<Dimension>::
~CRKSPHSpecificThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHSpecificThermalEnergyPolicy<Dimension>::
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
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, Scalar());
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, Scalar());
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, Scalar());
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, Scalar());
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Vector> acceleration = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> DepsDt0 = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);

  const FieldList<Dimension, Scalar> eps0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const FieldList<Dimension, vector<Vector> > pairAccelerations = derivs.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  const ConnectivityMap<Dimension>& connectivityMap = mDataBasePtr->connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  CHECK(nodeLists.size() == numFields);

  // Prepare a counter to keep track of how we go through the pair-accelerations.
  FieldList<Dimension, Scalar> DepsDt(FieldStorageType::CopyFields);
  FieldList<Dimension, int> offset(FieldStorageType::CopyFields);
  for (size_t nodeListi = 0; nodeListi != numFields; ++nodeListi) {
    DepsDt.appendNewField("delta E", *nodeLists[nodeListi], 0.0);
    offset.appendNewField("offset", *nodeLists[nodeListi], 0);
  }

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
      const Vector& ri = position(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Scalar& voli = volume(nodeListi, i);
      const Scalar& Pi = pressure(nodeListi, i);
      const Scalar& rhoi = rho(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& ui = eps0(nodeListi, i);
      const Vector& ai = acceleration(nodeListi, i);
      const Scalar& DepsDt0i = DepsDt0(nodeListi, i);
      const Scalar& A0i = A0(nodeListi, i);
      const Vector vi12 = vi + ai*hdt;
      const vector<Vector>& pacci = pairAccelerations(nodeListi, i);
      CHECK(pacci.size() == connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i) + 1);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Get the connectivity (neighbor set) for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);

      // The total energy change due to work on point i.
      const Scalar dEi = -mi*vi12.dot(ai);

      // Find the discrepency from the time evolution form.
      const Scalar deltai = dEi/dt - DepsDt0i;
      DepsDti += DepsDt0i + A0i*voli*mWT.kernelValue(0.0, Hdeti)*deltai;

      // Distribute the discrepancy across our neighbors.
      // Iterate over the neighbor NodeLists.
      for (size_t nodeListj = 0; nodeListj != numFields; ++nodeListj) {

        // The set of neighbors from this NodeList.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors in this NodeList.
          for (vector<int>::const_iterator jitr = connectivity.begin();
               jitr != connectivity.end();
               ++jitr) {
            const int j = *jitr;

            Scalar& DepsDtj = DepsDt(nodeListj, j);
            const Vector& rj = position(nodeListj, j);
            const Scalar& mj = mass(nodeListj, j);
            const Scalar& volj = volume(nodeListj, j);
            const Scalar& Pj = pressure(nodeListj, j);
            const Scalar& rhoj = rho(nodeListj, j);
            const Vector& vj = velocity(nodeListj, j);
            const Scalar& uj = eps0(nodeListj, j);
            const Vector& aj = acceleration(nodeListj, j);
            const Scalar& DepsDt0j = DepsDt0(nodeListj, j);
            const Scalar& A0j = A0(nodeListj, j);
            const Vector vj12 = vj + aj*hdt;
            const vector<Vector>& paccj = pairAccelerations(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();

            const Vector rij = ri - rj;
            // const Scalar etai = (Hi*rij).magnitude();
            const Scalar etaj = (Hj*rij).magnitude();
            // const Scalar Wi = A0j*mWT.kernelValue(etai, Hdeti);
            const Scalar Wj = A0i*mWT.kernelValue(etaj, Hdetj);

            DepsDtj += A0i*volj*Wj*deltai;

              // DepsDtj += dEi*Wj/mj;
              // DepsDti += dEj*Wi/mi;

              // const Scalar dEij = dEi*Wj + dEj*Wi;

              // CHECK(offset(nodeListi, i) < pacci.size());
              // const Vector& pai = pacci[offset(nodeListi, i)];
              // ++offset(nodeListi, i);

              // CHECK(offset(nodeListj, j) < paccj.size());
              // const Vector& paj = paccj[offset(nodeListj, j)];
              // ++offset(nodeListj, j);

              // const Scalar dEij = -(mi*vi12.dot(pai) + mj*vj12.dot(paj));
              // const Scalar duij = dEij/mi;
              // // const Scalar wi = monotonicWeighting(ui, uj, mi, mj, dEij);
              // // const Scalar wi = standardWeighting(ui, uj, mi, mj, duij)
              // // const Scalar wi = PoverRho2Weighting(Pi, rhoi, Pj, rhoj);
              // const Scalar wi = weighting(ui, uj, mi, mj, duij, dt);

              // // const Scalar wi = DepsDt0i*safeInv(DepsDt0i + DepsDt0j);

              // CHECK(wi >= 0.0 and wi <= 1.0);
              // CHECK(fuzzyEqual(wi + weighting(uj, ui, mj, mi, dEij/mj, dt), 1.0, 1.0e-10));
              // DepsDti += wi*duij;
              // DepsDtj += (1.0 - wi)*dEij/mj;
            // }
          }
        }
      }

      // // Add the self-contribution.
      // const Scalar duii = -vi12.dot(pacci.back());
      // DepsDti += duii;

      // // Grab the self-interaction term.  We will distribute this amongst our neighbors.
      // const unsigned numNeighbors = connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i);
      // const Scalar dEii = -mi*vi12.dot(pacci.back())/numNeighbors;
      // for (size_t nodeListj = 0; nodeListj != numFields; ++nodeListj) {
      //   const vector<int>& connectivity = fullConnectivity[nodeListj];
      //   if (connectivity.size() > 0) {
      // 	  const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
      //     for (vector<int>::const_iterator jitr = connectivity.begin();
      //          jitr != connectivity.end();
      //          ++jitr) {
      //       const int j = *jitr;
      // 	    Scalar& DepsDtj = DepsDt(nodeListj, j);
      // 	    const Scalar& mj = mass(nodeListj, j);
      // 	    DepsDtj += 0.5*dEii;
      // 	    DepsDti += 0.5*dEii;
      // 	  }
      // 	}
      // }
      //      DepsDti += dEii;

      // Now we can update the energy.
      CHECK(offset(nodeListi, i) == pacci.size());
      eps(nodeListi, i) += DepsDti*multiplier;
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CRKSPHSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const CRKSPHSpecificThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const CRKSPHSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

