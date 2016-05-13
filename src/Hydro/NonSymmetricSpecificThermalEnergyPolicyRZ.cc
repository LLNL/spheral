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

namespace {

//------------------------------------------------------------------------------
// The entropy weighted energy form.
//------------------------------------------------------------------------------
inline
double
entropyWeighting(const double si,
                 const double sj,
                 const double ri,
                 const double rj,
                 const double zi,
                 const double zj,
                 const double duij) {
  if (ri <= 0.0) { //  and fuzzyEqual(zi, zj, 1.0e-8)) {
    return 0.0;
  } else if (rj <= 0.0) { // and fuzzyEqual(zi, zj, 1.0e-8)) {
    return 1.0;
  } else {
    double result = 0.5;
    const double smin = min(abs(si), abs(sj));
    const double smax = max(abs(si), abs(sj));
    if (smax > 1.0e-15) {
      CHECK(smin + smax > 1.0e-15);
      if (duij > 0.0) {    // Heating
        if (si > sj) {
          result = smin/(smin + smax);
        } else {
          result = smax/(smin + smax);
        }
      } else {             // Cooling
        if (si > sj) {
          result = smax/(smin + smax);
        } else {
          result = smin/(smin + smax);
        }
      }
    }
    CHECK(result >= 0.0 and result <= 1.0);
    return result;
  }
}

}

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
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, Scalar());
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Vector> acceleration = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> eps0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", Scalar());
  const FieldList<Dimension, Scalar> entropy = state.fields(HydroFieldNames::entropy, Scalar());
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
    for (ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // State for node i.
      Scalar& DepsDti = DepsDt(nodeListi, i);
      const Scalar ri = abs(position(nodeListi, i).y());
      const Scalar zi = position(nodeListi, i).x();
      const Scalar mi = mass(nodeListi, i);
      const Scalar mRZi = mi/(2.0*M_PI*ri);
      const Scalar si = entropy(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar ui = eps0(nodeListi, i);
      const Vector& ai = acceleration(nodeListi, i);
      const Vector vi12 = vi + ai*hdt;
      // vi12.x(vi12.x()/(2.0*M_PI*ri));
      const vector<Vector>& pacci = pairAccelerations(nodeListi, i);
      CHECK(ri > 0.0);
      CHECK(pacci.size() == connectivityMap.numNeighborsForNode(nodeLists[nodeListi], i) + 1);

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
              const Scalar rj = abs(position(nodeListj, j).y());
              const Scalar zj = position(nodeListj, j).x();
              const Scalar mj = mass(nodeListj, j);
              const Scalar mRZj = mj/(2.0*M_PI*ri);
              const Scalar sj = entropy(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar uj = eps0(nodeListj, j);
              const Vector& aj = acceleration(nodeListj, j);
              const Vector vj12 = vj + aj*hdt;
              // vj12.x(vj12.x()/(2.0*M_PI*rj));
              const vector<Vector>& paccj = pairAccelerations(nodeListj, j);
              CHECK(rj > 0.0);
              CHECK(j >= firstGhostNodej or paccj.size() == (connectivityMap.numNeighborsForNode(nodeLists[nodeListj], j) + 1));

              CHECK(offset(nodeListi, i) < pacci.size());
              const Vector& pai = pacci[offset(nodeListi, i)];
              ++offset(nodeListi, i);

              CHECK(offset(nodeListj, j) < paccj.size());
              const Vector& paj = paccj[offset(nodeListj, j)];
              ++offset(nodeListj, j);

              const Scalar dEij = -(mi*vi12.dot(pai) + mj*vj12.dot(paj));
              const Scalar duij = dEij/mi;
              const Scalar wi = entropyWeighting(si, sj, ri, rj, zi, zj, duij);

              // const Scalar dERZij = -(mRZi*vi12.dot(pai) + mRZj*vj12.dot(paj));
              // const Scalar duRZij = dERZij/mRZi;
              // const Scalar dEij = -(mi*vi12.dot(pai) + mj*vj12.dot(paj));
              // const Scalar duij = dEij/mi;
              // const Scalar wi = entropyWeighting(si, sj, duRZij);

              CHECK(wi >= 0.0 and wi <= 1.0);
              CHECK(fuzzyEqual(wi + entropyWeighting(sj, si, zj, zi, dEij/mj), 1.0, 1.0e-10));
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

