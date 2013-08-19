//---------------------------------Spheral++----------------------------------//
// CompatibleFaceSpecificThermalEnergyPolicy
//
// Created by JMO, Fri Aug 16 14:47:57 PDT 2013
//----------------------------------------------------------------------------//
#include <vector>

#include "CompatibleFaceSpecificThermalEnergyPolicy.hh"
#include "computeSVPHCorrectionsOnFaces.hh"
#include "Hydro/HydroFieldNames.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Kernel/TableKernel.hh"
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
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

using namespace std;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using MeshSpace::Mesh;
using NeighborSpace::Neighbor;
using ArtificialViscositySpace::ArtificialViscosity;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
CompatibleFaceSpecificThermalEnergyPolicy(const TableKernel<Dimension>& W,
                                          const DataBase<Dimension>& dataBase,
                                          const ArtificialViscosity<Dimension>& Q,
                                          const bool linearConsistent):
  IncrementState<Dimension, typename Dimension::Scalar>(HydroFieldNames::mesh,
                                                        HydroFieldNames::volume,
                                                        HydroFieldNames::mass,
                                                        HydroFieldNames::position,
                                                        HydroFieldNames::velocity,
                                                        HydroFieldNames::H),
  mW(W),
  mDataBase(dataBase),
  mQ(Q),
  mLinearConsistent(linearConsistent) {
  mFired = false;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
~CompatibleFaceSpecificThermalEnergyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  typedef typename Mesh<Dimension>::Zone Zone;
  typedef typename Mesh<Dimension>::Face Face;

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy);

  if (not mFired) {
    mFired = true;

    // Get the state fields.
    FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
    const Mesh<Dimension>& mesh = state.mesh();
    const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
    const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
    const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> DepsDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
    const ConnectivityMap<Dimension>& connectivityMap = mDataBase.connectivityMap();
    const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
    const size_t numNodeLists = nodeLists.size();

    const unsigned numFaces = mesh.numFaces();
    vector<Scalar> Pface(numFaces, 0.0);
    vector<Vector> dAface(numFaces, Vector::zero), posFace(numFaces, Vector::zero), velFace(numFaces, Vector::zero);
    vector<Tensor> Qface(numFaces, Tensor::zero);

    // Compute the SVPH corrections.
    vector<Scalar> A;
    vector<Vector> B;
    SVPHSpace::computeSVPHCorrectionsOnFaces(mesh, mW,
                                             volume, position, H,
                                             A, B);
    if (not mLinearConsistent) B = vector<Vector>(numFaces, Vector::zero);

    // Walk the faces and interpolate the current velocity and old pressure.
    const SymTensor Hface = 1.0e100*SymTensor::one;
    unsigned i, j, k, z1id, z2id, nodeListi, nodeListj;
    Scalar DepsDti;
    for (k = 0; k != numFaces; ++k) {
      const Face& face = mesh.face(k);
      posFace[k] = face.position();
      dAface[k] = face.area() * face.unitNormal();

      const Scalar& Ai = A[k];
      const Vector& Bi = B[k];

      // Set the neighbors for this face.
      Neighbor<Dimension>::setMasterNeighborGroup(posFace[k], Hface,
                                                  nodeLists.begin(), nodeLists.end(),
                                                  mW.kernelExtent());

      // Iterate over the NodeLists.
      for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const NodeList<Dimension>& nodeList = *nodeLists[nodeListj];
        Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
        neighbor.setRefineNeighborList(posFace[k], Hface);
        for (typename Neighbor<Dimension>::const_iterator neighborItr = neighbor.refineNeighborBegin();
             neighborItr != neighbor.refineNeighborEnd();
             ++neighborItr) {
          j = *neighborItr;
      
          // Get the state for node j
          const Vector& rj = position(nodeListj, j);
          const Vector& vj = velocity(nodeListj, j);
          const Scalar& Pj = pressure(nodeListj, j);
          const SymTensor& Hj = H(nodeListj, j);
          const Scalar& Vj = volume(nodeListj, j);
          const Scalar Hdetj = Hj.Determinant();
          CHECK(Vj > 0.0);
          CHECK(Hdetj > 0.0);

          // Pair-wise kernel type stuff.
          const Vector rij = posFace[k] - rj;
          const Vector etaj = Hj*rij;
          const Scalar Wj = mW.kernelValue(etaj.magnitude(), Hdetj);

          // Increment the face fluid properties.
          const Scalar VWRj = Vj*(1.0 + Bi.dot(rij))*Wj;
          velFace[k] += VWRj*vj;
          Pface[k] += VWRj*Pj;
        }
      }

      // Finish the face state.
      CHECK2(Ai >= 0.0, i << " " << Ai);
      velFace[k] *= Ai;
      Pface[k] *= Ai;

      // Find the SVPH nodes on either side of this face.
      z1id = Mesh<Dimension>::positiveID(face.zone1ID());
      z2id = Mesh<Dimension>::positiveID(face.zone2ID());
      if (z1id != Mesh<Dimension>::UNSETID and
          z2id != Mesh<Dimension>::UNSETID) {
        mesh.lookupNodeListID(z1id, nodeListi, i);
        mesh.lookupNodeListID(z2id, nodeListj, j);

        // Get the node properties.
        const Vector& ri = position(nodeListi, i);
        const Vector& vi = velocity(nodeListi, i);
        const Scalar& rhoi = massDensity(nodeListi, i);
        const Scalar& ci = soundSpeed(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);

        const Vector& rj = position(nodeListj, j);
        const Vector& vj = velocity(nodeListj, j);
        const Scalar& rhoj = massDensity(nodeListj, j);
        const Scalar& cj = soundSpeed(nodeListj, j);
        const SymTensor& Hj = H(nodeListj, j);

        const Vector rij = ri - rj;
        const Vector etai = Hi*rij;
        const Vector etaj = Hj*rij;

        // Get the face Q values (in this case P/rho^2).
        const pair<Tensor, Tensor> QPiij = mQ.Piij(nodeListi, i, nodeListj, j,
                                                   ri, etai, vi, rhoi, ci, Hi,
                                                   rj, etaj, vj, rhoj, cj, Hj);
        Qface[k] = 0.5*(rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second);
      }
    }

    // Start our big loop over all FluidNodeLists.
    nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = mDataBase.fluidNodeListBegin();
         itr != mDataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      const NodeList<Dimension>& nodeList = **itr;
      const int firstGhostNodei = nodeList.firstGhostNode();

      // Iterate over the internal nodes in this NodeList.
      for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
           iItr != connectivityMap.end(nodeListi);
           ++iItr) {
        const int i = *iItr;
        const Zone& zonei = mesh.zone(nodeListi, i);
        const vector<int>& faceIDs = zonei.faceIDs();
        const unsigned nfaces = faceIDs.size();

        // Walk the faces and sum the face work.
        DepsDti = 0.0;
        for (k = 0; k != nfaces; ++k) {  
          const unsigned fid = Mesh<Dimension>::positiveID(faceIDs[k]);
          const Vector dA = dAface[fid] * sgn(faceIDs[k]);
          const Vector fforce = Pface[fid]*dA + Qface[fid]*dA;
          DepsDti += fforce.dot(velFace[fid]);
        }
        DepsDti /= mass(nodeListi, i);

        // Now the final DepsDt is the average of the initial and final states.
        DepsDti = 0.5*(DepsDti + DepsDt(nodeListi, i));
        eps(nodeListi, i) += multiplier*DepsDti;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also a CompatibleSpecificThermalEnergy operator.
  const CompatibleFaceSpecificThermalEnergyPolicy<Dimension>* rhsPtr = dynamic_cast<const CompatibleFaceSpecificThermalEnergyPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

//------------------------------------------------------------------------------
// Static initializations.
//------------------------------------------------------------------------------
template<typename Dimension> bool CompatibleFaceSpecificThermalEnergyPolicy<Dimension>::mFired = false;

}

