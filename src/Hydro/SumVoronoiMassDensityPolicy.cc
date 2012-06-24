//---------------------------------Spheral++----------------------------------//
// SumVoronoiMassDensityPolicy -- An implementation of UpdatePolicyBase 
// specialized for the updating the mass density according to the specific 
// volume from the Voronoi tesselation.  This version uses a weighted sum of the
// local mass and volume to do the deed.
//
// Created by JMO, Mon Aug  1 10:48:03 PDT 2011
//----------------------------------------------------------------------------//

#include "SumVoronoiMassDensityPolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using FieldSpace::Field;
using NeighborSpace::ConnectivityMap;
using KernelSpace::TableKernel;
using PhysicsSpace::Physics;
using BoundarySpace::Boundary;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SumVoronoiMassDensityPolicy<Dimension>::
SumVoronoiMassDensityPolicy(const TableKernel<Dimension>& W, 
                            const Physics<Dimension>& package,
                            const double rhoMin, const double rhoMax):
  ReplaceState<Dimension, typename Dimension::Scalar>(HydroFieldNames::mass,
                                                      HydroFieldNames::volume),
  mW(W),
  mPackage(package),
  mRhoMin(rhoMin),
  mRhoMax(rhoMax) {
  REQUIRE(rhoMin <= rhoMax);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SumVoronoiMassDensityPolicy<Dimension>::
~SumVoronoiMassDensityPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SumVoronoiMassDensityPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity);
  Field<Dimension, Scalar>& massDensity = state.field(key, 0.0);

  // Intitialize the mass density.
  massDensity = 0.0;

  // We'll also need a place to keep the volume as we accumulate it.
  const NodeList<Dimension>& nodeList = massDensity.nodeList();
  Field<Dimension, Scalar> volEff("effective volume", nodeList);

  // Grab the state we need to do our job.
  const KeyType massKey = State<Dimension>::buildFieldKey(HydroFieldNames::mass, nodeListKey);
  const KeyType volKey = State<Dimension>::buildFieldKey(HydroFieldNames::volume, nodeListKey);
  const KeyType posKey = State<Dimension>::buildFieldKey(HydroFieldNames::position, nodeListKey);
  const KeyType HKey = State<Dimension>::buildFieldKey(HydroFieldNames::H, nodeListKey);
  CHECK(state.registered(massKey));
  CHECK(state.registered(volKey));
  CHECK(state.registered(posKey));
  CHECK(state.registered(HKey));
  const Field<Dimension, Scalar>& mass = state.field(massKey, 0.0);
  const Field<Dimension, Scalar>& volume = state.field(volKey, 0.0);
  const Field<Dimension, Vector>& pos = state.field(posKey, Vector::zero);
  const Field<Dimension, SymTensor>& H = state.field(HKey, SymTensor::zero);
  const ConnectivityMap<Dimension>& cm = state.connectivityMap();

  // // The volume still needs to have boundary conditions enforced.
  // typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  // for (ConstBoundaryIterator boundaryItr = mPackage.boundaryBegin();
  //      boundaryItr != mPackage.boundaryEnd();
  //      ++boundaryItr) {
  //   (*boundaryItr)->applyGhostBoundary(mass);
  //   (*boundaryItr)->applyGhostBoundary(volume);
  // }
  // for (ConstBoundaryIterator boundaryItr = mPackage.boundaryBegin();
  //      boundaryItr != mPackage.boundaryEnd();
  //      ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // Which NodeList are we doing in the ConnectivityMap?
  const unsigned nodeListi = cm.nodeListIndex(&nodeList);

  // Walk the local nodes.
  typename ConnectivityMap<Dimension>::const_iterator iItr;
  typename vector<int>::const_iterator jItr;
  unsigned i, j;
  Scalar Hdeti, Hdetj, etai, etaj, Wi, Wj, W0;
  Vector xij;
  const unsigned firstGhostNodei = nodeList.firstGhostNode();
  for (iItr = cm.begin(nodeListi); iItr != cm.end(nodeListi); ++iItr) {
    i = *iItr;
    const Scalar& Vi = volume(i);
    const Scalar& mi = mass(i);
    const Vector& xi = pos(i);
    const SymTensor& Hi = H(i);
    Hdeti = Hi.Determinant();
    Scalar& rhoi = massDensity(i);
    Scalar& Veffi = volEff(i);

    // Walk the neighbors of this node.
    const vector<int>& connectivity = cm.connectivityForNode(nodeListi, i)[nodeListi];
    for (jItr = connectivity.begin(); jItr != connectivity.end(); ++jItr) {
      j = *jItr;
      if (cm.calculatePairInteraction(nodeListi, i, 
                                      nodeListi, j,
                                      firstGhostNodei)) {
        const Scalar& Vj = volume(j);
        const Scalar& mj = mass(j);
        const Vector& xj = pos(j);
        const SymTensor& Hj = H(j);
        Hdetj = Hj.Determinant();
        Scalar& rhoj = massDensity(j);
        Scalar& Veffj = volEff(j);

        xij = xi - xj;
        etai = (Hi*xij).magnitude();
        etaj = (Hj*xij).magnitude();
        Wi = mW.kernelValue(etai, Hdeti);
        Wj = mW.kernelValue(etaj, Hdetj);

        Veffi += Vj*Wi;
        rhoi +=  mj*Wi;

        Veffj += Vi*Wj;
        rhoj +=  mi*Wj;
      }
    }

    // Finalize the density for node i.
    W0 = mW.kernelValue(0.0, Hdeti);
    rhoi = max(mRhoMin, min(mRhoMax, (rhoi + mi*W0) * safeInv(Veffi + Vi*W0)));
  }
}

//------------------------------------------------------------------------------
// Update the density as an increment.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SumVoronoiMassDensityPolicy<Dimension>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {

  // Find the matching derivative field from the StateDerivatives.
  KeyType incrementKey = IncrementState<Dimension, Scalar>::prefix() + key;
  FieldSpace::Field<Dimension, Scalar>& f = state.field(key, 0.0);
  const FieldSpace::Field<Dimension, Scalar>& df = derivs.field(incrementKey, 0.0);

  // Loop over the internal values of the field.
  for (unsigned i = 0; i != f.nodeList().numInternalNodes(); ++i) {
    f(i) = max(mRhoMin, min(mRhoMax, f(i) + multiplier*(df(i))));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SumVoronoiMassDensityPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const SumVoronoiMassDensityPolicy<Dimension>* rhsPtr = dynamic_cast<const SumVoronoiMassDensityPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class SumVoronoiMassDensityPolicy<Dim<1> >;
  template class SumVoronoiMassDensityPolicy<Dim<2> >;
  template class SumVoronoiMassDensityPolicy<Dim<3> >;
}

