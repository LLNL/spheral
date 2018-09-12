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
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/DBC.hh"

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
SumVoronoiMassDensityPolicy<Dimension>::
SumVoronoiMassDensityPolicy(const TableKernel<Dimension>& W, 
                            const Physics<Dimension>& package,
                            const double rhoMin, const double rhoMax):
  ReplaceFieldList<Dimension, typename Dimension::Scalar>(HydroFieldNames::mass,
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
  REQUIRE(fieldKey == HydroFieldNames::massDensity and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> massDensity = state.fields(fieldKey, Scalar());
  const unsigned numFields = massDensity.numFields();

  // Intitialize the mass density.
  massDensity = 0.0;

  // Grab the state we need to do our job.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
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

  // Walk the fields.
  for (unsigned nodeListi = 0; nodeListi != numFields; ++nodeListi) {

    // We need a place to keep the volume as we accumulate it.
    const NodeList<Dimension>& nodeList = massDensity[nodeListi]->nodeList();
    Field<Dimension, Scalar> volEff("effective volume", nodeList);

    // Walk the local nodes.
    typename ConnectivityMap<Dimension>::const_iterator iItr;
    typename vector<int>::const_iterator jItr;
    unsigned i, j;
    Scalar Hdeti, Hdetj, etai, etaj, Wi, Wj, W0;
    Vector xij;
    const unsigned firstGhostNodei = nodeList.firstGhostNode();
    for (iItr = cm.begin(nodeListi); iItr != cm.end(nodeListi); ++iItr) {
      i = *iItr;
      const Scalar& Vi = volume(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Vector& xi = pos(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      Hdeti = Hi.Determinant();
      Scalar& rhoi = massDensity(nodeListi, i);
      Scalar& Veffi = volEff(i);

      // Walk the neighbors of this node.
      const vector<int>& connectivity = cm.connectivityForNode(nodeListi, i)[nodeListi];
      for (jItr = connectivity.begin(); jItr != connectivity.end(); ++jItr) {
        j = *jItr;
        if (cm.calculatePairInteraction(nodeListi, i, 
                                        nodeListi, j,
                                        firstGhostNodei)) {
          const Scalar& Vj = volume(nodeListi, j);
          const Scalar& mj = mass(nodeListi, j);
          const Vector& xj = pos(nodeListi, j);
          const SymTensor& Hj = H(nodeListi, j);
          Hdetj = Hj.Determinant();
          Scalar& rhoj = massDensity(nodeListi, j);
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

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::massDensity and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the matching derivative field from the StateDerivatives.
  KeyType incrementKey = IncrementFieldList<Dimension, Scalar>::prefix() + fieldKey;
  FieldList<Dimension, Scalar> f = state.fields(fieldKey, 0.0);
  const FieldList<Dimension, Scalar> df = derivs.fields(incrementKey, 0.0);

  // Loop over the internal values of the field.
  const unsigned numFields = f.numFields();
  for (unsigned i = 0; i != numFields; ++i) {
    for (unsigned j = 0; j != f[i]->numInternalElements(); ++j) {
      f(i,j) = max(mRhoMin, min(mRhoMax, f(i,j) + multiplier*(df(i,j))));
    }
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

