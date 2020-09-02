//---------------------------------Spheral++----------------------------------//
// An experimental hour glass control algorithm based on an estimate of the
// local third moment of the node distribution.
//
// Created by JMO, Thu Apr  2 09:02:00 PDT 2009
//----------------------------------------------------------------------------//
#include "ThirdMomentHourglassControl.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/FastMath.hh"
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

using namespace FastMath;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ThirdMomentHourglassControl<Dimension>::
ThirdMomentHourglassControl(const DataBase<Dimension>& dataBase,
                            const TableKernel<Dimension>& W,
                            const double multiplier,
                            const double maxAccelerationFactor):
  Physics<Dimension>(),
  mW(W),
  mMultiplier(multiplier),
  mMaxAccelerationFactor(maxAccelerationFactor),
  mThirdMoment(dataBase.newFluidFieldList(ThirdRankTensor(), "Third moment")) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ThirdMomentHourglassControl<Dimension>::
~ThirdMomentHourglassControl() {
}

//------------------------------------------------------------------------------
// Determine the principle derivatives for the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ThirdMomentHourglassControl<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  const double tiny = 1.0e-30;

  // Get the state fields.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, Scalar());
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity,  Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, Scalar());
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, Scalar());
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);

  // Derivative fields.
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, Scalar());
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());

  // Get the connectivity map.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Compute the magnitudes of the accelerations before we tweak them.
  // Note we also need the ghost node values, though we can defer finalizing
  // the boundary conditions here 'cause that will be done after computing
  // the third moment.
  FieldList<Dimension, Scalar> DvDtmag = dataBase.newFluidFieldList(Scalar(0.0), "starting acceleration magnitude");
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    for (size_t i = 0; i != nodeLists[nodeListi]->numInternalNodes(); ++i) {
      DvDtmag(nodeListi, i) = DvDt(nodeListi, i).magnitude();
    }
  }
  for (typename Physics<Dimension>::ConstBoundaryIterator itr = this->boundaryBegin();
       itr != this->boundaryEnd();
       ++itr) (*itr)->applyFieldListGhostBoundary(DvDtmag);

  // Compute the third moment of the node distribution.
  mThirdMoment.Zero();
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
    const Scalar W0 = mW.kernelValue(0.0, 1.0);
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {

      // State of node I.
      const int i = *iItr;
      const Vector& ri = position(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Iterate over the neighboring NodeLists.
      for (auto nodeListj = 0u; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors in this NodeList.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            CHECK(j < nodeLists[nodeListj]->numNodes());

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {

              // State for node J.
              const Vector& rj = position(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);

              // Kernel weighting and gradient.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Scalar Wi = mW.kernelValue(etai.magnitude(), 1.0)/W0;
              const Scalar Wj = mW.kernelValue(etaj.magnitude(), 1.0)/W0;

              //const Scalar rij2 = rij.magnitude2();
              //const Scalar hi2 = rij2/(etai.magnitude2() + tiny);
              //const Scalar hj2 = rij2/(etaj.magnitude2() + tiny);

              // Pair-wise contribution.
              const Vector rijUnit = rij.unitVector();
              const ThirdRankTensor thpt = outerProduct<Dimension>(rijUnit, outerProduct<Dimension>(rijUnit, rijUnit));
              mThirdMoment(nodeListi, i) += pow3(Wi)*thpt;
              mThirdMoment(nodeListj, j) -= pow3(Wj)*thpt;
            }
          }
        }
      }
    }
  }

  // Apply boundary conditions to the third moment.
  for (typename Physics<Dimension>::ConstBoundaryIterator itr = this->boundaryBegin();
       itr != this->boundaryEnd();
       ++itr) (*itr)->applyFieldListGhostBoundary(mThirdMoment);
  for (typename Physics<Dimension>::ConstBoundaryIterator itr = this->boundaryBegin();
       itr != this->boundaryEnd();
       ++itr) (*itr)->finalizeGhostBoundary();

  // Prepare an empty FieldList to tell us which node pair in the pair accelerations
  // we're incrementing as we loop over the nodes.
  FieldList<Dimension, int> pairAccelerationsOffset = dataBase.newFluidFieldList(0, "offset");

  const bool compatibleEnergyEvolution = derivatives.registered(HydroFieldNames::pairAccelerations);   // BLAGO!!

  // Now generate the corrective accelerations based on the third moment.
  for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
    const int firstGhostNodei = nodeLists[nodeListi]->firstGhostNode();
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {

      // State of node I.
      const int i = *iItr;
      const Scalar& mi = mass(nodeListi, i);
      const Vector& ri = position(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const ThirdRankTensor& Ti = mThirdMoment(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      CHECK(rhoi > 0.0);
      CHECK(fullConnectivity.size() == numNodeLists);

      // Iterate over the neighboring NodeLists.
      for (auto nodeListj = 0u; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Iterate over the neighbors in this NodeList.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            CHECK(j < nodeLists[nodeListj]->numNodes());

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {

              // State for node J.
              const Scalar& mj = mass(nodeListj, j);
              const Vector& rj = position(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar& rhoj = massDensity(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const ThirdRankTensor& Tj = mThirdMoment(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              CHECK(rhoj > 0.0);
  
              // Kernel weighting and gradient.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Vector gradWi = Hi*etai.unitVector()*mW.gradValue(etai.magnitude(), Hdeti);
              const Vector gradWj = Hj*etaj.unitVector()*mW.gradValue(etaj.magnitude(), Hdetj);
              const Vector gradWij = 0.5*(gradWi + gradWj);

              // Smoothing scales and other symmetrized properties.
              const Scalar rij2 = rij.magnitude2();
              const Scalar hi2 = rij2/(etai.magnitude2() + tiny);
              const Scalar hj2 = rij2/(etaj.magnitude2() + tiny);
              const Scalar hij2 = 0.5*(hi2 + hj2);
              CHECK(hij2 > 0.0);

              // Compute the grad^2 T term.
              const ThirdRankTensor Tji = Tj - Ti;
              const Scalar safety = 0.01*hij2;
              CHECK(safety > 0.0);
              Vector thpt;
              for (size_t g = 0; g != Dimension::nDim; ++g) {
                for (size_t b = 0; b != Dimension::nDim; ++b) {
                  for (size_t a = 0; a != Dimension::nDim; ++a) {
                    thpt(a) += Tji(g,b,a)*rij(g)/(pow2(rij(g)) + safety)*gradWij(b);
                  }
                }
              }
              const Scalar rhoij = 0.5*(rhoi + rhoj);
              thpt *= 0.5*mMultiplier*hij2*rhoij*(1.0/(rhoi*rhoi) + 1.0/(rhoj*rhoj))*(DvDtmag(nodeListi, i) + DvDtmag(nodeListj, j));
              const Vector DvDtij =  mj*thpt;
              const Vector DvDtji = -mi*thpt;

              // Increment the accelerations.
              DvDt(nodeListi, i) += DvDtij;
              DvDt(nodeListj, j) += DvDtji;

              // Increment the work.
              const Vector vij = vi - vj;
              DepsDt(nodeListi, i) += 0.5*(vij.dot(DvDtij));
              DepsDt(nodeListj, j) -= 0.5*(vij.dot(DvDtji));

              // In compatible energy mode, we need to increment the pair-wise
              // accelerations.
              if (compatibleEnergyEvolution) {
                if (!(i >= firstGhostNodei or pairAccelerationsOffset(nodeListi, i) < (int)pairAccelerationsi.size())) 
                  cerr << i << " "
                       << firstGhostNodei << " "
                       << pairAccelerationsOffset(nodeListi, i) << " "
                       << pairAccelerationsi.size() << " "
                       << endl;
                CHECK(i >= firstGhostNodei or pairAccelerationsOffset(nodeListi, i) < pairAccelerationsi.size());
                CHECK(j >= firstGhostNodej or pairAccelerationsOffset(nodeListj, j) < pairAccelerationsj.size());
                if (i < firstGhostNodei) pairAccelerationsi[pairAccelerationsOffset(nodeListi, i)] += DvDtij;
                if (j < firstGhostNodej) pairAccelerationsj[pairAccelerationsOffset(nodeListj, j)] += DvDtji;
              }
              if (i < firstGhostNodei) ++pairAccelerationsOffset(nodeListi, i);
              if (j < firstGhostNodej) ++pairAccelerationsOffset(nodeListj, j);
            }
          }
        }
      }
    }
  }

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    int nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      if (compatibleEnergyEvolution) {
        for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
          if (!(pairAccelerations(nodeListi, i).size() == pairAccelerationsOffset(nodeListi, i))) {
            cerr << nodeListi << " "
                 << i << " "
                 << pairAccelerations(nodeListi, i).size() << " "
                 << pairAccelerationsOffset(nodeListi, i) << endl;
          }
          ENSURE(pairAccelerations(nodeListi, i).size() == pairAccelerationsOffset(nodeListi, i));
        }
      }
    }
  }
  END_CONTRACT_SCOPE

}

//------------------------------------------------------------------------------
// Calculate the timestep constraint.
//------------------------------------------------------------------------------
template<typename Dimension>
typename ThirdMomentHourglassControl<Dimension>::TimeStepType
ThirdMomentHourglassControl<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const typename Dimension::Scalar /*currentTime*/) const {
  return TimeStepType(FLT_MAX, "No vote.");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ThirdMomentHourglassControl<Dimension>::
registerState(DataBase<Dimension>& /*dataBase*/,
              State<Dimension>& /*state*/) {
//   REQUIRE(mThirdMoment.numFields() == dataBase.numFluidNodeLists);
//   for (size_t i = 0; i != dataBase.numFluidNodeLists(); ++i) {
//     derivs.registerField(*mThirdMoment[i]);
//   }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ThirdMomentHourglassControl<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& /*derivs*/) {
}

}

