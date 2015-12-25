//---------------------------------Spheral++----------------------------------//
// PSPHHydroBase -- The PSPH/APSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Wed Dec 16 20:52:02 PST 2015
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "PSPHHydroBase.hh"
#include "computeSPHSumMassDensity.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computePSPHCorrections.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/GammaPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "FileIO/FileIO.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"

namespace Spheral {
namespace SPHSpace {

using namespace std;
using PhysicsSpace::GenericHydro;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using MeshSpace::Mesh;

using PhysicsSpace::MassDensityType;
using PhysicsSpace::HEvolutionType;

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
PSPHHydroBase<Dimension>::
PSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
              ArtificialViscosity<Dimension>& Q,
              const TableKernel<Dimension>& W,
              const TableKernel<Dimension>& WPi,
              const double filter,
              const double cfl,
              const bool useVelocityMagnitudeForDt,
              const bool compatibleEnergyEvolution,
              const bool gradhCorrection,
              const bool XSPH,
              const bool HopkinsConductivity,
              const bool sumMassDensityOverAllNodeLists,
              const MassDensityType densityUpdate,
              const HEvolutionType HUpdate,
              const double epsTensile,
              const double nTensile,
              const Vector& xmin,
              const Vector& xmax):
  SPHHydroBase<Dimension>(smoothingScaleMethod,
                          Q,
                          W,
                          WPi,
                          filter,
                          cfl,
                          useVelocityMagnitudeForDt,
                          compatibleEnergyEvolution,
                          gradhCorrection,
                          XSPH,
                          false,
                          sumMassDensityOverAllNodeLists,
                          densityUpdate,
                          HUpdate,
                          epsTensile,
                          nTensile,
                          xmin,
                          xmax),
  mHopkinsConductivity(HopkinsConductivity),
  mGamma(FieldSpace::Copy),
  mPSPHpbar(FieldSpace::Copy),
  mPSPHcorrection(FieldSpace::Copy) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PSPHHydroBase<Dimension>::
~PSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // SPH does most of it.
  SPHHydroBase<Dimension>::initializeProblemStartup(dataBase);
  dataBase.fluidGamma(mGamma);

  // Create storage for our internal state.
  mGamma = dataBase.newFluidFieldList(0.0, HydroFieldNames::gamma);
  mPSPHpbar = dataBase.newFluidFieldList(0.0, HydroFieldNames::PSPHpbar);
  mPSPHcorrection = dataBase.newFluidFieldList(0.0, HydroFieldNames::PSPHcorrection);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // SPH does most of it.
  SPHHydroBase<Dimension>::registerState(dataBase, state);

  // We also require the fluid gamma.
  PolicyPointer gammaPolicy(new GammaPolicy<Dimension>());
  state.enroll(mGamma, gammaPolicy);

  // Register the PSPH correction terms
  // We deliberately make this non-dynamic here.  These corrections are computed
  // during PSPHHydroBase::initialize, not as part of our usual state update.
  state.enroll(mPSPHpbar);
  state.enroll(mPSPHcorrection);
}

//------------------------------------------------------------------------------
// Initialize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // SPH does most of it.
  SPHHydroBase<Dimension>::initialize(time, dt, dataBase, state, derivs);

  // Do the PSPH corrections.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> PSPHpbar = state.fields(HydroFieldNames::PSPHpbar, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  computePSPHCorrections(connectivityMap, W, mass, position, specificThermalEnergy, gamma, H, PSPHpbar, PSPHcorrection);

  // Replace the pressure in the state with the PSPH sum definition.
  FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);
  P.assignFields(PSPHpbar);

  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(PSPHpbar);
    (*boundItr)->applyFieldListGhostBoundary(PSPHcorrection);
    (*boundItr)->applyFieldListGhostBoundary(P);
  }
  // We depend on the caller knowing to finalize the ghost boundaries!
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const double tiny = 1.0e-10;
  const Scalar W0 = W(0.0, 1.0);

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  const FieldList<Dimension, Scalar> PSPHpbar = state.fields(HydroFieldNames::PSPHpbar, 0.0);
  const FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  const FieldList<Dimension, Scalar> reducingViscosityMultiplierQ = state.fields(HydroFieldNames::reducingViscosityMultiplierQ, 0.0);
  const FieldList<Dimension, Scalar> reducingViscosityMultiplierL = state.fields(HydroFieldNames::reducingViscosityMultiplierL, 0.0);

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(gamma.size() == numNodeLists);
  CHECK(PSPHpbar.size() == numNodeLists);
  CHECK(PSPHcorrection.size() == numNodeLists);
  CHECK((not mHopkinsConductivity) or (reducingViscosityMultiplierQ.size() == numNodeLists));
  CHECK((not mHopkinsConductivity) or (reducingViscosityMultiplierL.size() == numNodeLists));

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> M = derivatives.fields(HydroFieldNames::M_CRKSPH, Tensor::zero);
  FieldList<Dimension, Tensor> localM = derivatives.fields("local " + HydroFieldNames::M_CRKSPH, Tensor::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, Scalar> effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  FieldList<Dimension, Scalar> viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Scalar> XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (this->mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const NodeList<Dimension>& nodeList = **itr;
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar hmin = nodeList.hmin();
    const Scalar hmax = nodeList.hmax();
    const Scalar hminratio = nodeList.hminratio();
    const int maxNumNeighbors = nodeList.maxNumNeighbors();
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();

    // The scale for the tensile correction.
    const Scalar WnPerh = W(1.0/nPerh, 1.0);

    // Get the work field for this NodeList.
    Field<Dimension, Scalar>& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Scalar& epsi = specificThermalEnergy(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar& ci = soundSpeed(nodeListi, i);
      const Scalar& gammai = gamma(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      Scalar& rhoSumi = rhoSum(nodeListi, i);
      Scalar& normi = normalization(nodeListi, i);
      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      Tensor& Mi = M(nodeListi, i);
      Tensor& localMi = localM(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      Scalar& effViscousPressurei = effViscousPressure(nodeListi, i);
      Scalar& viscousWorki = viscousWork(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Scalar& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = 1.0; // (nodeListi == nodeListj ? 1.0 : 0.2);
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const Vector& rj = position(nodeListj, j);
              const Scalar& mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar& rhoj = massDensity(nodeListj, j);
              const Scalar& epsj = specificThermalEnergy(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar& cj = soundSpeed(nodeListj, j);
              const Scalar& gammaj = gamma(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);

              Scalar& rhoSumj = rhoSum(nodeListj, j);
              Scalar& normj = normalization(nodeListj, j);
              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Tensor& Mj = M(nodeListj, j);
              Tensor& localMj = localM(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              Scalar& effViscousPressurej = effViscousPressure(nodeListj, j);
              Scalar& viscousWorkj = viscousWork(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              Scalar& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
              Vector& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              Scalar& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              SymTensor& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Node displacement.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);

              // Symmetrized kernel weight and gradient.
              const Vector Hetai = Hi*etai.unitVector();
              const std::pair<double, double> WWi = W.kernelAndGradValue(etaMagi, Hdeti);
              const Scalar Wi = WWi.first;
              const Scalar gWi = WWi.second;
              const Vector gradWi = gWi*Hetai;
              const Vector gradWQi = WQ.gradValue(etaMagi, Hdeti) * Hetai;

              const Vector Hetaj = Hj*etaj.unitVector();
              const std::pair<double, double> WWj = W.kernelAndGradValue(etaMagj, Hdetj);
              const Scalar Wj = WWj.first;
              const Scalar gWj = WWj.second;
              const Vector gradWj = gWj*Hetaj;
              const Vector gradWQj = WQ.gradValue(etaMagj, Hdetj) * Hetaj;

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()*safeInv(rij2*FastMath::square(Dimension::pownu12(rij2)), tiny);
              weightedNeighborSumi += fweightij*std::abs(gWi);
              weightedNeighborSumj += fweightij*std::abs(gWj);
              massSecondMomenti += fweightij*gradWi.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWj.magnitude2()*thpt;

              // Contribution to the sum density.
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wi;
                rhoSumj += mi*Wj;
                normi += mi/rhoi*Wi;
                normj += mj/rhoj*Wj;
              }

              // Mass density evolution.
              const Vector vij = vi - vj;
              const double deltaDrhoDti = vij.dot(gradWi);
              const double deltaDrhoDtj = vij.dot(gradWj);
              DrhoDti += deltaDrhoDti;
              DrhoDtj += deltaDrhoDtj;

              // Compute the pair-wise artificial viscosity.
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Vector Qacci = 0.5*(QPiij.first *gradWQi);
              const Vector Qaccj = 0.5*(QPiij.second*gradWQj);
              // const Scalar workQi = 0.5*(QPiij.first *vij).dot(gradWQi);
              // const Scalar workQj = 0.5*(QPiij.second*vij).dot(gradWQj);
              const Scalar workQi = vij.dot(Qacci);
              const Scalar workQj = vij.dot(Qaccj);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += mj/rhoj * Qi * Wi;
              effViscousPressurej += mi/rhoi * Qj * Wj;
              viscousWorki += mj*workQi;
              viscousWorkj += mi*workQj;

              // Acceleration.
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const Scalar& Fcorri=PSPHcorrection(nodeListi, i);
              const Scalar& Fcorrj=PSPHcorrection(nodeListj, j);
              const Scalar& Pbari=PSPHpbar(nodeListi, i);
              const Scalar& Pbarj=PSPHpbar(nodeListj, j);
              const Scalar Fij=1.0-Fcorri*safeInv(mj*epsj, tiny);
              const Scalar Fji=1.0-Fcorrj*safeInv(mi*epsi, tiny);
              const double engCoef=(gammai-1)*(gammaj-1)*epsi*epsj;
              const Vector deltaDvDt = engCoef*(gradWi*Fij*safeInv(Pbari, tiny) + gradWj*Fji*safeInv(Pbarj, tiny)) + Qacci + Qaccj;

              DvDti -= mj*deltaDvDt;
              DvDtj += mi*deltaDvDt;

              // Specific thermal energy evolution.
              // DepsDti += mj*(engCoef*deltaDrhoDti*Fij/max(Pbari,tiny) + workQi);
              // DepsDtj += mi*(engCoef*deltaDrhoDtj*Fji/max(Pbarj,tiny) + workQj);
              DepsDti += mj*(engCoef*deltaDrhoDti*Fij*safeInv(Pbari, tiny) + workQi);
              DepsDtj += mi*(engCoef*deltaDrhoDtj*Fji*safeInv(Pbarj, tiny) + workQj);

              //ADD ARITIFICIAL CONDUCTIVITY IN HOPKINS 2014A
              if (mHopkinsConductivity) {
                const Scalar alph_c = 0.25;//Parameter = 0.25 in Hopkins 2014
                const Scalar Vs = ci+cj-3.0*vij.dot(rij.unitVector());
                const Scalar& Qalpha_i = reducingViscosityMultiplierL(nodeListi, i); //Both L and Q corrections are the same for Cullen Viscosity
                const Scalar& Qalpha_j = reducingViscosityMultiplierL(nodeListj, j); //Both L and Q corrections are the same for Cullen Viscosity
                //DepsDti += (Vs > 0.0)*alph_c*mi*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj+1e-30)*(rhoi+rhoj+1e-30));
                //DepsDtj += (Vs > 0.0)*alph_c*mi*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj+1e-30)*(rhoi+rhoj+1e-30));
                //DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pbari-Pbarj)*((gradWi).dot(rij.unitVector()))/max((Pbari+Pbarj)*0.5*(rhoi+rhoj),tiny) : 0.0;
                //DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pbari-Pbarj)*((gradWj).dot(rij.unitVector()))/max((Pbari+Pbarj)*0.5*(rhoi+rhoj),tiny) : 0.0;
                DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj),tiny) : 0.0;
                DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj),tiny) : 0.0;
                //const Scalar tmpi = (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj)) : 0.0;
                //const Scalar tmpj = (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj)) : 0.0;
                //DepsDti += tmpi;
                //DepsDtj += tmpj;
                //DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj)) : 0.0;
                //DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj)) : 0.0;
                //DepsDti += (Vs > 0.0 && Pbari > 1e-4 && Pbarj > 1e-4) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj),tiny) : 0.0;
                //DepsDtj += (Vs > 0.0 && Pbari > 1e-4 && Pbarj > 1e-4) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pbari-Pbarj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pbari+Pbarj)*(rhoi+rhoj),tiny) : 0.0;
              }

              if (this->mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(-mj*deltaDvDt);
                pairAccelerationsj.push_back( mi*deltaDvDt);
              }

              // Velocity gradient.
              const Tensor deltaDvDxi = vij.dyad(gradWi);
              const Tensor deltaDvDxj = vij.dyad(gradWj);
              DvDxi -= mj*deltaDvDxi;
              DvDxj -= mi*deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= mj*deltaDvDxi;
                localDvDxj -= mi*deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (this->mXSPH and (nodeListi == nodeListj)) {
                const double fXSPH = max(0.0, min(1.0, abs(vij.dot(rij)*safeInv(vij.magnitude()*rij.magnitude()))));
                CHECK(fXSPH >= 0.0 and fXSPH <= 1.0);
                XSPHWeightSumi += fXSPH*mj/rhoj*Wi;
                XSPHWeightSumj += fXSPH*mi/rhoi*Wj;
                XSPHDeltaVi -= fXSPH*mj/rhoj*Wi*vij;
                XSPHDeltaVj += fXSPH*mi/rhoi*Wj*vij;
              }

              // Linear gradient correction term.
              Mi -= mj/rhoj*rij.dyad(gradWi);
              Mj -= mi/rhoi*rij.dyad(gradWj);
              if (nodeListi == nodeListj) {
                localMi -= mj/rhoj*rij.dyad(gradWi);
                localMj -= mi/rhoi*rij.dyad(gradWj);
              }
            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not this->mCompatibleEnergyEvolution or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;
      normi += mi/rhoi*W0*Hdeti;

      // Finish the continuity equation.
      DrhoDti *= mi;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      DvDxi /= rhoi;
      localDvDxi /= rhoi;
      if (this->correctVelocityGradient()) {
        DvDxi = Mi*DvDxi;
        localDvDxi = localMi*DvDxi;
      }

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (this->mXSPH) {
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti = vi + XSPHDeltaVi*safeInvVar(XSPHWeightSumi, tiny);
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = this->mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                                   ri,
                                                                   DvDxi,
                                                                   hmin,
                                                                   hmax,
                                                                   hminratio,
                                                                   nPerh);
      Hideali = this->mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                              ri,
                                                              weightedNeighborSumi,
                                                              massSecondMomenti,
                                                              W,
                                                              hmin,
                                                              hmax,
                                                              hminratio,
                                                              nPerh,
                                                              connectivityMap,
                                                              nodeListi,
                                                              i);

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());

      // Now add the pairwise time for each neighbor we computed here.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
          Field<Dimension, Scalar>& workFieldj = nodeLists[nodeListj]->work();
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              workFieldj(j) += deltaTimePair;
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (this->compatibleEnergyEvolution()) {

    FieldList<Dimension, Vector> accelerations = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // SPH does most of it.
  SPHHydroBase<Dimension>::applyGhostBoundaries(state, derivs);

  // Apply boundary conditions to the PSPH corrections.
  FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  FieldList<Dimension, Scalar> PSPHpbar = state.fields(HydroFieldNames::PSPHpbar, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(gamma);
    (*boundaryItr)->applyFieldListGhostBoundary(PSPHpbar);
    (*boundaryItr)->applyFieldListGhostBoundary(PSPHcorrection);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // SPH does most of it.
  SPHHydroBase<Dimension>::enforceBoundaries(state, derivs);

  // Enforce boundary conditions on the PSPH corrections.
  FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  FieldList<Dimension, Scalar> PSPHpbar = state.fields(HydroFieldNames::PSPHpbar, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(gamma);
    (*boundaryItr)->enforceFieldListBoundary(PSPHpbar);
    (*boundaryItr)->enforceFieldListBoundary(PSPHcorrection);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
dumpState(FileIO& file, string pathName) const {

  // SPH does most of it.
  SPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mGamma, pathName + "/gamma");
  file.write(mPSPHpbar, pathName + "/PSPHpbar");
  file.write(mPSPHcorrection, pathName + "/PSPHcorrection");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
restoreState(const FileIO& file, string pathName) {
 
  // SPH does most of it.
  SPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mGamma, pathName + "/gamma");
  file.read(mPSPHpbar, pathName + "/PSPHpbar");
  file.read(mPSPHcorrection, pathName + "/PSPHcorrection");
}

}
}

