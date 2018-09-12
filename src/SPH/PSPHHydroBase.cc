//---------------------------------Spheral++----------------------------------//
// PSPHHydroBase -- The PSPH/APSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Wed Dec 16 20:52:02 PST 2015
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
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
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"

#include "PSPHHydroBase.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
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
              const bool evolveTotalEnergy,
              const bool XSPH,
              const bool correctVelocityGradient,
              const bool HopkinsConductivity,
              const bool sumMassDensityOverAllNodeLists,
              const MassDensityType densityUpdate,
              const HEvolutionType HUpdate,
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
                          evolveTotalEnergy,
                          true,
                          XSPH,
                          correctVelocityGradient,
                          sumMassDensityOverAllNodeLists,
                          densityUpdate,
                          HUpdate,
                          0.0,
                          1.0,
                          xmin,
                          xmax),
  mHopkinsConductivity(HopkinsConductivity),
  mGamma(FieldStorageType::CopyFields),
  mPSPHcorrection(FieldStorageType::CopyFields) {
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

  // Create storage for our internal state.
  mGamma = dataBase.newFluidFieldList(0.0, HydroFieldNames::gamma);
  mPSPHcorrection = dataBase.newFluidFieldList(0.0, HydroFieldNames::PSPHcorrection);
  dataBase.fluidGamma(mGamma);
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

  // Make sure we're the right size.
  dataBase.resizeFluidFieldList(mGamma, 0.0, HydroFieldNames::gamma, false);
  dataBase.resizeFluidFieldList(mPSPHcorrection, 0.0, HydroFieldNames::PSPHcorrection, false);

  // SPH does most of it.
  SPHHydroBase<Dimension>::registerState(dataBase, state);

  // We also require the fluid gamma.
  PolicyPointer gammaPolicy(new GammaPolicy<Dimension>());
  state.enroll(mGamma, gammaPolicy);

  // Override the default policies for pressure and sound speed.  We'll compute those
  // specially in the postStateUpdate.  Same goes for registering the PSPHcorrections.
  state.enroll(mPSPHcorrection);
  state.removePolicy(this->mPressure);
  state.removePolicy(this->mSoundSpeed);
}

//------------------------------------------------------------------------------
// Pre-step initializations.  Since the topology has just been changed we need
// to recompute the PSPH corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Do the PSPH corrections.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  computePSPHCorrections(connectivityMap, W, mass, position, specificThermalEnergy, gamma, H, 
                         (this->mDensityUpdate != MassDensityType::IntegrateDensity),
                         rho, P, cs, PSPHcorrection);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(rho);
    (*boundItr)->applyFieldListGhostBoundary(P);
    (*boundItr)->applyFieldListGhostBoundary(cs);
    (*boundItr)->applyFieldListGhostBoundary(PSPHcorrection);
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
}

//------------------------------------------------------------------------------
// Post-state update.  This is where we can recompute the PSPH pressure and
// corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PSPHHydroBase<Dimension>::
postStateUpdate(const Scalar time, 
                const Scalar dt,
                const DataBase<Dimension>& dataBase, 
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivatives) {

  // First we need out boundary conditions completed, which the time integrator hasn't 
  // verified yet.
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  
  // Do the PSPH corrections.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, Scalar> gamma = state.fields(HydroFieldNames::gamma, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  computePSPHCorrections(connectivityMap, W, mass, position, specificThermalEnergy, gamma, H,
                         (this->mDensityUpdate != MassDensityType::IntegrateDensity),
                         rho, P, cs, PSPHcorrection);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(rho);
    (*boundItr)->applyFieldListGhostBoundary(P);
    (*boundItr)->applyFieldListGhostBoundary(cs);
    (*boundItr)->applyFieldListGhostBoundary(PSPHcorrection);
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
  const FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);
  const FieldList<Dimension, Scalar> reducingViscosityMultiplierQ = state.fields(HydroFieldNames::ArtificialViscousCqMultiplier, 0.0);
  const FieldList<Dimension, Scalar> reducingViscosityMultiplierL = state.fields(HydroFieldNames::ArtificialViscousClMultiplier, 0.0);

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(gamma.size() == numNodeLists);
  CHECK(PSPHcorrection.size() == numNodeLists);
  CHECK((not mHopkinsConductivity) or (reducingViscosityMultiplierQ.size() == numNodeLists));
  CHECK((not mHopkinsConductivity) or (reducingViscosityMultiplierL.size() == numNodeLists));

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  FieldList<Dimension, Tensor> localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
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
      const Scalar& Pi = pressure(nodeListi, i);
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
              const Scalar& Pj = pressure(nodeListj, j);
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
              const double fweightij = nodeListi == nodeListj ? 1.0 : mj*rhoi/(mi*rhoj);
              const double rij2 = rij.magnitude2();
              const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
              weightedNeighborSumi +=     fweightij*std::abs(gWi);
              weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
              massSecondMomenti +=     fweightij*gradWi.magnitude2()*thpt;
              massSecondMomentj += 1.0/fweightij*gradWj.magnitude2()*thpt;

              // Contribution to the sum density.
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wj;
                rhoSumj += mi*Wi;
                normi += mi/rhoi*Wj;
                normj += mj/rhoj*Wi;
              }

              // Compute the pair-wise artificial viscosity.
              const Vector vij = vi - vj;
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Vector Qacci = 0.5*rhoi*QPiij.first  * (gradWi + gradWj)/(rhoi + rhoj);
              const Vector Qaccj = 0.5*rhoj*QPiij.second * (gradWi + gradWj)/(rhoi + rhoj);
              // const Vector Qacci = 0.5*(QPiij.first *gradWQi);
              // const Vector Qaccj = 0.5*(QPiij.second*gradWQj);
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
              const Scalar Fij=1.0-Fcorri*safeInv(mj*epsj, tiny);
              const Scalar Fji=1.0-Fcorrj*safeInv(mi*epsi, tiny);
              const double engCoef=(gammai-1)*(gammaj-1)*epsi*epsj;
              const Vector deltaDvDt = engCoef*(gradWi*Fij*safeInv(Pi, tiny) + gradWj*Fji*safeInv(Pj, tiny)) + Qacci + Qaccj;

              DvDti -= mj*deltaDvDt;
              DvDtj += mi*deltaDvDt;

              // Specific thermal energy evolution.
              DepsDti += mj*(engCoef*Fij*safeInv(Pi, tiny)*vij.dot(gradWi) + workQi);
              DepsDtj += mi*(engCoef*Fji*safeInv(Pj, tiny)*vij.dot(gradWj) + workQj);
              if (this->mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(-mj*deltaDvDt);
                pairAccelerationsj.push_back( mi*deltaDvDt);
              }

              //ADD ARITIFICIAL CONDUCTIVITY IN HOPKINS 2014A
              if (mHopkinsConductivity) {
                const Scalar alph_c = 0.25;//Parameter = 0.25 in Hopkins 2014
                const Scalar Vs = ci+cj-3.0*vij.dot(rij.unitVector());
                const Scalar& Qalpha_i = reducingViscosityMultiplierL(nodeListi, i); //Both L and Q corrections are the same for Cullen Viscosity
                const Scalar& Qalpha_j = reducingViscosityMultiplierL(nodeListj, j); //Both L and Q corrections are the same for Cullen Viscosity
                //DepsDti += (Vs > 0.0)*alph_c*mi*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj+1e-30)*(rhoi+rhoj+1e-30));
                //DepsDtj += (Vs > 0.0)*alph_c*mi*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj+1e-30)*(rhoi+rhoj+1e-30));
                //DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi).dot(rij.unitVector()))/max((Pi+Pj)*0.5*(rhoi+rhoj),tiny) : 0.0;
                //DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWj).dot(rij.unitVector()))/max((Pi+Pj)*0.5*(rhoi+rhoj),tiny) : 0.0;
                DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
                DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
                //const Scalar tmpi = (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
                //const Scalar tmpj = (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
                //DepsDti += tmpi;
                //DepsDtj += tmpj;
                //DepsDti += (Vs > 0.0) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
                //DepsDtj += (Vs > 0.0) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj)) : 0.0;
                //DepsDti += (Vs > 0.0 && Pi > 1e-4 && Pj > 1e-4) ? alph_c*mj*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsi-epsj)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
                //DepsDtj += (Vs > 0.0 && Pi > 1e-4 && Pj > 1e-4) ? alph_c*mi*(Qalpha_i+Qalpha_j)*0.5*Vs*(epsj-epsi)*abs(Pi-Pj)*((gradWi+gradWj).dot(rij.unitVector()))*safeInv((Pi+Pj)*(rhoi+rhoj),tiny) : 0.0;
              }

              // Velocity gradient.
              const Tensor deltaDvDxi = mj*vij.dyad(gradWi);
              const Tensor deltaDvDxj = mi*vij.dyad(gradWj);
              DvDxi -= deltaDvDxi;
              DvDxj -= deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= deltaDvDxi;
                localDvDxj -= deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (this->mXSPH and (nodeListi == nodeListj)) {
                const double wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
                XSPHWeightSumi += wXSPHij;
                XSPHWeightSumj += wXSPHij;
                XSPHDeltaVi -= wXSPHij*vij;
                XSPHDeltaVj += wXSPHij*vij;
              }

              // Linear gradient correction term.
              Mi -= mj*rij.dyad(gradWi);
              Mj -= mi*rij.dyad(gradWj);
              if (nodeListi == nodeListj) {
                localMi -= mj*rij.dyad(gradWi);
                localMj -= mi*rij.dyad(gradWj);
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

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        localMi = localMi.Inverse();
        DvDxi = DvDxi*Mi;
        localDvDxi = localDvDxi*localMi;
      } else {
        DvDxi /= rhoi;
        localDvDxi /= rhoi;
      }

      // Evaluate the continuity equation.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (this->mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

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
  if (this->mCompatibleEnergyEvolution) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    auto DepsDt = derivs.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
      (*boundaryItr)->applyFieldListGhostBoundary(DepsDt);
    }
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
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(gamma);
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
  FieldList<Dimension, Scalar> PSPHcorrection = state.fields(HydroFieldNames::PSPHcorrection, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(gamma);
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
  file.read(mPSPHcorrection, pathName + "/PSPHcorrection");
}

}
