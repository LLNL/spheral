//---------------------------------Spheral++----------------------------------//
// Hydro -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "CRKSPHUtilities.hh"
#include "computeVoronoiVolume.hh"
#include "computeHullVolumes.hh"
#include "computeCRKSPHSumVolume.hh"
#include "computeHVolumes.hh"
#include "editMultimaterialSurfaceTopology.hh"
#include "zerothOrderSurfaceCorrections.hh"
#include "SurfaceNodeCoupling.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "SPH/correctSPHSumMassDensity.hh"
#include "computeCRKSPHSumMassDensity.hh"
#include "computeCRKSPHMoments.hh"
#include "detectSurface.hh"
#include "computeCRKSPHCorrections.hh"
#include "computeCRKSPHIntegral.hh"
#include "gradientCRKSPH.hh"
#include "centerOfMass.hh"
#include "volumeSpacing.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/entropyWeightingFunction.hh"
#include "Strength/SolidFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "ContinuityVolumePolicy.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"

#include "CRKSPHHydroBase.hh"

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

namespace {

double fluxlimiterVL(const double x) {
  // if (x > 0.0) {
  //   return min(1.0, x);                           // minmod
  // } else {
  //   return 0.0;
  // }

  // if (x > 0.0) {
  //   return 2.0/(1.0 + x)*2.0*x/(1.0 + x);                       // van Leer
  // } else {
  //   return 0.0;
  // }
  return (x + abs(x))/(1.0 + abs(x));                       // van Leer
}

}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHHydroBase<Dimension>::
CRKSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                ArtificialViscosity<Dimension>& Q,
                const TableKernel<Dimension>& W,
                const TableKernel<Dimension>& WPi,
                const double filter,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool evolveTotalEnergy,
                const bool XSPH,
                const MassDensityType densityUpdate,
                const HEvolutionType HUpdate,
                const CRKOrder correctionOrder,
                const CRKVolumeType volumeType,
                const double epsTensile,
                const double nTensile,
                const bool limitMultimaterialTopology):
  GenericHydro<Dimension>(W, WPi, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCorrectionOrder(correctionOrder),
  mVolumeType(volumeType),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mLimitMultimaterialTopology(limitMultimaterialTopology),
  mfilter(filter),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mSpecificThermalEnergy0(FieldStorageType::CopyFields),
  mEntropy(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
  mMaxViscousPressure(FieldStorageType::CopyFields),
  mEffViscousPressure(FieldStorageType::CopyFields),
  mViscousWork(FieldStorageType::CopyFields),
  mVolume(FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldStorageType::CopyFields),
  mMassSecondMoment(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mPairAccelerations(),
  mA(FieldStorageType::CopyFields),
  mB(FieldStorageType::CopyFields),
  mC(FieldStorageType::CopyFields),
  mGradA(FieldStorageType::CopyFields),
  mGradB(FieldStorageType::CopyFields),
  mGradC(FieldStorageType::CopyFields),
  mM0(FieldStorageType::CopyFields),
  mM1(FieldStorageType::CopyFields),
  mM2(FieldStorageType::CopyFields),
  mM3(FieldStorageType::CopyFields),
  mM4(FieldStorageType::CopyFields),
  mGradm0(FieldStorageType::CopyFields),
  mGradm1(FieldStorageType::CopyFields),
  mGradm2(FieldStorageType::CopyFields),
  mGradm3(FieldStorageType::CopyFields),
  mGradm4(FieldStorageType::CopyFields),
  mSurfacePoint(FieldStorageType::CopyFields),
  mEtaVoidPoints(FieldStorageType::CopyFields),
  mCells(FieldStorageType::CopyFields),
  mCellFaceFlags(FieldStorageType::CopyFields),
  mVoidBoundary(mSurfacePoint, mEtaVoidPoints),
  mRestart(registerWithRestart(*this)) {
  // this->appendBoundary(mVoidBoundary);  // Suspend actually building the void points.
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHHydroBase<Dimension>::
~CRKSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  mEntropy = dataBase.newFluidFieldList(0.0, HydroFieldNames::entropy);
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
  mEffViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::effectiveViscousPressure);
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mViscousWork = dataBase.newFluidFieldList(0.0, HydroFieldNames::viscousWork);
  mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mPairAccelerations.clear();
  mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero, "delta centroid");

  mA = dataBase.newFluidFieldList(0.0,                        HydroFieldNames::A_CRKSPH);
  mB = dataBase.newFluidFieldList(Vector::zero,               HydroFieldNames::B_CRKSPH);
  mGradA = dataBase.newFluidFieldList(Vector::zero,           HydroFieldNames::gradA_CRKSPH);
  mGradB = dataBase.newFluidFieldList(Tensor::zero,           HydroFieldNames::gradB_CRKSPH);
  mM0 = dataBase.newFluidFieldList(0.0,                       HydroFieldNames::m0_CRKSPH);
  mM1 = dataBase.newFluidFieldList(Vector::zero,              HydroFieldNames::m1_CRKSPH);
  mM2 = dataBase.newFluidFieldList(SymTensor::zero,           HydroFieldNames::m2_CRKSPH);
  mGradm0 = dataBase.newFluidFieldList(Vector::zero,          HydroFieldNames::gradM0_CRKSPH);
  mGradm1 = dataBase.newFluidFieldList(Tensor::zero,          HydroFieldNames::gradM1_CRKSPH);
  mGradm2 = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::gradM2_CRKSPH);
  if (mCorrectionOrder == CRKOrder::QuadraticOrder) {
    mC = dataBase.newFluidFieldList(Tensor::zero,                HydroFieldNames::C_CRKSPH);
    mGradC = dataBase.newFluidFieldList(ThirdRankTensor::zero,   HydroFieldNames::gradC_CRKSPH);
    mM3 = dataBase.newFluidFieldList(ThirdRankTensor::zero,      HydroFieldNames::m3_CRKSPH);
    mM4 = dataBase.newFluidFieldList(FourthRankTensor::zero,     HydroFieldNames::m4_CRKSPH);
    mGradm3 = dataBase.newFluidFieldList(FourthRankTensor::zero, HydroFieldNames::gradM3_CRKSPH);
    mGradm4 = dataBase.newFluidFieldList(FifthRankTensor::zero,  HydroFieldNames::gradM4_CRKSPH);
  }

  // We need volumes in order to prepare the surface detection.
  mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::etaVoidPoints);
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    mCells = dataBase.newFluidFieldList(FacetedVolume(), HydroFieldNames::cells);
    mCellFaceFlags = dataBase.newFluidFieldList(vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags);
  }
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();

  // Compute the volumes for real.
  if (mVolumeType == CRKVolumeType::CRKMassOverDensity) {
    mVolume.assignFields(mass/massDensity);
  } else if (mVolumeType == CRKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, W, position, mass, H, mVolume);
  } else if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    mVolume.assignFields(mass/massDensity);
    const FieldList<Dimension, typename Dimension::SymTensor> damage = dataBase.solidEffectiveDamage();
    computeVoronoiVolume(position, H, connectivityMap, damage,
                         vector<typename Dimension::FacetedVolume>(),               // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),      // no holes
                         vector<Boundary<Dimension>*>(this->boundaryBegin(),        // boundaries
                                                      this->boundaryEnd()),
                         FieldList<Dimension, typename Dimension::Scalar>(),        // no weights
                         mSurfacePoint, mVolume, mDeltaCentroid, mEtaVoidPoints,    // return values
                         mCells,                                                    // return cells
                         mCellFaceFlags);                                           // node cell multimaterial faces
    if (mLimitMultimaterialTopology) {
      for (ConstBoundaryIterator boundItr = this->boundaryBegin();
           boundItr != this->boundaryEnd();
           ++boundItr) (*boundItr)->applyFieldListGhostBoundary(mSurfacePoint);
      for (ConstBoundaryIterator boundItr = this->boundaryBegin();
           boundItr != this->boundaryEnd();
           ++boundItr) (*boundItr)->finalizeGhostBoundary();
      // editMultimaterialSurfaceTopology(mSurfacePoint,
      //                                  const_cast<ConnectivityMap<Dimension>&>(connectivityMap));
    }
  } else if (mVolumeType == CRKVolumeType::CRKHullVolume) {
    computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, mVolume);
  } else if (mVolumeType == CRKVolumeType::HVolume) {
    const Scalar nPerh = mVolume.nodeListPtrs()[0]->nodesPerSmoothingScale();
    computeHVolumes(nPerh, H, mVolume);
  } else {
    VERIFY2(false, "Unknown CRK volume weighting.");
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mVolume);
    if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(mSurfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
  //   // flagSurfaceNeighbors(mSurfacePoint, connectivityMap);
  //   // mVolume = computeShepardsInterpolation(mVolume,
  //   //                                        connectivityMap,
  //   //                                        W,
  //   //                                        position,
  //   //                                        H,
  //   //                                        mVolume);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->applyFieldListGhostBoundary(mVolume);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // }

  // Compute the corrections.
  const NodeCoupling couple;
  computeCRKSPHMoments(connectivityMap, W, mVolume, position, H, correctionOrder(), couple, mM0, mM1, mM2, mM3, mM4, mGradm0, mGradm1, mGradm2, mGradm3, mGradm4);
  computeCRKSPHCorrections(mM0, mM1, mM2, mM3, mM4, mGradm0, mGradm1, mGradm2, mGradm3, mGradm4, H, mSurfacePoint, correctionOrder(), mA, mB, mC, mGradA, mGradB, mGradC);
  if (mLimitMultimaterialTopology) zerothOrderSurfaceCorrections(mA, mB, mC, mGradA, mGradB, mGradC, mM0, mGradm0, mSurfacePoint);

  // This breaks domain independence, so we'll try being inconsistent on the first step.
  // // We need to initialize the velocity gradient if we're using the CRKSPH artificial viscosity.
  // const FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  // mDvDx.assignFields(gradientCRKSPH(velocity, position, mVolume, H, mA, mB, mC, mGradA, mGradB, mGradC, connectivityMap, correctionOrder(), W, NodeCoupling()));

  // Initialize the pressure, sound speed, and entropy.
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
  dataBase.fluidEntropy(mEntropy);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Create the local storage for time step mask, pressure, sound speed, and correction fields.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  // dataBase.fluidPressure(mPressure);
  // dataBase.fluidSoundSpeed(mSoundSpeed);
  dataBase.resizeFluidFieldList(mEntropy,    0.0,                   HydroFieldNames::entropy, false);
  dataBase.resizeFluidFieldList(mPressure,   0.0,                   HydroFieldNames::pressure, false);
  dataBase.resizeFluidFieldList(mSoundSpeed, 0.0,                   HydroFieldNames::soundSpeed, false);
  dataBase.resizeFluidFieldList(mVolume,     0.0,                   HydroFieldNames::volume, false);
  dataBase.resizeFluidFieldList(mA,          0.0,                   HydroFieldNames::A_CRKSPH, false);
  dataBase.resizeFluidFieldList(mB,          Vector::zero,          HydroFieldNames::B_CRKSPH, false);
  dataBase.resizeFluidFieldList(mGradA,      Vector::zero,          HydroFieldNames::gradA_CRKSPH, false);
  dataBase.resizeFluidFieldList(mGradB,      Tensor::zero,          HydroFieldNames::gradB_CRKSPH, false);
  dataBase.resizeFluidFieldList(mM0,         0.0,                   HydroFieldNames::m0_CRKSPH, false);
  dataBase.resizeFluidFieldList(mM1,         Vector::zero,          HydroFieldNames::m1_CRKSPH, false);
  dataBase.resizeFluidFieldList(mM2,         SymTensor::zero,       HydroFieldNames::m2_CRKSPH, false);
  dataBase.resizeFluidFieldList(mGradm0,     Vector::zero,          HydroFieldNames::gradM0_CRKSPH, false);
  dataBase.resizeFluidFieldList(mGradm1,     Tensor::zero,          HydroFieldNames::gradM1_CRKSPH, false);
  dataBase.resizeFluidFieldList(mGradm2,     ThirdRankTensor::zero, HydroFieldNames::gradM2_CRKSPH, false);
  if (mCorrectionOrder == CRKOrder::QuadraticOrder) {
    dataBase.resizeFluidFieldList(mC,        Tensor::zero,          HydroFieldNames::C_CRKSPH, false);
    dataBase.resizeFluidFieldList(mGradC,    ThirdRankTensor::zero, HydroFieldNames::gradC_CRKSPH, false);
    dataBase.resizeFluidFieldList(mM3,       ThirdRankTensor::zero, HydroFieldNames::m3_CRKSPH, false);
    dataBase.resizeFluidFieldList(mM4,       FourthRankTensor::zero,HydroFieldNames::m4_CRKSPH, false);
    dataBase.resizeFluidFieldList(mGradm3,   FourthRankTensor::zero,HydroFieldNames::gradM3_CRKSPH, false);
    dataBase.resizeFluidFieldList(mGradm4,   FifthRankTensor::zero, HydroFieldNames::gradM4_CRKSPH, false);
  }
  dataBase.resizeFluidFieldList(mSurfacePoint, 0, HydroFieldNames::surfacePoint, false);
  dataBase.resizeFluidFieldList(mEtaVoidPoints, vector<Vector>(), HydroFieldNames::etaVoidPoints, false);
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    dataBase.resizeFluidFieldList(mCells, FacetedVolume(), HydroFieldNames::cells, false);
    dataBase.resizeFluidFieldList(mCellFaceFlags, vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags, false);
  }

  // We have to choose either compatible or total energy evolution.
  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "CRKSPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // If we're using the compatibile energy discretization, prepare to maintain a copy
  // of the thermal energy.
  dataBase.resizeFluidFieldList(mSpecificThermalEnergy0, 0.0);
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      *mSpecificThermalEnergy0[nodeListi] = (*itr)->specificThermalEnergy();
      (*mSpecificThermalEnergy0[nodeListi]).name(HydroFieldNames::specificThermalEnergy + "0");
    }
  }

  // Now register away.
  // Mass.
  FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  state.enroll(mass);

  // Volume.
  PolicyPointer volumePolicy(new ContinuityVolumePolicy<Dimension>());
  state.enroll(mVolume, volumePolicy);

  // We need to build up CompositeFieldListPolicies for the mass density and H fields
  // in order to enforce NodeList dependent limits.
  FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  std::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > rhoPolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  std::shared_ptr<CompositeFieldListPolicy<Dimension, SymTensor> > Hpolicy(new CompositeFieldListPolicy<Dimension, SymTensor>());
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    rhoPolicy->push_back(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                      (*itr)->rhoMax()));
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == HEvolutionType::IntegrateH) {
      Hpolicy->push_back(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    } else {
      CHECK(HEvolution() == HEvolutionType::IdealH);
      Hpolicy->push_back(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    }
  }
  state.enroll(massDensity, rhoPolicy);
  state.enroll(Hfield, Hpolicy);

  // Register the position update, which depends on whether we're using XSPH or not.
  FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  if (true) { // (mXSPH) {
    PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
    state.enroll(position, positionPolicy);
  } else {
    PolicyPointer positionPolicy(new PositionPolicy<Dimension>());
    state.enroll(position, positionPolicy);
  }

  // Register the entropy.
  PolicyPointer entropyPolicy(new EntropyPolicy<Dimension>());
  state.enroll(mEntropy, entropyPolicy);

  // Are we using the compatible energy evolution scheme?
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  if (compatibleEnergyEvolution()) {
    // The compatible energy update.
    PolicyPointer thermalEnergyPolicy(new SpecificThermalEnergyPolicy<Dimension>(dataBase));
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::specificThermalEnergy,
                                                                           true));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
    state.enroll(mSpecificThermalEnergy0);

  } else if (mEvolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    PolicyPointer thermalEnergyPolicy(new SpecificFromTotalThermalEnergyPolicy<Dimension>());
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::specificThermalEnergy,
                                                                           true));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);

  } else {
    // Otherwise we're just time-evolving the specific energy.
    PolicyPointer thermalEnergyPolicy(new IncrementFieldList<Dimension, Scalar>());
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(true));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
  }

  // Register the time step mask, initialized to 1 so that everything defaults to being
  // checked.
  state.enroll(mTimeStepMask);

  // Compute and register the pressure and sound speed.
  PolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
  PolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mSoundSpeed, csPolicy);

  // Register the CRKSPH correction fields.
  // We deliberately make these non-dynamic here.  This corrections are computed
  // during CRKSPHHydroBase::initialize, not as part of our usual state update.
  // This is necessary 'cause we need boundary conditions *and* the current set of
  // neighbors before we compute these suckers.
  state.enroll(mA);
  state.enroll(mB);
  state.enroll(mC);
  state.enroll(mGradA);
  state.enroll(mGradB);
  state.enroll(mGradC);
  state.enroll(mM0);
  state.enroll(mM1);
  state.enroll(mM2);
  state.enroll(mM3);
  state.enroll(mM4);
  state.enroll(mGradm0);
  state.enroll(mGradm1);
  state.enroll(mGradm2);
  state.enroll(mGradm3);
  state.enroll(mGradm4);

  // Multimaterial information
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    state.enroll(mCells);
    state.enroll(mCellFaceFlags);
  }
  state.enroll(mSurfacePoint);

  // We also register the delta centroid for visualiation purposes.
  if (mfilter > 0.0 and mVolumeType == CRKVolumeType::CRKVoronoiVolume) state.enroll(mDeltaCentroid);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  typedef typename StateDerivatives<Dimension>::KeyType Key;
  const string DxDtName = IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const string DvDtName = HydroFieldNames::hydroAcceleration;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  dataBase.resizeFluidFieldList(mEffViscousPressure, 0.0, HydroFieldNames::effectiveViscousPressure, false);
  dataBase.resizeFluidFieldList(mViscousWork, 0.0, HydroFieldNames::viscousWork, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);

  derivs.enroll(mHideal);
  derivs.enroll(mMaxViscousPressure);
  derivs.enroll(mEffViscousPressure);
  derivs.enroll(mViscousWork);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mXSPHDeltaV);

  // These two (the position and velocity updates) may be registered
  // by other physics packages as well, so we need to be careful
  // not to duplicate if so.
  if (not derivs.registered(mDxDt)) derivs.enroll(mDxDt);
  if (not derivs.registered(mDvDt)) derivs.enroll(mDvDt);

  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDHDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enrollAny(HydroFieldNames::pairAccelerations, mPairAccelerations);
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Initialize the volume.
  const auto& W = this->kernel();
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  auto        vol = state.fields(HydroFieldNames::volume, 0.0);
  auto        surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, FacetedVolume> cells;
  FieldList<Dimension, vector<CellFaceFlag>> cellFaceFlags;
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    cells = state.fields(HydroFieldNames::cells, FacetedVolume());
    cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, vector<CellFaceFlag>());
  }
  if (mVolumeType == CRKVolumeType::CRKMassOverDensity) {
    vol.assignFields(mass/massDensity);
  } else if (mVolumeType == CRKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, W, position, mass, H, vol);
  } else if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    vol.assignFields(mass/massDensity);
    computeVoronoiVolume(position, H, connectivityMap, damage,
                         vector<typename Dimension::FacetedVolume>(),                // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),       // no holes
                         vector<Boundary<Dimension>*>(this->boundaryBegin(),         // boundaries
                                                      this->boundaryEnd()),
                         FieldList<Dimension, typename Dimension::Scalar>(),         // no weights
                         surfacePoint, vol, mDeltaCentroid, mEtaVoidPoints,          // return values
                         cells,                                                      // return cells
                         cellFaceFlags);                                             // node cell multimaterial faces
    if (mLimitMultimaterialTopology) {
      for (ConstBoundaryIterator boundItr = this->boundaryBegin();
           boundItr != this->boundaryEnd();
           ++boundItr) (*boundItr)->applyFieldListGhostBoundary(surfacePoint);
      for (ConstBoundaryIterator boundItr = this->boundaryBegin();
           boundItr != this->boundaryEnd();
           ++boundItr) (*boundItr)->finalizeGhostBoundary();
      // editMultimaterialSurfaceTopology(surfacePoint,
      //                                  const_cast<ConnectivityMap<Dimension>&>(connectivityMap));
    }
  } else if (mVolumeType == CRKVolumeType::CRKHullVolume) {
    computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, vol);
  } else if (mVolumeType == CRKVolumeType::HVolume) {
    const Scalar nPerh = vol.nodeListPtrs()[0]->nodesPerSmoothingScale();
    computeHVolumes(nPerh, H, vol);
  } else {
    VERIFY2(false, "Unknown CRK volume weighting.");
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(vol);
    if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(cells);
      (*boundItr)->applyFieldListGhostBoundary(surfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
  //   // flagSurfaceNeighbors(surfacePoint, connectivityMap);
  //   vol = computeShepardsInterpolation(vol,
  //                                      connectivityMap,
  //                                      W,
  //                                      position,
  //                                      H,
  //                                      vol);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->applyFieldListGhostBoundary(vol);
  //   for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //        boundItr != this->boundaryEnd();
  //        ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // }

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (densityUpdate() == MassDensityType::RigorousSumDensity) {
    computeCRKSPHSumMassDensity(connectivityMap, W, position, mass, vol, H, massDensity);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr != this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr != this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  } else if (densityUpdate() == MassDensityType::VoronoiCellDensity) {
    massDensity.assignFields(mass/vol);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr != this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr != this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro before evaluating derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // Compute the kernel correction fields.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, int> surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  FieldList<Dimension, ThirdRankTensor> gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
  FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CRKSPH, 0.0);
  FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CRKSPH, Vector::zero);
  FieldList<Dimension, SymTensor> m2 = state.fields(HydroFieldNames::m2_CRKSPH, SymTensor::zero);
  FieldList<Dimension, ThirdRankTensor> m3 = state.fields(HydroFieldNames::m3_CRKSPH, ThirdRankTensor::zero);
  FieldList<Dimension, FourthRankTensor> m4 = state.fields(HydroFieldNames::m4_CRKSPH, FourthRankTensor::zero);
  FieldList<Dimension, Vector> gradm0 = state.fields(HydroFieldNames::gradM0_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradm1 = state.fields(HydroFieldNames::gradM1_CRKSPH, Tensor::zero);
  FieldList<Dimension, ThirdRankTensor> gradm2 = state.fields(HydroFieldNames::gradM2_CRKSPH, ThirdRankTensor::zero);
  FieldList<Dimension, FourthRankTensor> gradm3 = state.fields(HydroFieldNames::gradM3_CRKSPH, FourthRankTensor::zero);
  FieldList<Dimension, FifthRankTensor> gradm4 = state.fields(HydroFieldNames::gradM4_CRKSPH, FifthRankTensor::zero);

  // Change CRKSPH weights here if need be!
  const NodeCoupling couple;
  computeCRKSPHMoments(connectivityMap, W, vol, position, H, correctionOrder(), couple, m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
  computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, H, surfacePoint, correctionOrder(), A, B, C, gradA, gradB, gradC);
  if (mLimitMultimaterialTopology) zerothOrderSurfaceCorrections(A, B, C, gradA, gradB, gradC, m0, gradm0, surfacePoint);

  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(A);
    (*boundItr)->applyFieldListGhostBoundary(B);
    (*boundItr)->applyFieldListGhostBoundary(C);
    (*boundItr)->applyFieldListGhostBoundary(gradA);
    (*boundItr)->applyFieldListGhostBoundary(gradB);
    (*boundItr)->applyFieldListGhostBoundary(gradC);
  }

  // Get the artificial viscosity and initialize it.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();
  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               W);
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
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
// Provide a hook to be called after the state has been updated and 
// boundary conditions have been enforced.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
postStateUpdate(const typename Dimension::Scalar time,
                const typename Dimension::Scalar dt,
                const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Apply boundary conditions to the basic fluid state Fields.

  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto entropy = state.fields(HydroFieldNames::entropy, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) {
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
  }

  auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, vector<Vector>());

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(vol);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(entropy);
    if (compatibleEnergyEvolution()) {
      (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
    }
    (*boundaryItr)->applyFieldListGhostBoundary(A);
    (*boundaryItr)->applyFieldListGhostBoundary(B);
    (*boundaryItr)->applyFieldListGhostBoundary(C);
    (*boundaryItr)->applyFieldListGhostBoundary(gradA);
    (*boundaryItr)->applyFieldListGhostBoundary(gradB);
    (*boundaryItr)->applyFieldListGhostBoundary(gradC);
    (*boundaryItr)->applyFieldListGhostBoundary(surfacePoint);
    (*boundaryItr)->applyFieldListGhostBoundary(etaVoidPoints);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Enforce boundary conditions on the fluid state Fields.
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto entropy = state.fields(HydroFieldNames::entropy, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) {
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
  }

  auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(vol);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->enforceFieldListBoundary(entropy);
    if (compatibleEnergyEvolution()) {
      (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
    }
    (*boundaryItr)->enforceFieldListBoundary(A);
    (*boundaryItr)->enforceFieldListBoundary(B);
    (*boundaryItr)->enforceFieldListBoundary(C);
    (*boundaryItr)->enforceFieldListBoundary(gradA);
    (*boundaryItr)->enforceFieldListBoundary(gradB);
    (*boundaryItr)->enforceFieldListBoundary(gradC);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.write(mEntropy, pathName + "/entropy");
  file.write(mHideal, pathName + "/Hideal");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.write(mEffViscousPressure, pathName + "/effViscousPressure");
  file.write(mViscousWork, pathName + "/viscousWork");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mVolume, pathName + "/Volume");
  file.write(mA, pathName + "/A");
  file.write(mB, pathName + "/B");
  file.write(mC, pathName + "/C");
  file.write(mGradA, pathName + "/gradA");
  file.write(mGradB, pathName + "/gradB");
  file.write(mGradC, pathName + "/gradC");
  file.write(mSurfacePoint, pathName + "/surfacePoint");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.read(mEntropy, pathName + "/entropy");
  file.read(mHideal, pathName + "/Hideal");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.read(mEffViscousPressure, pathName + "/effViscousPressure");
  file.read(mViscousWork, pathName + "/viscousWork");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mVolume, pathName + "/Volume");
  file.read(mA, pathName + "/A");
  file.read(mB, pathName + "/B");
  file.read(mC, pathName + "/C");
  file.read(mGradA, pathName + "/gradA");
  file.read(mGradB, pathName + "/gradB");
  file.read(mGradC, pathName + "/gradC");
  file.read(mSurfacePoint, pathName + "/surfacePoint");
}

}

