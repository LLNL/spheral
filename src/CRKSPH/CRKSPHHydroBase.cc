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
#include "flagSurfaceNeighbors.hh"
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

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using PhysicsSpace::Physics;
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
using Geometry::innerProduct;
using Geometry::outerProduct;
using BoundarySpace::CRKSPHVoidBoundary;

using PhysicsSpace::MassDensityType;
using PhysicsSpace::HEvolutionType;

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
                const bool detectSurfaces,
                const double detectThreshold,
                const double sweepAngle,
                const double detectRange,
                const double epsTensile,
                const double nTensile):
  GenericHydro<Dimension>(W, WPi, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCorrectionOrder(correctionOrder),
  mVolumeType(volumeType),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mfilter(filter),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mDetectSurfaces(detectSurfaces),
  mDetectThreshold(detectThreshold),
  mSweepAngle(sweepAngle),
  mDetectRange(detectRange),
  mCorrectionMin(std::numeric_limits<Scalar>::lowest()),
  mCorrectionMax(std::numeric_limits<Scalar>::max()),
  mTimeStepMask(FieldSpace::FieldStorageType::CopyFields),
  mPressure(FieldSpace::FieldStorageType::CopyFields),
  mSoundSpeed(FieldSpace::FieldStorageType::CopyFields),
  mSpecificThermalEnergy0(FieldSpace::FieldStorageType::CopyFields),
  mEntropy(FieldSpace::FieldStorageType::CopyFields),
  mHideal(FieldSpace::FieldStorageType::CopyFields),
  mMaxViscousPressure(FieldSpace::FieldStorageType::CopyFields),
  mEffViscousPressure(FieldSpace::FieldStorageType::CopyFields),
  mViscousWork(FieldSpace::FieldStorageType::CopyFields),
  mVolume(FieldSpace::FieldStorageType::CopyFields),
  mMassDensityGradient(FieldSpace::FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldSpace::FieldStorageType::CopyFields),
  mMassSecondMoment(FieldSpace::FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldSpace::FieldStorageType::CopyFields),
  mDxDt(FieldSpace::FieldStorageType::CopyFields),
  mDvDt(FieldSpace::FieldStorageType::CopyFields),
  mDmassDensityDt(FieldSpace::FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldSpace::FieldStorageType::CopyFields),
  mDHDt(FieldSpace::FieldStorageType::CopyFields),
  mDvDx(FieldSpace::FieldStorageType::CopyFields),
  mInternalDvDx(FieldSpace::FieldStorageType::CopyFields),
  mPairAccelerations(FieldSpace::FieldStorageType::CopyFields),
  mA(FieldSpace::FieldStorageType::CopyFields),
  mB(FieldSpace::FieldStorageType::CopyFields),
  mC(FieldSpace::FieldStorageType::CopyFields),
  mGradA(FieldSpace::FieldStorageType::CopyFields),
  mGradB(FieldSpace::FieldStorageType::CopyFields),
  mGradC(FieldSpace::FieldStorageType::CopyFields),
  mM0(FieldSpace::FieldStorageType::CopyFields),
  mM1(FieldSpace::FieldStorageType::CopyFields),
  mM2(FieldSpace::FieldStorageType::CopyFields),
  mM3(FieldSpace::FieldStorageType::CopyFields),
  mM4(FieldSpace::FieldStorageType::CopyFields),
  mGradm0(FieldSpace::FieldStorageType::CopyFields),
  mGradm1(FieldSpace::FieldStorageType::CopyFields),
  mGradm2(FieldSpace::FieldStorageType::CopyFields),
  mGradm3(FieldSpace::FieldStorageType::CopyFields),
  mGradm4(FieldSpace::FieldStorageType::CopyFields),
  mSurfacePoint(FieldSpace::FieldStorageType::CopyFields),
  mVoidPoint(FieldSpace::FieldStorageType::CopyFields),
  mEtaVoidPoints(FieldSpace::FieldStorageType::CopyFields),
  mVoidBoundary(mSurfacePoint, mEtaVoidPoints),
  mRestart(DataOutput::registerWithRestart(*this)) {
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
  mMassDensityGradient = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::massDensityGradient);
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
  mPairAccelerations = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::pairAccelerations);
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
  mVoidPoint = dataBase.newFluidFieldList(0, HydroFieldNames::voidPoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::etaVoidPoints);
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  if (mDetectSurfaces) {
    mVolume.assignFields(mass/massDensity);
    for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr) (*boundItr)->applyFieldListGhostBoundary(mVolume);
    for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr) (*boundItr)->finalizeGhostBoundary();
    const NodeCoupling couple;
    computeCRKSPHMoments(connectivityMap, W, mVolume, position, H, correctionOrder(), couple, mM0, mM1, mM2, mM3, mM4, mGradm0, mGradm1, mGradm2, mGradm3, mGradm4);
    detectSurface(connectivityMap, mM0, mM1, position, H, mDetectThreshold, mDetectRange*W.kernelExtent(), mSweepAngle, mSurfacePoint);
  }

  // Compute the volumes for real.
  if (mVolumeType == CRKVolumeType::CRKMassOverDensity) {
    mVolume.assignFields(mass/massDensity);
  } else if (mVolumeType == CRKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, W, position, mass, H, mVolume);
  } else if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    mVolume.assignFields(mass/massDensity);
    FieldList<Dimension, typename Dimension::FacetedVolume> cells;
    computeVoronoiVolume(position, H, massDensity, mMassDensityGradient, connectivityMap, 
                         vector<typename Dimension::FacetedVolume>(),               // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),      // no holes
                         FieldList<Dimension, typename Dimension::Scalar>(),        // no weights
                         mVoidPoint,                                                // void point flags
                         mSurfacePoint, mVolume, mDeltaCentroid, mEtaVoidPoints,    // return values
                         cells);                                                    // no return cells
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
      (*boundItr)->applyFieldListGhostBoundary(mVolume);
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
  computeCRKSPHCorrections(mM0, mM1, mM2, mM3, mM4, mGradm0, mGradm1, mGradm2, mGradm3, mGradm4, H, correctionOrder(), mA, mB, mC, mGradA, mGradB, mGradC);

  // This breaks domain independence, so we'll try being inconsistent on the first step.
  // // We need to initialize the velocity gradient if we're using the CRKSPH artificial viscosity.
  // const FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  // mDvDx.assignFields(CRKSPHSpace::gradientCRKSPH(velocity, position, mVolume, H, mA, mB, mC, mGradA, mGradB, mGradC, connectivityMap, correctionOrder(), W, NodeCoupling()));

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
    dataBase.resizeFluidFieldList(mGradm3,   FourthRankTensor::zero,HydroFieldNames::m3_CRKSPH, false);
    dataBase.resizeFluidFieldList(mGradm4,   FifthRankTensor::zero, HydroFieldNames::m4_CRKSPH, false);
  }
  dataBase.resizeFluidFieldList(mSurfacePoint, 0, HydroFieldNames::surfacePoint, false);
  dataBase.resizeFluidFieldList(mVoidPoint, 0, HydroFieldNames::voidPoint, false);
  dataBase.resizeFluidFieldList(mEtaVoidPoints, vector<Vector>(), HydroFieldNames::etaVoidPoints, false);

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
    PolicyPointer thermalEnergyPolicy(new SpecificThermalEnergyPolicy<Dimension>(dataBase));   // Change back to non-symmetric if needed.
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
  state.enroll(mSurfacePoint);
  state.enroll(mVoidPoint);

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
  dataBase.resizeFluidFieldList(mMassDensityGradient, Vector::zero, HydroFieldNames::massDensityGradient, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mPairAccelerations, vector<Vector>(), HydroFieldNames::pairAccelerations, false);

  derivs.enroll(mHideal);
  derivs.enroll(mMaxViscousPressure);
  derivs.enroll(mEffViscousPressure);
  derivs.enroll(mViscousWork);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mMassDensityGradient);
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
  derivs.enroll(mPairAccelerations);
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
  const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  const NodeCoupling couple;
  computeCRKSPHMoments(connectivityMap, W, vol, position, H, correctionOrder(), couple, m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
  computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, H, correctionOrder(), A, B, C, gradA, gradB, gradC);

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
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  order = this->correctionOrder();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  const auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  const auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  const auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  const auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  const auto voidPoint = state.fields(HydroFieldNames::voidPoint, 0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(A.size() == numNodeLists);
  CHECK(B.size() == numNodeLists or order == CRKOrder::ZerothOrder);
  CHECK(C.size() == numNodeLists or order != CRKOrder::QuadraticOrder);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists or order == CRKOrder::ZerothOrder);
  CHECK(gradC.size() == numNodeLists or order != CRKOrder::QuadraticOrder);
  CHECK(surfacePoint.size() == numNodeLists);
  CHECK(voidPoint.size() == numNodeLists);

  // Derivative FieldLists.
  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  auto DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  auto DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto DepsDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto DHDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto gradRho = derivatives.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(gradRho.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (mCompatibleEnergyEvolution) {
    auto nodeListi = 0;
    for (auto itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        const size_t n = connectivityMap.numNeighborsForNode(*itr, i);
        pairAccelerations(nodeListi, i).reserve(n);
      }
    }
  }

  // Some scratch variables.
  Scalar Ai, Aj;
  Vector gradAi, gradAj, forceij, forceji;
  Vector Bi = Vector::zero, Bj = Vector::zero;
  Tensor Ci = Tensor::zero, Cj = Tensor::zero;
  Tensor gradBi = Tensor::zero, gradBj = Tensor::zero;
  ThirdRankTensor gradCi = ThirdRankTensor::zero, gradCj = ThirdRankTensor::zero;
  Scalar gWi, gWj, Wi, Wj, gW0i, gW0j, W0i, W0j;
  Vector gradWi, gradWj, gradW0i, gradW0j;
  Vector deltagrad;

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (auto itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const auto& nodeList = **itr;
    const auto firstGhostNodei = nodeList.firstGhostNode();
    const auto hmin = nodeList.hmin();
    const auto hmax = nodeList.hmax();
    const auto hminratio = nodeList.hminratio();
    const auto maxNumNeighbors = nodeList.maxNumNeighbors();
    const auto nPerh = nodeList.nodesPerSmoothingScale();

    // Get the work field for this NodeList.
    auto& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (auto iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Prepare to accumulate the time.
      const auto start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto rhoi = massDensity(nodeListi, i);
      const auto epsi = specificThermalEnergy(nodeListi, i);
      const auto Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto ci = soundSpeed(nodeListi, i);
      Ai = A(nodeListi, i);
      gradAi = gradA(nodeListi, i);
      if (order != CRKOrder::ZerothOrder) {
        Bi = B(nodeListi, i);
        gradBi = gradB(nodeListi, i);
      }
      if (order == CRKOrder::QuadraticOrder) {
        Ci = C(nodeListi, i);
        gradCi = gradC(nodeListi, i);
      }
      const auto Hdeti = Hi.Determinant();
      const auto weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);
      // CHECK2(Ai > 0.0, i << " " << Ai);
      CHECK2(Hdeti > 0.0, i << " " << Hdeti);
      CHECK2(weighti > 0.0, i << " " << weighti);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& viscousWorki = viscousWork(nodeListi, i);
      auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& gradRhoi = gradRho(nodeListi, i);
      auto& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const auto firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const auto& rj = position(nodeListj, j);
              const auto  mj = mass(nodeListj, j);
              const auto& vj = velocity(nodeListj, j);
              const auto  rhoj = massDensity(nodeListj, j);
              const auto  epsj = specificThermalEnergy(nodeListj, j);
              const auto  Pj = pressure(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              const auto  cj = soundSpeed(nodeListj, j);
              Aj = A(nodeListj, j);
              gradAj = gradA(nodeListj, j);
              if (order != CRKOrder::ZerothOrder) {
                Bj = B(nodeListj, j);
                gradBj = gradB(nodeListj, j);
              }
              if (order == CRKOrder::QuadraticOrder) {
                Cj = C(nodeListj, j);
                gradCj = gradC(nodeListj, j);
              }
              const auto Hdetj = Hj.Determinant();
              const auto weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Aj > 0.0 or j >= firstGhostNodej);
              CHECK(Hdetj > 0.0);
              CHECK(weightj > 0.0);

              auto& DxDtj = DxDt(nodeListj, j);
              auto& DrhoDtj = DrhoDt(nodeListj, j);
              auto& DvDtj = DvDt(nodeListj, j);
              auto& DepsDtj = DepsDt(nodeListj, j);
              auto& DvDxj = DvDx(nodeListj, j);
              auto& localDvDxj = localDvDx(nodeListj, j);
              auto& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              auto& effViscousPressurej = effViscousPressure(nodeListj, j);
              auto& viscousWorkj = viscousWork(nodeListj, j);
              auto& pairAccelerationsj = pairAccelerations(nodeListj, j);
              auto& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              auto& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              auto& massSecondMomentj = massSecondMoment(nodeListj, j);
              auto& gradRhoj = gradRho(nodeListj, j);

              // Find the effective weights of i->j and j->i.
              // const auto wi = 2.0*weighti*weightj/(weighti + weightj);
              // const auto wij = 0.5*(weighti + weightj);

              // Node displacement.
              const auto rij = ri - rj;
              const auto etai = Hi*rij;
              const auto etaj = Hj*rij;
              const auto etaMagi = etai.magnitude();
              const auto etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);
              const auto vij = vi - vj;

              // Symmetrized kernel weight and gradient.
              CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, order,  rij,  etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, mCorrectionMin, mCorrectionMax);
              CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, order, -rij, -etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj, mCorrectionMin, mCorrectionMax);
              deltagrad = gradWj - gradWi;
              const auto gradWSPHi = (Hi*etai.unitVector())*W.gradValue(etai.magnitude(), Hdeti);
              const auto gradWSPHj = (Hj*etaj.unitVector())*W.gradValue(etaj.magnitude(), Hdetj);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              if (voidPoint(nodeListi, i) == 0 and voidPoint(nodeListj, j) == 0) {
                const auto fweightij = nodeListi == nodeListj ? 1.0 : mj*rhoi/(mi*rhoj);
                const auto rij2 = rij.magnitude2();
                const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
                weightedNeighborSumi +=     fweightij*std::abs(gWi);
                weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
                massSecondMomenti +=     fweightij*gradWSPHi.magnitude2()*thpt;
                massSecondMomentj += 1.0/fweightij*gradWSPHj.magnitude2()*thpt;
              }

              // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
              const auto QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                        ri, etai, vi, rhoi, ci, Hi,
                                        rj, etaj, vj, rhoj, cj, Hj);
              const auto Qaccij = (rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second).dot(deltagrad);
              // const auto workQij = 0.5*(vij.dot(Qaccij));
              const auto workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);                // CRK
              const auto workQj = rhoi*rhoi*QPiij.first .dot(vij).dot(deltagrad);                // CRK
              // const auto workQVi =  vij.dot((rhoj*rhoj*QPiij.second).dot(gradWj));               //RK V and RK I Work
              // const auto workQVj =  vij.dot((rhoi*rhoi*QPiij.first).dot(gradWi));                //RK V and RK I Work
              const auto Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const auto Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                                 // We need tighter timestep controls on the Q with CRK
              maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
              //effViscousPressurei += wij * Qi * Wj;
              //effViscousPressurej += wij * Qj * Wi;
              effViscousPressurei += weightj * Qi * Wj;
              effViscousPressurej += weighti * Qj * Wi;
              //viscousWorki += 0.5*wij*wij/mi*workQi;
              //viscousWorkj += 0.5*wij*wij/mj*workQj;
              viscousWorki += 0.5*weighti*weightj/mi*workQi;
              viscousWorkj += 0.5*weighti*weightj/mj*workQj;

              // Velocity gradient.
              //DvDxi -= wij*vij.dyad(gradWj);
              //DvDxj += wij*vij.dyad(gradWi);
              //if (nodeListi == nodeListj) {
                //localDvDxi -= wij*vij.dyad(gradWj);
                //localDvDxj += wij*vij.dyad(gradWi);
              //}
              const Tensor deltaDvDxi = -weightj*vij.dyad(gradWj);
              const Tensor deltaDvDxj =  weighti*vij.dyad(gradWi);
              DvDxi += deltaDvDxi;
              DvDxj += deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi += deltaDvDxi;
                localDvDxj += deltaDvDxj;
              }

              // Mass density gradient.
              //gradRhoi += wij*(rhoj - rhoi)*gradWj;
              //gradRhoj += wij*(rhoi - rhoj)*gradWi;
              gradRhoi += weightj*(rhoj - rhoi)*gradWj;
              gradRhoj += weighti*(rhoi - rhoj)*gradWi;

              // We decide between RK and CRK for the momentum and energy equations based on the surface condition.
              // Momentum
              //forceij = (true ? // surfacePoint(nodeListi, i) <= 1 ? 
              //           0.5*wij*wij*((Pi + Pj)*deltagrad + Qaccij) :                    // Type III CRK interpoint force.
              //           mi*wij*((Pj - Pi)/rhoi*gradWj + rhoi*QPiij.first.dot(gradWj))); // RK
              //forceji = (true ? // surfacePoint(nodeListj, j) <= 1 ? 
              //           0.5*wij*wij*((Pi + Pj)*deltagrad + Qaccij) :                    // Type III CRK interpoint force.
              //           mj*wij*((Pj - Pi)/rhoj*gradWi - rhoj*QPiij.second.dot(gradWi)));// RK
              forceij = (true ? // surfacePoint(nodeListi, i) <= 1 ? 
                         0.5*weighti*weightj*((Pi + Pj)*deltagrad + Qaccij) :                    // Type III CRK interpoint force.
                         mi*weightj*((Pj - Pi)/rhoi*gradWj + rhoi*QPiij.first.dot(gradWj))); // RK
              forceji = (true ? // surfacePoint(nodeListj, j) <= 1 ? 
                         0.5*weighti*weightj*((Pi + Pj)*deltagrad + Qaccij) :                    // Type III CRK interpoint force.
                         mj*weighti*((Pj - Pi)/rhoj*gradWi - rhoj*QPiij.second.dot(gradWi)));// RK
              DvDti -= forceij/mi;
              DvDtj += forceji/mj; 
              if (mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(-forceij/mi);
                pairAccelerationsj.push_back( forceji/mj);
              }

              // Energy
              //DepsDti += (true ? // surfacePoint(nodeListi, i) <= 1 ? 
              //            0.5*wij*wij*(Pj*vij.dot(deltagrad) + workQi)/mi :              // CRK
              //            wij*rhoi*QPiij.first.dot(vij).dot(gradWj));                    // RK
              //DepsDtj += (true ? // surfacePoint(nodeListj, j) <= 1 ? 
              //            0.5*wij*wij*(Pi*vij.dot(deltagrad) + workQj)/mj :              // CRK
              //           -wij*rhoj*QPiij.second.dot(vij).dot(gradWi));                   // RK
              DepsDti += (true ? // surfacePoint(nodeListi, i) <= 1 ? 
                          0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + workQi)/mi :              // CRK
                          weightj*rhoi*QPiij.first.dot(vij).dot(gradWj));                    // RK
              DepsDtj += (true ? // surfacePoint(nodeListj, j) <= 1 ? 
                          0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + workQj)/mj :              // CRK
                         -weighti*rhoj*QPiij.second.dot(vij).dot(gradWi));                   // RK

              // Estimate of delta v (for XSPH).
              if (mXSPH and (nodeListi == nodeListj)) {
                //XSPHDeltaVi -= wij*Wj*vij;
                //XSPHDeltaVj += wij*Wi*vij;
                XSPHDeltaVi -= weightj*Wj*vij;
                XSPHDeltaVj += weighti*Wi*vij;
              }
                
            }
          }
        }
      }
      const auto numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not mCompatibleEnergyEvolution or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // // For a surface point, add the RK thermal energy evolution.
      // // DepsDti -= Pi/rhoi*DvDxi.Trace();
      // if (surfacePoint(nodeListi, i) > 1) DepsDti -= Pi/rhoi*DvDxi.Trace();

      // Get the time for pairwise interactions.
      const auto deltaTimePair = Timing::difference(start, Timing::currentTime())/max(size_t(1), ncalc);

      // Time evolution of the mass density.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             ri,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
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
      for (auto nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const auto& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const auto firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
          auto& workFieldj = nodeLists[nodeListj]->work();
#pragma vector always
          for (auto jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const auto j = *jItr;
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
// Finalize the state after state has been updated and boundary conditions 
// enforced.  For CRKSPH this is where we update the volumes and RPKM corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
postStateUpdate(const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                const StateDerivatives<Dimension>& derivs) const {

  // // Grab state we're going to use.
  //const TableKernel<Dimension>& W = this->kernel();
  //const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  //const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  //const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  // const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);

  // Compute the volume per node.
  //FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  //computeCRKSPHSumVolume(connectivityMap, W, position, H, vol);
  // computeHullVolumes(connectivityMap, position, vol);

  // // We need boundary conditions enforced on the volume before we can compute corrections.
  //for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //     boundItr != this->boundaryEnd();
  //     ++boundItr) (*boundItr)->applyFieldListGhostBoundary(vol);
  //for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //     boundItr != this->boundaryEnd();
  //     ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // // Compute the kernel correction fields.
  // FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  // FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  // FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CRKSPH, Vector::zero);
  // FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CRKSPH, Tensor::zero);
  // FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  // FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  // computeCRKSPHCorrections(connectivityMap, W, vol, position, H, A, B, C, D, gradA, gradB);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) {
  //   (*boundItr)->applyFieldListGhostBoundary(A);
  //   (*boundItr)->applyFieldListGhostBoundary(B);
  //   (*boundItr)->applyFieldListGhostBoundary(C);
  //   (*boundItr)->applyFieldListGhostBoundary(D);
  //   (*boundItr)->applyFieldListGhostBoundary(gradA);
  //   (*boundItr)->applyFieldListGhostBoundary(gradB);
  // }
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Base class finalization.
  GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Volume.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> gradRho = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  const FieldList<Dimension, int> voidPoint = state.fields(HydroFieldNames::voidPoint, 0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, int> surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  if (mVolumeType == CRKVolumeType::CRKMassOverDensity) {
    vol.assignFields(mass/massDensity);
  } else if (mVolumeType == CRKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, W, position, mass, H, vol);
  } else if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    vol.assignFields(mass/massDensity);
    FieldList<Dimension, typename Dimension::FacetedVolume> cells;
    computeVoronoiVolume(position, H, massDensity, gradRho, connectivityMap, 
                         vector<typename Dimension::FacetedVolume>(),                // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),       // no holes
                         FieldList<Dimension, typename Dimension::Scalar>(),         // no weights
                         voidPoint,                                                  // void point flags
                         surfacePoint, vol, mDeltaCentroid, mEtaVoidPoints,          // return values
                         cells);                                                     // no return cells
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
    computeCRKSPHSumMassDensity(connectivityMap, W, position, mass, vol, H, mVoidPoint, massDensity);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    }
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  } else if (densityUpdate() == MassDensityType::VoronoiCellDensity) {
    massDensity.assignFields(mass/vol);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    }
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }

  // Update the surface point flags.
  if (mDetectSurfaces) {
    const FieldList<Dimension, Scalar> m0 = state.fields(HydroFieldNames::m0_CRKSPH, 0.0);
    const FieldList<Dimension, Vector> m1 = state.fields(HydroFieldNames::m1_CRKSPH, Vector::zero);
    detectSurface(connectivityMap, m0, m1, position, H, mDetectThreshold, mDetectRange*W.kernelExtent(), mSweepAngle, surfacePoint);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(surfacePoint);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }

  // Add any filtering component to the node movement.
  // This form uses the deltaCentroid computed by computeVoronoiVolume, so only works if we're using that volume definition.
  if (mfilter > 0.0 and mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);  // Gotta get a non-const version now.
    const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
    const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
    const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
    const unsigned numNodeLists = position.numFields();
    Scalar minmag2, dcmag2, fi, mi, rhoi, Vi, V0i;
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        dcmag2 = mDeltaCentroid(nodeListi, i).magnitude2();
        // minmag2 = max(FastMath::square(soundSpeed(nodeListi, i)),
        //               velocity(nodeListi, i).magnitude2())*dt*dt;
        // minmag2 = min(FastMath::square(min(soundSpeed(nodeListi, i), abs(DvDx(nodeListi, i).Trace())/H(nodeListi, i).eigenValues().maxElement())),
        //               velocity(nodeListi, i).magnitude2())*dt*dt;
        minmag2 = FastMath::square(DvDx(nodeListi, i).eigenValues().maxAbsElement()/H(nodeListi, i).eigenValues().maxElement()*dt);
        // mi = mass(nodeListi, i);
        // rhoi = massDensity(nodeListi, i);
        // Vi = vol(nodeListi, i);
        fi = mfilter; // *max(0.0, min(1.0, max(V0i/Vi, Vi/V0i) - 1.0));
        position(nodeListi, i) += fi*sqrt(min(minmag2, dcmag2)*safeInvVar(dcmag2))*mDeltaCentroid(nodeListi, i);
      }
    }

    // Check for any boundary violations.
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
    this->enforceBoundaries(state, derivs);
  }
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
  auto voidPoint = state.fields(HydroFieldNames::voidPoint, 0);
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
    (*boundaryItr)->applyFieldListGhostBoundary(voidPoint);
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
dumpState(FileIO& file, string pathName) const {
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
  file.write(mMassDensityGradient, pathName + "/massDensityGradient");
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
restoreState(const FileIO& file, string pathName) {
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
  file.read(mMassDensityGradient, pathName + "/massDensityGradient");
  file.read(mA, pathName + "/A");
  file.read(mB, pathName + "/B");
  file.read(mC, pathName + "/C");
  file.read(mGradA, pathName + "/gradA");
  file.read(mGradB, pathName + "/gradB");
  file.read(mGradC, pathName + "/gradC");
  file.read(mSurfacePoint, pathName + "/surfacePoint");
}

}
}

