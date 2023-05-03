//---------------------------------Spheral++----------------------------------//
// SolidFSISPHHydroBase -- SolidSPHHydro modified for large density discontinuities
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

#include "Physics/GenericHydro.hh"                 

#include "GSPH/computeSPHVolume.hh"
#include "GSPH/Policies/PureReplaceFieldList.hh"
#include "GSPH/Policies/ReplaceWithRatioPolicy.hh"

#include "Hydro/HydroFieldNames.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"

#include "Strength/SolidFieldNames.hh"
#include "Strength/DeviatoricStressPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/PlasticStrainPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"
#include "Strength/StrengthSoundSpeedPolicy.hh"

#include "Damage/DamagedPressurePolicy.hh"

#include "NodeList/SolidNodeList.hh"
#include "NodeList/SmoothingScaleBase.hh"

#include "SolidMaterial/SolidEquationOfState.hh" 

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"

#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/Timer.hh"

#include "FSISPH/SolidFSISPHHydroBase.hh"
#include "FSISPH/FSIFieldNames.hh"
#include "FSISPH/computeFSISPHSumMassDensity.hh"
#include "FSISPH/computeHWeightedFSISPHSumMassDensity.hh"
#include "FSISPH/computeInterfacePressureCorrectedSumMassDensity.hh"
#include "FSISPH/SlideSurface.hh"


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
using std::make_shared;

namespace Spheral {


inline
Dim<1>::SymTensor
tensileStressCorrection(const Dim<1>::SymTensor& sigma) {
  if (sigma.xx() > 0.0) {
    return -sigma;
  } else {
    return Dim<1>::SymTensor::zero;
  }
}

inline
Dim<2>::SymTensor
tensileStressCorrection(const Dim<2>::SymTensor& sigma) {
  const EigenStruct<2> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  Dim<2>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

inline
Dim<3>::SymTensor
tensileStressCorrection(const Dim<3>::SymTensor& sigma) {
  const EigenStruct<3> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  const double lambdaz = eigen.eigenValues.z();
  Dim<3>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,                              0.0,
                           0.0,                              (lambday > 0.0 ? -lambday : 0.0), 0.0,
                           0.0,                              0.0,                              (lambdaz > 0.0 ? -lambdaz : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}


//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidFSISPHHydroBase<Dimension>::
SolidFSISPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                  DataBase<Dimension>& dataBase,
                  ArtificialViscosity<Dimension>& Q,
                  SlideSurface<Dimension>& slides,
                  const TableKernel<Dimension>& W,
                  const double filter,
                  const double cfl,
                  const double surfaceForceCoefficient,
                  const double densityStabilizationCoefficient,
                  const double specificThermalEnergyDiffusionCoefficient,
                  const double xsphCoefficient,
                  const InterfaceMethod interfaceMethod,
                  const KernelAveragingMethod kernelAveragingMethod,
                  const std::vector<int> sumDensityNodeLists,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool evolveTotalEnergy,
                  const bool gradhCorrection,
                  const bool XSPH,
                  const bool correctVelocityGradient,
                  const MassDensityType densityUpdate,
                  const HEvolutionType HUpdate,
                  const double epsTensile,
                  const double nTensile,
                  const double interfacePmin,
                  const bool planeStrain,
                  const bool damageRelieveRubble,
                  const bool strengthInDamage,
                  const Vector& xmin,
                  const Vector& xmax):
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mKernel(W),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mSlideSurface(slides),
  mPlaneStrain(planeStrain),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mCorrectVelocityGradient(correctVelocityGradient),
  mInterfacePmin(interfacePmin),
  mSurfaceForceCoefficient(surfaceForceCoefficient),
  mDensityStabilizationCoefficient(densityStabilizationCoefficient),
  mSpecificThermalEnergyDiffusionCoefficient(specificThermalEnergyDiffusionCoefficient),
  mXSPHCoefficient(xsphCoefficient),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mxmin(xmin),
  mxmax(xmax),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mInterfaceMethod(interfaceMethod),
  mKernelAveragingMethod(kernelAveragingMethod),
  mApplySelectDensitySum(false),
  mSumDensityNodeLists(sumDensityNodeLists),
  mTimeStepMask(FieldStorageType::CopyFields),
  mVolume(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mRawPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mBulkModulus(FieldStorageType::CopyFields),
  mShearModulus(FieldStorageType::CopyFields),
  mYieldStrength(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
  mMaxViscousPressure(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldStorageType::CopyFields),
  mMassSecondMoment(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDdeviatoricStressDt(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mDPDx(FieldStorageType::CopyFields),
  mDepsDx(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  mLocalM(FieldStorageType::CopyFields),
  mInterfaceFlags(FieldStorageType::CopyFields),
  mInterfaceAreaVectors(FieldStorageType::CopyFields),
  mInterfaceNormals(FieldStorageType::CopyFields),
  mInterfaceSmoothness(FieldStorageType::CopyFields),
  mNewInterfaceFlags(FieldStorageType::CopyFields),
  mNewInterfaceAreaVectors(FieldStorageType::CopyFields),
  mNewInterfaceNormals(FieldStorageType::CopyFields),
  mInterfaceSmoothnessNormalization(FieldStorageType::CopyFields),
  mInterfaceFraction(FieldStorageType::CopyFields),
  mNewInterfaceSmoothness(FieldStorageType::CopyFields),
  mInterfaceAngles(FieldStorageType::CopyFields),
  mPairDepsDt(),
  mPairAccelerations(){

    mPairDepsDt.clear();
    mPairAccelerations.clear();

    // see if we're summing density for any nodelist
    auto numNodeLists = dataBase.numNodeLists();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      if (sumDensityNodeLists[nodeListi]==1){
        mApplySelectDensitySum = true;
      } 
    }
    
    mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
    mVolume = dataBase.newFluidFieldList(0.0,HydroFieldNames::volume);
    mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
    mRawPressure = dataBase.newFluidFieldList(0.0, FSIFieldNames::rawPressure);
    mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
    mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
    mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
    mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
    mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
    mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
    mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
    mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
    mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
    mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
    mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
    mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position);
    mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
    mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
    mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
    mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
    mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H);
    mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
    mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
    mDPDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::pressureGradient);
    mDepsDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::specificThermalEnergyGradient);
    mM = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::M_SPHCorrection);
    mLocalM = dataBase.newFluidFieldList(Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection);
    mInterfaceFlags = dataBase.newFluidFieldList(int(0),  FSIFieldNames::interfaceFlags);
    mInterfaceAreaVectors = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::interfaceAreaVectors);
    mInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::interfaceNormals);
    mInterfaceSmoothness = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceSmoothness);
    mNewInterfaceFlags = dataBase.newFluidFieldList(int(0), PureReplaceFieldList<Dimension,int>::prefix() + FSIFieldNames::interfaceFlags);
    mNewInterfaceAreaVectors = dataBase.newFluidFieldList(Vector::one, PureReplaceFieldList<Dimension,Vector>::prefix() + FSIFieldNames::interfaceAreaVectors);
    mNewInterfaceNormals = dataBase.newFluidFieldList(Vector::one, PureReplaceFieldList<Dimension,Vector>::prefix() + FSIFieldNames::interfaceNormals);
    mInterfaceSmoothnessNormalization = dataBase.newFluidFieldList(0.0, FSIFieldNames::interfaceSmoothnessNormalization);
    mInterfaceFraction = dataBase.newFluidFieldList(0.0, FSIFieldNames::interfaceFraction);
    mNewInterfaceSmoothness = dataBase.newFluidFieldList(0.0, PureReplaceFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness);
    mInterfaceAngles = dataBase.newFluidFieldList(0.0, FSIFieldNames::interfaceAngles);
    
  }

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SolidFSISPHHydroBase<Dimension>::
~SolidFSISPHHydroBase() {
}

//------------------------------------------------------------------------------
// initialization on problem start up
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase){

  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
  mRawPressure+=this->pressure();

  auto mass = dataBase.fluidMass();
  auto massDensity = dataBase.fluidMassDensity();
  computeSPHVolume(mass,massDensity,mVolume);

  auto nodeListi = 0;
  for (auto itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    (*itr)->bulkModulus(*mBulkModulus[nodeListi]);
    (*itr)->shearModulus(*mShearModulus[nodeListi]);
    (*itr)->yieldStrength(*mYieldStrength[nodeListi]);
  }

}

//------------------------------------------------------------------------------
// Register states
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("SolidFSISPHregisterState");

  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "FSISPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  FieldList<Dimension, SymTensor> deviatoricStress = dataBase.solidDeviatoricStress();
  FieldList<Dimension, Scalar> plasticStrain = dataBase.solidPlasticStrain();
  FieldList<Dimension, SymTensor> damage = dataBase.solidDamage();
  FieldList<Dimension, int> fragIDs = dataBase.solidFragmentIDs();
  FieldList<Dimension, int> pTypes = dataBase.solidParticleTypes();

  CHECK(mass.numFields() == dataBase.numFluidNodeLists());
  CHECK(massDensity.numFields() == dataBase.numFluidNodeLists());
  CHECK(Hfield.numFields() == dataBase.numFluidNodeLists());
  CHECK(position.numFields() == dataBase.numFluidNodeLists());
  CHECK(specificThermalEnergy.numFields() == dataBase.numFluidNodeLists());
  CHECK(velocity.numFields() == dataBase.numFluidNodeLists());
  CHECK(deviatoricStress.numFields() == dataBase.numFluidNodeLists());
  CHECK(plasticStrain.numFields() == dataBase.numFluidNodeLists());
  CHECK(damage.numFields() == dataBase.numFluidNodeLists());
  CHECK(fragIDs.numFields() == dataBase.numFluidNodeLists());
  CHECK(pTypes.numFields() == dataBase.numFluidNodeLists());

  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeFluidFieldList(mRawPressure, 0.0, FSIFieldNames::rawPressure, false);
  dataBase.resizeFluidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeFluidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeFluidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);
  dataBase.resizeFluidFieldList(mInterfaceFlags, int(0), FSIFieldNames::interfaceFlags,false);
  dataBase.resizeFluidFieldList(mInterfaceAreaVectors, Vector::zero, FSIFieldNames::interfaceAreaVectors,false);
  dataBase.resizeFluidFieldList(mInterfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mInterfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness,false);
  dataBase.resizeFluidFieldList(mVolume, 0.0, HydroFieldNames::volume,false);

  auto positionPolicy = make_shared<IncrementFieldList<Dimension, Vector>>();
  auto plasticStrainPolicy = make_shared<PlasticStrainPolicy<Dimension>>();
  auto bulkModulusPolicy = make_shared<BulkModulusPolicy<Dimension>>();
  auto shearModulusPolicy = make_shared<ShearModulusPolicy<Dimension>>();
  auto yieldStrengthPolicy = make_shared<YieldStrengthPolicy<Dimension>>();
  auto soundSpeedPolicy = make_shared<StrengthSoundSpeedPolicy<Dimension>>();
  auto pressurePolicy = make_shared<DamagedPressurePolicy<Dimension>>();
  auto rawPressurePolicy = make_shared<PressurePolicy<Dimension>>();
  auto interfaceFlagsPolicy = make_shared<PureReplaceFieldList<Dimension,int>>(PureReplaceFieldList<Dimension,int>::prefix() + FSIFieldNames::interfaceFlags);
  auto interfaceAreaVectorsPolicy = make_shared<PureReplaceFieldList<Dimension,Vector>>(PureReplaceFieldList<Dimension,Vector>::prefix() + FSIFieldNames::interfaceAreaVectors);
  auto interfaceNormalsPolicy = make_shared<PureReplaceFieldList<Dimension,Vector>>(PureReplaceFieldList<Dimension,Vector>::prefix() + FSIFieldNames::interfaceNormals);
  auto interfaceSmoothnessPolicy = make_shared<PureReplaceFieldList<Dimension,Scalar>>(PureReplaceFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness);
  auto volumePolicy = make_shared<ReplaceWithRatioPolicy<Dimension,Scalar>>(HydroFieldNames::mass,HydroFieldNames::massDensity,HydroFieldNames::massDensity);
  
  std::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > rhoPolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  std::shared_ptr<CompositeFieldListPolicy<Dimension, SymTensor> > Hpolicy(new CompositeFieldListPolicy<Dimension, SymTensor>());
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    rhoPolicy->push_back(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),(*itr)->rhoMax()));
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == HEvolutionType::IntegrateH) {
      Hpolicy->push_back(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    } else {
      CHECK(HEvolution() == HEvolutionType::IdealH);
      Hpolicy->push_back(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    }
  }

  if(this->compatibleEnergyEvolution()){
    auto thermalEnergyPolicy = make_shared<CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>>(dataBase);
    auto velocityPolicy = make_shared<IncrementFieldList<Dimension, Vector>>(HydroFieldNames::position,HydroFieldNames::specificThermalEnergy,true);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
  }else if (this->evolveTotalEnergy()) {
    auto thermalEnergyPolicy = make_shared<SpecificFromTotalThermalEnergyPolicy<Dimension>>();
    auto velocityPolicy = make_shared<IncrementFieldList<Dimension, Vector>>(HydroFieldNames::position,HydroFieldNames::specificThermalEnergy,true);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
  } else {
    auto thermalEnergyPolicy = make_shared<IncrementFieldList<Dimension, Scalar>>();
    auto velocityPolicy = make_shared<IncrementFieldList<Dimension, Vector>>(HydroFieldNames::position,true);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
  }

  if(this->planeStrain()){
    auto deviatoricStressPolicy = make_shared<IncrementFieldList<Dimension, SymTensor>>();
    state.enroll(deviatoricStress, deviatoricStressPolicy);
  }else{
    auto deviatoricStressPolicy = make_shared<DeviatoricStressPolicy<Dimension>>();
    state.enroll(deviatoricStress, deviatoricStressPolicy);
  }

  state.enroll(mTimeStepMask);
  state.enroll(mass);
  state.enroll(mVolume,volumePolicy);
  state.enroll(massDensity, rhoPolicy);
  state.enroll(Hfield, Hpolicy);
  state.enroll(position, positionPolicy);
  state.enroll(mBulkModulus, bulkModulusPolicy);
  state.enroll(mShearModulus, shearModulusPolicy);
  state.enroll(mYieldStrength, yieldStrengthPolicy);
  state.enroll(plasticStrain, plasticStrainPolicy);
  state.enroll(mSoundSpeed, soundSpeedPolicy);
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mRawPressure,rawPressurePolicy);
  state.enroll(mInterfaceFlags,interfaceFlagsPolicy);
  state.enroll(mInterfaceAreaVectors,interfaceAreaVectorsPolicy); 
  state.enroll(mInterfaceNormals,interfaceNormalsPolicy); 
  state.enroll(mInterfaceSmoothness,interfaceSmoothnessPolicy); 
  state.enroll(damage);
  state.enroll(fragIDs);
  state.enroll(pTypes);


  TIME_END("SolidFSISPHregisterState");
}

//------------------------------------------------------------------------------
// Register Derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>&  dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidFSISPHregisterDerivs");

  FieldList<Dimension, Scalar> plasticStrainRate = dataBase.solidPlasticStrainRate();
  
  // make sure we're tracking the right number of node lists
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  dataBase.resizeFluidFieldList(mNormalization, 0.0, HydroFieldNames::normalization, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSPHWeightSum, 0.0, HydroFieldNames::XSPHWeightSum, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mDPDx, Vector::zero, FSIFieldNames::pressureGradient, false);
  dataBase.resizeFluidFieldList(mDepsDx, Vector::zero, FSIFieldNames::specificThermalEnergyGradient, false);
  dataBase.resizeFluidFieldList(mM, Tensor::zero, HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mLocalM, Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mNewInterfaceFlags, 0,  PureReplaceFieldList<Dimension,int>::prefix() + FSIFieldNames::interfaceFlags,false);
  dataBase.resizeFluidFieldList(mNewInterfaceAreaVectors, Vector::zero,  PureReplaceFieldList<Dimension,Vector>::prefix() + FSIFieldNames::interfaceAreaVectors,false);
  dataBase.resizeFluidFieldList(mNewInterfaceNormals, Vector::zero,  PureReplaceFieldList<Dimension,Vector>::prefix() + FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mInterfaceSmoothnessNormalization, 0.0, FSIFieldNames::interfaceSmoothnessNormalization,false); 
  dataBase.resizeFluidFieldList(mInterfaceFraction, 0.0, FSIFieldNames::interfaceFraction,false); 
  dataBase.resizeFluidFieldList(mNewInterfaceSmoothness, 0.0,  PureReplaceFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness,false);
  dataBase.resizeFluidFieldList(mInterfaceAngles, 0.0,  FSIFieldNames::interfaceAngles,false);

  // Check if someone already registered DxDt.
  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }

  // Check that no-one else is trying to control the hydro vote for DvDt.
  CHECK(not derivs.registered(mDvDt));
  derivs.enroll(mDvDt);

  // enroll 
  derivs.enroll(mHideal);
  derivs.enroll(mMaxViscousPressure);
  derivs.enroll(mNormalization);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mXSPHWeightSum);
  derivs.enroll(mXSPHDeltaV);
  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDdeviatoricStressDt);
  derivs.enroll(plasticStrainRate);
  derivs.enroll(mDHDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enroll(mDPDx);
  derivs.enroll(mDepsDx);
  derivs.enroll(mM);
  derivs.enroll(mLocalM);
  derivs.enroll(mNewInterfaceFlags);
  derivs.enroll(mNewInterfaceAreaVectors);
  derivs.enroll(mNewInterfaceNormals);
  derivs.enroll(mInterfaceSmoothnessNormalization);
  derivs.enroll(mInterfaceFraction);
  derivs.enroll(mNewInterfaceSmoothness);
  derivs.enroll(mInterfaceAngles);
  derivs.enrollAny(HydroFieldNames::pairAccelerations, mPairAccelerations);
  derivs.enrollAny(HydroFieldNames::pairWork,  mPairDepsDt);

  TIME_END("SolidFSISPHregisterDerivs");
}


//------------------------------------------------------------------------------
// evaluate the derivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                          StateDerivatives<Dimension>& derivatives) const {

  // bump out 
  const double interfaceNeighborAngleThreshold = 0.707; // (equiv to 90 deg opening angle threshold)

  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // Get the SlideSurfaces.
  auto& slides = this->slideSurface();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  const auto tiny = std::numeric_limits<double>::epsilon();
  const auto tinyScalarDamage = 1.0e-2;
  const auto tinyNonDimensional = 1.0e-6;

  const auto interfacePmin = this->interfacePmin();
  const auto W0 = W(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto totalEnergy = this->evolveTotalEnergy();
  const auto epsDiffusionCoeff = this->specificThermalEnergyDiffusionCoefficient();
  const auto rhoStabilizeCoeff = this->densityStabilizationCoefficient();
  const auto surfaceForceCoeff = this->surfaceForceCoefficient();
  const auto xsphCoeff = this->xsphCoefficient();
  const auto XSPH = xsphCoeff > tiny;
  const auto diffuseEnergy = epsDiffusionCoeff>tiny and compatibleEnergy;
  const auto stabilizeDensity = rhoStabilizeCoeff>tiny;
  const auto alwaysAverageKernels = (mKernelAveragingMethod==KernelAveragingMethod::AlwaysAverageKernels);
  const auto averageInterfaceKernels = (mKernelAveragingMethod==KernelAveragingMethod::AverageInterfaceKernels);
  const auto constructHLLC = (mInterfaceMethod == InterfaceMethod::HLLCInterface);
  const auto activateConstruction = !(mInterfaceMethod == InterfaceMethod::NoInterface);
  const auto strainTraceConstant =  (this->planeStrain() ? 1.0/3.0 : 1.0/Dimension::nDim);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  numPairs = pairs.size();
  const auto  nPerh = nodeLists[0]->nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);

  // Get the state and derivative FieldLists.
  const auto interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  const auto interfaceFlags = state.fields(FSIFieldNames::interfaceFlags, int(0));
  const auto interfaceAreaVectors = state.fields(FSIFieldNames::interfaceAreaVectors, Vector::zero);
  const auto interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  CHECK(fragIDs.size()==numNodeLists);
  CHECK(interfaceFlags.size() == numNodeLists);
  CHECK(interfaceAreaVectors.size() == numNodeLists);
  CHECK(interfaceNormals.size() == numNodeLists);
  CHECK(interfaceSmoothness.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(rawPressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(K.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);

  // Derivative FieldLists.
  const auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  const auto  DepsDx = derivatives.fields(FSIFieldNames::specificThermalEnergyGradient, Vector::zero);
  const auto  DPDx = derivatives.fields(FSIFieldNames::pressureGradient, Vector::zero);
  auto  newInterfaceNormals = derivatives.fields(PureReplaceFieldList<Dimension, Vector>::prefix() + FSIFieldNames::interfaceNormals, Vector::zero);
  auto  newInterfaceFlags = derivatives.fields(PureReplaceFieldList<Dimension, int>::prefix() + FSIFieldNames::interfaceFlags, int(0));
  auto  newInterfaceAreaVectors = derivatives.fields(PureReplaceFieldList<Dimension, Vector>::prefix() + FSIFieldNames::interfaceAreaVectors, Vector::zero);
  auto  interfaceSmoothnessNormalization = derivatives.fields(FSIFieldNames::interfaceSmoothnessNormalization, 0.0);
  auto  interfaceFraction = derivatives.fields(FSIFieldNames::interfaceFraction, 0.0);
  auto  newInterfaceSmoothness = derivatives.fields(PureReplaceFieldList<Dimension, Scalar>::prefix() + FSIFieldNames::interfaceSmoothness, 0.0);
  auto  interfaceAngles = derivatives.fields(FSIFieldNames::interfaceAngles, 0.0);
  auto  normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto  DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto& pairDepsDt = derivatives.getAny(HydroFieldNames::pairWork, vector<Scalar>());
  
  CHECK(newInterfaceFlags.size() == numNodeLists);
  CHECK(newInterfaceAreaVectors.size() == numNodeLists);
  CHECK(newInterfaceNormals.size() == numNodeLists);
  CHECK(interfaceSmoothnessNormalization.size() == numNodeLists);
  CHECK(interfaceFraction.size() == numNodeLists);
  CHECK(newInterfaceSmoothness.size() == numNodeLists);
  CHECK(interfaceAngles.size() == numNodeLists);
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
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if(compatibleEnergy){
    pairAccelerations.resize(numPairs);
    pairDepsDt.resize(2*numPairs);
  }

  this->computeMCorrection(time,dt,dataBase,state,derivatives);

// Now we calculate  the hydro deriviatives
// Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, Wj, gWj, PLineari, PLinearj, epsLineari, epsLinearj;
    Tensor QPiij, QPiji;
    SymTensor sigmai, sigmaj;
    Vector sigmarhoi, sigmarhoj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto newInterfaceFlags_thread = newInterfaceFlags.threadCopy(threadStack, ThreadReduction::MAX);
    auto newInterfaceAreaVectors_thread = newInterfaceAreaVectors.threadCopy(threadStack);
    auto newInterfaceNormals_thread = newInterfaceNormals.threadCopy(threadStack);
    auto newInterfaceSmoothness_thread = newInterfaceSmoothness.threadCopy(threadStack);
    auto interfaceSmoothnessNormalization_thread = interfaceSmoothnessNormalization.threadCopy(threadStack);
    auto interfaceFraction_thread = interfaceFraction.threadCopy(threadStack);
    auto interfaceAngles_thread = interfaceAngles.threadCopy(threadStack, ThreadReduction::MAX);
    auto normalization_thread = normalization.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DrhoDt_thread = DrhoDt.threadCopy(threadStack);
    auto DSDt_thread = DSDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    
#pragma omp for
    for (auto kk = 0u; kk < numPairs; ++kk) {
      
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& interfaceFlagsi = interfaceFlags(nodeListi,i);
      const auto& interfaceAreaVectorsi = interfaceAreaVectors(nodeListi,i);
      const auto& interfaceNormalsi = interfaceNormals(nodeListi,i);
      const auto& interfaceSmoothnessi = interfaceSmoothness(nodeListi,i);
      const auto& ri = position(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& voli = volume(nodeListi,i);
      const auto& Hi = H(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi,i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& rPi = rawPressure(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& pTypei = pTypes(nodeListi, i);
      const auto  mui = max(mu(nodeListi,i),tiny);
      const auto  Ki = max(tiny,K(nodeListi,i))+4.0/3.0*mui;
      const auto  Hdeti = Hi.Determinant();
      const auto fragIDi = fragIDs(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DepsDxi = DepsDx(nodeListi, i);
      const auto& DPDxi = DPDx(nodeListi, i);
      const auto& Mi = M(nodeListi, i);
      auto& normi = normalization_thread(nodeListi,i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& newInterfaceFlagsi = newInterfaceFlags_thread(nodeListi,i);
      auto& newInterfaceAreaVectorsi = newInterfaceAreaVectors_thread(nodeListi,i);
      auto& newInterfaceNormalsi = newInterfaceNormals_thread(nodeListi,i);
      auto& newInterfaceSmoothnessi = newInterfaceSmoothness_thread(nodeListi,i);
      auto& interfaceFractioni = interfaceFraction_thread(nodeListi,i);
      auto& interfaceSmoothnessNormalizationi = interfaceSmoothnessNormalization_thread(nodeListi,i);
      auto& minNeighborAnglei = interfaceAngles_thread(nodeListi,i);
      
      // Get the state for node j
      const auto& interfaceFlagsj = interfaceFlags(nodeListj,j);
      const auto& interfaceAreaVectorsj = interfaceAreaVectors(nodeListj,j);
      const auto& interfaceNormalsj = interfaceNormals(nodeListj,j);
      const auto& interfaceSmoothnessj = interfaceSmoothness(nodeListj,j);
      const auto& rj = position(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto  volj = volume(nodeListj,j);
      const auto& Hj = H(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj,j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& rPj = rawPressure(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto& pTypej = pTypes(nodeListj, j);
      const auto&  muj = max(mu(nodeListj,j),tiny);
      const auto  Kj = max(tiny,K(nodeListj,j))+4.0/3.0*muj;
      const auto  Hdetj = Hj.Determinant();
      const auto fragIDj = fragIDs(nodeListj, j);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      const auto& DepsDxj = DepsDx(nodeListj, j);
      const auto& DPDxj = DPDx(nodeListj, j);
      const auto& Mj = M(nodeListj,j);
      auto& normj = normalization_thread(nodeListj,j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& newInterfaceFlagsj = newInterfaceFlags_thread(nodeListj,j);
      auto& newInterfaceAreaVectorsj = newInterfaceAreaVectors_thread(nodeListj,j);
      auto& newInterfaceNormalsj = newInterfaceNormals_thread(nodeListj,j);
      auto& newInterfaceSmoothnessj = newInterfaceSmoothness_thread(nodeListj,j);
      auto& interfaceFractionj = interfaceFraction_thread(nodeListj,j);
      auto& interfaceSmoothnessNormalizationj = interfaceSmoothnessNormalization_thread(nodeListj,j);
      auto& minNeighborAnglej = interfaceAngles_thread(nodeListj,j);

      // line of action
      const auto rij = ri - rj;
      const auto rhatij = rij.unitVector();
      const auto rMagij = rij.magnitude();

      // decoupling and boolean switches
      //-------------------------------------------------------
      // Flag if this is a contiguous material pair or not.
      const auto sameMatij =  (nodeListi == nodeListj and fragIDi==fragIDj);
      const auto differentMatij = !sameMatij; 
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);
      
      // some thresholds
      const auto tinyPressure = tinyNonDimensional*min(Ki,Kj);

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);
      
      // pairwise damage and nodal damage
      //const auto fDij = (sameMatij ? pairs[kk].f_couple : 0.0);
      const auto Di = std::max(0.0, std::min(1.0, damage(nodeListi,i).eigenValues().maxElement()));
      const auto Dj = std::max(0.0, std::min(1.0, damage(nodeListj,j).eigenValues().maxElement()));
      const auto fDi =  (sameMatij ? max((1.0-Di)*(1.0-Di),tiny) : 0.0 );
      const auto fDj =  (sameMatij ? max((1.0-Dj)*(1.0-Dj),tiny) : 0.0 );
      const auto fDij = (sameMatij ? pow(1.0-std::abs(Di-Dj),2.0) : 0.0 );

      // is Pmin being activated (Pmin fsi has additional pmin for intefaces)
      const auto pMini = (sameMatij ? (Pi-tinyPressure) : interfacePmin);
      const auto pMinj = (sameMatij ? (Pj-tinyPressure) : interfacePmin);
      const auto pminActivei = (rPi < pMini);
      const auto pminActivej = (rPj < pMinj);
      
      // decoupling criteria 
      const auto isExpanding = (ri-rj).dot(vi-vj) > 0.0;
      const auto isFullyDamaged = (fDi<tinyScalarDamage) or (fDj<tinyScalarDamage);
      const auto isPastAdhesionThreshold = pminActivei or pminActivej;
      const auto decouple = isExpanding and (isFullyDamaged and isPastAdhesionThreshold);

      // do we need to construct our interface velocity?
      const auto constructInterface = (fDij < 1.0-tinyScalarDamage) and activateConstruction;
      const auto negligableShearWave = max(mui,muj) < tinyPressure;

      // do we reduce our deviatoric stress
      const auto isTensile = (((Si+Sj)-(Pi+Pj)*SymTensor::one).dot(rhatij)).dot(rhatij) > 0;
      const auto damageReduceStress = isTensile or differentMatij;

      // Kernels
      //--------------------------------------
      const auto Hij = 0.5*(Hi+Hj);
      const auto etaij = Hij*rij;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagij = etaij.magnitude();
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagij >= 0.0);
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      const auto Hetai = Hi*etai.unitVector();
      const auto Hetaj = Hj*etaj.unitVector();
      auto gradWi = gWi*Hetai;
      auto gradWj = gWj*Hetaj;
      auto gradWiMi = gradWi;
      auto gradWjMj = gradWj;

      // average our kernels
      const auto gradWij = 0.5*(gradWi+gradWj);
      const auto Wij = 0.5*(Wi+Wj);
      if(averageKernelij){
        const auto gWij = 0.5*(gWi+gWj);
        Wi = Wij;
        Wj = Wij;
        gWi = gWij;
        gWj = gWij;
        gradWi = gradWij;
        gradWj = gradWij;
      }

      if(this->correctVelocityGradient()){
        gradWiMi = Mi.Transpose()*gradWi;
        gradWjMj = Mj.Transpose()*gradWj;
      }

      // interface normals 
      //-----------------------------------------------------------
      const auto fSij = ( sameMatij ? 1.0 : 0.0);              // direction parameter
      const auto fSij3 = ( sameMatij ? 1.0 : -1.0);              // direction parameter
      const auto AijMij = fSij3*(voli*volj)*(gradWjMj+gradWiMi); // surface area vector

      const auto alignment = (interfaceFlagsi == 5 or interfaceFlagsj == 5 ? 
                              0.0:
                              std::max(fSij3*interfaceNormalsi.dot(interfaceNormalsj),0.0));

      newInterfaceAreaVectorsi -= AijMij;
      newInterfaceAreaVectorsj += AijMij;

      if (sameMatij){
        interfaceFractioni += volj * Wij;
        interfaceFractionj += voli * Wij;
        // we use the angle between normal and neighbors to check for interface control particles (type 2 or 4)
        if (interfaceFlagsj != 5) minNeighborAnglei = std::max(-interfaceAreaVectorsi.unitVector().dot(rhatij),minNeighborAnglei);
        if (interfaceFlagsi != 5) minNeighborAnglej = std::max( interfaceAreaVectorsj.unitVector().dot(rhatij),minNeighborAnglej);
        if (interfaceFlagsj == 2) newInterfaceFlagsi = std::max(newInterfaceFlagsi,1);
        if (interfaceFlagsi == 2) newInterfaceFlagsj = std::max(newInterfaceFlagsj,1);
      }

      if(differentMatij and (interfaceFlagsi!=5 and interfaceFlagsj!=5)){
        newInterfaceFlagsi=3;
        newInterfaceFlagsj=3;
      }
      if (interfaceFlagsj == 4) newInterfaceFlagsi = std::max(newInterfaceFlagsi,3);
      if (interfaceFlagsi == 4) newInterfaceFlagsj = std::max(newInterfaceFlagsj,3);

      if (interfaceFlagsi > 0 and interfaceFlagsj > 0){
        const double proxWeighti = 1.0 - (interfaceFlagsi % 2);
        const double proxWeightj = 1.0 - (interfaceFlagsj % 2);
        newInterfaceNormalsi += proxWeightj*fSij3*volj*interfaceAreaVectorsj*Wij;
        newInterfaceNormalsj += proxWeighti*fSij3*voli*interfaceAreaVectorsi*Wij;
        // interfaceSmoothnessNormalizationi -= volj*(fSij3*interfaceNormalsi-interfaceNormalsj).dot(gradWiMi);//volj*interfaceFlagsj*Wij;
        // interfaceSmoothnessNormalizationj -= voli*(fSij3*interfaceNormalsi-interfaceNormalsj).dot(gradWjMj);//voli*interfaceFlagsi*Wij;  
      }
      const auto ninja = interfaceNormalsi.magnitude() * interfaceNormalsj.magnitude();
      interfaceSmoothnessNormalizationi -= ninja*volj*(rij).dot(gradWiMi);//volj*interfaceFlagsj*Wij;
      interfaceSmoothnessNormalizationj -= ninja*voli*(rij).dot(gradWjMj);//voli*interfaceFlagsi*Wij;
      newInterfaceSmoothnessi += ninja*volj*(fSij3*interfaceNormalsi-interfaceNormalsj).dot(gradWiMi);//alignment*interfaceFlagsj*volj*Wij;
      newInterfaceSmoothnessj += ninja*voli*(fSij3*interfaceNormalsi-interfaceNormalsj).dot(gradWjMj);//alignment*interfaceFlagsi*voli*Wij;
      //newInterfaceSmoothnessi += 2.0*fSij*volj*rij.dot(gradWiMi)*safeInv(rMagij,tiny);//fSij*std::min(interfaceFlagsj,1)*std::min(interfaceFlagsi,1)*volj*(fSij3*interfaceNormalsi-interfaceNormalsj).dot(gradWiMi);//alignment*interfaceFlagsj*volj*Wij;
      //newInterfaceSmoothnessj += 2.0*fSij*voli*rij.dot(gradWjMj)*safeInv(rMagij,tiny);
      // normalization 
      //-----------------------------------------------------------
      normi += volj*Wi;
      normj += voli*Wj;

      if (!decouple){

        // Zero'th and second moment of the node distribution -- used for the
        // ideal H calculation.
        //---------------------------------------------------------------
        const auto rij2 = rij.magnitude2();
        const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
        weightedNeighborSumi += abs(gWi);
        weightedNeighborSumj += abs(gWj);
        massSecondMomenti += gradWi.magnitude2()*thpt;
        massSecondMomentj += gradWj.magnitude2()*thpt;

        // Stress state
        //---------------------------------------------------------------
        const auto rhoij = (rhoi*rhoi+rhoj*rhoj)/(rhoi+rhoj); 
        const auto cij = (rhoi*ci+rhoj*cj)/(rhoi+rhoj); 
        const auto vij = vi - vj;

        std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                        ri, etaij, vi, rhoij, cij, Hij,  
                                        rj, etaij, vj, rhoij, cij, Hij); 

        //slide correction
        if (slides.isSlideSurface(nodeListi,nodeListj)){
          
          const auto slideCorr = slides.slideCorrection(interfaceSmoothnessi,
                                                        interfaceSmoothnessj,
                                                        interfaceNormalsi,
                                                        interfaceNormalsj,
                                                        vi,
                                                        vj);
          QPiij *= slideCorr;
          QPiji *= slideCorr;
        }

        // save our max pressure from the AV for each node
        maxViscousPressurei = max(maxViscousPressurei, rhoi*rhoj * QPiij.diagonalElements().maxAbsElement());
        maxViscousPressurej = max(maxViscousPressurej, rhoi*rhoj * QPiji.diagonalElements().maxAbsElement());

        // stress tensor
        const auto Seffi = (damageReduceStress ? fDij : 1.0) * Si;
        const auto Seffj = (damageReduceStress ? fDij : 1.0) * Sj;
        const auto Peffi = (differentMatij ? max(Pi,interfacePmin) : Pi);
        const auto Peffj = (differentMatij ? max(Pj,interfacePmin) : Pj);
        sigmai = Seffi - Peffi * SymTensor::one;
        sigmaj = Seffj - Peffj * SymTensor::one;

        // Compute the tensile correction to add to the stress as described in 
        // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
        const auto fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
        const auto fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
        const auto Ri = fi*tensileStressCorrection(sigmai);
        const auto Rj = fj*tensileStressCorrection(sigmaj);
        sigmai += Ri;
        sigmaj += Rj;

        // accelerations
        //---------------------------------------------------------------
        const auto rhoirhoj = safeInv(rhoi*rhoj,tiny);
        const auto sf = (sameMatij ? 1.0 : 1.0 + surfaceForceCoeff*abs((rhoi-rhoj)*safeInv(rhoi+rhoj,tiny)));
        sigmarhoi = sf*(rhoirhoj*sigmai-0.5*QPiij)*gradWiMi;
        sigmarhoj = sf*(rhoirhoj*sigmaj-0.5*QPiji)*gradWjMj;

        if (averageKernelij){
          const auto sigmarhoij = 0.5*(sigmarhoi+sigmarhoj);
          sigmarhoi = sigmarhoij;
          sigmarhoj = sigmarhoij;
        }
      
        const auto deltaDvDt = sigmarhoi + sigmarhoj;

        if (freeParticle) {
          DvDti += mj*deltaDvDt;
          DvDtj -= mi*deltaDvDt;
        } 
      
        // Velocity Gradient
        //-----------------------------------------------------------

        // construct our interface velocity 
        auto vstar = 0.5*(vi+vj);

        linearReconstruction(ri,rj,Pi,Pj,DPDxi,DPDxj,PLineari,PLinearj);

        if (constructInterface){
          
          // components
          const auto ui = vi.dot(rhatij);
          const auto uj = vj.dot(rhatij);
          const auto wi = vi - ui*rhatij;
          const auto wj = vj - uj*rhatij;

          // acoustic wave speeds
          const auto Ci =  (constructHLLC ? std::sqrt(rhoi*Ki)  : Ki  ) + tiny;
          const auto Cj =  (constructHLLC ? std::sqrt(rhoj*Kj)  : Kj  ) + tiny;
          const auto Csi = (constructHLLC ? std::sqrt(rhoi*mui) : mui ) + tiny;
          const auto Csj = (constructHLLC ? std::sqrt(rhoj*muj) : muj ) + tiny;

          const auto CiCjInv = safeInv(Ci+Cj,tiny);
          const auto CsiCsjInv = safeInv(Csi+Csj,tiny);
          
          // weights
          const auto weightUi = max(0.0, min(1.0, Ci*CiCjInv));
          const auto weightUj = 1.0 - weightUi;
          const auto weightWi = (negligableShearWave ? weightUi : max(0.0, min(1.0, Csi*CsiCsjInv )) );
          const auto weightWj = 1.0 - weightWi;

          // interface velocity
          const auto ustar = weightUi*ui + weightUj*uj + (constructHLLC ? (PLinearj - PLineari)*CiCjInv : 0.0); 
          const auto wstar = weightWi*wi + weightWj*wj;
          vstar = fDij * vstar + (1.0-fDij)*(ustar*rhatij + wstar);
  
        }

        // local velocity gradient for DSDt
        if (sameMatij) {
          localDvDxi -=  2.0*volj*((vi-vstar).dyad(gradWi));
          localDvDxj -=  2.0*voli*((vstar-vj).dyad(gradWj)); 
        }

        // diffuse to stabilize things
        if (stabilizeDensity and (ci>tiny and cj>tiny)){
          const auto cFactor = 1.0 + max(min( (vi-vj).dot(rhatij)/max(cij,tiny), 0.0), -1.0);
          const auto effCoeff = (differentMatij ? 1.0 : rhoStabilizeCoeff*cFactor);
          vstar += (constructHLLC ? fDij : 1.0) * effCoeff * rhatij * cij * min(max((PLinearj-PLineari)/(Ki + Kj),-0.25),0.25);
        }

        // global velocity gradient
        DvDxi -= 2.0*volj*(vi-vstar).dyad(gradWiMi);
        DvDxj -= 2.0*voli*(vstar-vj).dyad(gradWjMj);

        // energy conservation
        // ----------------------------------------------------------
        const auto deltaDepsDti = 2.0*sigmarhoi.dot(vi-vstar);
        const auto deltaDepsDtj = 2.0*sigmarhoj.dot(vstar-vj);

        DepsDti -= mj*deltaDepsDti;
        DepsDtj -= mi*deltaDepsDtj;

        if(compatibleEnergy){
          pairAccelerations[kk] = - deltaDvDt;
          pairDepsDt[2*kk]   = - deltaDepsDti; 
          pairDepsDt[2*kk+1] = - deltaDepsDtj;
        }
        
        // thermal diffusion
        //-----------------------------------------------------------
        if (sameMatij and diffuseEnergy){
          linearReconstruction(ri,rj,epsi,epsj,DepsDxi,DepsDxj,epsLineari,epsLinearj);
          const auto cijEff = max(min(cij + (vi-vj).dot(rhatij), cij),0.0);
          const auto diffusion =  epsDiffusionCoeff*cijEff*(epsLineari-epsLinearj)*etaij.dot(gradWij)/(rhoij*etaMagij*etaMagij+tiny);
          pairDepsDt[2*kk]   += diffusion; 
          pairDepsDt[2*kk+1] -= diffusion;
        }

        // XSPH -- we use this to handle tensile instability here
        //-----------------------------------------------------------
        if (sameMatij and XSPH) {
          const auto fxsph  = (min(Pi,Pj) < 0.0 ? 1.0 : 0.0);
          XSPHWeightSumi += fxsph*volj*Wi;
          XSPHWeightSumj += fxsph*voli*Wj;
          XSPHDeltaVi -= volj*Wi*(vi-vstar);
          XSPHDeltaVj -= voli*Wj*(vj-vstar);
        }

      } // if damageDecouple 
    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
  } // OpenMP parallel region


  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& mui = mu(nodeListi, i);
      const auto& interfaceFlagsi = interfaceFlags(nodeListi,i);
      const auto& interfaceAreaVectorsi = interfaceAreaVectors(nodeListi,i);
      const auto  Hdeti = Hi.Determinant();
      const auto psi = Hdeti*mi/rhoi*W0;
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      const auto& DvDti = DvDt(nodeListi,i);
      const auto& localMi = localM(nodeListi, i);
      auto& normi = normalization(nodeListi,i);
      auto& DepsDti = DepsDt(nodeListi,i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);
      auto& newInterfaceNormalsi = newInterfaceNormals(nodeListi,i);
      auto& newInterfaceSmoothnessi = newInterfaceSmoothness(nodeListi,i);
      auto& interfaceSmoothnessNormalizationi = interfaceSmoothnessNormalization(nodeListi,i);
      auto& interfaceFractioni = interfaceFraction(nodeListi,i);
      auto& newInterfaceFlagsi = newInterfaceFlags(nodeListi,i);
      auto& minNeighborAnglei = interfaceAngles(nodeListi,i);

      // finish our normalization
      normi +=  psi;

      // finish our interface fields.
      if(minNeighborAnglei < interfaceNeighborAngleThreshold)  newInterfaceFlagsi = std::max(2,newInterfaceFlagsi+1);
      if(interfaceFractioni/(1.0-psi)<0.1) newInterfaceFlagsi = 5.0;

      interfaceFractioni += psi;

      if ( interfaceFlagsi > 0 ){
        const double proxWeighti = 100*(1.0 - interfaceFlagsi % 2);
        newInterfaceNormalsi = (newInterfaceNormalsi + proxWeighti * psi * interfaceAreaVectorsi).unitVector();
      } else {
        newInterfaceNormalsi = Vector::zero;
        //newInterfaceSmoothnessi = 0.0;
      }
      //interfaceSmoothnessNormalizationi += psi * interfaceFlagsi;
      //newInterfaceSmoothnessi = min(1.0,max(0.0,(newInterfaceSmoothnessi+psi * interfaceFlagsi)/max(interfaceSmoothnessNormalizationi,tiny)));
      newInterfaceSmoothnessi = newInterfaceSmoothnessi/max(interfaceSmoothnessNormalizationi,tiny);
      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;
 
      // continuity
      DrhoDti -=  rhoi*DvDxi.Trace();

      // energy
      if (totalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // XSPH
      DxDti = vi;
      if (XSPH) {
        XSPHWeightSumi += psi;
        DxDti += xsphCoeff*XSPHWeightSumi/max(normi*normi,tiny) * XSPHDeltaVi;
      }


      // H - Evolution
      DHDti = smoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                            ri,
                                                            DvDxi,
                                                            hmin,
                                                            hmax,
                                                            hminratio,
                                                            nPerh);
      
      Hideali = smoothingScaleMethod.newSmoothingScale(Hi,
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

      if(this->correctVelocityGradient()) localDvDxi = localDvDxi*localMi;

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()*strainTraceConstant)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti += spinCorrection + 2.0*mui*deviatoricDeformation;
      
    } //loop-nodes
  }   //loop-nodeLists
}     // evaluateDerivatives method


//------------------------------------------------------------------------------
// EvalDerivs subroutine for spatial derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
computeMCorrection(const typename Dimension::Scalar /*time*/,
                   const typename Dimension::Scalar /*dt*/,
                   const DataBase<Dimension>& dataBase,
                   const State<Dimension>& state,
                         StateDerivatives<Dimension>& derivatives) const {

  // The kernels and such.
  const auto& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  const auto alwaysAverageKernels = (mKernelAveragingMethod==KernelAveragingMethod::AlwaysAverageKernels);
  const auto averageInterfaceKernels = (mKernelAveragingMethod==KernelAveragingMethod::AverageInterfaceKernels);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  numPairs = pairs.size();

  // Get the state and derivative FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DepsDx = derivatives.fields(FSIFieldNames::specificThermalEnergyGradient, Vector::zero);
  auto  DPDx = derivatives.fields(FSIFieldNames::pressureGradient, Vector::zero);
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  
  CHECK(DepsDx.size() == numNodeLists);
  CHECK(DPDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);

#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    auto DPDx_thread = DPDx.threadCopy(threadStack);
    auto DepsDx_thread = DepsDx.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < numPairs; ++kk) {

      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& fragIDi = fragIDs(nodeListi,i);
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      // Get the state for node j
      const auto& fragIDj = fragIDs(nodeListj,j);
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DPDxi = DPDx_thread(nodeListi, i);
      auto& DPDxj = DPDx_thread(nodeListj, j);
      auto& DepsDxi = DepsDx_thread(nodeListi, i);
      auto& DepsDxj = DepsDx_thread(nodeListj, j);
      auto& localMi = localM_thread(nodeListi,i);
      auto& localMj = localM_thread(nodeListj,j);
      auto& Mi = M_thread(nodeListi,i);
      auto& Mj = M_thread(nodeListj,j);

      const auto rij = ri - rj;
      const auto Pij = Pi - Pj;
      const auto epsij = epsi - epsj;

      // logic
      //---------------------------------------
      const auto sameMatij = (nodeListi == nodeListj and fragIDi == fragIDj);
      const auto differentMatij = (nodeListi!=nodeListj);
      const auto averageKernelij = ( (differentMatij and averageInterfaceKernels) or alwaysAverageKernels);

      // Kernels
      //--------------------------------------
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      const auto gWi = W.gradValue(etaMagi, Hdeti);
      const auto gWj = W.gradValue(etaMagj, Hdetj);
      
      const auto Hetai = Hi*etai.unitVector();
      const auto Hetaj = Hj*etaj.unitVector();
      
      auto gradWi = gWi*Hetai;
      auto gradWj = gWj*Hetaj;
      
      //Wi & Wj --> Wij for interface better agreement DrhoDt and DepsDt
      if(averageKernelij){
        const auto gradWij = 0.5*(gradWi+gradWj);
        gradWi = gradWij;
        gradWj = gradWij;
      }

      gradWi *= mj/rhoj;
      gradWj *= mi/rhoi;

      // spatial gradients and correction
      //---------------------------------------------------------------
      const auto deltaRi = rij.dyad(gradWi);
      const auto deltaRj = rij.dyad(gradWj);

      Mi -= deltaRi;
      Mj -= deltaRj;

      DPDxi -= Pij*gradWi;
      DPDxj -= Pij*gradWj;

      if(sameMatij){
        localMi -=  deltaRi;
        localMj -=  deltaRj;
        DepsDxi -= epsij*gradWi;
        DepsDxj -= epsij*gradWj;
      }
    } // loop over pairs
      // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);
  }   // OpenMP parallel region

   
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto& nodeList = mass[nodeListi]->nodeList();
      const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
      for (auto i = 0u; i < ni; ++i) {
        const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
        auto& Mi = M(nodeListi, i);
        auto& localMi = localM(nodeListi, i);
        auto& DepsDxi = DepsDx(nodeListi, i);
        auto& DPDxi = DPDx(nodeListi, i);

        const auto Mdeti = Mi.Determinant();
        const auto goodM = ( Mdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
        Mi =  (goodM ? Mi.Inverse() : Tensor::one);
        
        const auto localMdeti = localMi.Determinant();
        const auto goodLocalM = ( localMdeti > 1.0e-2 and numNeighborsi > Dimension::pownu(2));
        localMi =  (goodLocalM ? localMi.Inverse() : Tensor::one);

        DPDxi = Mi.Transpose()*DPDxi;
        DepsDxi = localMi.Transpose()*DepsDxi;

      } // for each node
    }   // for each nodelist

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(M);
      (*boundaryItr)->applyFieldListGhostBoundary(DPDx);
      (*boundaryItr)->applyFieldListGhostBoundary(DepsDx);
    }
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

} // method 


//------------------------------------------------------------------------------
// FSI specialized density summmation
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_BEGIN("SolidFSISPHpreStepInitialize");
  //const auto dummySwitchVar = 2;
  //switch(dummySwitchVar){
  //  case 1:
      // const auto& connectivityMap = dataBase.connectivityMap();
      // const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
      // const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
      // const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
      // const auto& W = this->kernel();
      //       auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      // computeFSISPHSumMassDensity(connectivityMap, W, mSumDensityNodeLists, position, mass, H, massDensity);
      // for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      // for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  //    break;
  //   case 2:
      const auto& W = this->kernel();
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
            auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
            auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
            auto volume = state.fields(HydroFieldNames::volume, 0.0);
            auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeInterfacePressureCorrectedSumMassDensity(connectivityMap, 
                                                      W, 
                                                      mSumDensityNodeLists, 
                                                      position, 
                                                      mass, 
                                                      H,
                                                      volume,
                                                      pressure,
                                                      soundSpeed,
                                                      massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  // //     break;
  //   case 3:
  //     const auto& connectivityMap = dataBase.connectivityMap();
  //     const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
  //     const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
  //     const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
  //     const auto& W = this->kernel();
  //           auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //     computeHWeightedFSISPHSumMassDensity(connectivityMap, W, mSumDensityNodeLists, position, mass, H, massDensity);
  //     for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
  //     for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  //     break;
  //   default:
  //     break;

  // } // switch
  
  TIME_END("SolidFSISPHpreStepInitialize");
}

//------------------------------------------------------------------------------
// FSI specialized of the initialize method
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("SolidFSISPHinitialize");

  const TableKernel<Dimension>& W = this->kernel();

  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               W);
  
  // We depend on the caller knowing to finalize the ghost boundaries!
  TIME_END("SolidFSISPHinitialize");
}

//------------------------------------------------------------------------------
// For compatible energy we need to apply the bc conditions to acceleration
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
finalizeDerivatives(const Scalar /*time*/, 
                    const Scalar  /*dt*/,
                    const DataBase<Dimension>&  /*dataBase*/, 
                    const State<Dimension>& /*state*/,
                          StateDerivatives<Dimension>&  derivs) const {
                            
  if (this->compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) {
           (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    }
    
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }

} // finalize


//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  FieldList<Dimension, int> pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  FieldList<Dimension, int> interfaceFlags = state.fields(FSIFieldNames::interfaceFlags, int(0));
  FieldList<Dimension, Vector> interfaceAreaVectors = state.fields(FSIFieldNames::interfaceAreaVectors, Vector::zero);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(rawPressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(S);
    (*boundaryItr)->applyFieldListGhostBoundary(K);
    (*boundaryItr)->applyFieldListGhostBoundary(mu);
    (*boundaryItr)->applyFieldListGhostBoundary(Y);
    (*boundaryItr)->applyFieldListGhostBoundary(fragIDs);
    (*boundaryItr)->applyFieldListGhostBoundary(pTypes);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceFlags);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceAreaVectors);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceNormals);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceSmoothness);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  FieldList<Dimension, int> pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  FieldList<Dimension, int> interfaceFlags = state.fields(FSIFieldNames::interfaceFlags, int(0));
  FieldList<Dimension, Vector> interfaceAreaVectors = state.fields(FSIFieldNames::interfaceAreaVectors, Vector::zero);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(volume);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(rawPressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->enforceFieldListBoundary(S);
    (*boundaryItr)->enforceFieldListBoundary(K);
    (*boundaryItr)->enforceFieldListBoundary(mu);
    (*boundaryItr)->enforceFieldListBoundary(Y);
    (*boundaryItr)->enforceFieldListBoundary(fragIDs);
    (*boundaryItr)->enforceFieldListBoundary(pTypes);
    (*boundaryItr)->enforceFieldListBoundary(interfaceFlags);
    (*boundaryItr)->enforceFieldListBoundary(interfaceAreaVectors);
    (*boundaryItr)->enforceFieldListBoundary(interfaceNormals);
    (*boundaryItr)->enforceFieldListBoundary(interfaceSmoothness);
  }

}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mVolume, pathName + "/volume");
  file.write(mPressure, pathName + "/pressure");
  file.write(mRawPressure, pathName + "/rawEosPressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mBulkModulus, pathName + "/bulkModulus");
  file.write(mShearModulus, pathName + "/shearModulus");
  file.write(mYieldStrength, pathName + "/yieldStrength");

  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mDPDx, pathName + "/DpDx");
  file.write(mDepsDx, pathName + "/DepsDx");
  file.write(mM, pathName + "/M");
  file.write(mLocalM, pathName + "/localM");

  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.write(mDHDt, pathName + "/DHDt");

  file.write(mInterfaceFlags, pathName + "/interfaceFlags");
  file.write(mInterfaceAreaVectors, pathName + "/interfaceAreaVectors");
  file.write(mInterfaceNormals, pathName + "/interfaceNormals");
  file.write(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  
  file.write(mNewInterfaceFlags, pathName + "/newInterfaceFlags");
  file.write(mNewInterfaceAreaVectors, pathName + "/newInterfaceAreaVectors");
  file.write(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.write(mInterfaceSmoothnessNormalization, pathName + "/interfaceSmoothnessNormalization");
  file.write(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");
  file.write(mInterfaceAngles, pathName + "/interfaceAngles");

  file.write(mHideal, pathName + "/Hideal");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.write(mNormalization, pathName + "/normalization");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");

}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mVolume, pathName + "/volume");
  file.read(mPressure, pathName + "/pressure");
  file.read(mRawPressure, pathName + "/rawEosPressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mBulkModulus, pathName + "/bulkModulus");
  file.read(mShearModulus, pathName + "/shearModulus");
  file.read(mYieldStrength, pathName + "/yieldStrength");

  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mDPDx, pathName + "/DpDx");
  file.read(mDepsDx, pathName + "/DepsDx");
  file.read(mM, pathName + "/M");
  file.read(mLocalM, pathName + "/localM");

  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.read(mDHDt, pathName + "/DHDt");

  file.read(mInterfaceFlags, pathName + "/interfaceFlags");
  file.read(mInterfaceAreaVectors, pathName + "/interfaceAreaVectors");
  file.read(mInterfaceNormals, pathName + "/interfaceNormals");
  file.read(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  
  file.read(mNewInterfaceFlags, pathName + "/newInterfaceFlags");
  file.read(mNewInterfaceAreaVectors, pathName + "/newInterfaceAreaVectors");
  file.read(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.read(mInterfaceSmoothnessNormalization, pathName + "/interfaceSmoothnessNormalization");
  file.read(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");
  file.read(mInterfaceAngles, pathName + "/interfaceAngles");

  file.read(mHideal, pathName + "/Hideal");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.read(mNormalization, pathName + "/normalization");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");
}

//------------------------------------------------------------------------------
// method for limited linear reconstruction between nodes
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
linearReconstruction(const typename Dimension::Vector& ri,
                     const typename Dimension::Vector& rj,
                     const typename Dimension::Scalar& yi,
                     const typename Dimension::Scalar& yj,
                     const typename Dimension::Vector& DyDxi,
                     const typename Dimension::Vector& DyDxj,
                           typename Dimension::Scalar& ytildei,
                           typename Dimension::Scalar& ytildej) const {
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  const auto rij = (ri-rj);

  // relavant deltas in field value
  const auto Dy0 = (yi-yj);
  const auto Dyi = 0.5*DyDxi.dot(rij);
  const auto Dyj = 0.5*DyDxj.dot(rij);

  // ratios of SPH derivs to ij particle difference
  const auto denom = 2.0 / (sgn(Dy0) * std::max(tiny,abs(Dy0)));
  const auto xi = Dyi * denom;
  const auto xj = Dyj * denom;

  // limiter function - vanleer 1979
  const auto phii = ( xi > 0.0 ?  min(4.0*xi/((1.0 + xi)*(1.0 + xi)),1.0) : 0.0 );
  const auto phij = ( xj > 0.0 ?  min(4.0*xj/((1.0 + xj)*(1.0 + xj)),1.0) : 0.0 );        
  const auto phi = 0.5*(phii+phij);

  // linear constructed inteface values
  ytildei = yi - phi * Dyi;
  ytildej = yj + phi * Dyj;
}

} // Spheral namespace



