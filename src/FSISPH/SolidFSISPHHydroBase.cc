//---------------------------------Spheral++----------------------------------//
// SolidFSISPHHydroBase -- SolidSPHHydro modified for large density discontinuities
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

#include "Physics/GenericHydro.hh"

#include "GSPH/Policies/PureReplaceFieldList.hh"
#include "NodeList/SolidNodeList.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "SolidMaterial/SolidEquationOfState.hh" 

#include "Hydro/HydroFieldNames.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/PositionPolicy.hh"

#include "Strength/SolidFieldNames.hh"
#include "Strength/DeviatoricStressPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/PlasticStrainPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"
#include "Strength/StrengthSoundSpeedPolicy.hh"
#include "Damage/DamagedPressurePolicy.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
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
#include "Utilities/globalBoundingVolumes.hh"

#include "FSISPH/SolidFSISPHHydroBase.hh"
#include "FSISPH/FSIFieldNames.hh"
#include "FSISPH/computeFSISPHSumMassDensity.hh"
#include "FSISPH/computeHWeightedFSISPHSumMassDensity.hh"
#include "FSISPH/SlideSurface.hh"
#include "FSISPH/Policies/PairwisePlasticStrainPolicy.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

using std::make_shared;
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
                  const bool damageRelieveRubble,
                  const bool strengthInDamage,
                  const Vector& xmin,
                  const Vector& xmax):
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mKernel(W),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mSlideSurface(slides),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mInterfaceMethod(interfaceMethod),
  mKernelAveragingMethod(kernelAveragingMethod),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mCorrectVelocityGradient(correctVelocityGradient),
  mApplySelectDensitySum(false),
  mSumDensityNodeLists(sumDensityNodeLists),
  mSurfaceForceCoefficient(surfaceForceCoefficient),
  mDensityStabilizationCoefficient(densityStabilizationCoefficient),
  mSpecificThermalEnergyDiffusionCoefficient(specificThermalEnergyDiffusionCoefficient),
  mXSPHCoefficient(xsphCoefficient),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mxmin(xmin),
  mxmax(xmax),
  mPairAccelerations(),
  mPairDepsDt(),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mRawPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mBulkModulus(FieldStorageType::CopyFields),
  mShearModulus(FieldStorageType::CopyFields),
  mYieldStrength(FieldStorageType::CopyFields),
  mPlasticStrain0(FieldStorageType::CopyFields),
  mInverseEquivalentDeviatoricStress(FieldStorageType::CopyFields),
  mVolume(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDdeviatoricStressDt(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
  mDPDx(FieldStorageType::CopyFields),
  mDepsDx(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  mLocalM(FieldStorageType::CopyFields),
  mMaxViscousPressure(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldStorageType::CopyFields),
  mMassSecondMoment(FieldStorageType::CopyFields),
  mInterfaceNormals(FieldStorageType::CopyFields),
  mInterfaceFraction(FieldStorageType::CopyFields),
  mInterfaceSmoothness(FieldStorageType::CopyFields),
  mNewInterfaceNormals(FieldStorageType::CopyFields),
  mSmoothedInterfaceNormals(FieldStorageType::CopyFields),
  mNewInterfaceFraction(FieldStorageType::CopyFields),
  mNewInterfaceSmoothness(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {

    // see if we're summing density for any nodelist
    auto numNodeLists = dataBase.numNodeLists();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      if (sumDensityNodeLists[nodeListi]==1){
        mApplySelectDensitySum = true;
      } 
    }

    mPairDepsDt.clear();
    mPairAccelerations.clear();

    mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
    mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
    mRawPressure = dataBase.newFluidFieldList(0.0, FSIFieldNames::rawPressure);
    mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
    mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
    mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
    mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
    mPlasticStrain0 = dataBase.newSolidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
    mInverseEquivalentDeviatoricStress = dataBase.newFluidFieldList(0.0, FSIFieldNames::inverseEquivalentDeviatoricStress);
    mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position);
    mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
    mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
    mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
    mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
    mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
    mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
    mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H);
    mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
    mDPDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::pressureGradient);
    mDepsDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::specificThermalEnergyGradient);
    mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
    mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
    mM = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::M_SPHCorrection);
    mLocalM = dataBase.newFluidFieldList(Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection);
    mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
    mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
    mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
    mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
    mInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::interfaceNormals);
    mInterfaceFraction = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceFraction);
    mInterfaceSmoothness = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceSmoothness);
    mNewInterfaceNormals = dataBase.newFluidFieldList(Vector::one, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals);
    mSmoothedInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::smoothedInterfaceNormals);
    mNewInterfaceFraction = dataBase.newFluidFieldList(0.0, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceFraction);
    mNewInterfaceSmoothness = dataBase.newFluidFieldList(0.0, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness); 
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

  // Set the moduli.
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

  FieldList<Dimension, Vector> position = dataBase.solidPosition();
  FieldList<Dimension, Scalar> mass = dataBase.solidMass();
  FieldList<Dimension, Scalar> massDensity = dataBase.solidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.solidHfield();
  FieldList<Dimension, Vector> velocity = dataBase.solidVelocity();
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.solidSpecificThermalEnergy();
  FieldList<Dimension, SymTensor> deviatoricStress = dataBase.solidDeviatoricStress();
  FieldList<Dimension, Scalar> plasticStrain = dataBase.solidPlasticStrain();
  FieldList<Dimension, SymTensor> damage = dataBase.solidDamage();
  FieldList<Dimension, int> fragIDs = dataBase.solidFragmentIDs();
  FieldList<Dimension, int> pTypes = dataBase.solidParticleTypes();

  CHECK(position.numFields() == dataBase.numSolidNodeLists());
  CHECK(mass.numFields() == dataBase.numSolidNodeLists());
  CHECK(massDensity.numFields() == dataBase.numSolidNodeLists());
  CHECK(Hfield.numFields() == dataBase.numSolidNodeLists());
  CHECK(velocity.numFields() == dataBase.numSolidNodeLists());
  CHECK(specificThermalEnergy.numFields() == dataBase.numSolidNodeLists());
  CHECK(deviatoricStress.numFields() == dataBase.numSolidNodeLists());
  CHECK(plasticStrain.numFields() == dataBase.numSolidNodeLists());
  CHECK(damage.numFields() == dataBase.numSolidNodeLists());
  CHECK(fragIDs.numFields() == dataBase.numSolidNodeLists());
  CHECK(pTypes.numFields() == dataBase.numSolidNodeLists());
  
  auto nodeListi = 0;
  for (auto itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    *mPlasticStrain0[nodeListi] = (*itr)->plasticStrain();
    (*mPlasticStrain0[nodeListi]).name(SolidFieldNames::plasticStrain + "0");
  }

  // Create the local storage.
  dataBase.resizeSolidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeSolidFieldList(mVolume, 0.0, HydroFieldNames::volume, false);
  dataBase.resizeSolidFieldList(mPressure, 0.0, HydroFieldNames::pressure, false);
  dataBase.resizeSolidFieldList(mSoundSpeed, 0.0, HydroFieldNames::soundSpeed, false);
  dataBase.resizeSolidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeSolidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeSolidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);
  dataBase.resizeSolidFieldList(mPlasticStrain0, 0.0, SolidFieldNames::plasticStrain + "0", false);
  dataBase.resizeSolidFieldList(mRawPressure, 0.0, FSIFieldNames::rawPressure, false);
  dataBase.resizeSolidFieldList(mInterfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals, false);
  dataBase.resizeSolidFieldList(mInterfaceFraction, 0.0, FSIFieldNames::interfaceFraction, false); 
  dataBase.resizeSolidFieldList(mInterfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness, false);
  dataBase.resizeSolidFieldList(mInverseEquivalentDeviatoricStress, 0.0, FSIFieldNames::inverseEquivalentDeviatoricStress, false);

  // Register the deviatoric stress and plastic strain to be evolved.
  auto positionPolicy = make_shared<IncrementFieldList<Dimension, Vector>>();
  auto deviatoricStressPolicy = make_shared<DeviatoricStressPolicy<Dimension>>();
  auto bulkModulusPolicy = make_shared<BulkModulusPolicy<Dimension>>();
  auto shearModulusPolicy = make_shared<ShearModulusPolicy<Dimension>>();
  auto yieldStrengthPolicy = make_shared<YieldStrengthPolicy<Dimension>>();
  auto csPolicy = make_shared<StrengthSoundSpeedPolicy<Dimension>>();
  auto Ppolicy = make_shared<DamagedPressurePolicy<Dimension>>();
  auto plasticStrainPolicy = make_shared<PairwisePlasticStrainPolicy<Dimension>>();
  auto rawPressurePolicy = make_shared<PressurePolicy<Dimension>>();
  auto interfaceNormalsPolicy = make_shared<PureReplaceFieldList<Dimension,Vector>>(ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals);
  auto interfaceFractionPolicy = make_shared<PureReplaceFieldList<Dimension,Scalar>>(ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceFraction);
  auto interfaceSmoothnessPolicy = make_shared<PureReplaceFieldList<Dimension,Scalar>>(ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness);
  auto velocityPolicy = make_shared<IncrementFieldList<Dimension, Vector>>(HydroFieldNames::position,
                                                                           HydroFieldNames::specificThermalEnergy,
                                                                           true);

  if(this->compatibleEnergyEvolution()){
    auto  thermalEnergyPolicy = make_shared<CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>>(dataBase);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }else if (mEvolveTotalEnergy) {
    auto  thermalEnergyPolicy = make_shared<SpecificFromTotalThermalEnergyPolicy<Dimension>>();
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  } else {
    auto  thermalEnergyPolicy = make_shared<IncrementFieldList<Dimension, Scalar>>();
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }

  std::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > rhoPolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  std::shared_ptr<CompositeFieldListPolicy<Dimension, SymTensor> > Hpolicy(new CompositeFieldListPolicy<Dimension, SymTensor>());
  for (typename DataBase<Dimension>::SolidNodeListIterator itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
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

  state.enroll(position,             positionPolicy);
  state.enroll(velocity,             velocityPolicy);
  state.enroll(massDensity,          rhoPolicy);
  state.enroll(deviatoricStress,     deviatoricStressPolicy);
  state.enroll(Hfield,               Hpolicy);
  state.enroll(mSoundSpeed,          csPolicy);
  state.enroll(mPressure,            Ppolicy);
  state.enroll(mRawPressure,         rawPressurePolicy);
  state.enroll(mBulkModulus,         bulkModulusPolicy);
  state.enroll(mShearModulus,        shearModulusPolicy);
  state.enroll(mYieldStrength,       yieldStrengthPolicy);
  state.enroll(plasticStrain,        plasticStrainPolicy);
  state.enroll(mInterfaceNormals,    interfaceNormalsPolicy); 
  state.enroll(mInterfaceFraction,   interfaceFractionPolicy);
  state.enroll(mInterfaceSmoothness, interfaceSmoothnessPolicy); 

  state.enroll(mTimeStepMask);
  state.enroll(mass);
  state.enroll(mInverseEquivalentDeviatoricStress);
  state.enroll(damage);
  state.enroll(fragIDs);
  state.enroll(pTypes);
  state.enroll(mPlasticStrain0);

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

  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mXSPHWeightSum, 0.0, HydroFieldNames::XSPHWeightSum, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDPDx, Vector::zero, FSIFieldNames::pressureGradient, false);
  dataBase.resizeFluidFieldList(mDepsDx, Vector::zero, FSIFieldNames::specificThermalEnergyGradient, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mM, Tensor::zero, HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mLocalM, Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  dataBase.resizeFluidFieldList(mNormalization, 0.0, HydroFieldNames::normalization, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mNewInterfaceNormals, Vector::zero,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mSmoothedInterfaceNormals, Vector::zero,  FSIFieldNames::smoothedInterfaceNormals,false);
  dataBase.resizeFluidFieldList(mNewInterfaceFraction, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() +  FSIFieldNames::interfaceFraction,false); 
  dataBase.resizeFluidFieldList(mNewInterfaceSmoothness, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness,false);

  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }

  CHECK(not derivs.registered(mDvDt));

  derivs.enrollAny(HydroFieldNames::pairAccelerations, mPairAccelerations);
  derivs.enrollAny(HydroFieldNames::pairWork,          mPairDepsDt);

  derivs.enroll(plasticStrainRate);
  derivs.enroll(mXSPHDeltaV);
  derivs.enroll(mXSPHWeightSum);
  derivs.enroll(mDvDt);
  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDdeviatoricStressDt);
  derivs.enroll(mDHDt);
  derivs.enroll(mHideal);
  derivs.enroll(mDPDx);
  derivs.enroll(mDepsDx);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enroll(mM);
  derivs.enroll(mLocalM);
  derivs.enroll(mMaxViscousPressure);
  derivs.enroll(mNormalization);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mNewInterfaceNormals);
  derivs.enroll(mSmoothedInterfaceNormals);
  derivs.enroll(mNewInterfaceFraction);
  derivs.enroll(mNewInterfaceSmoothness);

  TIME_END("SolidFSISPHregisterDerivs");
}

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
  if (mApplySelectDensitySum){
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
      const auto& W = this->kernel();
            auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeHWeightedFSISPHSumMassDensity(connectivityMap, W, mSumDensityNodeLists, position, mass, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
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
  FieldList<Dimension, Scalar> invSqrtJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);
  FieldList<Dimension, Scalar> interfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
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
    (*boundaryItr)->applyFieldListGhostBoundary(invSqrtJ2);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceFraction);
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
  FieldList<Dimension, Scalar> invSqrtJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);
  FieldList<Dimension, Scalar> interfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
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
    (*boundaryItr)->enforceFieldListBoundary(invSqrtJ2);
    (*boundaryItr)->enforceFieldListBoundary(interfaceFraction);
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
  file.write(mPressure, pathName + "/pressure");
  file.write(mRawPressure, pathName + "/rawEosPressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mVolume, pathName + "/volume");
  file.write(mBulkModulus, pathName + "/bulkModulus");
  file.write(mShearModulus, pathName + "/shearModulus");
  file.write(mYieldStrength, pathName + "/yieldStrength");
  file.write(mPlasticStrain0, pathName + "/plasticStrain0");
  file.write(mInverseEquivalentDeviatoricStress, pathName + "/inverseEquivalentDeviatoricStress");
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mHideal, pathName + "/Hideal");
  file.write(mDPDx, pathName + "/DpDx");
  file.write(mDepsDx, pathName + "/DepsDx");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mM, pathName + "/M");
  file.write(mLocalM, pathName + "/localM");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.write(mNormalization, pathName + "/normalization");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mInterfaceNormals, pathName + "/interfaceNormals");
  file.write(mInterfaceFraction, pathName + "/interfaceFraction");
  file.write(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.write(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.write(mSmoothedInterfaceNormals, pathName + "/smoothedInterfaceNormals");
  file.write(mNewInterfaceFraction, pathName + "/newInterfaceFraction");
  file.write(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mRawPressure, pathName + "/rawEosPressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mVolume, pathName + "/volume");
  file.read(mBulkModulus, pathName + "/bulkModulus");
  file.read(mShearModulus, pathName + "/shearModulus");
  file.read(mYieldStrength, pathName + "/yieldStrength");
  file.read(mPlasticStrain0, pathName + "/plasticStrain0");
  file.read(mInverseEquivalentDeviatoricStress, pathName + "/inverseEquivalentDeviatoricStress");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mHideal, pathName + "/Hideal");
  file.read(mDPDx, pathName + "/DpDx");
  file.read(mDepsDx, pathName + "/DepsDx");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mM, pathName + "/M");
  file.read(mLocalM, pathName + "/localM");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.read(mNormalization, pathName + "/normalization");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mInterfaceNormals, pathName + "/interfaceNormals");
  file.read(mInterfaceFraction, pathName + "/interfaceFraction");
  file.read(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.read(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.read(mSmoothedInterfaceNormals, pathName + "/smoothedInterfaceNormals");
  file.read(mNewInterfaceFraction, pathName + "/newInterfaceFraction");
  file.read(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");

  // For backwards compatibility on change 3597 -- drop in the near future.
  for (auto DvDtPtr: mDvDt) DvDtPtr->name(HydroFieldNames::hydroAcceleration);

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



