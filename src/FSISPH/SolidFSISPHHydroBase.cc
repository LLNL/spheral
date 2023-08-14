//---------------------------------Spheral++----------------------------------//
// SolidFSISPHHydroBase -- SolidSPHHydro modified for large density discontinuities
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

//#include "SPH/SPHHydroBase.hh"                 
#include "SPH/SPHHydroBase.hh"

#include "GSPH/Policies/PureReplaceFieldList.hh"
#include "NodeList/SolidNodeList.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "SolidMaterial/SolidEquationOfState.hh" 

#include "Hydro/HydroFieldNames.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"

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
  SPHHydroBase<Dimension>(smoothingScaleMethod, 
                          dataBase,
                          Q,
                          W,
                          W,
                          filter,
                          cfl,
                          useVelocityMagnitudeForDt,
                          compatibleEnergyEvolution,
                          evolveTotalEnergy,
                          false,
                          XSPH,
                          correctVelocityGradient,
                          true,
                          densityUpdate,
                          HUpdate,
                          epsTensile,
                          nTensile,
                          xmin,
                          xmax),
  mSlideSurface(slides),
  mSurfaceForceCoefficient(surfaceForceCoefficient),
  mDensityStabilizationCoefficient(densityStabilizationCoefficient),
  mSpecificThermalEnergyDiffusionCoefficient(specificThermalEnergyDiffusionCoefficient),
  mXSPHCoefficient(xsphCoefficient),
  mInterfaceMethod(interfaceMethod),
  mKernelAveragingMethod(kernelAveragingMethod),
  mApplySelectDensitySum(false),
  mSumDensityNodeLists(sumDensityNodeLists),
  mPairDepsDt(),
  mRawPressure(FieldStorageType::CopyFields),
  mDPDx(FieldStorageType::CopyFields),
  mDepsDx(FieldStorageType::CopyFields),
  mInterfaceNormals(FieldStorageType::CopyFields),
  mInterfaceFraction(FieldStorageType::CopyFields),
  mInterfaceSmoothness(FieldStorageType::CopyFields),
  mNewInterfaceNormals(FieldStorageType::CopyFields),
  mSmoothedInterfaceNormals(FieldStorageType::CopyFields),
  mNewInterfaceFraction(FieldStorageType::CopyFields),
  mNewInterfaceSmoothness(FieldStorageType::CopyFields),
  mDdeviatoricStressDt(FieldStorageType::CopyFields),
  mBulkModulus(FieldStorageType::CopyFields),
  mShearModulus(FieldStorageType::CopyFields),
  mYieldStrength(FieldStorageType::CopyFields),
  mPlasticStrain0(FieldStorageType::CopyFields),
  mHfield0(FieldStorageType::CopyFields),
  mInverseEquivalentDeviatoricStress(FieldStorageType::CopyFields){

    mPairDepsDt.clear();

    // see if we're summing density for any nodelist
    auto numNodeLists = dataBase.numNodeLists();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      if (sumDensityNodeLists[nodeListi]==1){
        mApplySelectDensitySum = true;
      } 
    }
    
    mRawPressure = dataBase.newFluidFieldList(0.0, FSIFieldNames::rawPressure);
    mDPDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::pressureGradient);
    mDepsDx = dataBase.newFluidFieldList(Vector::zero, FSIFieldNames::specificThermalEnergyGradient);
    mInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::interfaceNormals);
    mInterfaceFraction = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceFraction);
    mInterfaceSmoothness = dataBase.newFluidFieldList(0.0,  FSIFieldNames::interfaceSmoothness);
    mNewInterfaceNormals = dataBase.newFluidFieldList(Vector::one, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals);
    mSmoothedInterfaceNormals = dataBase.newFluidFieldList(Vector::one,  FSIFieldNames::smoothedInterfaceNormals);
    mNewInterfaceFraction = dataBase.newFluidFieldList(0.0, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceFraction);
    mNewInterfaceSmoothness = dataBase.newFluidFieldList(0.0, ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness);
    
    mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
    mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
    mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
    mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
    mPlasticStrain0 = dataBase.newSolidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
    mHfield0 = dataBase.newSolidFieldList(SymTensor::zero, HydroFieldNames::H + "0");

    mInverseEquivalentDeviatoricStress = dataBase.newFluidFieldList(0.0, FSIFieldNames::inverseEquivalentDeviatoricStress);
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

  SPHHydroBase<Dimension>::initializeProblemStartup(dataBase);

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

  // Copy the initial H field to apply to nodes as they become damaged.
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  mHfield0.assignFields(H);
  
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

  SPHHydroBase<Dimension>::registerState(dataBase,state);
  
  //typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  FieldList<Dimension, SymTensor> deviatoricStress = dataBase.solidDeviatoricStress();
  FieldList<Dimension, Scalar> plasticStrain = dataBase.solidPlasticStrain();
  FieldList<Dimension, SymTensor> damage = dataBase.solidDamage();
  FieldList<Dimension, int> fragIDs = dataBase.solidFragmentIDs();
  FieldList<Dimension, int> pTypes = dataBase.solidParticleTypes();

  FieldList<Dimension, Scalar> cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);

  CHECK(cs.numFields() == dataBase.numFluidNodeLists());
  CHECK(P.numFields() == dataBase.numFluidNodeLists());
  CHECK(deviatoricStress.numFields() == dataBase.numFluidNodeLists());
  CHECK(plasticStrain.numFields() == dataBase.numFluidNodeLists());
  CHECK(damage.numFields() == dataBase.numFluidNodeLists());
  CHECK(fragIDs.numFields() == dataBase.numFluidNodeLists());
  CHECK(pTypes.numFields() == dataBase.numFluidNodeLists());
  CHECK(plasticStrain.numFields() == dataBase.numFluidNodeLists());

  auto nodeListi = 0;
  for (auto itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    *mPlasticStrain0[nodeListi] = (*itr)->plasticStrain();
    (*mPlasticStrain0[nodeListi]).name(SolidFieldNames::plasticStrain + "0");
  }

  // Create the local storage.
  dataBase.resizeFluidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeFluidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeFluidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);
  dataBase.resizeFluidFieldList(mPlasticStrain0, 0.0, SolidFieldNames::plasticStrain + "0", false);
  dataBase.resizeFluidFieldList(mRawPressure, 0.0, FSIFieldNames::rawPressure, false);
  dataBase.resizeFluidFieldList(mInterfaceNormals, Vector::zero, FSIFieldNames::interfaceNormals, false);
  dataBase.resizeFluidFieldList(mInterfaceFraction, 0.0, FSIFieldNames::interfaceFraction, false); 
  dataBase.resizeFluidFieldList(mInterfaceSmoothness, 0.0, FSIFieldNames::interfaceSmoothness, false);
  dataBase.resizeFluidFieldList(mInverseEquivalentDeviatoricStress, 0.0, FSIFieldNames::inverseEquivalentDeviatoricStress, false);

  // Register the deviatoric stress and plastic strain to be evolved.
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
  
  // Override the specific thermal energy policy if compatible
  if(this->compatibleEnergyEvolution()){
    auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
    CHECK(specificThermalEnergy.numFields() == dataBase.numFluidNodeLists());
    auto  epsPolicy = make_shared<CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>>(dataBase);
    state.enroll(specificThermalEnergy, epsPolicy);
  }

  state.enroll(mInverseEquivalentDeviatoricStress);
  state.enroll(plasticStrain, plasticStrainPolicy);
  state.enroll(mInterfaceNormals,interfaceNormalsPolicy); 
  state.enroll(mInterfaceFraction,interfaceFractionPolicy);
  state.enroll(mInterfaceSmoothness,interfaceSmoothnessPolicy); 
  state.enroll(deviatoricStress, deviatoricStressPolicy);
  state.enroll(mBulkModulus, bulkModulusPolicy);
  state.enroll(mShearModulus, shearModulusPolicy);
  state.enroll(mYieldStrength, yieldStrengthPolicy);
  state.enroll(cs, csPolicy);
  state.enroll(P, Ppolicy);
  state.enroll(mRawPressure,rawPressurePolicy);
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

  // Call the ancestor method.
  SPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);  
  
  // database fields
  FieldList<Dimension, Scalar> plasticStrainRate = dataBase.solidPlasticStrainRate();

  // make sure we're tracking the right number of node lists
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress, false);
  dataBase.resizeFluidFieldList(mDPDx, Vector::zero, FSIFieldNames::pressureGradient, false);
  dataBase.resizeFluidFieldList(mDepsDx, Vector::zero, FSIFieldNames::specificThermalEnergyGradient, false);
  dataBase.resizeFluidFieldList(mNewInterfaceNormals, Vector::zero,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceNormals,false);
  dataBase.resizeFluidFieldList(mSmoothedInterfaceNormals, Vector::zero,  FSIFieldNames::smoothedInterfaceNormals,false);
  dataBase.resizeFluidFieldList(mNewInterfaceFraction, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() +  FSIFieldNames::interfaceFraction,false); 
  dataBase.resizeFluidFieldList(mNewInterfaceSmoothness, 0.0,  ReplaceBoundedFieldList<Dimension,Scalar>::prefix() + FSIFieldNames::interfaceSmoothness,false);

  // enroll 
  derivs.enrollAny(HydroFieldNames::pairWork,  mPairDepsDt);
  derivs.enroll(mDdeviatoricStressDt);
  derivs.enroll(plasticStrainRate);
  derivs.enroll(mDPDx);
  derivs.enroll(mDepsDx);
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

  SPHHydroBase<Dimension>::applyGhostBoundaries(state,derivs);
  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  FieldList<Dimension, Scalar> interfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);
  FieldList<Dimension, Scalar> rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);
  //FieldList<Dimension, Scalar> yield = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, Scalar> invSqrtJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(S);
    (*boundaryItr)->applyFieldListGhostBoundary(K);
    (*boundaryItr)->applyFieldListGhostBoundary(mu);
    (*boundaryItr)->applyFieldListGhostBoundary(Y);
    (*boundaryItr)->applyFieldListGhostBoundary(fragIDs);
    (*boundaryItr)->applyFieldListGhostBoundary(pTypes);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceFraction);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceNormals);
    (*boundaryItr)->applyFieldListGhostBoundary(interfaceSmoothness);
    (*boundaryItr)->applyFieldListGhostBoundary(rawPressure);
    //(*boundaryItr)->applyFieldListGhostBoundary(yield);
    (*boundaryItr)->applyFieldListGhostBoundary(invSqrtJ2);
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

  SPHHydroBase<Dimension>::enforceBoundaries(state,derivs);

  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  FieldList<Dimension, Scalar> interfaceFraction = state.fields(FSIFieldNames::interfaceFraction, 0.0);
  FieldList<Dimension, Vector> interfaceNormals = state.fields(FSIFieldNames::interfaceNormals, Vector::zero);
  FieldList<Dimension, Scalar> interfaceSmoothness = state.fields(FSIFieldNames::interfaceSmoothness, 0.0);
  FieldList<Dimension, Scalar> rawPressure = state.fields(FSIFieldNames::rawPressure, 0.0);
  //FieldList<Dimension, Scalar> yield = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, Scalar> invSqrtJ2 = state.fields(FSIFieldNames::inverseEquivalentDeviatoricStress, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(S);
    (*boundaryItr)->enforceFieldListBoundary(K);
    (*boundaryItr)->enforceFieldListBoundary(mu);
    (*boundaryItr)->enforceFieldListBoundary(Y);
    (*boundaryItr)->enforceFieldListBoundary(fragIDs);
    (*boundaryItr)->enforceFieldListBoundary(pTypes);
    (*boundaryItr)->enforceFieldListBoundary(interfaceFraction);
    (*boundaryItr)->enforceFieldListBoundary(interfaceNormals);
    (*boundaryItr)->enforceFieldListBoundary(interfaceSmoothness);
    (*boundaryItr)->enforceFieldListBoundary(rawPressure);
    //(*boundaryItr)->enforceFieldListBoundary(yield);
    (*boundaryItr)->enforceFieldListBoundary(invSqrtJ2);
  }

}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  SPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.write(mBulkModulus, pathName + "/bulkModulus");
  file.write(mShearModulus, pathName + "/shearModulus");
  file.write(mYieldStrength, pathName + "/yieldStrength");
  file.write(mPlasticStrain0, pathName + "/plasticStrain0");
  file.write(mHfield0, pathName + "/Hfield0");

  file.write(mRawPressure, pathName + "/rawEosPressure");
  file.write(mDPDx, pathName + "/DpDx");
  file.write(mDepsDx, pathName + "/DepsDx");
  file.write(mInterfaceNormals, pathName + "/interfaceNormals");
  file.write(mInterfaceFraction, pathName + "/interfaceFraction");
  file.write(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.write(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.write(mSmoothedInterfaceNormals, pathName + "/smoothedInterfaceNormals");
  file.write(mNewInterfaceFraction, pathName + "/newInterfaceFraction");
  file.write(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");

  file.write(mInverseEquivalentDeviatoricStress, pathName + "/inverseEquivalentDeviatoricStress");
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  SPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.read(mBulkModulus, pathName + "/bulkModulus");
  file.read(mShearModulus, pathName + "/shearModulus");
  file.read(mYieldStrength, pathName + "/yieldStrength");
  file.read(mPlasticStrain0, pathName + "/plasticStrain0");
  file.read(mHfield0, pathName + "/Hfield0");

  file.read(mRawPressure, pathName + "/rawEosPressure");
  file.read(mDPDx, pathName + "/DpDx");
  file.read(mDepsDx, pathName + "/DepsDx");
  file.read(mInterfaceNormals, pathName + "/interfaceNormals");
  file.read(mInterfaceFraction, pathName + "/interfaceFraction");
  file.read(mInterfaceSmoothness, pathName + "/interfaceSmoothness");
  file.read(mNewInterfaceNormals, pathName + "/newInterfaceNormals");
  file.read(mSmoothedInterfaceNormals, pathName + "/smoothedInterfaceNormals");
  file.read(mNewInterfaceFraction, pathName + "/newInterfaceFraction");
  file.read(mNewInterfaceSmoothness, pathName + "/newInterfaceSmoothness");

  file.read(mInverseEquivalentDeviatoricStress, pathName + "/inverseEquivalentDeviatoricStress");
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



