//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "SolidCRKSPHHydroBase.hh"
#include "CRKSPHHydroBase.hh"
#include "CRKSPHUtilities.hh"
#include "volumeSpacing.hh"
#include "computeCRKSPHSumVolume.hh"
#include "computeCRKSPHMoments.hh"
#include "computeCRKSPHCorrections.hh"
#include "computeCRKSPHSumMassDensity.hh"
#include "computeSolidCRKSPHSumMassDensity.hh"
#include "gradientCRKSPH.hh"
#include "Physics/GenericHydro.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "Strength/SolidNodeList.hh"
#include "Strength/DeviatoricStressPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
#include "Strength/PlasticStrainPolicy.hh"
#include "Strength/ShearModulusPolicy.hh"
#include "Strength/YieldStrengthPolicy.hh"
#include "Strength/StrengthSoundSpeedPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "FileIO/FileIO.hh"
#include "SolidSPH/DamagedNodeCouplingWithFrags.hh"
#include "SolidMaterial/SolidEquationOfState.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using SolidMaterial::SolidNodeList;
using SolidMaterial::SolidEquationOfState;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Compute the artificial tensile stress correction tensor for the given 
// stress tensor
//------------------------------------------------------------------------------
inline
Dim<1>::SymTensor
tensileStressCorrection(const Dim<1>::SymTensor& sigma) {
  if (sigma.xx() > 0.0) {
    return -sigma;
  } else {
    return Dim<1>::SymTensor();
  }
}

inline
Dim<2>::SymTensor
tensileStressCorrection(const Dim<2>::SymTensor& sigma) {
  const EigenStruct<2> eigen = sigma.eigenVectors();
  const double lambdax = eigen.eigenValues.x();
  const double lambday = eigen.eigenValues.y();
  Dim<2>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
                           0.0, (lambday > 0.0 ? -lambday : 0.0));
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
  Dim<3>::SymTensor result((lambdax > 0.0 ? -lambdax : 0.0), 0.0, 0.0,
                           0.0, (lambday > 0.0 ? -lambday : 0.0), 0.0,
                           0.0, 0.0, (lambdaz > 0.0 ? -lambdaz : 0.0));
  result.rotationalTransform(eigen.eigenVectors);
  return result;
}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidCRKSPHHydroBase<Dimension>::
SolidCRKSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                     ArtificialViscosity<Dimension>& Q,
                     const TableKernel<Dimension>& W,
                     const TableKernel<Dimension>& WPi,
                     const double filter,
                     const double cfl,
                     const bool useVelocityMagnitudeForDt,
                     const bool compatibleEnergyEvolution,
                     const bool XSPH,
                     const PhysicsSpace::MassDensityType densityUpdate,
                     const PhysicsSpace::HEvolutionType HUpdate,
                     const CRKSPHSpace::CRKOrder correctionOrder,
                     const CRKSPHSpace::CRKVolumeType volumeType,
                     const double epsTensile,
                     const double nTensile):
  CRKSPHHydroBase<Dimension>(smoothingScaleMethod, 
                             Q,
                             W,
                             WPi,
                             filter,
                             cfl,
                             useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution,
                             XSPH,
                             densityUpdate,
                             HUpdate,
                             correctionOrder,
                             volumeType,
                             epsTensile,
                             nTensile),
  mDdeviatoricStressDt(FieldSpace::Copy),
  mBulkModulus(FieldSpace::Copy),
  mShearModulus(FieldSpace::Copy),
  mYieldStrength(FieldSpace::Copy),
  mPlasticStrain0(FieldSpace::Copy),
  mFragIDs(FieldSpace::Reference),
  mAdamage(FieldSpace::Copy),
  mBdamage(FieldSpace::Copy),
  mCdamage(FieldSpace::Copy),
  mGradAdamage(FieldSpace::Copy),
  mGradBdamage(FieldSpace::Copy),
  mGradCdamage(FieldSpace::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SolidCRKSPHHydroBase<Dimension>::
~SolidCRKSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Call the ancestor.
  CRKSPHHydroBase<Dimension>::initializeProblemStartup(dataBase);

  // Create storage for the state we're holding.
  mDdeviatoricStressDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
  mBulkModulus = dataBase.newFluidFieldList(0.0, SolidFieldNames::bulkModulus);
  mShearModulus = dataBase.newFluidFieldList(0.0, SolidFieldNames::shearModulus);
  mYieldStrength = dataBase.newFluidFieldList(0.0, SolidFieldNames::yieldStrength);
  mPlasticStrain0 = dataBase.newFluidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
  mAdamage = dataBase.newFluidFieldList(0.0,              HydroFieldNames::A_CRKSPH + " damage");
  mBdamage = dataBase.newFluidFieldList(Vector::zero,     HydroFieldNames::B_CRKSPH + " damage");
  mCdamage = dataBase.newFluidFieldList(Tensor::zero,     HydroFieldNames::C_CRKSPH + " damage");
  mGradAdamage = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::gradA_CRKSPH + " damage");
  mGradBdamage = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::gradB_CRKSPH + " damage");
  mGradCdamage = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::gradC_CRKSPH + " damage");

  // Copy the standard CRK corrections to the damage forms.
  mAdamage.assignFields(this->A());
  mBdamage.assignFields(this->B());
  mGradAdamage.assignFields(this->gradA());
  mGradBdamage.assignFields(this->gradB());
  if (this->correctionOrder() == QuadraticOrder) {
    mCdamage.assignFields(this->C());
    mGradCdamage.assignFields(this->gradC());
  }

  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::SolidNodeListIterator itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {

    // Add the NodeList fragment IDs to our local FieldList.
    mFragIDs.appendField((*itr)->fragmentIDs());

    // Set the moduli.
    (*itr)->bulkModulus(*mBulkModulus[nodeListi]);
    (*itr)->shearModulus(*mShearModulus[nodeListi]);
    (*itr)->yieldStrength(*mYieldStrength[nodeListi]);
  }
}


//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Invoke CRKSPHHydro's state.
  CRKSPHHydroBase<Dimension>::registerState(dataBase, state);

  // Create the local storage.
  dataBase.resizeFluidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeFluidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeFluidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);
  dataBase.resizeFluidFieldList(mPlasticStrain0, 0.0, SolidFieldNames::plasticStrain + "0", false);

  dataBase.resizeFluidFieldList(mAdamage,     0.0,             HydroFieldNames::A_CRKSPH + " damage", false);
  dataBase.resizeFluidFieldList(mBdamage,     Vector::zero,    HydroFieldNames::B_CRKSPH + " damage", false);
  dataBase.resizeFluidFieldList(mGradAdamage, Vector::zero,    HydroFieldNames::gradA_CRKSPH + " damage", false);
  dataBase.resizeFluidFieldList(mGradBdamage, Tensor::zero,    HydroFieldNames::gradB_CRKSPH + " damage", false);
  if (this->correctionOrder() == QuadraticOrder) {
    dataBase.resizeFluidFieldList(mCdamage,     Tensor::zero,    HydroFieldNames::C_CRKSPH + " damage", false);
    dataBase.resizeFluidFieldList(mGradCdamage, ThirdRankTensor::zero, HydroFieldNames::gradC_CRKSPH + " damage", false);
  }

  // Grab the normal Hydro's registered version of the sound speed.
  FieldList<Dimension, Scalar> cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  CHECK(cs.numFields() == dataBase.numFluidNodeLists());

  // Build the FieldList versions of our state.
  FieldList<Dimension, SymTensor> S, D;
  FieldList<Dimension, Scalar> ps;
  FieldList<Dimension, Vector> gradD;
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::SolidNodeListIterator itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    S.appendField((*itr)->deviatoricStress());
    ps.appendField((*itr)->plasticStrain());
    D.appendField((*itr)->effectiveDamage());
    gradD.appendField((*itr)->damageGradient());

    // Make a copy of the beginning plastic strain.
    *mPlasticStrain0[nodeListi] = (*itr)->plasticStrain();
    (*mPlasticStrain0[nodeListi]).name(SolidFieldNames::plasticStrain + "0");
  }

  // Register the deviatoric stress and plastic strain to be evolved.
  PolicyPointer deviatoricStressPolicy(new DeviatoricStressPolicy<Dimension>());
  PolicyPointer plasticStrainPolicy(new PlasticStrainPolicy<Dimension>());
  state.enroll(S, deviatoricStressPolicy);
  state.enroll(ps, plasticStrainPolicy);

  // Register the bulk modulus, shear modulus, and yield strength.
  PolicyPointer bulkModulusPolicy(new BulkModulusPolicy<Dimension>());
  PolicyPointer shearModulusPolicy(new ShearModulusPolicy<Dimension>());
  PolicyPointer yieldStrengthPolicy(new YieldStrengthPolicy<Dimension>());
  state.enroll(mBulkModulus, bulkModulusPolicy);
  state.enroll(mShearModulus, shearModulusPolicy);
  state.enroll(mYieldStrength, yieldStrengthPolicy);

  // Override the policy for the sound speed.
  PolicyPointer csPolicy(new StrengthSoundSpeedPolicy<Dimension>());
  state.enroll(cs, csPolicy);

  // Register the effective damage and damage gradient with default no-op updates.
  // If there are any damage models running they can override these choices.
  state.enroll(D);
  state.enroll(gradD);

  // Register the fragment IDs.
  state.enroll(mFragIDs);

  // And finally the intial plastic strain.
  state.enroll(mPlasticStrain0);

  // Register the CRKSPH correction fields.
  // We deliberately make these non-dynamic here.  This corrections are computed
  // during CRKSPHHydroBase::initialize, not as part of our usual state update.
  // This is necessary 'cause we need boundary conditions *and* the current set of
  // neighbors before we compute these suckers.
  state.enroll(mAdamage);
  state.enroll(mBdamage);
  state.enroll(mCdamage);
  state.enroll(mGradAdamage);
  state.enroll(mGradBdamage);
  state.enroll(mGradCdamage);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  // Call the ancestor method.
  CRKSPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  const string DSDtName = IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress;
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, DSDtName, false);

  derivs.enroll(mDdeviatoricStressDt);

  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    SolidNodeList<Dimension>* solidNodeListPtr = dynamic_cast<SolidNodeList<Dimension>*>(*itr);
    CHECK(solidNodeListPtr != 0);
    derivs.enroll(solidNodeListPtr->plasticStrainRate());
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro before evaluating derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // The fluid CRK bit.
  // Compute the kernel correction fields.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, SymTensor> D = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const FieldList<Dimension, Vector> gradD = state.fields(SolidFieldNames::damageGradient, Vector::zero);
  const FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
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

  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  FieldList<Dimension, ThirdRankTensor> gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
  FieldList<Dimension, Scalar> Adamage = state.fields(HydroFieldNames::A_CRKSPH + " damage", 0.0);
  FieldList<Dimension, Vector> Bdamage = state.fields(HydroFieldNames::B_CRKSPH + " damage", Vector::zero);
  FieldList<Dimension, Tensor> Cdamage = state.fields(HydroFieldNames::C_CRKSPH + " damage", Tensor::zero);
  FieldList<Dimension, Vector> gradAdamage = state.fields(HydroFieldNames::gradA_CRKSPH + " damage", Vector::zero);
  FieldList<Dimension, Tensor> gradBdamage = state.fields(HydroFieldNames::gradB_CRKSPH + " damage", Tensor::zero);
  FieldList<Dimension, ThirdRankTensor> gradCdamage = state.fields(HydroFieldNames::gradC_CRKSPH + " damage", ThirdRankTensor::zero);
  DamagedNodeCouplingWithFrags<Dimension> nodeCoupling(D, gradD, H, fragIDs);
  
  // Compute the volume per node.
  // Change CRKSPH weights here if need be!
  FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  if (this->volumeType() == MassOverDensity) {
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    vol.assignFields(mass/massDensity);
  } else if (this->volumeType() == SumVolume) {
    computeCRKSPHSumVolume(connectivityMap, W, position, H, vol);
  } else {
    VERIFY2(false, "Unknown CRK volume weighting.");
  }

  // We need boundary conditions enforced on the volume before we can compute corrections.
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->applyFieldListGhostBoundary(vol);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  computeCRKSPHMoments(connectivityMap, W, vol, position, H, this->correctionOrder(), NodeCoupling(), m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
  computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, this->correctionOrder(), A, B, C, gradA, gradB, gradC);
  computeCRKSPHMoments(connectivityMap, W, vol, position, H, this->correctionOrder(), nodeCoupling, m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4);
  computeCRKSPHCorrections(m0, m1, m2, m3, m4, gradm0, gradm1, gradm2, gradm3, gradm4, this->correctionOrder(), Adamage, Bdamage, Cdamage, gradAdamage, gradBdamage, gradCdamage);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(A);
    (*boundItr)->applyFieldListGhostBoundary(B);
    (*boundItr)->applyFieldListGhostBoundary(C);
    (*boundItr)->applyFieldListGhostBoundary(gradA);
    (*boundItr)->applyFieldListGhostBoundary(gradB);
    (*boundItr)->applyFieldListGhostBoundary(gradC);
    (*boundItr)->applyFieldListGhostBoundary(Adamage);
    (*boundItr)->applyFieldListGhostBoundary(Bdamage);
    (*boundItr)->applyFieldListGhostBoundary(Cdamage);
    (*boundItr)->applyFieldListGhostBoundary(gradAdamage);
    (*boundItr)->applyFieldListGhostBoundary(gradBdamage);
    (*boundItr)->applyFieldListGhostBoundary(gradCdamage);
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
SolidCRKSPHHydroBase<Dimension>::
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
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const double tiny = 1.0e-30;
  const bool compatibleEnergy = this->compatibleEnergyEvolution();
  const bool XSPH = this->XSPH();
  const Scalar epsTensile = this->epsilonTensile();
  const CRKOrder order = this->correctionOrder();

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const FieldList<Dimension, SymTensor> damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const FieldList<Dimension, Vector> gradDamage = state.fields(SolidFieldNames::damageGradient, Vector::zero);
  const FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  const FieldList<Dimension, Tensor> C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  const FieldList<Dimension, ThirdRankTensor> gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
  const FieldList<Dimension, Scalar> Adamage = state.fields(HydroFieldNames::A_CRKSPH + " damage", 0.0);
  const FieldList<Dimension, Vector> Bdamage = state.fields(HydroFieldNames::B_CRKSPH + " damage", Vector::zero);
  const FieldList<Dimension, Tensor> Cdamage = state.fields(HydroFieldNames::C_CRKSPH + " damage", Tensor::zero);
  const FieldList<Dimension, Vector> gradAdamage = state.fields(HydroFieldNames::gradA_CRKSPH + " damage", Vector::zero);
  const FieldList<Dimension, Tensor> gradBdamage = state.fields(HydroFieldNames::gradB_CRKSPH + " damage", Tensor::zero);
  const FieldList<Dimension, ThirdRankTensor> gradCdamage = state.fields(HydroFieldNames::gradC_CRKSPH + " damage", ThirdRankTensor::zero);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(gradDamage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(A.size() == numNodeLists);
  CHECK(B.size() == numNodeLists);
  CHECK(C.size() == numNodeLists or order != QuadraticOrder);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists);
  CHECK(gradC.size() == numNodeLists or order != QuadraticOrder);
  CHECK(Adamage.size() == numNodeLists);
  CHECK(Bdamage.size() == numNodeLists);
  CHECK(Cdamage.size() == numNodeLists or order != QuadraticOrder);
  CHECK(gradAdamage.size() == numNodeLists);
  CHECK(gradBdamage.size() == numNodeLists);
  CHECK(gradCdamage.size() == numNodeLists or order != QuadraticOrder);

  // Derivative FieldLists.
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, Scalar> effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  FieldList<Dimension, Scalar> viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  FieldList<Dimension, SymTensor> DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
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
  CHECK(DSDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Some scratch variables.
  Vector Bi = Vector::zero, Bj = Vector::zero, Bdami = Vector::zero, Bdamj = Vector::zero;
  Tensor Ci = Tensor::zero, Cj = Tensor::zero, Cdami = Tensor::zero, Cdamj = Tensor::zero;
  Tensor gradBi = Tensor::zero, gradBj = Tensor::zero, gradBdami = Tensor::zero, gradBdamj = Tensor::zero;
  ThirdRankTensor gradCi = ThirdRankTensor::zero, gradCj = ThirdRankTensor::zero, gradCdami = ThirdRankTensor::zero, gradCdamj = ThirdRankTensor::zero;

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::ConstSolidNodeListIterator itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    const SolidNodeList<Dimension>& nodeList = **itr;
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

    // Build the functor we use to compute the effective coupling between nodes.
    DamagedNodeCouplingWithFrags<Dimension> coupling(damage, gradDamage, H, fragIDs);

    // Check if we can identify a reference density.
    Scalar rho0 = 0.0;
    try {
      rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(nodeList.equationOfState()).referenceDensity();
      // cerr << "Setting reference density to " << rho0 << endl;
    } catch(...) {
      // cerr << "BLAGO!" << endl;
    }

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
      const SymTensor& Si = S(nodeListi, i);
      const Scalar& mui = mu(nodeListi, i);
      const Scalar Ai = A(nodeListi, i);
      const Vector& gradAi = gradA(nodeListi, i);
      const Scalar Adami = Adamage(nodeListi, i);
      const Vector& gradAdami = gradAdamage(nodeListi, i);
      if (order != ZerothOrder) {
        Bi = B(nodeListi, i);
        gradBi = gradB(nodeListi, i);
        Bdami = Bdamage(nodeListi, i);
        gradBdami = gradBdamage(nodeListi, i);
      }
      if (order == QuadraticOrder) {
        Ci = C(nodeListi, i);
        gradCi = gradC(nodeListi, i);
        Cdami = Cdamage(nodeListi, i);
        gradCdami = gradCdamage(nodeListi, i);
      }
      const Scalar Hdeti = Hi.Determinant();
      const Scalar weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Ai > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(weighti > 0.0);

      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      Scalar& effViscousPressurei = effViscousPressure(nodeListi, i);
      Scalar& viscousWorki = viscousWork(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      SymTensor& DSDti = DSDt(nodeListi, i);
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
              const Scalar& Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar& cj = soundSpeed(nodeListj, j);
              const Scalar Aj = A(nodeListj, j);
              const Vector& gradAj = gradA(nodeListj, j);
              const Scalar Adamj = Adamage(nodeListj, j);
              const Vector& gradAdamj = gradAdamage(nodeListj, j);
              if (order != ZerothOrder) {
                Bj = B(nodeListj, j);
                gradBj = gradB(nodeListj, j);
                Bdamj = Bdamage(nodeListj, j);
                gradBdamj = gradBdamage(nodeListj, j);
              }
              if (order == QuadraticOrder) {
                Cj = C(nodeListj, j);
                gradCj = gradC(nodeListj, j);
                Cdamj = Cdamage(nodeListj, j);
                gradCdamj = gradCdamage(nodeListj, j);
              }
              const SymTensor& Sj = S(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);
              CHECK(weightj > 0.0);

              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              Scalar& effViscousPressurej = effViscousPressure(nodeListj, j);
              Scalar& viscousWorkj = viscousWork(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
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
              const Vector vij = vi - vj;

              // Symmetrized kernel weight and gradient.
              Scalar gWi, gWj, Wi, Wj, gWdami, gWdamj, Wdami, Wdamj;
              Vector gradWi, gradWj, gradWdami, gradWdamj;
              CRKSPHKernelAndGradient(W, CRKSPHHydroBase<Dimension>::correctionOrder(),  rij, -etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, Wj, gWj, gradWj);
              CRKSPHKernelAndGradient(W, CRKSPHHydroBase<Dimension>::correctionOrder(), -rij,  etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj, Wi, gWi, gradWi);
              CRKSPHKernelAndGradient(W, CRKSPHHydroBase<Dimension>::correctionOrder(),  rij, -etai, Hi, Hdeti,  etaj, Hj, Hdetj, Adami, Bdami, Cdami, gradAdami, gradBdami, gradCdami, Wdamj, gWdamj, gradWdamj);//Replace with Solid form of quadratic Cdami and gradCdami when implemented
              CRKSPHKernelAndGradient(W, CRKSPHHydroBase<Dimension>::correctionOrder(), -rij,  etaj, Hj, Hdetj, -etai, Hi, Hdeti, Adamj, Bdamj, Cdamj, gradAdamj, gradBdamj, gradCdamj, Wdami, gWdami, gradWdami);
              const Vector deltagrad = gradWj - gradWi;
              const Vector deltagraddam = gradWdamj - gradWdami;
              const Vector gradWSPHi = (Hi*etai.unitVector())*WQ.gradValue(etai.magnitude(), Hdeti);
              const Vector gradWSPHj = (Hj*etaj.unitVector())*WQ.gradValue(etaj.magnitude(), Hdetj);

              // Determine how we're applying damage.
              const Scalar fDeffij = coupling(nodeListi, i, nodeListj, j);
              
              // How different is the damage of these two points?
              // const Scalar Di = damage(nodeListi, i).Trace() / Dimension::nDim;
              // const Scalar Dj = damage(nodeListj, j).Trace() / Dimension::nDim;
              // CHECK(Di >= 0.0 and Di <= 1.0);
              // CHECK(Dj >= 0.0 and Dj <= 1.0);
              // const Scalar Ddisc = 1.0 + 8.0*abs(Di - Dj)*safeInv(Di + Dj);
              const Scalar Ddisc = 1.0; // + damage(nodeListi, i).Trace() + damage(nodeListj, j).Trace();

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()/max(tiny, rij2*FastMath::square(Dimension::pownu12(rij2)));
              weightedNeighborSumi += fweightij*std::abs(gWi) * Ddisc;
              weightedNeighborSumj += fweightij*std::abs(gWj) * Ddisc;
              massSecondMomenti += fweightij*gradWSPHi.magnitude2()*thpt;
              massSecondMomentj += gradWSPHj.magnitude2()*thpt;

              // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Scalar workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);
              const Scalar workQj = rhoi*rhoi*QPiij.first .dot(vij).dot(deltagrad);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += weightj * Qi * Wj;
              effViscousPressurej += weighti * Qj * Wi;
              viscousWorki += 0.5*weighti*weightj/mi*workQi;
              viscousWorkj += 0.5*weighti*weightj/mj*workQj;

              // Mass density evolution.
              DrhoDti += fDeffij*rhoi*weightj*vij.dot(gradWdamj);
              DrhoDtj -= fDeffij*rhoj*weighti*vij.dot(gradWdami);

              // Local (damaged) velocity gradient.
              localDvDxi -= fDeffij*weightj*vij*gradWdamj;
              localDvDxj += fDeffij*weighti*vij*gradWdami;

              // We treat positive and negative pressures distinctly, so split 'em up.
              const Scalar Pposi = max(0.0, Pi),
                           Pnegi = min(0.0, Pi),
                           Pposj = max(0.0, Pj),
                           Pnegj = min(0.0, Pj);

              // Compute the stress tensors.
              SymTensor sigmai = -Pnegi*SymTensor::one + Si;
              SymTensor sigmaj = -Pnegj*SymTensor::one + Sj;

              // // Compute the tensile correction to add to the stress as described in 
              // // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
              // const Scalar fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
              // const Scalar fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
              // const SymTensor Ri = fi*tensileStressCorrection(sigmai);
              // const SymTensor Rj = fj*tensileStressCorrection(sigmaj);
              // sigmai += Ri;
              // sigmaj += Rj;

              // Acceleration (CRKSPH form).
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              Vector deltaDvDti, deltaDvDtj;
              const Vector forceij = 0.5*weighti*weightj*((Pposi + Pposj)*deltagrad + ((rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second)*deltagrad) -
                                                          fDeffij*fDeffij*(sigmai + sigmaj)*deltagraddam);                // <- Type III, with CRKSPH Q forces
              // const Vector forceij = 0.5*weighti*weightj*(((rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second)*deltagrad) -
              //                                             fDeffij*fDeffij*(sigmai + sigmaj)*deltagraddam);                // <- Type III, with CRKSPH Q forces
              deltaDvDti = -forceij/mi;
              deltaDvDtj =  forceij/mj;
              DvDti += deltaDvDti;
              DvDtj += deltaDvDtj;
              if (compatibleEnergy) {
                pairAccelerationsi.push_back(deltaDvDti);
                pairAccelerationsj.push_back(deltaDvDtj);
              }

              // Specific thermal energy evolution.
              DepsDti += 0.5*weighti*weightj*(Pposj*vij.dot(deltagrad) + fDeffij*fDeffij*sigmaj.dot(vij).dot(deltagraddam) + workQi)/mi;
              DepsDtj += 0.5*weighti*weightj*(Pposi*vij.dot(deltagrad) + fDeffij*fDeffij*sigmai.dot(vij).dot(deltagraddam) + workQj)/mj;

              // DepsDti += 0.5*weighti*weightj*(fDeffij*fDeffij*sigmaj.dot(vij).dot(deltagraddam) + workQi)/mi;
              // DepsDtj += 0.5*weighti*weightj*(fDeffij*fDeffij*sigmai.dot(vij).dot(deltagraddam) + workQj)/mj;

              // Estimate of delta v (for XSPH).
              if (XSPH and (nodeListi == nodeListj)) {
                XSPHDeltaVi -= fDeffij*weightj*Wdamj*vij;
		XSPHDeltaVj += fDeffij*weighti*Wdami*vij;
              }

            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not this->compatibleEnergyEvolution() or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/max(size_t(1), ncalc);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
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

      // Determine the deviatoric stress evolution.
      const SymTensor deformation = localDvDxi.Symmetric();
      const Tensor spin = localDvDxi.SkewSymmetric();
      const SymTensor deviatoricDeformation = deformation - (deformation.Trace()/3.0)*SymTensor::one;
      const SymTensor spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;

      // In the presence of damage, add a term to reduce the stress on this point.
      const Scalar Di = max(0.0, min(1.0, damage(nodeListi, i).eigenValues().maxElement()));
      DSDti = (1.0 - Di)*DSDti - 0.25/dt*Di*Si;

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.25/dt*Di*(rhoi - rho0);

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
  // // For now make the localDvDx a simple copy of the total.
  // localDvDx.assignFields(DvDx);
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Base class finalization.
  PhysicsSpace::GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (this->densityUpdate() == PhysicsSpace::RigorousSumDensity) {
    const TableKernel<Dimension>& W = this->kernel();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, SymTensor> damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
    const FieldList<Dimension, Vector> gradDamage = state.fields(SolidFieldNames::damageGradient, Vector::zero);
    const FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
    DamagedNodeCouplingWithFrags<Dimension> coupling(damage, gradDamage, H, fragIDs);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    computeCRKSPHSumMassDensity(connectivityMap, this->kernel(), position, mass, H, this->boundaryBegin(), this->boundaryEnd(), massDensity);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // } else if (this->densityUpdate() == PhysicsSpace::SumDensity) {
  //   FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   FieldList<Dimension, Scalar> massDensitySum = derivs.fields(ReplaceFieldList<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
  //                                                               HydroFieldNames::massDensity, 0.0);
  //   massDensity.assignFields(massDensitySum);
  }

  // // Add any filtering component to the node movement.
  // // Note that the FacetedVolumes are in coordinates with the node at the origin already!
  // if (mfilter > 0.0) {
  //   const TableKernel<Dimension>& W = this->kernel();
  //   const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  //   const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  //   const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  //   FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);

  //   // Find the local hulls.
  //   FieldList<Dimension, FacetedVolume> polyvol = dataBase.newFluidFieldList(FacetedVolume(), "faceted volumes");
  //   FieldList<Dimension, Scalar> vol = dataBase.newFluidFieldList(0.0, "volume");
  //   computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, polyvol, vol);

  //   // Displace everyone toward their hull centroids.
  //   const unsigned numNodeLists = position.size();
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const unsigned n = position[nodeListi]->numInternalElements();
  //     for (unsigned i = 0; i != n; ++i) {
  //       const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
  //       if (mag0 > 0.0) {
  //         const Vector deltai = mfilter*polyvol(nodeListi, i).centroid();
  //         const Scalar deltamag = deltai.magnitude();
  //         const Scalar effmag = min(mag0, deltamag);
  //         position(nodeListi, i) += effmag*deltai.unitVector();
  //       }
  //     }
  //   }

  //   // Check for any boundary violations.
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
  //   this->enforceBoundaries(state, derivs);
  // }

  // // Move toward the hull neighbor centroids.
  // if (mfilter > 0.0) {
  //   const TableKernel<Dimension>& W = this->kernel();
  //   const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  //   FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  //   const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  //   const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  //   const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   const FieldList<Dimension, Vector> DrhoDx = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  //   const unsigned numNodeLists = mass.size();
  //   const Scalar W0 = W.kernelValue(0.0, 1.0);
  //   FieldList<Dimension, Vector> delta = dataBase.newFluidFieldList(Vector::zero, "delta position");
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
  //          iItr != connectivityMap.end(nodeListi);
  //          ++iItr) {
  //       const int i = *iItr;
  //       const Vector& ri = position(nodeListi, i);
  //       const Scalar mi = mass(nodeListi, i);
  //       const Scalar rhoi = massDensity(nodeListi, i);
  //       const Vector DrhoDxi = DrhoDx(nodeListi, i);
  //       const SymTensor& Hi = H(nodeListi, i);
  //       const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
  //       const FacetedVolume polyvoli = computeNeighborHull(fullConnectivity, 1.0, ri, Hi, position);
  //       const Vector com = centerOfMass(polyvoli, DrhoDxi);
  //       delta(nodeListi, i) = mfilter*com;
  //     }
  //   }

  //   // Apply the filtering.
  //   const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const unsigned n = position[nodeListi]->numInternalElements();
  //     for (unsigned i = 0; i != n; ++i) {
  //       const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
  //       if (mag0 > 0.0) {
  //         const Scalar deltamag = delta(nodeListi, i).magnitude();
  //         const Scalar effmag = min(mfilter*mag0, deltamag);
  //         position(nodeListi, i) += effmag*delta(nodeListi, i).unitVector();
  //       }
  //     }
  //   }

  //   // Check for any boundary violations.
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
  //   this->enforceBoundaries(state, derivs);
  // }

  // This form looks for points that are too close based on specific volume.
  const double filter = this->filter();
  if (filter > 0.0) {
    const TableKernel<Dimension>& W = this->kernel();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const unsigned numNodeLists = mass.size();
    const Scalar W0 = W.kernelValue(0.0, 1.0);
    FieldList<Dimension, Vector> deltar = dataBase.newFluidFieldList(Vector::zero, "delta position");
    FieldList<Dimension, Scalar> deltav = dataBase.newFluidFieldList(0.0, "delta velocity");
    FieldList<Dimension, Scalar> weightsum = dataBase.newFluidFieldList(0.0, "delta velocity weight sum");
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
           iItr != connectivityMap.end(nodeListi);
           ++iItr) {
        const int i = *iItr;
        const Vector& ri = position(nodeListi, i);
        const Vector& vi = velocity(nodeListi, i);
        const Scalar mi = mass(nodeListi, i);
        const Scalar rhoi = massDensity(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          for (typename vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
               jItr != fullConnectivity[nodeListj].end();
               ++jItr) {
            const unsigned j = *jItr;
            const Vector& rj = position(nodeListj, j);
            const Vector& vj = velocity(nodeListj, j);
            const Scalar mj = mass(nodeListj, j);
            const Scalar rhoj = massDensity(nodeListj, j);
            const Vector rji = rj - ri;
            const Vector rjihat = rji.unitVector();
            const Scalar deltai = 2.0*max(0.0, volumeSpacing<Dimension>(mi/rhoi) + volumeSpacing<Dimension>(mj/rhoj) - rji.magnitude());
            // const Scalar deltai = max(0.0, 2.0*volumeSpacing<Dimension>((mi + mj)/(rhoi + rhoj)) - rji.magnitude());
            // deltar(nodeListi, i) -= deltai*rjihat;
            const Scalar etai = (Hi*rji).magnitude();
            const Scalar weight = W.kernelValue(etai, 1.0)/W0 * (vj - vi).magnitude();
            deltar(nodeListi, i) -= weight*deltai*rjihat;
            weightsum(nodeListi, i) += weight;
            deltav(nodeListi, i) += weight*(vj - vi).magnitude();
          }
        }
      }
    }

    // Apply the filtering.
    const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        // const Scalar hi = 1.0/(H(nodeListi, i).eigenValues().maxElement());
        // const Scalar mag0 = DvDx(nodeListi, i).eigenValues().maxAbsElement()*hi*dt;
        // const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
        const Scalar mag0 = deltav(nodeListi, i)*safeInv(weightsum(nodeListi, i))*dt;
        if (mag0 > 0.0) {
          const Scalar deltamag = deltar(nodeListi, i).magnitude();
          const Scalar effmag = filter*deltamag;
          // const Scalar effmag = filter*min(mag0, deltamag);
          position(nodeListi, i) += effmag*deltar(nodeListi, i).unitVector();
        }
      }
    }

    // Check for any boundary violations.
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
    this->enforceBoundaries(state, derivs);
  }

  // // This form looks uses Voronoi centroids.
  // if (mfilter > 0.0) {
  //   FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  //   const FieldList<Dimension, Vector> centroids = computeVoronoiCentroids(position);
  //   const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  //   const unsigned numNodeLists = position.size();
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const unsigned n = position[nodeListi]->numInternalElements();
  //     for (unsigned i = 0; i != n; ++i) {
  //       const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
  //       if (mag0 > 0.0) {
  //         const Vector delta = centroids(nodeListi, i) - position(nodeListi, i);
  //         const Scalar deltamag = delta.magnitude();
  //         const Scalar effmag = mfilter*min(mfilter*mag0, deltamag);
  //         position(nodeListi, i) += effmag*delta.unitVector();
  //       }
  //     }
  //   }

  //   // Check for any boundary violations.
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
  //   this->enforceBoundaries(state, derivs);
  // }

}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Ancestor method.
  CRKSPHHydroBase<Dimension>::applyGhostBoundaries(state, derivs);

  // Apply boundary conditions to our extra strength variables.
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));

  FieldList<Dimension, Scalar> Adamage = state.fields(HydroFieldNames::A_CRKSPH + " damage", 0.0);
  FieldList<Dimension, Vector> Bdamage = state.fields(HydroFieldNames::B_CRKSPH + " damage", Vector::zero);
  FieldList<Dimension, Tensor> Cdamage = state.fields(HydroFieldNames::C_CRKSPH + " damage", Tensor::zero);
  FieldList<Dimension, Vector> gradAdamage = state.fields(HydroFieldNames::gradA_CRKSPH + " damage", Vector::zero);
  FieldList<Dimension, Tensor> gradBdamage = state.fields(HydroFieldNames::gradB_CRKSPH + " damage", Tensor::zero);
  FieldList<Dimension, ThirdRankTensor> gradCdamage = state.fields(HydroFieldNames::gradC_CRKSPH + " damage", ThirdRankTensor::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(S);
    (*boundaryItr)->applyFieldListGhostBoundary(K);
    (*boundaryItr)->applyFieldListGhostBoundary(mu);
    (*boundaryItr)->applyFieldListGhostBoundary(Y);
    (*boundaryItr)->applyFieldListGhostBoundary(fragIDs);
    (*boundaryItr)->applyFieldListGhostBoundary(Adamage);
    (*boundaryItr)->applyFieldListGhostBoundary(Bdamage);
    (*boundaryItr)->applyFieldListGhostBoundary(Cdamage);
    (*boundaryItr)->applyFieldListGhostBoundary(gradAdamage);
    (*boundaryItr)->applyFieldListGhostBoundary(gradBdamage);
    (*boundaryItr)->applyFieldListGhostBoundary(gradCdamage);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Ancestor method.
  CRKSPHHydroBase<Dimension>::enforceBoundaries(state, derivs);

  // Enforce boundary conditions on the extra strength variable.s
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));

  FieldList<Dimension, Scalar> Adamage = state.fields(HydroFieldNames::A_CRKSPH + " damage", 0.0);
  FieldList<Dimension, Vector> Bdamage = state.fields(HydroFieldNames::B_CRKSPH + " damage", Vector::zero);
  FieldList<Dimension, Tensor> Cdamage = state.fields(HydroFieldNames::C_CRKSPH + " damage", Tensor::zero);
  FieldList<Dimension, Vector> gradAdamage = state.fields(HydroFieldNames::gradA_CRKSPH + " damage", Vector::zero);
  FieldList<Dimension, Tensor> gradBdamage = state.fields(HydroFieldNames::gradB_CRKSPH + " damage", Tensor::zero);
  FieldList<Dimension, ThirdRankTensor> gradCdamage = state.fields(HydroFieldNames::gradC_CRKSPH + " damage", ThirdRankTensor::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(S);
    (*boundaryItr)->enforceFieldListBoundary(K);
    (*boundaryItr)->enforceFieldListBoundary(mu);
    (*boundaryItr)->enforceFieldListBoundary(Y);
    (*boundaryItr)->enforceFieldListBoundary(fragIDs);
    (*boundaryItr)->enforceFieldListBoundary(Adamage);
    (*boundaryItr)->enforceFieldListBoundary(Bdamage);
    (*boundaryItr)->enforceFieldListBoundary(Cdamage);
    (*boundaryItr)->enforceFieldListBoundary(gradAdamage);
    (*boundaryItr)->enforceFieldListBoundary(gradBdamage);
    (*boundaryItr)->enforceFieldListBoundary(gradCdamage);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Ancestor method.
  CRKSPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.write(mBulkModulus, pathName + "/bulkModulus");
  file.write(mShearModulus, pathName + "/shearModulus");
  file.write(mYieldStrength, pathName + "/yieldStrength");
  file.write(mPlasticStrain0, pathName + "/plasticStrain0");
  file.write(mFragIDs, pathName + "/fragIDs");
  file.write(mAdamage, pathName + "/Adamage");
  file.write(mBdamage, pathName + "/Bdamage");
  file.write(mBdamage, pathName + "/Cdamage");
  file.write(mGradAdamage, pathName + "/gradAdamage");
  file.write(mGradBdamage, pathName + "/gradBdamage");
  file.write(mGradBdamage, pathName + "/gradCdamage");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidCRKSPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  // Ancestor method.
  CRKSPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.read(mBulkModulus, pathName + "/bulkModulus");
  file.read(mShearModulus, pathName + "/shearModulus");
  file.read(mYieldStrength, pathName + "/yieldStrength");
  file.read(mPlasticStrain0, pathName + "/plasticStrain0");
  file.read(mFragIDs, pathName + "/fragIDs");
  file.read(mAdamage, pathName + "/Adamage");
  file.read(mBdamage, pathName + "/Bdamage");
  file.read(mBdamage, pathName + "/Cdamage");
  file.read(mGradAdamage, pathName + "/gradAdamage");
  file.read(mGradBdamage, pathName + "/gradBdamage");
  file.read(mGradBdamage, pathName + "/gradCdamage");
}

}
}

