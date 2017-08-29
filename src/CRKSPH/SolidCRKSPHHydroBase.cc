//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The solid CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "CRKSPHHydroBase.hh"
#include "CRKSPHUtilities.hh"
#include "volumeSpacing.hh"
#include "computeHullVolumes.hh"
#include "computeCRKSPHSumVolume.hh"
#include "computeCRKSPHMoments.hh"
#include "computeCRKSPHCorrections.hh"
#include "computeCRKSPHSumMassDensity.hh"
#include "computeSolidCRKSPHSumMassDensity.hh"
#include "gradientCRKSPH.hh"
#include "Physics/GenericHydro.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/entropyWeightingFunction.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
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
#include "SPH/DamagedNodeCouplingWithFrags.hh"
#include "SolidMaterial/SolidEquationOfState.hh"

#include "SolidCRKSPHHydroBase.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using NodeSpace::SolidNodeList;
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
                     const bool evolveTotalEnergy,
                     const bool XSPH,
                     const PhysicsSpace::MassDensityType densityUpdate,
                     const PhysicsSpace::HEvolutionType HUpdate,
                     const CRKSPHSpace::CRKOrder correctionOrder,
                     const CRKSPHSpace::CRKVolumeType volumeType,
                     const bool detectSurfaces,
                     const double detectThreshold,
                     const double sweepAngle,
                     const double detectRange,
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
                             evolveTotalEnergy,
                             XSPH,
                             densityUpdate,
                             HUpdate,
                             correctionOrder,
                             volumeType,
                             detectSurfaces,
                             detectThreshold,
                             sweepAngle,
                             detectRange,
                             epsTensile,
                             nTensile),
  mDdeviatoricStressDt(FieldSpace::FieldStorageType::CopyFields),
  mBulkModulus(FieldSpace::FieldStorageType::CopyFields),
  mShearModulus(FieldSpace::FieldStorageType::CopyFields),
  mYieldStrength(FieldSpace::FieldStorageType::CopyFields),
  mPlasticStrain0(FieldSpace::FieldStorageType::CopyFields),
  mHfield0(FieldSpace::FieldStorageType::CopyFields),
  mFragIDs(FieldSpace::FieldStorageType::ReferenceFields),
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
  mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
  mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
  mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
  mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
  mPlasticStrain0 = dataBase.newSolidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
  mHfield0 = dataBase.newSolidFieldList(SymTensor::zero, HydroFieldNames::H + "0");

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

  // Copy the initial H field to apply to nodes as they become damaged.
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  mHfield0.assignFields(H);
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

  // Grab the normal Hydro's registered version of the sound speed.
  auto cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  CHECK(cs.numFields() == dataBase.numFluidNodeLists());

  // Register the deviatoric stress and plastic strain to be evolved.
  auto ps = dataBase.solidPlasticStrain();
  auto S = dataBase.solidDeviatoricStress();
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
  auto D = dataBase.solidEffectiveDamage();
  auto gradD = dataBase.solidDamageGradient();
  state.enroll(D);
  state.enroll(gradD);

  // Register the fragment IDs.
  state.enroll(mFragIDs);

  // And finally the intial plastic strain.
  mPlasticStrain0 = ps;
  mPlasticStrain0.copyFields();
  state.enroll(mPlasticStrain0);
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

  auto psr = dataBase.solidPlasticStrainRate();
  derivs.enroll(psr);
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
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto XSPH = this->XSPH();
  const auto epsTensile = this->epsilonTensile();
  const auto order = this->correctionOrder();
  const auto correctionMin = this->correctionMin();
  const auto correctionMax = this->correctionMax();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto entropy = state.fields(HydroFieldNames::entropy, Scalar());
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const auto gradDamage = state.fields(SolidFieldNames::damageGradient, Vector::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  const auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  const auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  const auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  const auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(entropy.size() == numNodeLists);
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
  CHECK(C.size() == numNodeLists or order != CRKOrder::QuadraticOrder);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists);
  CHECK(gradC.size() == numNodeLists or order != CRKOrder::QuadraticOrder);
  CHECK(surfacePoint.size() == numNodeLists);

  // Derivative FieldLists.
  auto DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
  auto DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
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
  CHECK(DSDt.size() == numNodeLists);
  CHECK(gradRho.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy) {
    auto nodeListi = 0;
    for (auto itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
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
  for (auto itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    const auto& nodeList = **itr;
    const auto  firstGhostNodei = nodeList.firstGhostNode();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  maxNumNeighbors = nodeList.maxNumNeighbors();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    // Get the work field for this NodeList.
    auto& workFieldi = nodeList.work();

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
    for (auto iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const auto i = *iItr;

      // Prepare to accumulate the time.
      const auto start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  si = entropy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto  mui = mu(nodeListi, i);
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
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Ai > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(weighti > 0.0);

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
      auto& DSDti = DSDt(nodeListi, i);
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
              const auto  sj = entropy(nodeListj, j);
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
              const auto& Sj = S(nodeListj, j);
              const auto  Hdetj = Hj.Determinant();
              const auto  weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
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
              CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, CRKSPHHydroBase<Dimension>::correctionOrder(),  rij,  etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, correctionMin, correctionMax);
              CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, CRKSPHHydroBase<Dimension>::correctionOrder(), -rij, -etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj, correctionMin, correctionMax);
              deltagrad = gradWj - gradWi;
              const auto gradWSPHi = (Hi*etai.unitVector())*W.gradValue(etai.magnitude(), Hdeti);
              const auto gradWSPHj = (Hj*etaj.unitVector())*W.gradValue(etaj.magnitude(), Hdetj);

              // Find the damaged pair weighting scaling.
              const auto fij = coupling(nodeListi, i, nodeListj, j);
              CHECK(fij >= 0.0 and fij <= 1.0);

              // Find the effective weights of i->j and j->i.
              // const auto wi = 2.0*weighti*weightj/(weighti + weightj);
              const auto wij = 0.5*(weighti + weightj);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const auto fweightij = nodeListi == nodeListj ? 1.0 : mj*rhoi/(mi*rhoj);
              const auto rij2 = rij.magnitude2();
              const auto thpt = rij.selfdyad()/max(tiny, rij2*FastMath::square(Dimension::pownu12(rij2)));
              weightedNeighborSumi +=     fweightij*std::abs(gWi);
              weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
              massSecondMomenti +=     fweightij*gradWSPHi.magnitude2()*thpt;
              massSecondMomentj += 1.0/fweightij*gradWSPHj.magnitude2()*thpt;

              // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
              const auto QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                        ri, etai, vi, rhoi, ci, Hi,
                                        rj, etaj, vj, rhoj, cj, Hj);
              const auto Qaccij = (rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second).dot(deltagrad);
              // const auto workQij = 0.5*(vij.dot(Qaccij));
              const auto workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);                // CRK
              const auto workQj = rhoi*rhoi*QPiij.first .dot(vij).dot(deltagrad);                // CRK
              const auto Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const auto Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                                 // We need tighter timestep controls on the Q with CRK
              maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
              effViscousPressurei += wij * Qi * Wj;
              effViscousPressurej += wij * Qj * Wi;
              viscousWorki += 0.5*wij*wij/mi*workQi;
              viscousWorkj += 0.5*wij*wij/mj*workQj;

              // Velocity gradient.
              DvDxi -= wij*vij.dyad(gradWj);
              DvDxj += wij*vij.dyad(gradWi);
              localDvDxi -= fij*wij*vij.dyad(gradWj);
              localDvDxj += fij*wij*vij.dyad(gradWi);

              // Mass density gradient.
              gradRhoi += wij*(rhoj - rhoi)*gradWj;
              gradRhoj += wij*(rhoi - rhoj)*gradWi;

              // We treat positive and negative pressures distinctly, so split 'em up.
              const auto Pposi = max(0.0, Pi),
                         Pnegi = min(0.0, Pi),
                         Pposj = max(0.0, Pj),
                         Pnegj = min(0.0, Pj);

              // Compute the stress tensors.
              SymTensor sigmai, sigmaj;
              if (nodeListi == nodeListj) {
                sigmai = -Pnegi*SymTensor::one + Si;
                sigmaj = -Pnegj*SymTensor::one + Sj;
              }

              // We decide between RK and CRK for the momentum and energy equations based on the surface condition.
              if (surfacePoint(nodeListi, i) == 0) {
                // CRK
                forceij  = 0.5*wij*wij*((Pposi + Pposj)*deltagrad - fij*(sigmai + sigmaj)*deltagrad + Qaccij);   // Type III CRK interpoint force.
                DepsDti += 0.5*wij*wij*(Pposj*vij.dot(deltagrad) + fij*sigmaj.dot(vij).dot(deltagrad) + workQi)/mi;
              } else {
                // RK
                forceij = mi*wij*(((Pposj - Pposi)*gradWj - fij*(sigmaj - sigmai)*gradWj)/rhoi + rhoi*QPiij.first.dot(gradWj));
                DepsDti += wij*rhoi*QPiij.first.dot(vij).dot(gradWj);     // Q term only -- adiabatic portion added later
              }

              if (surfacePoint(nodeListj, j) == 0) {
                // CRK
                forceji  = 0.5*wij*wij*((Pposi + Pposj)*deltagrad - fij*(sigmai + sigmaj)*deltagrad + Qaccij);   // Type III CRK interpoint force.
                DepsDtj += 0.5*wij*wij*(Pposj*vij.dot(deltagrad) + fij*sigmaj.dot(vij).dot(deltagrad) + workQj)/mj;
              } else {
                // RK
                forceji = mj*wij*(((Pposj - Pposi)*gradWi - fij*(sigmaj - sigmai)*gradWi)/rhoj - rhoj*QPiij.second.dot(gradWi));
                DepsDtj -= wij*rhoj*QPiij.second.dot(vij).dot(gradWi);     // Q term only -- adiabatic portion added later
              }

              DvDti -= forceij/mi;
              DvDtj += forceji/mj;
              if (compatibleEnergy) {
                pairAccelerationsi.push_back(-forceij/mi);
                pairAccelerationsj.push_back( forceji/mj);
              }

              // Estimate of delta v (for XSPH).
              XSPHDeltaVi -= fij*wij*Wj*vij;
              XSPHDeltaVj += fij*wij*Wi*vij;
            }
          }
        }
      }
      const auto numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not this->compatibleEnergyEvolution() or NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent() or
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // For a surface point, add the RK thermal energy evolution.
      if (surfacePoint(nodeListi, i) != 0) DepsDti += (Si - Pi*SymTensor::one).doubledot(DvDxi)/rhoi;

      // Get the time for pairwise interactions.
      const auto deltaTimePair = Timing::difference(start, Timing::currentTime())/max(size_t(1), ncalc);

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

      // If this node is damaged we begin to force it back to it's original H.
      const auto Di = max(0.0, min(1.0, damage(nodeListi, i).eigenValues().maxElement()));
      Hideali = (1.0 - Di)*Hideali + Di*mHfield0(nodeListi, i);

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()/3.0)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;

      // In the presence of damage, add a term to reduce the stress on this point.
      // const auto Di = max(0.0, min(1.0, damage(nodeListi, i).eigenValues().maxElement()));
      DSDti = (1.0 - Di)*DSDti - 0.25/dt*Di*Si;

      // Time evolution of the mass density.
      DrhoDti = -rhoi*localDvDxi.Trace();

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.25/dt*Di*(rhoi - rho0);

      // If needed finish the total energy derivative.
      if (this->evolveTotalEnergy()) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

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

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(S);
    (*boundaryItr)->applyFieldListGhostBoundary(K);
    (*boundaryItr)->applyFieldListGhostBoundary(mu);
    (*boundaryItr)->applyFieldListGhostBoundary(Y);
    (*boundaryItr)->applyFieldListGhostBoundary(fragIDs);
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

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(S);
    (*boundaryItr)->enforceFieldListBoundary(K);
    (*boundaryItr)->enforceFieldListBoundary(mu);
    (*boundaryItr)->enforceFieldListBoundary(Y);
    (*boundaryItr)->enforceFieldListBoundary(fragIDs);
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
  file.write(mHfield0, pathName + "/Hfield0");
  file.write(mFragIDs, pathName + "/fragIDs");
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
  file.read(mHfield0, pathName + "/Hfield0");
  file.read(mFragIDs, pathName + "/fragIDs");
}

}
}

