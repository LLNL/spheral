//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The solid CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "CRKSPH/CRKSPHHydroBase.hh"
#include "RK/ReproducingKernel.hh"
#include "RK/RKFieldNames.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.hh"
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
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/updateStateFields.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "SolidMaterial/SolidEquationOfState.hh"

#include "SolidCRKSPHHydroBase.hh"

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
                     DataBase<Dimension>& dataBase,
                     ArtificialViscosity<Dimension>& Q,
                     const RKOrder order,
                     const double filter,
                     const double cfl,
                     const bool useVelocityMagnitudeForDt,
                     const bool compatibleEnergyEvolution,
                     const bool evolveTotalEnergy,
                     const bool XSPH,
                     const MassDensityType densityUpdate,
                     const HEvolutionType HUpdate,
                     const double epsTensile,
                     const double nTensile,
                     const bool damageRelieveRubble):
  CRKSPHHydroBase<Dimension>(smoothingScaleMethod, 
                             dataBase,
                             Q,
                             order,
                             filter,
                             cfl,
                             useVelocityMagnitudeForDt,
                             compatibleEnergyEvolution,
                             evolveTotalEnergy,
                             XSPH,
                             densityUpdate,
                             HUpdate,
                             epsTensile,
                             nTensile),
  mDamageRelieveRubble(damageRelieveRubble),
  mDdeviatoricStressDt(FieldStorageType::CopyFields),
  mBulkModulus(FieldStorageType::CopyFields),
  mShearModulus(FieldStorageType::CopyFields),
  mYieldStrength(FieldStorageType::CopyFields),
  mPlasticStrain0(FieldStorageType::CopyFields),
  mHfield0(FieldStorageType::CopyFields) {

  // Create storage for the state we're holding.
  mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
  mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
  mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
  mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
  mPlasticStrain0 = dataBase.newSolidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
  mHfield0 = dataBase.newSolidFieldList(SymTensor::zero, HydroFieldNames::H + "0");
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
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // Call the ancestor.
  CRKSPHHydroBase<Dimension>::initializeProblemStartupDependencies(dataBase, state, derivs);

  // Set the moduli.
  updateStateFields(SolidFieldNames::bulkModulus, state, derivs);
  updateStateFields(SolidFieldNames::shearModulus, state, derivs);
  updateStateFields(SolidFieldNames::yieldStrength, state, derivs);

  // Copy the initial H field to apply to nodes as they become damaged.
  const auto H = dataBase.fluidHfield();
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

  // Invoke CRKSPHHydro's state.
  CRKSPHHydroBase<Dimension>::registerState(dataBase, state);

  // Create the local storage.
  dataBase.resizeFluidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeFluidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeFluidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);

  // Register the deviatoric stress and plastic strain to be evolved.
  auto ps = dataBase.solidPlasticStrain();
  auto S = dataBase.solidDeviatoricStress();
  state.enroll(S, make_policy<DeviatoricStressPolicy<Dimension>>());
  state.enroll(ps, make_policy<PlasticStrainPolicy<Dimension>>());

  // Register the bulk modulus, shear modulus, and yield strength.
  state.enroll(mBulkModulus, make_policy<BulkModulusPolicy<Dimension>>());
  state.enroll(mShearModulus, make_policy<ShearModulusPolicy<Dimension>>());
  state.enroll(mYieldStrength, make_policy<YieldStrengthPolicy<Dimension>>());

  // Register the damage with a default no-op update.
  // If there are any damage models running they can override this choice.
  auto D = dataBase.solidDamage();
  state.enroll(D);

  // Register the fragment IDs.
  auto fragIDs = dataBase.solidFragmentIDs();
  state.enroll(fragIDs);

  // Register the particle types.
  auto pTypes = dataBase.solidParticleTypes();
  state.enroll(pTypes);

  // And finally the intial plastic strain.
  mPlasticStrain0 = ps;
  mPlasticStrain0.copyFields();
  for (auto itr = mPlasticStrain0.begin(); itr != mPlasticStrain0.end(); ++itr) (**itr).name(SolidFieldNames::plasticStrain + "0");
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
  const auto DSDtName = IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress;
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
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto  order = this->correctionOrder();
  const auto& WR = state.template getAny<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(order));
  const auto& WT = WR.kernel();      // Base TableKernel

  // A few useful constants we'll use in the following loop.
  //const double tiny = 1.0e-30;
  const auto  compatibleEnergy = this->compatibleEnergyEvolution();
  const auto  evolveTotalEnergy = this->evolveTotalEnergy();
  const auto  XSPH = this->XSPH();
  const auto  damageRelieveRubble = this->damageRelieveRubble();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

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
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  const auto corrections = state.fields(RKFieldNames::rkCorrections(order), RKCoefficients<Dimension>());
  const auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
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
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);
  CHECK(corrections.size() == numNodeLists);
  CHECK(surfacePoint.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DxDt = derivatives.fields(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massFirstMoment = derivatives.fields(HydroFieldNames::massFirstMoment, Vector::zero);
  auto  massSecondMomentEta = derivatives.fields(HydroFieldNames::massSecondMomentEta, SymTensor::zero);
  auto  massSecondMomentLab = derivatives.fields(HydroFieldNames::massSecondMomentLab, SymTensor::zero);
  auto  DSDt = derivatives.fields(IncrementState<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
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
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massFirstMoment.size() == numNodeLists);
  CHECK(massSecondMomentEta.size() == numNodeLists);
  CHECK(massSecondMomentLab.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy) pairAccelerations.resize(npairs);

  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar etaMagi, etaMagj, fweightij;
    Scalar Wi, Wj;
    Scalar WSPHi, WSPHj, WASPHi, WASPHj;
    Tensor QPiij, QPiji;
    Vector gradWi, gradWj;
    Vector deltagrad, forceij, forceji;
    Vector rij, vij, etai, etaj;
    SymTensor sigmai, sigmaj, rijdyad;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    auto viscousWork_thread = viscousWork.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massFirstMoment_thread = massFirstMoment.threadCopy(threadStack);
    auto massSecondMomentEta_thread = massSecondMomentEta.threadCopy(threadStack);
    auto massSecondMomentLab_thread = massSecondMomentLab.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      //const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto  pTypei = pTypes(nodeListi, i);
      const auto& correctionsi = corrections(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      CONTRACT_VAR(Hdeti);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(weighti > 0.0);

      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& viscousWorki = viscousWork_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massFirstMomenti = massFirstMoment_thread(nodeListi, i);
      auto& massSecondMomentEtai = massSecondMomentEta_thread(nodeListi, i);
      auto& massSecondMomentLabi = massSecondMomentLab_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto  mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      //const auto  epsj = specificThermalEnergy(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto& Sj = S(nodeListj, j);
      const auto  pTypej = pTypes(nodeListj, j);
      const auto& correctionsj = corrections(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      const auto  weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
      CONTRACT_VAR(Hdetj);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);
      CHECK(weightj > 0.0);

      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& viscousWorkj = viscousWork_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massFirstMomentj = massFirstMoment_thread(nodeListj, j);
      auto& massSecondMomentEtaj = massSecondMomentEta_thread(nodeListj, j);
      auto& massSecondMomentLabj = massSecondMomentLab_thread(nodeListj, j);

      // Node displacement.
      rij = ri - rj;
      etai = Hi*rij;
      etaj = Hj*rij;
      vij = vi - vj;
      etaMagi = etai.magnitude();
      etaMagj = etaj.magnitude();

      // Flag if this is a contiguous material pair or not.
      const auto sameMatij = true; // nodeListi == nodeListj; // (nodeListi == nodeListj and fragIDi == fragIDj);

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);

      // Symmetrized kernel weight and gradient.
      std::tie(Wj, gradWj) = WR.evaluateKernelAndGradient( rij, Hj, correctionsi);  // Hj because we compute RK using scatter formalism
      std::tie(Wi, gradWi) = WR.evaluateKernelAndGradient(-rij, Hi, correctionsj);
      deltagrad = gradWj - gradWi;

      // Find the damaged pair weighting scaling.
      const auto fDij = pairs[kk].f_couple;
      CHECK(fDij >= 0.0 and fDij <= 1.0);

      // Moments of the node distribution -- used for the ideal H calculation.
      WSPHi = WT.kernelValueSPH(etaMagi);
      WSPHj = WT.kernelValueSPH(etaMagj);
      WASPHi = WT.kernelValueASPH(etaMagi, nPerh);
      WASPHj = WT.kernelValueASPH(etaMagj, nPerh);
      fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      rijdyad = rij.selfdyad();
      weightedNeighborSumi +=     fweightij*WSPHi;
      weightedNeighborSumj += 1.0/fweightij*WSPHj;
      massFirstMomenti -=     fweightij*WSPHi*etai;
      massFirstMomentj += 1.0/fweightij*WSPHj*etaj;
      massSecondMomentEtai +=     fweightij*WASPHi*etai.selfdyad();
      massSecondMomentEtaj += 1.0/fweightij*WASPHj*etaj.selfdyad();
      massSecondMomentLabi +=     fweightij*WASPHi*rijdyad;
      massSecondMomentLabj += 1.0/fweightij*WASPHj*rijdyad;

      // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
      std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                      ri, etai, vi, rhoi, ci, Hi,
                                      rj, etaj, vj, rhoj, cj, Hj);
      const auto Qaccij = (rhoi*rhoi*QPiij + rhoj*rhoj*QPiji).dot(deltagrad);
      // const auto workQij = 0.5*(vij.dot(Qaccij));
      const auto workQi = rhoj*rhoj*QPiji.dot(vij).dot(deltagrad);                // CRK
      const auto workQj = rhoi*rhoi*QPiij.dot(vij).dot(deltagrad);                // CRK
      // const auto workQVi =  vij.dot((rhoj*rhoj*QPiji).dot(gradWj));               //RK V and RK I Work
      // const auto workQVj =  vij.dot((rhoi*rhoi*QPiij).dot(gradWi));               //RK V and RK I Work
      const auto Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
      const auto Qj = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
      maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                     // We need tighter timestep controls on the Q with CRK
      maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
      effViscousPressurei += weightj * Qi * Wj;
      effViscousPressurej += weighti * Qj * Wi;
      viscousWorki += 0.5*weighti*weightj/mi*workQi;
      viscousWorkj += 0.5*weighti*weightj/mj*workQj;

      // Velocity gradient.
      DvDxi -= fDij * weightj*vij.dyad(gradWj);
      DvDxj += fDij * weighti*vij.dyad(gradWi);
      if (sameMatij) {
        localDvDxi -= fDij * weightj*vij.dyad(gradWj);
        localDvDxj += fDij * weighti*vij.dyad(gradWi);
      }

      // // Mass density gradient.
      // gradRhoi += weightj*(rhoj - rhoi)*gradWj;
      // gradRhoj += weighti*(rhoi - rhoj)*gradWi;

      // Compute the stress tensors.
      if (sameMatij) {
        sigmai = fDij*Si - Pi * SymTensor::one;
        sigmaj = fDij*Sj - Pj * SymTensor::one;
      } else {
        sigmai = -Pi * SymTensor::one;
        sigmaj = -Pj * SymTensor::one;
      }

      // We decide between RK and CRK for the momentum and energy equations based on the surface condition.
      // Momentum
      forceij = (true ? // surfacePoint(nodeListi, i) <= 1 ? 
                 0.5*weighti*weightj*(-(sigmai + sigmaj)*deltagrad + Qaccij) :                                    // Type III CRK interpoint force.
                 -mi*weightj*((sigmaj - sigmai)*gradWj/rhoi + rhoi*QPiij.dot(gradWj)));                           // RK
      forceji = (true ? // surfacePoint(nodeListj, j) <= 1 ?
                 0.5*weighti*weightj*(-(sigmai + sigmaj)*deltagrad + Qaccij) :                                    // Type III CRK interpoint force.
                 -mj*weighti*((sigmaj - sigmai)*gradWi/rhoj - rhoj*QPiji.dot(gradWi)));                           // RK
      if (freeParticle) {
        DvDti -= forceij/mi;
        DvDtj += forceji/mj;
      }
      if (compatibleEnergy) pairAccelerations[kk] = -forceij/mi;                                                  // Acceleration for i (j anti-symmetric)

      // Energy
      DepsDti += (true ? // surfacePoint(nodeListi, i) <= 1 ?
                  0.5*weighti*weightj*(-sigmaj.dot(vij).dot(deltagrad) + workQi)/mi :                             // CRK
                  weightj*rhoi*QPiij.dot(vij).dot(gradWj));                                                       // RK, Q term only -- adiabatic portion added later
      DepsDtj += (true ? // surfacePoint(nodeListj, j) <= 1 ?
                  0.5*weighti*weightj*(-sigmai.dot(vij).dot(deltagrad) + workQj)/mj :                             // CRK
                  -weighti*rhoj*QPiji.dot(vij).dot(gradWi));                                                      // RK, Q term only -- adiabatic portion added later

      // Estimate of delta v (for XSPH).
      XSPHDeltaVi -= fDij*weightj*Wj*vij;
      XSPHDeltaVj += fDij*weighti*Wi*vij;
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OMP parallel

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
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto  mui = mu(nodeListi, i);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massFirstMomenti = massFirstMoment(nodeListi, i);
      auto& massSecondMomentEtai = massSecondMomentEta(nodeListi, i);
      auto& massSecondMomentLabi = massSecondMomentLab(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // Time evolution of the mass density.
      DrhoDti = -rhoi*localDvDxi.Trace();

      // For a surface point, add the RK thermal energy evolution.
      // if (surfacePoint(nodeListi, i) > 1) DepsDti += (Si - Pi*SymTensor::one).doubledot(DvDxi)/rhoi;

      // If needed finish the total energy derivative.
      if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi));

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
                                                       massFirstMomenti,
                                                       massSecondMomentEtai,
                                                       massSecondMomentLabi,
                                                       WR.kernel(),
                                                       hmin,
                                                       hmax,
                                                       hminratio,
                                                       nPerh,
                                                       connectivityMap,
                                                       nodeListi,
                                                       i);

      // Optionally use damage to ramp down stress on damaged material.
      const auto Di = (damageRelieveRubble ? 
                       max(0.0, min(1.0, damage(nodeListi, i).Trace() - 1.0)) :
                       0.0);
      // Hideali = (1.0 - Di)*Hideali + Di*mHfield0(nodeListi, i);
      // DHDti = (1.0 - Di)*DHDti + Di*(mHfield0(nodeListi, i) - Hi)*0.25/dt;

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()/3.0)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;

      // In the presence of damage, add a term to reduce the stress on this point.
      DSDti = (1.0 - Di)*DSDti - Di*Si*0.25/dt;
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
  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  for (auto boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(S);
    (*boundaryItr)->applyFieldListGhostBoundary(K);
    (*boundaryItr)->applyFieldListGhostBoundary(mu);
    (*boundaryItr)->applyFieldListGhostBoundary(Y);
    (*boundaryItr)->applyFieldListGhostBoundary(fragIDs);
    (*boundaryItr)->applyFieldListGhostBoundary(pTypes);
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

  // Enforce boundary conditions on the extra strength variables.
  auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  auto Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  for (auto boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(S);
    (*boundaryItr)->enforceFieldListBoundary(K);
    (*boundaryItr)->enforceFieldListBoundary(mu);
    (*boundaryItr)->enforceFieldListBoundary(Y);
    (*boundaryItr)->enforceFieldListBoundary(fragIDs);
    (*boundaryItr)->enforceFieldListBoundary(pTypes);
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
}

}

