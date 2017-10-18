//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "DamagedNodeCouplingWithFrags.hh"
#include "SPH/SPHHydroBase.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
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
#include "SolidMaterial/SolidEquationOfState.hh"

#include "SolidSPHHydroBase.hh"

#include "Eigen/Dense"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

namespace Spheral {
namespace SPHSpace {

using namespace std;
using SPHSpace::SPHHydroBase;
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
// inline
// Dim<1>::SymTensor::EigenType
// tensileStressCorrection(const Dim<1>::SymTensor::EigenType& sigma) {
//   if (sigma(0,0) > 0.0) {
//     return -sigma;
//   } else {
//     return Dim<1>::SymTensor::EigenType::Zero();
//   }
// }

// inline
// Dim<2>::SymTensor::EigenType
// tensileStressCorrection(const Dim<2>::SymTensor::EigenType& sigma) {
//   Eigen::SelfAdjointEigenSolver<Dim<2>::SymTensor::EigenType> eigensolver(sigma);
//   CHECK(eigensolver.info() == Eigen::Success);
//   const auto lambdax = eigensolver.eigenvalues().x();
//   const auto lambday = eigensolver.eigenvalues().y();
//   const auto vecs = eigensolver.eigenvectors();
//   Dim<2>::SymTensor::EigenType result;
//   result << ((lambdax > 0.0 ? -lambdax : 0.0), 0.0,
//              0.0,                              (lambday > 0.0 ? -lambday : 0.0));
//   result = (vecs * result * vecs.transpose()).eval();
//   return result;
// }

// inline
// Dim<3>::SymTensor::EigenType
// tensileStressCorrection(const Dim<3>::SymTensor::EigenType& sigma) {
//   Eigen::SelfAdjointEigenSolver<Dim<3>::SymTensor::EigenType> eigensolver(sigma);
//   CHECK(eigensolver.info() == Eigen::Success);
//   const auto lambdax = eigensolver.eigenvalues().x();
//   const auto lambday = eigensolver.eigenvalues().y();
//   const auto lambdaz = eigensolver.eigenvalues().z();
//   const auto vecs = eigensolver.eigenvectors();
//   Dim<3>::SymTensor::EigenType result;
//   result << ((lambdax > 0.0 ? -lambdax : 0.0), 0.0,                              0.0,
//              0.0,                              (lambday > 0.0 ? -lambday : 0.0), 0.0,
//              0.0,                              0.0,                              (lambdaz > 0.0 ? -lambdaz : 0.0));
//   result = (vecs * result * vecs.transpose()).eval();
//   return result;
// }

template<typename Dimension>
inline
typename Dimension::SymTensor::EigenType
tensileStressCorrection(const typename Dimension::SymTensor::EigenType& sigma) {
  Eigen::SelfAdjointEigenSolver<typename Dimension::SymTensor::EigenType> eigensolver(sigma);
  CHECK(eigensolver.info() == Eigen::Success);
  const auto lambdas = -(eigensolver.eigenvalues().cwiseMax(0.0));
  const auto vecs = eigensolver.eigenvectors();
  return vecs * lambdas.asDiagonal() * vecs.transpose();
}

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SolidSPHHydroBase<Dimension>::
SolidSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                  ArtificialViscosity<Dimension>& Q,
                  const TableKernel<Dimension>& W,
                  const TableKernel<Dimension>& WPi,
                  const TableKernel<Dimension>& WGrad,
                  const double filter,
                  const double cfl,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool evolveTotalEnergy,
                  const bool gradhCorrection,
                  const bool XSPH,
                  const bool correctVelocityGradient,
                  const bool sumMassDensityOverAllNodeLists,
                  const PhysicsSpace::MassDensityType densityUpdate,
                  const PhysicsSpace::HEvolutionType HUpdate,
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
                          evolveTotalEnergy,
                          gradhCorrection,
                          XSPH,
                          correctVelocityGradient,
                          sumMassDensityOverAllNodeLists,
                          densityUpdate,
                          HUpdate,
                          epsTensile,
                          nTensile,
                          xmin,
                          xmax),
  mGradKernel(WGrad),
  mDdeviatoricStressDt(FieldSpace::FieldStorageType::CopyFields),
  mBulkModulus(FieldSpace::FieldStorageType::CopyFields),
  mShearModulus(FieldSpace::FieldStorageType::CopyFields),
  mYieldStrength(FieldSpace::FieldStorageType::CopyFields),
  mPlasticStrain0(FieldSpace::FieldStorageType::CopyFields),
  mHfield0(FieldSpace::FieldStorageType::CopyFields),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SolidSPHHydroBase<Dimension>::
~SolidSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Call the ancestor.
  SPHHydroBase<Dimension>::initializeProblemStartup(dataBase);

  // Create storage for the state we're holding.
  mDdeviatoricStressDt = dataBase.newSolidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress);
  mBulkModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::bulkModulus);
  mShearModulus = dataBase.newSolidFieldList(0.0, SolidFieldNames::shearModulus);
  mYieldStrength = dataBase.newSolidFieldList(0.0, SolidFieldNames::yieldStrength);
  mPlasticStrain0 = dataBase.newSolidFieldList(0.0, SolidFieldNames::plasticStrain + "0");
  mHfield0 = dataBase.newSolidFieldList(SymTensor::zero, HydroFieldNames::H + "0");

  // Set the moduli.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::SolidNodeListIterator itr = dataBase.solidNodeListBegin();
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
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Invoke SPHHydro's state.
  SPHHydroBase<Dimension>::registerState(dataBase, state);

  // Create the local storage.
  dataBase.resizeFluidFieldList(mBulkModulus, 0.0, SolidFieldNames::bulkModulus, false);
  dataBase.resizeFluidFieldList(mShearModulus, 0.0, SolidFieldNames::shearModulus, false);
  dataBase.resizeFluidFieldList(mYieldStrength, 0.0, SolidFieldNames::yieldStrength, false);
  dataBase.resizeFluidFieldList(mPlasticStrain0, 0.0, SolidFieldNames::plasticStrain + "0", false);

  // Grab the normal Hydro's registered version of the sound speed.
  FieldList<Dimension, Scalar> cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  CHECK(cs.numFields() == dataBase.numFluidNodeLists());

  // Build the FieldList versions of our state.
  FieldList<Dimension, SymTensor> S, D;
  FieldList<Dimension, Scalar> ps;
  FieldList<Dimension, Vector> gradD;
  FieldList<Dimension, int> fragIDs;
  FieldList<Dimension, int> pTypes;
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::SolidNodeListIterator itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    S.appendField((*itr)->deviatoricStress());
    ps.appendField((*itr)->plasticStrain());
    D.appendField((*itr)->effectiveDamage());
    gradD.appendField((*itr)->damageGradient());
    fragIDs.appendField((*itr)->fragmentIDs());
    pTypes.appendField((*itr)->particleTypes());

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
  state.enroll(fragIDs);

  // Register the particle types.
  state.enroll(pTypes);

  // And finally the intial plastic strain.
  state.enroll(mPlasticStrain0);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  // Call the ancestor method.
  SPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  const string DSDtName = IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress;
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, DSDtName, false);

  derivs.enroll(mDdeviatoricStressDt);

  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::SolidNodeListIterator itr = dataBase.solidNodeListBegin();
       itr != dataBase.solidNodeListEnd();
       ++itr, ++nodeListi) {
    CHECK((*itr) != 0);
    derivs.enroll((*itr)->plasticStrainRate());
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  typedef typename Vector::EigenType EigenVector;
  typedef typename Tensor::EigenType EigenTensor;
  typedef typename SymTensor::EigenType EigenSymTensor;

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();
  const TableKernel<Dimension>& WG = this->GradKernel();
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const double tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);
  const Scalar WQ0 = WQ(0.0, 1.0);
  const Scalar epsTensile = this->epsilonTensile();
  const bool compatibleEnergy = this->compatibleEnergyEvolution();
  const bool XSPH = this->XSPH();

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
  const FieldList<Dimension, Scalar> omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  const FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const FieldList<Dimension, SymTensor> damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const FieldList<Dimension, Vector> gradDamage = state.fields(SolidFieldNames::damageGradient, Vector::zero);
  const FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const FieldList<Dimension, int> pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(gradDamage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
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
  FieldList<Dimension, Scalar> rhoSumCorrection = derivatives.fields(HydroFieldNames::massDensityCorrection, 0.0);
  FieldList<Dimension, Scalar> viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Scalar> XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  FieldList<Dimension, SymTensor> DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  CHECK(rhoSum.size() == numNodeLists);
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
  CHECK(rhoSumCorrection.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
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
  Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
  Tensor QPiij0, QPiji0;
  EigenSymTensor sigmai, sigmaj;

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const FluidNodeList<Dimension>& nodeList = **itr;
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
      const auto ri = position(nodeListi, i).eigen();
      const auto mi = mass(nodeListi, i);
      const auto vi = velocity(nodeListi, i).eigen();
      const auto rhoi = massDensity(nodeListi, i);
      const auto epsi = specificThermalEnergy(nodeListi, i);
      const auto Pi = pressure(nodeListi, i);
      const auto Hi = H(nodeListi, i).eigen();
      const auto ci = soundSpeed(nodeListi, i);
      const auto omegai = omega(nodeListi, i);
      const auto Si = S(nodeListi, i).eigen();
      const auto mui = mu(nodeListi, i);
      const auto Hdeti = Hi.determinant();
      const auto safeOmegai = safeInv(omegai, tiny);
      const auto fragIDi = fragIDs(nodeListi, i);
      const auto pTypei = pTypes(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure(nodeListi, i);
      auto& rhoSumCorrectioni = rhoSumCorrection(nodeListi, i);
      auto& viscousWorki = viscousWork(nodeListi, i);
      auto& pairAccelerationsi = pairAccelerations(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);
      auto& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const auto& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

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
            const auto j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const auto rj = position(nodeListj, j).eigen();
              const auto mj = mass(nodeListj, j);
              const auto vj = velocity(nodeListj, j).eigen();
              const auto rhoj = massDensity(nodeListj, j);
              const auto epsj = specificThermalEnergy(nodeListj, j);
              const auto Pj = pressure(nodeListj, j);
              const auto Hj = H(nodeListj, j).eigen();
              const auto cj = soundSpeed(nodeListj, j);
              const auto omegaj = omega(nodeListj, j);
              const auto Sj = S(nodeListj, j).eigen();
              const auto Hdetj = Hj.determinant();
              const auto safeOmegaj = safeInv(omegaj, tiny);
              const auto fragIDj = fragIDs(nodeListj, j);
              const auto pTypej = pTypes(nodeListj, j);
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);

              auto& rhoSumj = rhoSum(nodeListj, j);
              auto& DxDtj = DxDt(nodeListj, j);
              auto& DvDtj = DvDt(nodeListj, j);
              auto& DepsDtj = DepsDt(nodeListj, j);
              auto& DvDxj = DvDx(nodeListj, j);
              auto& localDvDxj = localDvDx(nodeListj, j);
              auto& Mj = M(nodeListj, j);
              auto& localMj = localM(nodeListj, j);
              auto& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              auto& effViscousPressurej = effViscousPressure(nodeListj, j);
              auto& rhoSumCorrectionj = rhoSumCorrection(nodeListj, j);
              auto& viscousWorkj = viscousWork(nodeListj, j);
              auto& pairAccelerationsj = pairAccelerations(nodeListj, j);
              auto& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
              auto& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              auto& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              auto& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Flag if this is a contiguous material pair or not.
              const bool sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

              // Flag if at least one particle is free (0).
              const bool freeParticle = (pTypei == 0 or pTypej == 0);

              // Node displacement.
              const auto rij = (ri - rj).eval();
              const auto etai = (Hi*rij).eval();
              const auto etaj = (Hj*rij).eval();
              const auto etaMagi = etai.norm();
              const auto etaMagj = etaj.norm();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);

              // Symmetrized kernel weight and gradient.
              std::tie(Wi, gWi) = W.kernelAndGradValue(etaMagi, Hdeti);
              std::tie(WQi, gWQi) = WQ.kernelAndGradValue(etaMagi, Hdeti);
              const auto Hetai = (Hi*etai.normalized()).eval();
              const auto gradWi = (gWi*Hetai).eval();
              const auto gradWQi = (gWQi*Hetai).eval();
              const auto gradWGi = (WG.gradValue(etaMagi, Hdeti) * Hetai).eval();

              std::tie(Wj, gWj) = W.kernelAndGradValue(etaMagj, Hdetj);
              std::tie(WQj, gWQj) = WQ.kernelAndGradValue(etaMagj, Hdetj);
              const auto Hetaj = (Hj*etaj.normalized()).eval();
              const auto gradWj = (gWj*Hetaj).eval();
              const auto gradWQj = (gWQj*Hetaj).eval();
              const auto gradWGj = (WG.gradValue(etaMagj, Hdetj) * Hetaj).eval();

              // Determine how we're applying damage.
              const auto fDeffij = coupling(nodeListi, i, nodeListj, j);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const auto fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
              const auto rij2 = rij.squaredNorm();
              const auto thpt = (rij*(rij.transpose())*safeInvVar(rij2*rij2*rij2)).eval();
              weightedNeighborSumi +=     fweightij*abs(gWi);
              weightedNeighborSumj += 1.0/fweightij*abs(gWj);
              massSecondMomenti +=     fweightij*gradWi.squaredNorm()*thpt;
              massSecondMomentj += 1.0/fweightij*gradWj.squaredNorm()*thpt;

              // Contribution to the sum density (only if the same material).
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wi;
                rhoSumj += mi*Wj;
              }

              // Contribution to the sum density correction
              rhoSumCorrectioni += mj * WQi / rhoj ;
              rhoSumCorrectionj += mi * WQj / rhoi ;

              // Compute the pair-wise artificial viscosity.
              const auto vij = vi - vj;
              std::tie(QPiij0, QPiji0) = Q.Piij(nodeListi, i, nodeListj, j,
                                                position(nodeListi, i), Vector(etai), velocity(nodeListi, i), rhoi, ci, H(nodeListi, i),
                                                position(nodeListj, j), Vector(etaj), velocity(nodeListj, j), rhoj, cj, H(nodeListj, j));
              const auto QPiij = QPiij0.eigen();
              const auto QPiji = QPiji0.eigen();
              const auto Qacci = 0.5*(QPiij*gradWQi);
              const auto Qaccj = 0.5*(QPiji*gradWQj);
              const auto workQi = vij.dot(Qacci);
              const auto workQj = vij.dot(Qaccj);
              const auto Qi = rhoi*rhoi*(QPiji.diagonal().cwiseAbs().maxCoeff());
              const auto Qj = rhoj*rhoj*(QPiji.diagonal().cwiseAbs().maxCoeff());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += mj*Qi*WQi/rhoj;
              effViscousPressurej += mi*Qj*WQj/rhoi;
              viscousWorki += mj*workQi;
              viscousWorkj += mi*workQj;

              // Damage scaling of negative pressures.
              const auto Peffi = (Pi > 0.0 ? Pi : fDeffij*Pi);
              const auto Peffj = (Pj > 0.0 ? Pj : fDeffij*Pj);

              // Compute the stress tensors.
              if (sameMatij) {
                sigmai.noalias() = fDeffij*Si - Peffi*EigenSymTensor::Identity();
                sigmaj.noalias() = fDeffij*Sj - Peffj*EigenSymTensor::Identity();
              } else {
                sigmai.noalias() = -Peffi*EigenSymTensor::Identity();
                sigmaj.noalias() = -Peffj*EigenSymTensor::Identity();
              }

              // Compute the tensile correction to add to the stress as described in 
              // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
              // const auto fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
              // const auto fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
              // const auto Ri = fi*tensileStressCorrection<Dimension>(sigmai);
              // const auto Rj = fj*tensileStressCorrection<Dimension>(sigmaj);
              // sigmai += Ri;
              // sigmaj += Rj;

              // Acceleration.
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const auto sigmarhoi = safeOmegai*sigmai/(rhoi*rhoi);
              const auto sigmarhoj = safeOmegaj*sigmaj/(rhoj*rhoj);
              const auto deltaDvDt = sigmarhoi*gradWi + sigmarhoj*gradWj - Qacci - Qaccj;
              if (freeParticle) {
                DvDti += mj*deltaDvDt;
                DvDtj -= mi*deltaDvDt;
              }

              // Pair-wise portion of grad velocity.
              const auto deltaDvDxi = fDeffij*vij*(gradWGi.transpose());
              const auto deltaDvDxj = fDeffij*vij*(gradWGj.transpose());

              // Specific thermal energy evolution.
              DepsDti -= mj*(fDeffij*(sigmarhoi*deltaDvDxi).trace() - workQi);
              DepsDtj -= mi*(fDeffij*(sigmarhoj*deltaDvDxj).trace() - workQj);
              if (compatibleEnergy) {
                pairAccelerationsi.push_back( mj*Vector(deltaDvDt));
                pairAccelerationsj.push_back(-mi*Vector(deltaDvDt));
              }

              // Velocity gradient.
              DvDxi -= mj*deltaDvDxi;
              DvDxj -= mi*deltaDvDxj;
              if (sameMatij) {
                localDvDxi -= mj*deltaDvDxi;
                localDvDxj -= mi*deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (XSPH and sameMatij) {
                const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
                XSPHWeightSumi += wXSPHij;
                XSPHWeightSumj += wXSPHij;
                XSPHDeltaVi -= wXSPHij*vij;
                XSPHDeltaVj += wXSPHij*vij;
              }

              // Linear gradient correction term.
              Mi -= mj*rij*(gradWGi.transpose());
              Mj -= mi*rij*(gradWGj.transpose());
              if (sameMatij) {
                localMi -= mj*rij*(gradWGi.transpose());
                localMj -= mi*rij*(gradWGj.transpose());
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
      const auto deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;

      // Add the self-contribution to density sum correction.
      rhoSumCorrectioni += mi*WQ0*Hdeti/rhoi ;

      // Correct the effective viscous pressure.
      effViscousPressurei /= rhoSumCorrectioni ;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
      } else {
        DvDxi /= rhoi;
      }
      if (this->mCorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      } else {
        localDvDxi /= rhoi;
      }

      // Evaluate the continuity equation.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (this->mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti.eigen()) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      DxDti = vi;
      if (XSPH) {
        CHECK(XSPHWeightSumi >= 0.0);
        XSPHWeightSumi += Hdeti*mi/rhoi*W0 + 1.0e-30;
        DxDti += XSPHDeltaVi/XSPHWeightSumi;
      }

      // The H tensor evolution.
      DHDti = smoothingScaleMethod.smoothingScaleDerivative(H(nodeListi, i),
                                                            position(nodeListi, i),
                                                            DvDxi,
                                                            hmin,
                                                            hmax,
                                                            hminratio,
                                                            nPerh);
      Hideali = smoothingScaleMethod.newSmoothingScale(H(nodeListi, i),
                                                       position(nodeListi, i),
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

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.25/dt*Di*(rhoi - rho0);

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()/3.0)*SymTensor::one;
      const auto spinCorrection = (spin*S(nodeListi, i) + S(nodeListi, i)*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;

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
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Ancestor method.
  SPHHydroBase<Dimension>::applyGhostBoundaries(state, derivs);

  // Apply boundary conditions to our extra strength variables.
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  FieldList<Dimension, int> pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
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
SolidSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Ancestor method.
  SPHHydroBase<Dimension>::enforceBoundaries(state, derivs);

  // Enforce boundary conditions on the extra strength variable.s
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);
  FieldList<Dimension, int> fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  FieldList<Dimension, int> pTypes = state.fields(SolidFieldNames::particleTypes, int(0));

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
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
SolidSPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // Ancestor method.
  SPHHydroBase<Dimension>::dumpState(file, pathName);

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
SolidSPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  // Ancestor method.
  SPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mDdeviatoricStressDt, pathName + "/DdeviatoricStressDt");
  file.read(mBulkModulus, pathName + "/bulkModulus");
  file.read(mShearModulus, pathName + "/shearModulus");
  file.read(mYieldStrength, pathName + "/yieldStrength");
  file.read(mPlasticStrain0, pathName + "/plasticStrain0");
  file.read(mHfield0, pathName + "/Hfield0");
}

}
}

