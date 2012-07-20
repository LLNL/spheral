//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBase -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "SolidSPHHydroBase.hh"
#include "SPH/SPHHydroBase.hh"
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
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "FileIO/FileIO.hh"

#include "TAU.h"

namespace Spheral {
namespace SolidSPHSpace {

using namespace std;
using SPHSpace::SPHHydroBase;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using SolidMaterial::SolidNodeList;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Construct a unit vector of the argument, going to zero as teh magnitude
// falls below a given "fuzz".
//------------------------------------------------------------------------------
template<typename Vector>
inline
Vector
unitVectorWithZero(const Vector& x,
                   const double fuzz = 0.01) {
  if (x.magnitude2() < fuzz) {
    return Vector::zero;
  } else {
    return x.unitVector();
  }
}

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
SolidSPHHydroBase<Dimension>::
SolidSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                  const TableKernel<Dimension>& W,
                  const TableKernel<Dimension>& WPi,
                  ArtificialViscosity<Dimension>& Q,
                  const double cfl,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool gradhCorrection,
                  const bool XSPH,
                  const PhysicsSpace::MassDensityType densityUpdate,
                  const PhysicsSpace::HEvolutionType HUpdate,
                  const double epsTensile,
                  const double nTensile,
                  const Vector& xmin,
                  const Vector& xmax):
  SPHHydroBase<Dimension>(smoothingScaleMethod, 
                          W,
                          WPi,
                          Q,
                          cfl,
                          useVelocityMagnitudeForDt,
                          compatibleEnergyEvolution,
                          gradhCorrection,
                          XSPH,
                          densityUpdate,
                          HUpdate,
                          epsTensile,
                          nTensile,
                          xmin,
                          xmax),
  mDdeviatoricStressDt(FieldList<Dimension, SymTensor>::Copy),
  mBulkModulus(FieldList<Dimension, Scalar>::Copy),
  mShearModulus(FieldList<Dimension, Scalar>::Copy),
  mYieldStrength(FieldList<Dimension, Scalar>::Copy),
  mPlasticStrain0(FieldList<Dimension, Scalar>::Copy),
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

  // TAU timers.
  TAU_PROFILE("SolidSPHHydroBase", "::initializeProblemStartup", TAU_USER);

  // Call the ancestor.
  SPHHydroBase<Dimension>::initializeProblemStartup(dataBase);

  // Create storage for the bulk and shear moduli, and yield strength.
  mBulkModulus = dataBase.newFluidFieldList(0.0, SolidFieldNames::bulkModulus);
  mShearModulus = dataBase.newFluidFieldList(0.0, SolidFieldNames::shearModulus);
  mYieldStrength = dataBase.newFluidFieldList(0.0, SolidFieldNames::yieldStrength);

  // Set those suckers.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    SolidNodeList<Dimension>* solidNodeListPtr = dynamic_cast<SolidNodeList<Dimension>*>(*itr);
    CHECK(solidNodeListPtr != 0);
    solidNodeListPtr->bulkModulus(*mBulkModulus[nodeListi]);
    solidNodeListPtr->shearModulus(*mShearModulus[nodeListi]);
    solidNodeListPtr->yieldStrength(*mYieldStrength[nodeListi]);
  }
}


//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // TAU timers.
  TAU_PROFILE("SolidSPHHydroBase", "::registerState", TAU_USER);

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

  // Now register away.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    SolidNodeList<Dimension>* solidNodeListPtr = dynamic_cast<SolidNodeList<Dimension>*>(*itr);
    CHECK(solidNodeListPtr != 0);

    // Register the deviatoric stress and plastic strain to be evolved.
    PolicyPointer deviatoricStressPolicy(new DeviatoricStressPolicy<Dimension>());
    PolicyPointer plasticStrainPolicy(new PlasticStrainPolicy<Dimension>());
    state.enroll(solidNodeListPtr->deviatoricStress(), deviatoricStressPolicy);
    state.enroll(solidNodeListPtr->plasticStrain(), plasticStrainPolicy);

    // Register the bulk modulus, shear modulus, and yield strength.
    PolicyPointer bulkModulusPolicy(new BulkModulusPolicy<Dimension>());
    PolicyPointer shearModulusPolicy(new ShearModulusPolicy<Dimension>());
    PolicyPointer yieldStrengthPolicy(new YieldStrengthPolicy<Dimension>());
    state.enroll(*mBulkModulus[nodeListi], bulkModulusPolicy);
    state.enroll(*mShearModulus[nodeListi], shearModulusPolicy);
    state.enroll(*mYieldStrength[nodeListi], yieldStrengthPolicy);

    // Override the policy for the sound speed.
    PolicyPointer csPolicy(new StrengthSoundSpeedPolicy<Dimension>());
    state.enroll(*cs[nodeListi], csPolicy);

    // Register the effective damage and damage gradient with default no-op updates.
    // If there are any damage models running they can override these choices.
    state.enroll(solidNodeListPtr->effectiveDamage());
    state.enroll(solidNodeListPtr->damageGradient());

    // Make a copy of the beginning plastic strain.
    *mPlasticStrain0[nodeListi] = solidNodeListPtr->plasticStrain();
    (*mPlasticStrain0[nodeListi]).name(SolidFieldNames::plasticStrain + "0");
    state.enroll(*mPlasticStrain0[nodeListi]);
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  // TAU timers.
  TAU_PROFILE("SolidSPHHydroBase", "::registerDerivatives", TAU_USER);

  // Call the ancestor method.
  SPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  const string DSDtName = IncrementState<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStress;
  dataBase.resizeFluidFieldList(mDdeviatoricStressDt, SymTensor::zero, DSDtName, false);

  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    SolidNodeList<Dimension>* solidNodeListPtr = dynamic_cast<SolidNodeList<Dimension>*>(*itr);
    CHECK(solidNodeListPtr != 0);
    derivs.enroll(*mDdeviatoricStressDt[nodeListi]);
    derivs.enroll(solidNodeListPtr->plasticStrainRate());
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

  // TAU timers.
  TAU_PROFILE("SolidSPHHydroBase", "::evaluateDerivatives", TAU_USER);
  TAU_PROFILE_TIMER(TimeSolidSPHHydroBaseCalcDerivs, "SolidSPHHydroBase", "::evaluateDerivatives : calc derivatives", TAU_USER);
  TAU_PROFILE_TIMER(TimeSolidSPHHydroBaseGetFieldLists, "SPHNode", "::evaluateDerivatives : get FieldLists", TAU_USER);
  TAU_PROFILE_TIMER(TimeSolidSPHHydroBaseNodeIState, "SPHNode", "::evaluateDerivatives : node i state", TAU_USER);

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const Scalar W0 = W(0.0, 1.0);
  const Scalar epsTensile = this->epsilonTensile();
  const bool compatibleEnergy = this->compatibleEnergyEvolution();
  const bool XSPH = this->XSPH();

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  TAU_PROFILE_START(TimeSolidSPHHydroBaseGetFieldLists);
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

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Scalar> XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  FieldList<Dimension, SymTensor> DSDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, SymTensor> >::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);
  TAU_PROFILE_STOP(TimeSolidSPHHydroBaseGetFieldLists);

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

  // Start our big loop over all FluidNodeLists.
  TAU_PROFILE_START(TimeSolidSPHHydroBaseCalcDerivs);
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

    // The scale for the tensile correction.
    const Scalar WnPerh = W(1.0/nPerh, 1.0);

    // Get the work field for this NodeList.
    Field<Dimension, Scalar>& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      TAU_PROFILE_START(TimeSolidSPHHydroBaseNodeIState);
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
      const Scalar& omegai = omega(nodeListi, i);
      const SymTensor& Si = S(nodeListi, i);
      const Scalar& mui = mu(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar safeOmegai = omegai/(omegai*omegai + 1.0e-4);
      const Scalar sDi = max(0.0, min(1.0, damage(nodeListi, i).eigenValues().maxElement()));
      const Vector gradDi = unitVectorWithZero(gradDamage(nodeListi, i)*Dimension::nDim/Hi.Trace());
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(omegai > 0.0);
      CHECK(Hdeti > 0.0);

      Scalar& rhoSumi = rhoSum(nodeListi, i);
      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Scalar& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      SymTensor& DSDti = DSDt(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);
      TAU_PROFILE_STOP(TimeSolidSPHHydroBaseNodeIState);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = (nodeListi == nodeListj ? 1.0 : 0.2);
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
              const Scalar& omegaj = omega(nodeListj, j);
              const SymTensor& Sj = S(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar safeOmegaj = omegaj/(omegaj*omegaj + 1.0e-4);
              const Scalar sDj = max(0.0, min(1.0, damage(nodeListj, j).eigenValues().maxElement()));
              const Vector gradDj = unitVectorWithZero(gradDamage(nodeListj, j)*Dimension::nDim/Hj.Trace());
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);

              Scalar& rhoSumj = rhoSum(nodeListj, j);
              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
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

              // Determine how we're applying damage.
              const Scalar gradDdot = gradDi.dot(gradDj);
              const Scalar phi = ((abs(gradDdot) < 0.1 or min(sDi, sDj) < 0.05) ? 
                                  1.0 :
                                  max(0.0, min(1.0, -gradDdot)));
              CHECK(phi >= 0.0 and phi <= 1.0);
              const Scalar fDeffij = FastMath::pow4(max(0.0, min(1.0, 1.0 - phi*max(sDi, sDj))));

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()/(rij2 + 1.0e-10) / FastMath::square(Dimension::pownu12(rij2 + 1.0e-10));
              weightedNeighborSumi += fweightij*abs(gWi);
              weightedNeighborSumj += fweightij*abs(gWj);
              massSecondMomenti += fweightij*gradWi.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWj.magnitude2()*thpt;

              // Contribution to the sum density (only if the same material).
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wi;
                rhoSumj += mi*Wj;
              }

              // Mass density evolution.
              const Vector vij = vi - vj;
              const double deltaDrhoDti = fDeffij*(vij.dot(gradWi));
              const double deltaDrhoDtj = fDeffij*(vij.dot(gradWj));
              DrhoDti += deltaDrhoDti;
              DrhoDtj += deltaDrhoDtj;

              // Compute the pair-wise artificial viscosity.
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Vector Qacci = 0.5*(QPiij.first *gradWQi);
              const Vector Qaccj = 0.5*(QPiij.second*gradWQj);
              const Scalar workQi = vij.dot(Qacci);
              const Scalar workQj = vij.dot(Qaccj);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);

              // Damage scaling of negative pressures.
              const Scalar Peffi = (Pi > 0.0 ? Pi : fDeffij*Pi);
              const Scalar Peffj = (Pj > 0.0 ? Pj : fDeffij*Pj);

              // Compute the stress tensors.
              SymTensor sigmai = -Peffi*SymTensor::one;
              SymTensor sigmaj = -Peffj*SymTensor::one;
              if (nodeListi == nodeListj) {
                sigmai += fDeffij*Si;
                sigmaj += fDeffij*Sj;
              }

              // Compute the tensile correction to add to the stress as described in 
              // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
              const Scalar fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
              const Scalar fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
              const SymTensor Ri = fi*tensileStressCorrection(sigmai);
              const SymTensor Rj = fj*tensileStressCorrection(sigmaj);
              sigmai += Ri;
              sigmaj += Rj;

              // Acceleration.
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const SymTensor sigmarhoi = sigmai/(rhoi*rhoi);
              const SymTensor sigmarhoj = sigmaj/(rhoj*rhoj);
              const Vector deltaDvDt = fDeffij*(sigmarhoi*safeOmegai*gradWi + sigmarhoj*safeOmegaj*gradWj) - Qacci - Qaccj;
              DvDti += mj*deltaDvDt;
              DvDtj -= mi*deltaDvDt;

              // Pair-wise portion of grad velocity.
              const Tensor deltaDvDxi = fDeffij*vij.dyad(gradWi);
              const Tensor deltaDvDxj = fDeffij*vij.dyad(gradWj);

              // Specific thermal energy evolution.
              DepsDti -= mj*(fDeffij*sigmarhoi.doubledot(deltaDvDxi.Symmetric()) - workQi);
              DepsDtj -= mi*(fDeffij*sigmarhoj.doubledot(deltaDvDxj.Symmetric()) - workQj);
              if (compatibleEnergy) {
                if (i < firstGhostNodei) pairAccelerationsi.push_back( mj*deltaDvDt);
                if (j < firstGhostNodej) pairAccelerationsj.push_back(-mi*deltaDvDt);
              }

              // Velocity gradient.
              DvDxi -= mj*deltaDvDxi;
              DvDxj -= mi*deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= mj*deltaDvDxi;
                localDvDxj -= mi*deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (XSPH and (nodeListi == nodeListj)) {
                const double fXSPH = fDeffij*max(0.0, min(1.0, abs(vij.dot(rij)*safeInv(vij.magnitude()*rij.magnitude()))));
                CHECK(fXSPH >= 0.0 and fXSPH <= 1.0);
                XSPHWeightSumi += fXSPH*mj/rhoj*Wi;
                XSPHWeightSumj += fXSPH*mi/rhoi*Wj;
                XSPHDeltaVi -= fXSPH*mj/rhoj*Wi*vij;
                XSPHDeltaVj += fXSPH*mi/rhoi*Wj*vij;
              }

            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not compatibleEnergy or 
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;

      // Finish the continuity equation.
      DrhoDti *= mi*safeOmegai;

      // Finish the thermal energy derivative.
      DepsDti *= safeOmegai;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      DvDxi *= safeOmegai/rhoi;
      localDvDxi *= safeOmegai/rhoi;

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        CHECK(XSPHWeightSumi >= 0.0);
        XSPHWeightSumi += Hdeti*mi/rhoi*W0 + 1.0e-30;
        DxDti = vi + XSPHDeltaVi/XSPHWeightSumi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = smoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                            DvDxi,
                                                            hmin,
                                                            hmax,
                                                            hminratio,
                                                            nPerh,
                                                            maxNumNeighbors);
      Hideali = smoothingScaleMethod.newSmoothingScale(Hi,
                                                       weightedNeighborSumi,
                                                       massSecondMomenti,
                                                       numNeighborsi,
                                                       W,
                                                       hmin,
                                                       hmax,
                                                       hminratio,
                                                       nPerh,
                                                       maxNumNeighbors);

      // Determine the deviatoric stress evolution.
      const SymTensor deformation = localDvDxi.Symmetric();
      const Tensor spin = localDvDxi.SkewSymmetric();
      const SymTensor deviatoricDeformation = deformation - (deformation.Trace()/3.0)*SymTensor::one;
      const SymTensor spinCorrection = (spin*Si + Si*spin).Symmetric();
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
  TAU_PROFILE_STOP(TimeSolidSPHHydroBaseCalcDerivs);
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // TAU timers.
  TAU_PROFILE("SolidSPHHydroBase", "::applyGhostBoundaries", TAU_USER);

  // Ancestor method.
  SPHHydroBase<Dimension>::applyGhostBoundaries(state, derivs);

  // Apply boundary conditions to our extra strength variables.
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(S);
    (*boundaryItr)->applyFieldListGhostBoundary(K);
    (*boundaryItr)->applyFieldListGhostBoundary(mu);
    (*boundaryItr)->applyFieldListGhostBoundary(Y);
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

  // TAU timers.
  TAU_PROFILE("SolidSPHHydroBase", "::enforceBoundaries", TAU_USER);

  // Ancestor method.
  SPHHydroBase<Dimension>::enforceBoundaries(state, derivs);

  // Enforce boundary conditions on the extra strength variable.s
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> K = state.fields(SolidFieldNames::bulkModulus, 0.0);
  FieldList<Dimension, Scalar> mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  FieldList<Dimension, Scalar> Y = state.fields(SolidFieldNames::yieldStrength, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(S);
    (*boundaryItr)->enforceFieldListBoundary(K);
    (*boundaryItr)->enforceFieldListBoundary(mu);
    (*boundaryItr)->enforceFieldListBoundary(Y);
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
}

}
}

//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidSPHSpace {
    template class SolidSPHHydroBase< Dim<1> >;
    template class SolidSPHHydroBase< Dim<2> >;
    template class SolidSPHHydroBase< Dim<3> >;
  }
}
