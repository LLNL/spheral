//---------------------------------Spheral++----------------------------------//
// CRKSPHHydroBaseRZ -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Thu May 12 15:25:24 PDT 2016
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "CRKSPHUtilities.hh"
#include "computeHullVolumes.hh"
#include "computeCRKSPHSumVolume.hh"
#include "computeHVolumes.hh"
#include "SPH/computeSPHSumMassDensity.hh"
#include "SPH/correctSPHSumMassDensity.hh"
#include "computeCRKSPHSumMassDensity.hh"
#include "computeCRKSPHMoments.hh"
#include "computeCRKSPHCorrections.hh"
#include "computeCRKSPHIntegral.hh"
#include "centerOfMass.hh"
#include "volumeSpacing.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
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
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "ContinuityVolumePolicyRZ.hh"
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
#include "SPH/computeSPHSumMassDensity.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"
#include "Kernel/NBSplineKernel.hh"

#include "CRKSPHHydroBaseRZ.hh"

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
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
CRKSPHHydroBaseRZ::
CRKSPHHydroBaseRZ(const SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
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
  CRKSPHHydroBase<Dim<2> >(smoothingScaleMethod,
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
                           epsTensile,
                           nTensile,
                           limitMultimaterialTopology) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
CRKSPHHydroBaseRZ::
~CRKSPHHydroBaseRZ() {
}

//------------------------------------------------------------------------------
// On problem start up, we set the RZ flag on the DataBase.
//------------------------------------------------------------------------------
void
CRKSPHHydroBaseRZ::
initializeProblemStartup(DataBase<Dim<2> >& dataBase) {
  dataBase.isRZ = true;

  // Correct the mass to mass/r.
  auto mass = dataBase.fluidMass();
  const auto pos = dataBase.fluidPosition();
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }

  // Do general initializations.
  CRKSPHHydroBase<Dim<2> >::initializeProblemStartup(dataBase);

  // Convert back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
CRKSPHHydroBaseRZ::
registerState(DataBase<Dim<2> >& dataBase,
              State<Dim<2> >& state) {

  typedef State<Dimension>::PolicyPointer PolicyPointer;

  // The base class does most of it.
  CRKSPHHydroBase<Dim<2> >::registerState(dataBase, state);

  // Reregister the volume update
  PolicyPointer volumePolicy(new ContinuityVolumePolicyRZ());
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  state.enroll(vol, volumePolicy);

  // Are we using the compatible energy evolution scheme?
  // If so we need to override the ordinary energy registration with a specialized version.
  if (mCompatibleEnergyEvolution) {
    FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
    PolicyPointer thermalEnergyPolicy(new NonSymmetricSpecificThermalEnergyPolicy<Dimension>(dataBase));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);

    // Get the policy for the position, and add the specific energy as a dependency.
    PolicyPointer positionPolicy = state.policy(state.buildFieldKey(HydroFieldNames::position, UpdatePolicyBase<Dimension>::wildcard()));
    positionPolicy->addDependency(HydroFieldNames::specificThermalEnergy);
  }
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
void
CRKSPHHydroBaseRZ::
preStepInitialize(const DataBase<Dim<2>>& dataBase, 
                  State<Dim<2>>& state,
                  StateDerivatives<Dim<2>>& derivs) {

  // Convert the mass to mass per unit length first.
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      mass(nodeListi, i) /= circi;
    }
  }

  // Base class does most of the work.
  CRKSPHHydroBase<Dimension>::preStepInitialize(dataBase, state, derivs);

  // Now convert back to true masses.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto& xi = pos(nodeListi, i);
      const auto circi = 2.0*M_PI*abs(xi.y());
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
CRKSPHHydroBaseRZ::
evaluateDerivatives(const Dim<2>::Scalar time,
                    const Dim<2>::Scalar dt,
                    const DataBase<Dim<2> >& dataBase,
                    const State<Dim<2> >& state,
                    StateDerivatives<Dim<2> >& derivatives) const {

  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  typedef Timing::Time Time;
  const auto tiny = 1.0e-30;
  const auto W0 = W.kernelValue(0.0, 1.0);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto  numNodeLists = nodeLists.size();
  const auto  order = this->correctionOrder();
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
  const auto A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const auto B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  const auto C = state.fields(HydroFieldNames::C_CRKSPH, Tensor::zero);
  const auto gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  const auto gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  const auto gradC = state.fields(HydroFieldNames::gradC_CRKSPH, ThirdRankTensor::zero);
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
  // CHECK(surfNorm.size() == numNodeLists);

  // Derivative FieldLists.
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
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
  CHECK(massSecondMoment.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (mCompatibleEnergyEvolution) pairAccelerations.resize(2*npairs);

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
    Tensor QPiij, QPiji;
    Scalar Ai, Aj;
    Vector gradAi, gradAj, forceij, forceji;
    Vector Bi = Vector::zero, Bj = Vector::zero;
    Tensor Ci = Tensor::zero, Cj = Tensor::zero;
    Tensor gradBi = Tensor::zero, gradBj = Tensor::zero;
    ThirdRankTensor gradCi = ThirdRankTensor::zero, gradCj = ThirdRankTensor::zero;
    Vector gradWi, gradWj;
    Vector deltagrad;

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
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0; kk < npairs; ++kk) {
      const auto start = Timing::currentTime();
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& posi = position(nodeListi, i);
      const auto  ri = abs(posi.y());
      const auto  circi = 2.0*M_PI*ri;
      const auto  mi = mass(nodeListi, i);
      const auto  mRZi = mi/circi;
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ci = soundSpeed(nodeListi, i);
      const auto  Ai = A(nodeListi, i);
      const auto& gradAi = gradA(nodeListi, i);
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
      const auto zetai = abs((Hi*posi).y());
      const auto hri = ri*safeInv(zetai);
      CHECK2(ri > 0.0, i << " " << ri);
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);
      // CHECK2(Ai > 0.0, i << " " << Ai);
      CHECK2(Hdeti > 0.0, i << " " << Hdeti);
      CHECK2(weighti > 0.0, i << " " << weighti);

      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& viscousWorki = viscousWork_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);

      // Get the state for node j
      const auto& posj = position(nodeListj, j);
      const auto  rj = abs(posj.y());
      const auto  circj = 2.0*M_PI*rj;
      const auto  mj = mass(nodeListj, j);
      const auto  mRZj = mj/circj;
      const auto& vj = velocity(nodeListj, j);
      const auto  rhoj = massDensity(nodeListj, j);
      const auto  epsj = specificThermalEnergy(nodeListj, j);
      const auto  Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  cj = soundSpeed(nodeListj, j);
      const auto  Aj = A(nodeListj, j);
      const auto& gradAj = gradA(nodeListj, j);
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
      const auto zetaj = abs((Hj*posj).y());
      CHECK2(rj > 0.0, j << " " << rj);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      // CHECK(Aj > 0.0 or j >= firstGhostNodej);
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
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);

      // Node displacement.
      const auto xij = posi - posj;
      const auto etai = Hi*xij;
      const auto etaj = Hj*xij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);
      const auto vij = vi - vj;

      // Symmetrized kernel weight and gradient.
      CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, order,  xij,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi);
      CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, order, -xij, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj);
      deltagrad = gradWj - gradWi;
      const auto gradWSPHi = (Hi*etai.unitVector())*W.gradValue(etai.magnitude(), Hdeti);
      const auto gradWSPHj = (Hj*etaj.unitVector())*W.gradValue(etaj.magnitude(), Hdetj);

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      const auto fweightij = nodeListi == nodeListj ? 1.0 : mRZj*rhoi/(mRZi*rhoj);
      const auto xij2 = xij.magnitude2();
      const auto thpt = xij.selfdyad()*safeInvVar(xij2*xij2*xij2);
      weightedNeighborSumi +=     fweightij*std::abs(gWi);
      weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
      massSecondMomenti +=     fweightij*gradWSPHi.magnitude2()*thpt;
      massSecondMomentj += 1.0/fweightij*gradWSPHj.magnitude2()*thpt;

      // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
      std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                      posi, etai, vi, rhoi, ci, Hi,
                                      posj, etaj, vj, rhoj, cj, Hj);
      const auto Qaccij = (rhoi*rhoi*QPiij + rhoj*rhoj*QPiji).dot(deltagrad);
      const auto workQi = rhoj*rhoj*QPiji.dot(vij).dot(deltagrad);                // CRK
      const auto workQj = rhoi*rhoi*QPiij.dot(vij).dot(deltagrad);                // CRK
      const auto Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
      const auto Qj = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
      maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                     // We need tighter timestep controls on the Q with CRK
      maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
      effViscousPressurei += weightj * Qi * Wj;
      effViscousPressurej += weighti * Qj * Wi;
      viscousWorki += 0.5*weighti*weightj/mi*workQi;
      viscousWorkj += 0.5*weighti*weightj/mj*workQj;

      // Velocity gradient.
      DvDxi -= weightj*vij.dyad(gradWj);
      DvDxj += weighti*vij.dyad(gradWi);
      if (nodeListi == nodeListj) {
        localDvDxi -= weightj*vij.dyad(gradWj);
        localDvDxj += weighti*vij.dyad(gradWi);
      }

      // Acceleration (CRKSPH form).
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      const auto forceij  = 0.5*weighti*weightj*((Pi + Pj)*deltagrad + Qaccij); // <- Type III, with CRKSPH Q forces
      DvDti -= forceij/mRZi; //CRK Acceleration
      DvDtj += forceij/mRZj; //CRK Acceleration
      if (mCompatibleEnergyEvolution) {
        pairAccelerations[2*kk]   = -forceij/mRZi;
        pairAccelerations[2*kk+1] =  forceij/mRZj;
      }

      DepsDti += 0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + workQi)/mRZi;    // CRK Q
      DepsDtj += 0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + workQj)/mRZj;    // CRK Q

      // Estimate of delta v (for XSPH).
      if ((mXSPH and (nodeListi == nodeListj)) or min(zetai, zetaj) < 1.0) {
        XSPHDeltaVi -= weightj*Wj*vij;
        XSPHDeltaVj += weighti*Wi*vij;
      }

      // Add timing info for work
      const auto deltaTimePair = 0.5*Timing::difference(start, Timing::currentTime());
#pragma omp atomic
      nodeLists[nodeListi]->work()(i) += deltaTimePair;
#pragma omp atomic
      nodeLists[nodeListj]->work()(j) += deltaTimePair;
    }

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dim<2>>(threadStack);

  }   // OMP parallel

  // Finish up the derivatives for each point.
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0; i < ni; ++i) {

      // Get the state for node i.
      const auto& posi = position(nodeListi, i);
      const auto  ri = abs(posi.y());
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto  epsi = specificThermalEnergy(nodeListi, i);
      const auto  Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  zetai = abs((Hi*posi).y());
      const auto  hri = ri*safeInv(zetai);
      const auto  riInv = safeInv(ri, 0.25*hri);
      CHECK2(ri > 0.0, i << " " << ri);
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);

      // Time evolution of the mass density.
      const auto vri = vi.y(); // + XSPHDeltaVi.y();
      DrhoDti = -rhoi*(DvDxi.Trace() + vri*riInv);

      // Finish the specific thermal energy evolution.
      DepsDti -= Pi/rhoi*vri*riInv;

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
                                                             posi,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                        posi,
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
    }
  }
}

// //------------------------------------------------------------------------------
// // Finalize the hydro.
// //------------------------------------------------------------------------------
// void
// CRKSPHHydroBaseRZ::
// finalize(const Dim<2>::Scalar time,
//          const Dim<2>::Scalar dt,
//          DataBase<Dim<2> >& dataBase,
//          State<Dim<2> >& state,
//          StateDerivatives<Dim<2> >& derivs) {

//   // Convert the mass to mass per unit length first.
//   auto mass = state.fields(HydroFieldNames::mass, 0.0);
//   const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
//   const unsigned numNodeLists = mass.numFields();
//   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
//     const unsigned n = mass[nodeListi]->numElements();
//     for (unsigned i = 0; i != n; ++i) {
//       const auto circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
//       mass(nodeListi, i) /= circi;
//     }
//   }

//   // Base class finalization does most of the work.
//   CRKSPHHydroBase<Dimension>::finalize(time, dt, dataBase, state, derivs);

//   // Now convert back to true masses.
//   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
//     const unsigned n = mass[nodeListi]->numElements();
//     for (unsigned i = 0; i != n; ++i) {
//       const auto& xi = pos(nodeListi, i);
//       const auto circi = 2.0*M_PI*abs(xi.y());
//       mass(nodeListi, i) *= circi;
//     }
//   }
// }

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
CRKSPHHydroBaseRZ::
applyGhostBoundaries(State<Dim<2> >& state,
                     StateDerivatives<Dim<2> >& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary CRKSPH BCs.
  CRKSPHHydroBase<Dim<2> >::applyGhostBoundaries(state, derivs);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Scale back to mass.
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) *= circi;
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
CRKSPHHydroBaseRZ::
enforceBoundaries(State<Dim<2> >& state,
                  StateDerivatives<Dim<2> >& derivs) {

  // Convert the mass to mass/length before BCs are applied.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
  const unsigned numNodeLists = mass.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
      CHECK(circi > 0.0);
      mass(nodeListi, i) /= circi;
    }
  }

  // Apply ordinary CRKSPH BCs.
  CRKSPHHydroBase<Dim<2> >::enforceBoundaries(state, derivs);

  // Scale back to mass.
  // We also ensure no point approaches the z-axis too closely.
  FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = mass[nodeListi]->numInternalElements();
    const Scalar nPerh = mass[nodeListi]->nodeList().nodesPerSmoothingScale();
    for (unsigned i = 0; i != n; ++i) {
      Vector& posi = pos(nodeListi, i);
      const Scalar circi = 2.0*M_PI*abs(posi.y());
      mass(nodeListi, i) *= circi;
    }
  }
}

}
