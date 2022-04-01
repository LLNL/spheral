//---------------------------------Spheral++----------------------------------//
// SphericalSPHHydroBase -- An SPH/ASPH hydrodynamic package for Spheral++,
//                          specialized for 1D Spherical (r) geometry.
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Tue Dec 22 10:04:21 PST 2020
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "computeSPHSumMassDensity.hh"
#include "correctSPHSumMassDensity.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "Hydro/SphericalPositionPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/Timer.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"

#include "SphericalSPHHydroBase.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
#include <memory>
using std::vector;
using std::map;
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

// Declare timers
extern Timer TIME_SPH;
extern Timer TIME_SPHinitializeStartup;
extern Timer TIME_SPHregister;
extern Timer TIME_SPHregisterDerivs;
extern Timer TIME_SPHpreStepInitialize;
extern Timer TIME_SPHinitialize;
extern Timer TIME_SPHfinalizeDerivs;
extern Timer TIME_SPHghostBounds;
extern Timer TIME_SPHupdateVol;
extern Timer TIME_SPHenforceBounds;
extern Timer TIME_SPHevalDerivs;
extern Timer TIME_SPHevalDerivs_initial;
extern Timer TIME_SPHevalDerivs_pairs;
extern Timer TIME_SPHevalDerivs_final;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SphericalSPHHydroBase::
SphericalSPHHydroBase(const SmoothingScaleBase<Dim<1>>& smoothingScaleMethod,
                      DataBase<Dimension>& dataBase,
                      ArtificialViscosity<Dim<1>>& Q,
                      const SphericalKernel& W,
                      const SphericalKernel& WPi,
                      const double filter,
                      const double cfl,
                      const bool useVelocityMagnitudeForDt,
                      const bool compatibleEnergyEvolution,
                      const bool evolveTotalEnergy,
                      const bool gradhCorrection,
                      const bool XSPH,
                      const bool correctVelocityGradient,
                      const bool sumMassDensityOverAllNodeLists,
                      const MassDensityType densityUpdate,
                      const HEvolutionType HUpdate,
                      const double epsTensile,
                      const double nTensile,
                      const Vector& xmin,
                      const Vector& xmax):
  SPHHydroBase<Dim<1>>(smoothingScaleMethod,
                       dataBase,
                       Q,
                       W.baseKernel1d(),
                       WPi.baseKernel1d(),
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
  mKernel(W),
  mPiKernel(WPi) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalSPHHydroBase::
~SphericalSPHHydroBase() {
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
registerState(DataBase<Dim<1>>& dataBase,
              State<Dim<1>>& state) {

  // The base class does most of it.
  SPHHydroBase<Dim<1>>::registerState(dataBase, state);

  // Re-regsiter the position update to prevent things going through the origin
    auto position = dataBase.fluidPosition();
  auto positionPolicy = make_shared<SphericalPositionPolicy>();
  state.enroll(position, positionPolicy);

  // Are we using the compatible energy evolution scheme?
  // If so we need to override the ordinary energy registration with a specialized version.
  if (mCompatibleEnergyEvolution) {
    auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
    auto thermalEnergyPolicy = make_shared<NonSymmetricSpecificThermalEnergyPolicy<Dim<1>>>(dataBase);
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // We have to override the density summation to use the correct kernel
  switch(densityUpdate()) {

  case MassDensityType::IntegrateDensity:
    break;

  case MassDensityType::RigorousSumDensity:
  case MassDensityType::CorrectedSumDensity:
    {
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
      auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      if (densityUpdate() == MassDensityType::CorrectedSumDensity) {
        correctSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
        for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
        for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      }
    }
    break;

  case MassDensityType::SumDensity:
    {
      auto       massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      const auto massDensitySum = derivs.fields(ReplaceFieldList<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                HydroFieldNames::massDensity, 0.0);
      massDensity.assignFields(massDensitySum);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
    }
    break;

  default:
    VERIFY2(false, "Unsupported mass density definition for Spherical SPH");
    break;
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
evaluateDerivatives(const Dim<1>::Scalar time,
                    const Dim<1>::Scalar dt,
                    const DataBase<Dim<1>>& dataBase,
                    const State<Dim<1>>& state,
                    StateDerivatives<Dim<1>>& derivs) const {
  TIME_SPHevalDerivs.start();
  TIME_SPHevalDerivs_initial.start();

  //static double totalLoopTime = 0.0;
  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = mKernel;
  const auto& WQ = mPiKernel;
  const auto  oneKernel = (W == WQ);
  const auto  etaMax = W.etamax();

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);

  // Derivative FieldLists.
  auto  rhoSum = derivs.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  normalization = derivs.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivs.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivs.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  M = derivs.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivs.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  DHDt = derivs.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivs.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivs.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto& pairAccelerations = derivs.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto  XSPHWeightSum = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivs.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  CHECK(rhoSum.size() == numNodeLists);
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
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  const auto nnodes = dataBase.numFluidInternalNodes();
  if (mCompatibleEnergyEvolution) pairAccelerations.resize(2u*npairs + nnodes);
  TIME_SPHevalDerivs_initial.stop();

  Vector gradSum_check;

  // Walk all the interacting pairs.
  TIME_SPHevalDerivs_pairs.start();
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Vector gradWii, gradWij, gradWQii, gradWQij;
    Vector gradWji, gradWjj, gradWQji, gradWQjj;
    Scalar Wii, gWii, WQii, gWQii, Wji, gWji, WQji, gWQji;
    Scalar Wij, gWij, WQij, gWQij, Wjj, gWjj, WQjj, gWQjj;
    Tensor QPiij, QPiji;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto rhoSum_thread = rhoSum.threadCopy(threadStack);
    auto normalization_thread = normalization.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& omegai = 1.0; //omega(nodeListi, i);
      const auto  safeOmegai = safeInv(omegai, tiny);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);

      auto& rhoSumi = rhoSum_thread(nodeListi, i);
      auto& normi = normalization_thread(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& Mi = M_thread(nodeListi, i);
      auto& localMi = localM_thread(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& omegaj = 1.0; // omega(nodeListj, j);
      const auto  safeOmegaj = safeInv(omegaj, tiny);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);

      auto& rhoSumj = rhoSum_thread(nodeListj, j);
      auto& normj = normalization_thread(nodeListj, j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& Mj = M_thread(nodeListj, j);
      auto& localMj = localM_thread(nodeListj, j);
      auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

      // Normalized node coordinates
      const auto etaii = Hi*ri;
      const auto etaji = Hi*rj;
      const auto etaij = Hj*ri;
      const auto etajj = Hj*rj;

      // Symmetrized kernel weight and gradient.
      // Note our subscript rules here:
      //   first subscript -> first node position in call
      //  second subscript -> smoothing scale used
      W.kernelAndGrad(etaji, etaii, Hi, Wji, gradWji, gWji);
      W.kernelAndGrad(etajj, etaij, Hj, Wjj, gradWjj, gWjj);
      W.kernelAndGrad(etaii, etaji, Hi, Wii, gradWii, gWii);
      W.kernelAndGrad(etaij, etajj, Hj, Wij, gradWij, gWij);
      if (oneKernel) {
        WQii = Wii;
        WQij = Wij;
        WQji = Wji;
        WQjj = Wjj;
        gradWQii = gradWii;
        gradWQij = gradWij;
        gradWQji = gradWji;
        gradWQjj = gradWjj;
      } else {
        WQ.kernelAndGrad(etaji, etaii, Hi, WQji, gradWQji, gWQji);
        WQ.kernelAndGrad(etajj, etaij, Hj, WQjj, gradWQjj, gWQjj);
        WQ.kernelAndGrad(etaii, etaji, Hi, WQii, gradWQii, gWQii);
        WQ.kernelAndGrad(etaij, etajj, Hj, WQij, gradWQij, gWQij);
      }

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.  Note in 1D there is no second moment.
      const auto fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      const auto rij = ri - rj;
      const auto rij2 = rij.magnitude2();
      const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      weightedNeighborSumi +=     fweightij*std::abs(gWii) * (etaii > etaMax ? 1.0 :
                                                              etaii < etaji ? 2.0 :
                                                              0.0);
      weightedNeighborSumj += 1.0/fweightij*std::abs(gWij) * (etajj > etaMax ? 1.0 :
                                                              etajj < etaij ? 2.0 :
                                                              0.0);

      // Contribution to the sum density.
      if (nodeListi == nodeListj) {
        rhoSumi += mj*Wii;
        rhoSumj += mi*Wij;
        normi += mi/rhoi*Wii;
        normj += mj/rhoj*Wij;
      }

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                      ri, etaii - etaji, vi, rhoi, ci, Hi,
                                      rj, etaij - etajj, vj, rhoj, cj, Hj);
      const auto Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
      const auto Qj = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
      maxViscousPressurei = max(maxViscousPressurei, Qi);
      maxViscousPressurej = max(maxViscousPressurej, Qj);
      effViscousPressurei += mj*Qi*WQii/rhoj;
      effViscousPressurej += mi*Qj*WQij/rhoi;

      //Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      const auto Prhoi = safeOmegai*Pi/(rhoi*rhoi);
      const auto Prhoj = safeOmegaj*Pj/(rhoj*rhoj);
      const auto fQi = ri.x()*safeInv(ri.x() + rj.x());
      const auto deltaDvDti = -mj*(Prhoi*gradWji + Prhoj*gradWjj + fQi        *(QPiij*gradWQji + QPiji*gradWQjj));
      const auto deltaDvDtj = -mi*(Prhoj*gradWij + Prhoi*gradWii + (1.0 - fQi)*(QPiji*gradWQij + QPiij*gradWQii));
      DvDti += deltaDvDti;
      DvDtj += deltaDvDtj;
      if (mCompatibleEnergyEvolution) {
        pairAccelerations[2*kk] = deltaDvDti;
        pairAccelerations[2*kk+1] = deltaDvDtj;
      }

      // Specific thermal energy evolution
      // DepsDti -= vi.dot(deltaDvDti);
      // DepsDtj -= vj.dot(deltaDvDtj);
      // DepsDti += mj*(Prhoi + fQi        *QPiij.xx())*vij.dot(gradWjj);
      // DepsDtj -= mi*(Prhoj + (1.0 - fQi)*QPiji.xx())*vij.dot(gradWii);
      DepsDti += mj*(Prhoi + 0.5*QPiij.xx())*(vj.dot(gradWij) + vi.dot(gradWjj));
      DepsDtj += mi*(Prhoj + 0.5*QPiji.xx())*(vi.dot(gradWji) + vj.dot(gradWii));

      // Velocity gradient.
      const auto deltaDvDxi = -mj*vij.dyad(gradWjj);
      const auto deltaDvDxj =  mi*vij.dyad(gradWii);
      DvDxi += deltaDvDxi; 
      DvDxj += deltaDvDxj;
      if (sameMatij) {
        localDvDxi += deltaDvDxi; 
        localDvDxj += deltaDvDxj;
      }

      // Estimate of delta v (for XSPH).
      if (mXSPH and (sameMatij)) {
        const auto wXSPHij = 0.5*(mi/rhoi*Wji + mj/rhoj*Wij);
        XSPHWeightSumi += wXSPHij;
        XSPHWeightSumj += wXSPHij;
        XSPHDeltaVi -= wXSPHij*vij;
        XSPHDeltaVj += wXSPHij*vij;
      }

      // Linear gradient correction term.
      const auto deltaMi = -mj*rij.dyad(gradWjj);
      const auto deltaMj =  mi*rij.dyad(gradWii);
      Mi += deltaMi;
      Mj += deltaMj;
      if (sameMatij) {
        localMi += deltaMi;
        localMj += deltaMj;
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  TIME_SPHevalDerivs_pairs.stop();

  // Finish up the derivatives for each point.
  TIME_SPHevalDerivs_final.start();
  size_t offset = 2u*npairs;
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
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& omegai = 1.0; // omega(nodeListi, i);
      const auto  safeOmegai = safeInv(omegai, tiny);
      const auto  Hdeti = Hi.Determinant();
      const auto  hi = 1.0/Hdeti;
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& rhoSumi = rhoSum(nodeListi, i);
      auto& normi = normalization(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);

      // Symmetrized kernel weight and gradient.
      const auto etaii = Hi*ri;
      const Vector etaQii = std::max(0.01, etaii[0]);
      double Wii, gWii, WQii, gWQii;
      Vector gradWii, gradWQii;
      W.kernelAndGrad(etaii, etaii, Hi, Wii, gradWii, gWii);
      WQ.kernelAndGrad(etaQii, etaQii, Hi, WQii, gradWQii, gWQii);

      // Add the self-contribution to density sum.
      rhoSumi += mi*Wii;
      normi += mi/rhoi*Wii;

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
      DrhoDti = -rhoi*(DvDxi.xx() + 2.0*vi.x()*safeInvVar(ri.x(), 0.1*hi));

      // If we're in range of the origin, compute an effective Q.
      Scalar Qi = 0.0;
      if (etaii.x() < etaMax) {
        const auto divv = -std::min(0.0, DvDxi.Trace() + 2.0*vi.x()*safeInvVar(ri.x(), 0.01*hi));
        Qi = 2.0*rhoi*hi*divv*(hi*divv + ci);
        maxViscousPressurei = std::max(maxViscousPressurei, Qi);
      }

      // Self-interaction for momentum (cause curvilinear coordinates are weird)
      const auto deltaDvDti = -mi*safeOmegai/(rhoi*rhoi)*(2.0*Pi*gradWii + Qi*gradWQii);
      DvDti += deltaDvDti;
      if (mCompatibleEnergyEvolution) pairAccelerations[offset + i] = deltaDvDti;
      // if (i == 0) {
      //   const auto veli = velocity(nodeListi, i);
      //   std::cerr << "   " << ri << " " << veli << " " << Pi << " " << Qi << " " << rhoi << " " << gradWii << " " << gradWQii << " : " << deltaDvDti << " " << DvDti << std::endl;
      // }

      // Specific thermal energy
      DepsDti -= vi.dot(deltaDvDti) + 2.0*Pi*vi.x()*safeInvVar(ri.x(), 0.01*hi);
      // DepsDti += vi.dot(deltaDvDti) + 0.5*deltaDvDti.dot(deltaDvDti)*dt - 2.0*Pi*vi.x()*safeInvVar(ri.x(), 0.1*hi);
      // DepsDti += 0.5*mi*QPiij.xx()*vi.dot(gradWii) - 2.0*Pi*vi.x()*safeInvVar(ri.x(), 0.1*hi);

      // If needed finish the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        XSPHWeightSumi += mi/rhoi*Wii;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
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
                                                        SymTensor::zero,
                                                        W.baseKernel1d(),
                                                        hmin,
                                                        hmax,
                                                        hminratio,
                                                        nPerh,
                                                        connectivityMap,
                                                        nodeListi,
                                                        i);
    }
    offset += ni;
  }
  CHECK(offset == 2u*npairs + nnodes);
  TIME_SPHevalDerivs_final.stop();
  TIME_SPHevalDerivs.stop();
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
applyGhostBoundaries(State<Dim<1>>& state,
                     StateDerivatives<Dim<1>>& derivs) {

  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) /= dA;
    }
  }

  // Apply ordinary SPH BCs.
  SPHHydroBase<Dim<1>>::applyGhostBoundaries(state, derivs);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Scale back to mass.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
enforceBoundaries(State<Dim<1>>& state,
                  StateDerivatives<Dim<1>>& derivs) {

  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) /= dA;
    }
  }

  // Apply ordinary SPH BCs.
  SPHHydroBase<Dim<1>>::enforceBoundaries(state, derivs);

  // Scale back to mass.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
const SphericalKernel&
SphericalSPHHydroBase::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// Access the kernel used for artificial viscosity gradients.
//------------------------------------------------------------------------------
const SphericalKernel&
SphericalSPHHydroBase::
PiKernel() const {
  return mPiKernel;
}

}
