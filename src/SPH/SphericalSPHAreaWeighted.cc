//---------------------------------Spheral++----------------------------------//
// SphericalSPHAreaWeighted -- An SPH/ASPH hydrodynamic package for Spheral++,
//                             specialized for 1D Spherical (r) geometry.
//
// This implementation uses linear-weighting (the spherical version of
// area-weighting in RZ cylindrical coordinates).n
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Thu Jun 29 13:59:18 PDT 2023
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
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"
#include "Utilities/Timer.hh"

#include "SphericalSPHAreaWeighted.hh"

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

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SphericalSPHAreaWeighted::
SphericalSPHAreaWeighted(const SmoothingScaleBase<Dim<1>>& smoothingScaleMethod,
                         DataBase<Dimension>& dataBase,
                         ArtificialViscosity<Dim<1>>& Q,
                         const TableKernel<Dim<1>>& W,
                         const TableKernel<Dim<1>>& WPi,
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
  mQself(0.0) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalSPHAreaWeighted::
~SphericalSPHAreaWeighted() {
}

//------------------------------------------------------------------------------
// Determine the timestep requirements for a hydro step.
//------------------------------------------------------------------------------
typename SphericalSPHAreaWeighted::TimeStepType
SphericalSPHAreaWeighted::
dt(const DataBase<Dim<1>>& dataBase,
   const State<Dim<1>>& state,
   const StateDerivatives<Dim<1>>& derivs,
   Dim<1>::Scalar currentTime) const {

  // Base timestep choice
  auto result = GenericHydro<Dim<1>>::dt(dataBase, state, derivs, currentTime);

  // Now enforce a constraint that no point can cross the axis in a single timestep
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto vel = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto numNodeLists = pos.numFields();
  CHECK(vel.numFields() == numNodeLists);
  for (auto k = 0u; k < numNodeLists; ++k) {
    const auto& fluidNodeList = **(dataBase.fluidNodeListBegin() + k);
    const auto  n = pos[k]->numInternalElements();
    const auto  rank = Process::getRank();
    for (auto i = 0u; i < n; ++i) {
      const auto& posi = pos(k,i);
      const auto& veli = vel(k,i);
      const auto  ri = std::abs(posi.x());
      const auto  vri = -std::min(0.0, veli.x());
      const auto  dti = 0.5*ri*safeInvVar(vri);
      if (dti < result.first) {
        result = std::make_pair(dti,
                                ("Axis crossing limit: dt = " + std::to_string(dti) + "\n" +
                                 "                     ri = " + std::to_string(ri) + "\n" +
                                 "                    vri = " + std::to_string(veli.x()) + "\n" +
                                 "               material = " + fluidNodeList.name() + "\n" +
                                 "  (nodeListID, i, rank) = (" + std::to_string(k) + " " + std::to_string(i) + " " + std::to_string(rank) + ")\n" +
                                 "             @ position = " + std::to_string(posi.x())));
      }
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SphericalSPHAreaWeighted::
registerState(DataBase<Dim<1>>& dataBase,
              State<Dim<1>>& state) {

  // The base class does most of it.
  SPHHydroBase<Dim<1>>::registerState(dataBase, state);

  // Re-register the position update to prevent things going through the origin
  auto position = dataBase.fluidPosition();
  auto positionPolicy = make_shared<SphericalPositionPolicy>();
  state.enroll(position, positionPolicy);

  // Are we using the compatible energy evolution scheme?
  // If so we need to override the ordinary energy registration with a specialized version.
  if (mCompatibleEnergyEvolution) {
    auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
    auto thermalEnergyPolicy = make_shared<NonSymmetricSpecificThermalEnergyPolicy<Dim<1>>>(dataBase,
                                                                                            [](const Scalar& mi, const Vector& posi) { return mi/FastMath::square(posi.x()); });
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  }
}

//------------------------------------------------------------------------------
// Stuff that occurs the beginning of a timestep
//------------------------------------------------------------------------------
void
SphericalSPHAreaWeighted::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // If we're going to do the SPH summation density, we need to convert the mass
  // to mass per unit area first.
  if (densityUpdate() == MassDensityType::RigorousSumDensity or
      densityUpdate() == MassDensityType::CorrectedSumDensity) {
    auto       mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
    const auto numNodeLists = mass.numFields();
    for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
      const auto n = mass[nodeListi]->numElements();
      for (auto i = 0u; i < n; ++i) {
        const auto areai = FastMath::square(pos(nodeListi, i).x());
        mass(nodeListi, i) /= areai;
      }
    }
  }

  // Base class finalization does most of the work.
  SPHHydroBase<Dimension>::preStepInitialize(dataBase, state, derivs);

  // Now convert back to true masses and mass densities.  We also apply the RZ
  // correction factor to the mass density.
  if (densityUpdate() == MassDensityType::RigorousSumDensity or
      densityUpdate() == MassDensityType::CorrectedSumDensity) {
    auto       mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
    const auto numNodeLists = mass.numFields();
    for (auto nodeListi = 0u; nodeListi != numNodeLists; ++nodeListi) {
      const auto n = mass[nodeListi]->numElements();
      for (auto i = 0u; i < n; ++i) {
        const auto areai = FastMath::square(pos(nodeListi, i).x());
        mass(nodeListi, i) *= areai;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SphericalSPHAreaWeighted::
evaluateDerivatives(const Dim<1>::Scalar time,
                    const Dim<1>::Scalar dt,
                    const DataBase<Dim<1>>& dataBase,
                    const State<Dim<1>>& state,
                    StateDerivatives<Dim<1>>& derivs) const {
  TIME_BEGIN("SphericalSPHevalDerivs");
  TIME_BEGIN("SphericalSPHevalDerivs_initial");

  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();
  const auto  oneKernel = (W == WQ);

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  const auto W0 = W(0.0, 1.0);

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
  // const auto effectiveRadius = state.fields(HydroFieldNames::reff, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);
  // CHECK(effectiveRadius.size() == numNodeLists);

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
  auto  viscousWork = derivs.fields(HydroFieldNames::viscousWork, 0.0);
  auto& pairAccelerations = derivs.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto  XSPHWeightSum = derivs.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivs.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivs.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivs.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
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
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  const auto nnodes = dataBase.numFluidInternalNodes();
  if (mCompatibleEnergyEvolution) pairAccelerations.resize(2u*npairs + nnodes);

  // // Assume all NodeLists have the same nPerh;
  // const auto& nodeList = mass[0]->nodeList();
  // const auto  nPerh = nodeList.nodesPerSmoothingScale();

  // Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
    Vector gradWi, gradWj, gradWQi, gradWQj;
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
    auto viscousWork_thread = viscousWork.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& posi = position(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  ri = std::abs(posi.x());
      const auto  retai = abs((Hi*posi).x());
      // const auto  hri = abs(posi.x())*safeInvVar(retai);
      const auto  mi = mass(nodeListi, i);
      const auto  mRi = mi*FastMath::square(safeInvVar(ri));
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& omegai = omega(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  safeOmegai = safeInv(omegai, tiny);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(retai > 0.0);
      // CHECK(hri > 0.0);
      CHECK(ri > 0.0);

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
      auto& viscousWorki = viscousWork_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);

      // Get the state for node j
      const auto& posj = position(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  rj = std::abs(posj.x());
      const auto  retaj = abs((Hj*posj).x());
      // const auto  hrj = abs(posj.x())*safeInvVar(retaj);
      const auto  mj = mass(nodeListj, j);
      const auto  mRj = mj*FastMath::square(safeInvVar(rj));
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& omegaj = omega(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      const auto  safeOmegaj = safeInv(omegaj, tiny);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);
      CHECK(retaj > 0.0);
      // CHECK(hrj > 0.0);
      CHECK(rj > 0.0);

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
      auto& viscousWorkj = viscousWork_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij = true; // (nodeListi == nodeListj and fragIDi == fragIDj);

      // Node displacement.
      const auto xij = posi - posj;
      const auto etai = Hi*xij;
      const auto etaj = Hj*xij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      const auto etaUnit = etai*safeInvVar(etaMagi);
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      W.kernelAndGradValue(etaMagi, Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaMagj, Hdetj, Wj, gWj);
      gradWi = gWi*Hi*etaUnit;
      gradWj = gWj*Hj*etaUnit;
      if (oneKernel) {
        WQi = Wi;
        WQj = Wj;
        gradWQi = gradWi;
        gradWQj = gradWj;
      } else {
        WQ.kernelAndGradValue(etaMagi, Hdeti, WQi, gWQi);
        WQ.kernelAndGradValue(etaMagj, Hdetj, WQj, gWQj);
        gradWQi = gWQi*Hi*etaUnit;
        gradWQj = gWQj*Hj*etaUnit;
      }

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      const auto fweightij = sameMatij ? 1.0 : mRj*rhoi/(mRi*rhoj);
      const auto xij2 = xij.magnitude2();
      const auto thpt = xij.selfdyad()*safeInvVar(xij2*xij2*xij2);
      weightedNeighborSumi +=     fweightij*std::abs(gWi);
      weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
      massSecondMomenti +=     fweightij*gradWi.magnitude2()*thpt;
      massSecondMomentj += 1.0/fweightij*gradWj.magnitude2()*thpt;

      // Contribution to the sum density.
      if (nodeListi == nodeListj) {
        rhoSumi += mRj*Wi;
        rhoSumj += mRi*Wj;
        normi += mRi/rhoi*Wi;
        normj += mRj/rhoj*Wj;
      }

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                      posi, etai, vi, rhoi, ci, Hi,
                                      posj, etaj, vj, rhoj, cj, Hj);
      const auto Qacci = 0.5*(QPiij*gradWQi);
      const auto Qaccj = 0.5*(QPiji*gradWQj);
      // const auto workQi = 0.5*(QPiij*vij).dot(gradWQi);
      // const auto workQj = 0.5*(QPiji*vij).dot(gradWQj);
      const auto workQi = vij.dot(Qacci);
      const auto workQj = vij.dot(Qaccj);
      const auto Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
      const auto Qj = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
      maxViscousPressurei = max(maxViscousPressurei, Qi);
      maxViscousPressurej = max(maxViscousPressurej, Qj);
      effViscousPressurei += mRj*Qi*WQi/rhoj;
      effViscousPressurej += mRi*Qj*WQj/rhoi;
      viscousWorki += mRj*workQi;
      viscousWorkj += mRi*workQj;

      // Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      const auto Prhoi = safeOmegai*Pi/(rhoi*rhoi);
      const auto Prhoj = safeOmegaj*Pj/(rhoj*rhoj);
      const auto deltaDvDt = Prhoi*gradWi + Prhoj*gradWj + Qacci + Qaccj;
      DvDti -= mRj*deltaDvDt;
      DvDtj += mRi*deltaDvDt;
      if (mCompatibleEnergyEvolution) {
        pairAccelerations[2*kk]   = -mRj*deltaDvDt;
        pairAccelerations[2*kk+1] =  mRi*deltaDvDt;
      }

      // Specific thermal energy evolution.
      DepsDti += mRj*(Prhoi*vij.dot(gradWi) + workQi);
      DepsDtj += mRi*(Prhoj*vij.dot(gradWj) + workQj);

      // Velocity gradient.
      const auto deltaDvDxi = mRj*vij.dyad(gradWi);
      const auto deltaDvDxj = mRi*vij.dyad(gradWj);
      DvDxi -= deltaDvDxi;
      DvDxj -= deltaDvDxj;
      if (sameMatij) {
        localDvDxi -= deltaDvDxi;
        localDvDxj -= deltaDvDxj;
      }

      // Estimate of delta v (for XSPH).
      if (sameMatij or min(retai, retaj) < 1.0) {
        const auto wXSPHij = 0.5*(mRi/rhoi*Wi + mRj/rhoj*Wj);
        XSPHWeightSumi += wXSPHij;
        XSPHWeightSumj += wXSPHij;
        XSPHDeltaVi -= wXSPHij*vij;
        XSPHDeltaVj += wXSPHij*vij;
      }

      // Linear gradient correction term.
      Mi -= mRj*xij.dyad(gradWi);
      Mj -= mRi*xij.dyad(gradWj);
      if (sameMatij) {
        localMi -= mRj*xij.dyad(gradWi);
        localMj -= mRi*xij.dyad(gradWj);
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region


  // Finish up the derivatives for each point.
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
      const auto& posi = position(nodeListi, i);
      const auto  ri = std::abs(posi.x());
      const auto& Hi = H(nodeListi, i);
      const auto  retai = abs((Hi*posi).x());
      const auto  hri = abs(posi.x())*safeInvVar(retai);
      const auto  mi = mass(nodeListi, i);
      const auto  areaInvi = FastMath::square(safeInvVar(ri));
      const auto  mRi = mi*areaInvi;
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  riInv = safeInvVar(ri);
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);
      CHECK(retai > 0.0);
      CHECK(hri > 0.0);

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
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);

      // Add the self-contribution to density sum.
      rhoSumi += mRi*W0*Hdeti;
      rhoSumi *= areaInvi;
      normi += mRi/rhoi*W0*Hdeti;

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

      // Finish the continuity equation.
      XSPHWeightSumi += Hdeti*mRi/rhoi*W0;
      CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      XSPHDeltaVi /= XSPHWeightSumi;
      const auto vri = vi.x(); // + XSPHDeltaVi.x();
      DrhoDti = -rhoi*(DvDxi.Trace() + 2.0*vri*riInv);

      // Extra hoop terms for the artificial viscosity
      const auto mu = std::max(0.0, -vri*riInv);
      CHECK(mu >= 0.0);
      const auto Qi = mQself*rhoi*hri*mu*(hri*mu + ci) * std::max(0.0,  W(retai, 1.0)/W0);
      VERIFY(Qi >= 0.0);
      maxViscousPressurei = std::max(maxViscousPressurei, Qi);
      const auto deltaDvDti = Qi/rhoi*riInv;
      DvDti[0] += deltaDvDti;
      if (mCompatibleEnergyEvolution) pairAccelerations[offset + i] = deltaDvDti;

      // Finish the specific thermal energy evolution.
      DepsDti -= (Pi + Qi)/rhoi*2.0*vri*riInv;

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

  TIME_END("SphericalSPHevalDerivs_final");
  TIME_END("SphericalSPHevalDerivs");
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SphericalSPHAreaWeighted::
applyGhostBoundaries(State<Dim<1>>& state,
                     StateDerivatives<Dim<1>>& derivs) {

  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::square(pos(nodeListi, i).x());
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
      const auto dA = FastMath::square(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SphericalSPHAreaWeighted::
enforceBoundaries(State<Dim<1>>& state,
                  StateDerivatives<Dim<1>>& derivs) {

  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::square(pos(nodeListi, i).x());
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
      const auto dA = FastMath::square(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

//------------------------------------------------------------------------------
// The self-Q multiplier
//------------------------------------------------------------------------------
double
SphericalSPHAreaWeighted::
Qself() const {
  return mQself;
}

void
SphericalSPHAreaWeighted::
Qself(const double x) {
  mQself = x;
}

}
