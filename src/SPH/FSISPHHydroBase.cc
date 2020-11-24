//---------------------------------Spheral++----------------------------------//
// SPHHydroBase -- SPHHydro modified for large density discontinuities
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "computeSPHSumMassDensity.hh"
#include "correctSPHSumMassDensity.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
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
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
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
#include "SPHHydroBase.hh"

#include "FSISPHHydroBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

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
template<typename Dimension>
FSISPHHydroBase<Dimension>::
FSISPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                  DataBase<Dimension>& dataBase,
                  ArtificialViscosity<Dimension>& Q,
                  const TableKernel<Dimension>& W,
                  const TableKernel<Dimension>& WPi,
                  const double alpha,
                  const std::vector<int> sumDensityNodeListSwitch,
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
  SPHHydroBase<Dimension>(smoothingScaleMethod,
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
  mAlpha(alpha),
  mSumDensityNodeListSwitch(sumDensityNodeListSwitch),
  mSurfaceNormal(FieldStorageType::CopyFields) {
  mSurfaceNormal = dataBase.newFluidFieldList(Vector::zero,"surface_Normal_between_materials");
  }

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
FSISPHHydroBase<Dimension>::
~FSISPHHydroBase() {
}


//------------------------------------------------------------------------------
// FSI specialized density summmation
//------------------------------------------------------------------------------
template<typename Dimension>
void
FSISPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  SPHHydroBase<Dimension>::preStepInitialize(dataBase,state,derivs);

  // test to see if any of the nodeLists require their own density sum
  bool applySelectDensitySum = false;
  auto numNodeLists = dataBase.numNodeLists();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    if (mSumDensityNodeListSwitch[nodeListi]==1){
      applySelectDensitySum = true;
    }
  }
  if (applySelectDensitySum){
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
      const auto& W = this->kernel();
      auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      this->computeFSISPHSumMassDensity(connectivityMap, W, position, mass, H, massDensity);
  }

}


//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
  TIME_SPHevalDerivs.start();
  TIME_SPHevalDerivs_initial.start();

  //static double totalLoopTime = 0.0;
  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);

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
  auto  rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  auto  effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  auto  viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto  XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
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
  if (mCompatibleEnergyEvolution) pairAccelerations.resize(npairs);

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  TIME_SPHevalDerivs_initial.stop();

  // Walk all the interacting pairs.
  TIME_SPHevalDerivs_pairs.start();


    // Walk all the interacting pairs and calculated linear correction tensor
  // based on taylor series expansion of gradient estimate
#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    //auto rSumDiff_thread = rSumDiff.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      const auto start = Timing::currentTime();
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto  mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto  rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto Ki = K(nodeListi,i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& Mi = M_thread(nodeListi, i);
      auto& localMi = localM_thread(nodeListi, i);
      //auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      //auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      //auto& rSumDiffi = rSumDiff(nodeListi, i);

      // Get the state for node j
      const auto rj = position(nodeListj, j);
      const auto mj = mass(nodeListj, j);
      const auto vj = velocity(nodeListj, j);
      const auto rhoj = massDensity(nodeListj, j);
      const auto Hj = H(nodeListj, j);
      const auto Hdetj = Hj.Determinant();
      const auto Kj = K(nodeListj,j);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);


      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& Mj = M_thread(nodeListj, j);
      auto& localMj = localM_thread(nodeListj, j);
      //auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      //auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      //auto& rSumDiffj = rSumDiff(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const auto sameMatij =  (nodeListi == nodeListj);

      // Node displacement.
      const auto rij = ri - rj;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = W.gradValue(etaMagi, Hdeti)*Hetai;

      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = W.gradValue(etaMagj, Hdetj)*Hetaj;

      // Compute the pair-wise artificial viscosity.
      const auto vij = vi - vj;
      const auto voli = mi/rhoi;
      const auto volj = mj/rhoj;

      // blending parameter or FSI based on bulk modulus
      auto kappai = 1.0;
      auto kappaj = 1.0;
      //auto kappaMax = 1.0;

      // treatment of FSI
      if(!sameMatij){
        kappai = 2.0*(Kj)/(Ki+Kj);
        kappaj = max(0.0,min(2.0,(2.0-kappai)));
        //rSumDiffi -= volj*rij;
        //rSumDiffj += voli*rij;
      }
      const auto gradRi = volj*rij.dyad(gradWi);
      const auto gradRj = voli*rij.dyad(gradWj);
      // Linear gradient correction term.
      Mi -= min(kappai,1.0)*gradRi;
      Mj -= min(kappaj,1.0)*gradRj;
      if (sameMatij) {
        localMi -= gradRi;
        localMj -= gradRj;
      }
      
      // not functional yet
      //if (XSPH) {
          //const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
      //    XSPHWeightSumi += kappai*volj*Wi;
      //    XSPHWeightSumj += kappaj*voli*Wj;
      //    XSPHDeltaVi -= kappai*vij*volj*Wi;//wXSPHij*vij;
      //    XSPHDeltaVj += kappaj*vij*voli*Wj;//wXSPHij*vij;
      //}

      // Add timing info for work
      const auto deltaTimePair = 0.5*Timing::difference(start, Timing::currentTime());
#pragma omp atomic
      nodeLists[nodeListi]->work()(i) += deltaTimePair;
#pragma omp atomic
      nodeLists[nodeListj]->work()(j) += deltaTimePair;

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region

 
 // Loop through nodes and invert the linear correction matrix
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {


      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);

      //auto& DvDxi = DvDx(nodeListi, i);
      //auto& localDvDxi = localDvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);

      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
      } else { 
        Mi = Tensor::one; 
      } 
      if (this->mCorrectVelocityGradient and
          std::abs(localMi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
          localMi = localMi.Inverse();
      } else { 
        localMi = Tensor::one;  
      }

    }
  }




#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj;
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
      const auto start = Timing::currentTime();
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
      const auto& omegai = omega(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  safeOmegai = safeInv(omegai, tiny);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

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
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& omegaj = omega(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      const auto  safeOmegaj = safeInv(omegaj, tiny);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

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
      const bool sameMatij =  (nodeListi == nodeListj); // and fragIDi == fragIDj);
      const auto kappai = (sameMatij ? 1.0 : max(0.0,min(2.0, 2.0*(Kj)/(Ki+Kj))));
      const auto kappaj = (sameMatij ? 1.0 : max(0.0,min(2.0, (2.0-kappai))));
      const auto kappaMax = max(max(kappai,kappaj),1.0);

      // Node displacement.
      const auto rij = ri - rj;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      std::tie(Wi, gWi) = W.kernelAndGradValue(etaMagi, Hdeti);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      std::tie(Wj, gWj) = W.kernelAndGradValue(etaMagj, Hdetj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      //const auto fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      const auto rij2 = rij.magnitude2();
      const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      weightedNeighborSumi += std::abs(gWi);
      weightedNeighborSumj += std::abs(gWj);
      massSecondMomenti += gradWi.magnitude2()*thpt;
      massSecondMomentj += gradWj.magnitude2()*thpt;

      // Contribution to the sum density.
      //if (nodeListi == nodeListj) {
      //  rhoSumi += mj*Wi;
      //  rhoSumj += mi*Wj;
      //  normi += mi/rhoi*Wi;
      //  normj += mj/rhoj*Wj;
      //}

      // Compute the pair-wise artificial viscosity.
      const auto rhoij = (rhoi+rhoj)/2.0;
      const auto cij = (ci+cj)/2.0;
      const auto vij = vi - vj;

      std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                      ri, etai, vi, rhoij, cij, Hi,
                                      rj, etaj, vj, rhoij, cij, Hj);

      //const auto Qacci = 0.5*(QPiij*gradWi);
      //const auto Qaccj = 0.5*(QPiji*gradWj);
      // const auto workQi = 0.5*(QPiij*vij).dot(gradWQi);
      // const auto workQj = 0.5*(QPiji*vij).dot(gradWQj);
      //const auto workQi = vij.dot(Qacci);
      //const auto workQj = vij.dot(Qaccj);
      //const auto Qi = rhoi*rhoi*(QPiij.diagonalElements().maxAbsElement());
      //const auto Qj = rhoj*rhoj*(QPiji.diagonalElements().maxAbsElement());
      //maxViscousPressurei = max(maxViscousPressurei, Qi);
      //maxViscousPressurej = max(maxViscousPressurej, Qj);
      //effViscousPressurei += mj*Qi*WQi/rhoj;
      //effViscousPressurej += mi*Qj*WQj/rhoi;
      //viscousWorki += mj*workQi;
      //viscousWorkj += mi*workQj;

      // Determine an effective pressure including a term to fight the tensile instability.
//             const auto fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
      const auto fij = mEpsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
      const auto Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
      const auto Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
      const auto Peffi = Pi + Ri;
      const auto Peffj = Pj + Rj;

      // Acceleration.
      const auto Prhoi = Peffi/(pow(rhoi,alpha)*pow(rhoj,2.0-alpha))*Tensor::one+0.5*(QPiij);
      const auto Prhoj = Peffj/(pow(rhoj,alpha)*pow(rhoi,2.0-alpha))*Tensor::one+0.5*(QPiji);
      const auto deltaDvDt = safeOmegai*Prhoi*gradWi + safeOmegaj*Prhoj*gradWj;
      DvDti -= mj*deltaDvDt;
      DvDtj += mi*deltaDvDt;

      if (mCompatibleEnergyEvolution) pairAccelerations[kk] = -mj*deltaDvDt;  // Acceleration for i (j anti-symmetric)

      // velocity gradient with FSI treatment
      const auto voli = mi/rhoi;
      const auto volj = mj/rhoj;

      auto deltaDvDxi = kappai*vij.dyad(gradWi);
      auto deltaDvDxj = kappaj*vij.dyad(gradWj);

      if (sameMatij) {
        localDvDxi -= volj*(deltaDvDxi)*localMi;
        localDvDxj -= voli*(deltaDvDxj)*localMj;
      }

      deltaDvDxi = deltaDvDxi*Mi;
      deltaDvDxj = deltaDvDxj*Mj;

      DvDxi -= volj*(deltaDvDxi);
      DvDxj -= voli*(deltaDvDxj);

      // DepsDt eqn -- possible remove symmetric?
      DepsDti -= mj*(sigmarhoi.doubledot(deltaDvDxi.Symmetric()));
      DepsDtj -= mi*(sigmarhoj.doubledot(deltaDvDxj.Symmetric()));
      

      // Estimate of delta v (for XSPH).
      if (mXSPH and (sameMatij)) {
        const auto wXSPHij = 0.5*(mi/rhoi*Wi + mj/rhoj*Wj);
        XSPHWeightSumi += wXSPHij;
        XSPHWeightSumj += wXSPHij;
        XSPHDeltaVi -= wXSPHij*vij;
        XSPHDeltaVj += wXSPHij*vij;
      }

      // Add timing info for work
      const auto deltaTimePair = 0.5*Timing::difference(start, Timing::currentTime());
#pragma omp atomic
      nodeLists[nodeListi]->work()(i) += deltaTimePair;
#pragma omp atomic
      nodeLists[nodeListj]->work()(j) += deltaTimePair;

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  TIME_SPHevalDerivs_pairs.stop();

  // Finish up the derivatives for each point.
  TIME_SPHevalDerivs_final.start();
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
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
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
      //auto& Mi = M(nodeListi, i);
      //auto& localMi = localM(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;
      normi += mi/rhoi*W0*Hdeti;

      // Evaluate the continuity equation.
      DrhoDti = -rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
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
  TIME_SPHevalDerivs_final.stop();
  TIME_SPHevalDerivs.stop();
}


//===================================================================================
// Allows density sum to only be applied to some nodeLists
//===================================================================================
template<typename Dimension>
void
FSISolidSPHHydroBase<Dimension>::
computeFSISPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // First the self contribution.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = massDensity[nodeListi]->numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      if (mSumDensityNodeListSwitch[nodeListi]==1){
        const auto  mi = mass(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        massDensity(nodeListi, i) = mi*Hdeti*W0;
      }
    }
  }

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto massDensity_thread = massDensity.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      if(mSumDensityNodeListSwitch[nodeListi]==1 || mSumDensityNodeListSwitch[nodeListj]==1){
        // State for node i
        const auto& ri = position(nodeListi, i);
        const auto  mi = mass(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        //const auto rhoi = massDensity(nodeListi, i)
      
        // State for node j
        const auto& rj = position(nodeListj, j);
        const auto  mj = mass(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  Hdetj = Hj.Determinant();
      //const auto rhoj = massDensity(nodeListj, j)
      
        // Kernel weighting and gradient.
        const auto rij = ri - rj;
        const auto etai = (Hi*rij).magnitude();
        const auto etaj = (Hj*rij).magnitude();
        const auto Wi = W.kernelValue(etai, Hdeti);
        const auto Wj = W.kernelValue(etaj, Hdetj);

        // Sum the pair-wise contributions.
        //massDensity_thread(nodeListi, i) += (nodeListi == nodeListj ? mj : mi)*Wj ;
        //massDensity_thread(nodeListj, j) += (nodeListi == nodeListj ? mi : mj)*Wi ;
    
        massDensity_thread(nodeListi, i) += (mSumDensityNodeListSwitch[nodeListi]==1 ? (nodeListi == nodeListj ? mj : mi)*Wi : 0.0);
        massDensity_thread(nodeListj, j) += (mSumDensityNodeListSwitch[nodeListj]==1 ? (nodeListi == nodeListj ? mi : mj)*Wj : 0.0);
      }
    }

#pragma omp critical
    {
      massDensity_thread.threadReduce();
    }
  }
}

}
