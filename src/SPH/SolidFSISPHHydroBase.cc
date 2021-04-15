//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBase -- SolidSPHHydro modified for large density discontinuities
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Utilities/DamagedNodeCoupling.hh"
#include "SPH/SPHHydroBase.hh"
#include "SPH/SolidSPHHydroBase.hh"
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
#include "Utilities/Timer.hh"

#include "SolidFSISPHHydroBase.hh"

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

extern Timer TIME_SolidFSISPHregisterDerivs;

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
                  const TableKernel<Dimension>& W,
                  const double filter,
                  const double cfl,
                  const double surfaceForceCoefficient,
                  const double densityStabilizationCoefficient,
                  const double densityDiffusionCoefficient,
                  const double specificThermalEnergyDiffusionCoefficient,
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
                  const bool negativePressureInDamage,
                  const bool strengthInDamage,
                  const Vector& xmin,
                  const Vector& xmax):
  SolidSPHHydroBase<Dimension>(smoothingScaleMethod,
                               dataBase,
                               Q,
                               W,
                               W,  //WPi
                               W,  //WGrad
                               filter,
                               cfl,
                               useVelocityMagnitudeForDt,
                               compatibleEnergyEvolution,
                               evolveTotalEnergy,
                               gradhCorrection, 
                               XSPH,
                               correctVelocityGradient,
                               true, // sumMassDensityOverAllNodeLists
                               densityUpdate,
                               HUpdate,
                               epsTensile,
                               nTensile,
                               damageRelieveRubble,
                               negativePressureInDamage,
                               strengthInDamage,
                               xmin,
                               xmax),
  mSurfaceForceCoefficient(surfaceForceCoefficient),
  mDensityStabilizationCoefficient(densityStabilizationCoefficient),
  mDensityDiffusionCoefficient(densityDiffusionCoefficient),
  mSpecificThermalEnergyDiffusionCoefficient(specificThermalEnergyDiffusionCoefficient),
  mSumDensityNodeLists(sumDensityNodeLists),
  mTensorDepsDt(FieldStorageType::CopyFields){
     mTensorDepsDt = dataBase.newFluidFieldList(Tensor::zero,  IncrementFieldList<Dimension, Scalar>::prefix() + "Tensor" + HydroFieldNames::specificThermalEnergy);
  }

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SolidFSISPHHydroBase<Dimension>::
~SolidFSISPHHydroBase() {
}

//------------------------------------------------------------------------------
// Register Derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_SolidFSISPHregisterDerivs.start();

  // Call the ancestor method.
  SolidSPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);

  const auto tensorDepsDtName = IncrementFieldList<Dimension, Scalar>::prefix() + "Tensor" + HydroFieldNames::specificThermalEnergy;
  dataBase.resizeFluidFieldList(mTensorDepsDt, Tensor::zero, tensorDepsDtName , false);
  derivs.enroll(mTensorDepsDt);

  TIME_SolidFSISPHregisterDerivs.stop();
}
//------------------------------------------------------------------------------
// FSI specialized density summmation
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  SolidSPHHydroBase<Dimension>::preStepInitialize(dataBase,state,derivs);

  // test to see if any of the nodeLists require their own density sum
  bool applySelectDensitySum = false;
  auto numNodeLists = dataBase.numNodeLists();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    if (mSumDensityNodeLists[nodeListi]==1){
      applySelectDensitySum = true;
    }
  }
  //std::cout << applySelectDensitySum;
  if (applySelectDensitySum){
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto& mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto& K = state.fields(SolidFieldNames::bulkModulus, 0.0);
      const auto& p = state.fields(HydroFieldNames::pressure, 0.0);
      const auto& H = state.fields(HydroFieldNames::H, SymTensor::zero);
      const auto& W = this->kernel();
      auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      this->computeFSISPHSumMassDensity(connectivityMap, W, position, mass, K, p, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }

}


//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {


  // Get the ArtificialViscosity.
  auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  const auto tiny = 1.0e-30;
  //const auto W0 = W(0.0, 1.0);
  const auto epsTensile = this->epsilonTensile();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto damageRelieveRubble = this->damageRelieveRubble();
  const auto rhoDiffusionCoeff = this->densityDiffusionCoefficient();
  const auto epsDiffusionCoeff = this->specificThermalEnergyDiffusionCoefficient();
  const auto rhoStabilizeCoeff = this->densityStabilizationCoefficient();
  const auto surfaceForceCoeff = this->surfaceForceCoefficient();
  //const auto alpha = this->alpha();
  //const auto diffCoeff = this->diffusionCoefficient();
  //const auto HUpdatePolicy = this-> HEvolution();
  //const auto strengthInDamage = this->strengthInDamage();
  //const auto negativePressureInDamage = this->negativePressureInDamage();
  const auto XSPH = this->XSPH();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  //auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  const auto S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  const auto mu = state.fields(SolidFieldNames::shearModulus, 0.0);
  const auto damage = state.fields(SolidFieldNames::tensorDamage, SymTensor::zero);
  //const auto gradDamage = state.fields(SolidFieldNames::damageGradient, Vector::zero);
  const auto fragIDs = state.fields(SolidFieldNames::fragmentIDs, int(1));
  const auto pTypes = state.fields(SolidFieldNames::particleTypes, int(0));
  //const auto K = state.fields(SolidFieldNames::bulkModulus, 0.0);
    
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  //CHECK(omega.size() == numNodeLists);
  CHECK(S.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  //CHECK(gradDamage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);
  //CHECK(K.size() == numNodeLists);

  // Derivative FieldLists.
  //auto  rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  tensorDepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() +"Tensor"+HydroFieldNames::specificThermalEnergy, Tensor::zero);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  //auto  maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  //auto  effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  //auto  rhoSumCorrection = derivatives.fields(HydroFieldNames::massDensityCorrection, 0.0);
  //auto  viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  //auto  XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  //auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto  DSDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStress, SymTensor::zero);
  //CHECK(rhoSum.size() == numNodeLists);
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
  CHECK(tensorDepsDt.size() == numNodeLists);
  //CHECK(maxViscousPressure.size() == numNodeLists);
  //CHECK(effViscousPressure.size() == numNodeLists);
  //CHECK(rhoSumCorrection.size() == numNodeLists);
  //CHECK(viscousWork.size() == numNodeLists);
  //CHECK(XSPHWeightSum.size() == numNodeLists);
  //CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DSDt.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy) pairAccelerations.resize(npairs);

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);

  //Walk all the interacting pairs and calculated linear correction tensor
  //based on taylor series expansion of gradient estimate
//   if (this->mCorrectVelocityGradient){
        
// #pragma omp parallel
//   {
//     // Thread private  scratch variables.
//     int i, j, nodeListi, nodeListj;

//     typename SpheralThreads<Dimension>::FieldListStack threadStack;
//     //auto DvDx_thread = DvDx.threadCopy(threadStack);
//     auto M_thread = M.threadCopy(threadStack);
//     //auto localM_thread = localM.threadCopy(threadStack);

// #pragma omp for
//     for (auto kk = 0u; kk < npairs; ++kk) {
//       const auto start = Timing::currentTime();
//       i = pairs[kk].i_node;
//       j = pairs[kk].j_node;
//       nodeListi = pairs[kk].i_list;
//       nodeListj = pairs[kk].j_list;

//       // Get the state for node i.
//       const auto& ri = position(nodeListi, i);
//       const auto  mi = mass(nodeListi, i);
//       const auto  rhoi = massDensity(nodeListi, i);
//       const auto& Hi = H(nodeListi, i);
//       const auto  Hdeti = Hi.Determinant();
//       CHECK(mi > 0.0);
//       CHECK(rhoi > 0.0);
//       CHECK(Hdeti > 0.0);

//       // Get the state for node j
//       const auto& rj = position(nodeListj, j);
//       const auto mj = mass(nodeListj, j);
//       const auto rhoj = massDensity(nodeListj, j);
//       const auto& Hj = H(nodeListj, j);
//       const auto Hdetj = Hj.Determinant();
//       CHECK(mj > 0.0);
//       CHECK(rhoj > 0.0);
//       CHECK(Hdetj > 0.0);

//       auto& Mi = M_thread(nodeListi, i);
//       auto& Mj = M_thread(nodeListj, j);

//       const auto rij = ri - rj;

//       const auto etai = Hi*rij;
//       const auto etaj = Hj*rij;
//       const auto etaMagi = etai.magnitude();
//       const auto etaMagj = etaj.magnitude();
//       CHECK(etaMagi >= 0.0);
//       CHECK(etaMagj >= 0.0);

//       const auto Hetai = Hi*etai.unitVector();
//       const auto gradWi = W.gradValue(etaMagi, Hdeti)*Hetai;

//       const auto Hetaj = Hj*etaj.unitVector();
//       const auto gradWj = W.gradValue(etaMagj, Hdetj)*Hetaj;

//       // Linear gradient correction term.
//       Mi -= mj/rhoj*rij.dyad(gradWi);
//       Mj -= mi/rhoi*rij.dyad(gradWj);

//       // Add timing info for work
//       const auto deltaTimePair = 0.5*Timing::difference(start, Timing::currentTime());
// #pragma omp atomic
//       nodeLists[nodeListi]->work()(i) += deltaTimePair;
// #pragma omp atomic
//       nodeLists[nodeListj]->work()(j) += deltaTimePair;

//     } // loop over pairs

//     // Reduce the thread values to the master.
//     threadReduceFieldLists<Dimension>(threadStack);

//   }   // OpenMP parallel region

 
//  // Loop through nodes and invert the linear correction matrix
//   for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
//     const auto& nodeList = mass[nodeListi]->nodeList();
//     const auto ni = nodeList.numInternalNodes();
// #pragma omp parallel for
//     for (auto i = 0u; i < ni; ++i) {
//       const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
//       auto& Mi = M(nodeListi, i);

//       if (std::abs(Mi.Determinant()) > 1.0e-10 and
//           numNeighborsi > Dimension::pownu(2)) {
//         Mi = Mi.Inverse();
//       } else { 
//         Mi = Tensor::one; 
//       } 

//     }
//   }
  
  
//    } // if statement for vel grad correction


// Now we calculate  the hydro deriviatives
// Walk all the interacting pairs.
#pragma omp parallel
  {
    // Thread private  scratch variables.
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, Wj, gWj, fDeffi, fDeffj;
    Tensor QPiij, QPiji;
    SymTensor sigmai, sigmaj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto tensorDepsDt_thread = tensorDepsDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);
    auto DSDt_thread = DSDt.threadCopy(threadStack);
    auto DrhoDt_thread = DrhoDt.threadCopy(threadStack);

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
      const auto& epsi = specificThermalEnergy(nodeListi,i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto  Di = damage(nodeListi,i).eigenValues().maxElement();
      const auto  Hdeti = Hi.Determinant();
      const auto& pTypei = pTypes(nodeListi, i);
      //const auto fragIDi = fragIDs(nodeListi, i);
      const auto& mui = mu(nodeListi,i);
      const auto Ki = max(tiny,rhoi*ci*ci)+4/3*mui;
      CHECK(mi > 0.0);
      CHECK(ci > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& Mi = M_thread(nodeListi, i);
      auto& localMi = localM_thread(nodeListi, i);
      auto& DrhoDti = DrhoDt_thread(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& tensorDepsDti = tensorDepsDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
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
      const auto& epsj = specificThermalEnergy(nodeListj,j);
      const auto& Sj = S(nodeListj, j);
      const auto  Dj = damage(nodeListj,j).eigenValues().maxElement();
      const auto  Hdetj = Hj.Determinant();
      //const auto fragIDj = fragIDs(nodeListj, j);
      const auto& pTypej = pTypes(nodeListj, j);
      const auto& muj = mu(nodeListj,j);
      const auto Kj = max(tiny,rhoj*cj*cj)+4/3*muj; // P - wave modulus
      CHECK(mj > 0.0);
      CHECK(cj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      //auto& dwdhj = dwdh_thread(nodeListj,j);
      auto& Mj = M_thread(nodeListj, j);
      auto& localMj = localM_thread(nodeListj, j);
      auto& DrhoDtj = DrhoDt_thread(nodeListj, j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& tensorDepsDtj = tensorDepsDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const auto sameMatij =  (nodeListi == nodeListj);// and fragIDi == fragIDj);  //11/05/2020

      // Flag if at least one particle is free (0).
      const auto freeParticle = (pTypei == 0 or pTypej == 0);

      // Node displacement.
      const auto rij = ri - rj;
      const auto Hij = 0.5*(Hi+Hj);
      const auto etaij = Hij*rij;
      const auto Hdetij = Hij.Determinant();
      const auto etaMagij = etaij.magnitude();
      CHECK(etaMagij >= 0.0);

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
      
      // Determine how we're applying damage.
      const auto fDij = pairs[kk].f_couple;

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      const auto rij2 = rij.magnitude2();
      const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      weightedNeighborSumi += abs(gWi);
      weightedNeighborSumj += abs(gWj);
      massSecondMomenti += gradWi.magnitude2()*thpt;
      massSecondMomentj += gradWj.magnitude2()*thpt;

      // averaged things.
      const auto rhoij = 0.5*(rhoi+rhoj); 
      const auto cij = 0.5*(ci+cj);  
      const auto gradWij = 0.5*(gradWi+gradWj);
       
      // volumes
      const auto voli = mi/rhoi;
      const auto volj = mj/rhoj;

      std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
                                      ri, etaij, vi, rhoij, cij, Hij,  
                                      rj, etaij, vj, rhoij, cij, Hij);      
            
      // material interface pressure.
      if (sameMatij) {
        sigmai = fDij*Si - Pi * SymTensor::one;
        sigmaj = fDij*Sj - Pj * SymTensor::one;
      }else {
        const auto PSi = rij.dot(Si.dot(rij))/rij2;
        const auto PSj = rij.dot(Sj.dot(rij))/rij2;
        const auto Pstar = ((Pi)*rhoj+(Pj)*rhoi)/(rhoi+rhoj);
        //const auto Pstar = (Pi+Pj)*0.5;
        sigmai = -max(Pstar,0.0)*SymTensor::one;
        sigmaj = -max(Pstar,0.0)*SymTensor::one;
      }
      //const auto sigmaij = 0.5*(sigmai + sigmaj)
      // Compute the tensile correction to add to the stress as described in 
      // Gray, Monaghan, & Swift (Comput. Methods Appl. Mech. Eng., 190, 2001)
      const auto fi = epsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
      const auto fj = epsTensile*FastMath::pow4(Wj/(Hdetj*WnPerh));
      const auto Ri = fi*tensileStressCorrection(sigmai);
      const auto Rj = fj*tensileStressCorrection(sigmaj);
      sigmai += Ri;
      sigmaj += Rj;

      
      // Force Calc
      //=====================================================
      // (Monaghan 2013) - interface force
      const auto sf = (sameMatij ? 1.0 : 1.0 + surfaceForceCoeff*abs((rhoi-rhoj)/(rhoi+rhoj+tiny)));
      const auto sigmarhoi = sf*sigmai/(rhoi*rhoj)-0.5*QPiij;
      const auto sigmarhoj = sf*sigmaj/(rhoj*rhoi)-0.5*QPiji;
      const auto deltaDvDt = sigmarhoi*gradWi + sigmarhoj*gradWj;
      if (freeParticle) {
        DvDti += mj*deltaDvDt;
        DvDtj -= mi*deltaDvDt;
      } 
      if(compatibleEnergy) pairAccelerations[kk] = mj*deltaDvDt; 
      

     // dwdhi = -(Dimension::nDim*Wi + etaMagi*gWi)*Dimension::rootnu(Hdeti);
     // dwdhj = -(Dimension::nDim*Wj + etaMagj*gWj)*Dimension::rootnu(Hdetj);
      const auto rhatij = rij.unitVector();
      const auto ui = vi.dot(rhatij);
      const auto uj = vj.dot(rhatij);
      const auto umax = max(ui,uj);
      const auto umin = min(ui,uj);
      const auto uij = ui-uj;
      
      //const auto wi = vi-ui*rhatij;
      //const auto wj = vi-ui*rhatij;
      //const auto isExpanding   = uij>0.0;
      //const auto isCompressing = uij<0.0;
      // if (sameMatij){ 
      //   const auto ustarStabilizer =  (rhoj-rhoi) * safeInv(max(tiny,max(rhoi,rhoj)*etaMagij));
      //   const auto deltaUstar = rhoStabilizeCoeff * min(0.1,max(-0.1, ustarStabilizer)) * cij
      //   const auto ustar =  max(umin, min(umax, 0.5*(ui+uj) + deltaUstar));
      //   auto vstari = (ustar - 0.5*(ui+uj))*rhatij + 0.5*(vi+vj);
      //   auto vstarj = (ustar - 0.5*(ui+uj))*rhatij + 0.5*(vi+vj);
      // }else{
      //   // adjust velocity gradient for diff compressibilities
      

      // //auto wstar =  (sameMatij    ?
      // //               0.5*(wi+wj)  :
      // //              (mui*volj*gWi*wi + muj*voli*gWj*wj)/
      // //              (mui*volj*gWi    + muj*voli*gWj));

      // // put our diffusion stabilizer here
      //   const auto denom = safeInv(max(tiny,max(rhoi*ci*ci,rhoj*cj*cj)*etaMagij));

      //   const auto ustarStabilizer =  (Pj-Pi) * safeInv( max(tiny, max(rhoi*ci*ci,rhoj*cj*cj) *etaMagij ) );
      
      //   auto ustar =  (Ki*volj*gWi*ui + Kj*voli*gWj*uj)/
      //                 (Ki*volj*gWi    + Kj*voli*gWj);

      //   ustar += rhoStabilizeCoeff * min(0.1, max(-0.1, ustarStabilizer)) * cij;
      
      //   ustar = max(umin, min(umax, ustar) );
     
      //   vstari =  (true               ?
      //             (ustar - ui)*rhatij + vi :
      //             (ustar - 0.5*(ui+uj))*rhatij + 0.5*(vi+vj));

      //   vstarj = (true                ?
      //                     (ustar - uj)*rhatij + vj  :
      //                     (ustar - 0.5*(ui+uj))*rhatij + 0.5*(vi+vj));
      // }
      // adjust velocity gradient for diff compressibilities

      auto ustar =  (sameMatij    ?
                     0.5*(ui+uj)  :
                    (Ki*volj*gWi*ui + Kj*voli*gWj*uj)/
                    (Ki*volj*gWi    + Kj*voli*gWj));

      // put our diffusion stabilizer here
      const auto denom = (sameMatij                                   ?
                          safeInv(max(tiny,max(rhoi,rhoj)*etaMagij))  :
                          safeInv(max(tiny,max(rhoi*ci*ci,rhoj*cj*cj)*etaMagij)) );

      const auto ustarStabilizer =  (sameMatij  ?
                                    (rhoj-rhoi) :
                                    (Pj-Pi)     ) *denom;

      ustar += rhoStabilizeCoeff * min(0.1, max(-0.1, ustarStabilizer)) * cij;
      ustar =  max(umin, min(umax, ustar) );
     
      const auto vstari =  (!sameMatij               ?
                           (ustar - ui)*rhatij + vi :
                           (ustar - 0.5*(ui+uj))*rhatij + 0.5*(vi+vj));

      const auto vstarj = (!sameMatij                ?
                          (ustar - uj)*rhatij + vj  :
                          (ustar - 0.5*(ui+uj))*rhatij + 0.5*(vi+vj));

      //const auto kappai = (sameMatij ?
      //                     0.5       :
      //                     max(0.0, min(1.0, Kj*voli*gWj/(Ki*volj*gWi + Kj*voli*gWj)))
      //                     );

      //const auto kappaj = max(0.0, min(1.0, 1.0-kappai));
      //const auto vij = vi-vj;
      //const auto deltaEps = mi*mj*deltaDvDt.dot(vij);
      //auto deltaDvDxi = 2.0*kappai*(vij).dyad(gradWi);
      //auto deltaDvDxj = 2.0*kappaj*(vij).dyad(gradWj);

      auto deltaDvDxi = 2.0*(vi-vstari).dyad(gradWi);
      auto deltaDvDxj = 2.0*(vstarj-vj).dyad(gradWj);
      if (sameMatij) {
        localMi -= volj*rij.dyad(gradWi);
        localMj -= voli*rij.dyad(gradWj);
        localDvDxi -= volj*(deltaDvDxi);
        localDvDxj -= voli*(deltaDvDxj); 
      }
      if (this->mCorrectVelocityGradient){
        Mi -=  volj*rij.dyad(gradWi);
        Mj -=  voli*rij.dyad(gradWj);
      }
      DvDxi -= volj*deltaDvDxi;
      DvDxj -= voli*deltaDvDxj;

      // Conservation of Energy
      //--------------------------
      tensorDepsDti -= mj*(deltaDvDxi.Transpose()*sigmarhoi);
      tensorDepsDtj -= mi*(deltaDvDxj.Transpose()*sigmarhoj);

      // Dissipation
      //-------------------------
      if (sameMatij and rhoDiffusionCoeff>tiny){
        const auto diffusion = rhoDiffusionCoeff*(rhoi-rhoj)*cij*etaij.dot(gradWij)/(etaMagij*etaMagij+tiny);
        DrhoDti += volj*diffusion;
        DrhoDtj -= voli*diffusion;
      }

      if (sameMatij and epsDiffusionCoeff>tiny){
        const auto diffusion = epsDiffusionCoeff*(epsi-epsj)*cij/rhoij*etaij.dot(gradWij)/(etaMagij*etaMagij+tiny);
        DepsDti += mj*diffusion;
        DepsDtj -= mi*diffusion;
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



  // Finish up the derivatives for each point.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    // Check if we can identify a reference density.
    auto rho0 = 0.0;
    try {
      rho0 = dynamic_cast<const SolidEquationOfState<Dimension>&>(dynamic_cast<const FluidNodeList<Dimension>&>(nodeList).equationOfState()).referenceDensity();
      // cerr << "Setting reference density to " << rho0 << endl;
    } catch(...) {
      // cerr << "BLAGO!" << endl;
    }

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& Si = S(nodeListi, i);
      const auto& mui = mu(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& tensorDepsDti = tensorDepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& localMi = localM(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      auto& DSDti = DSDt(nodeListi, i);
      
      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;
      
      if (this->mCorrectVelocityGradient and 
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi=DvDxi*Mi;
        DepsDti += tensorDepsDti.doubledot(Mi);
      } else { 
        DepsDti += tensorDepsDti.Trace(); 
      } 

      DrhoDti -=  rhoi*DvDxi.Trace();

      DxDti = vi;

      if (this->mEvolveTotalEnergy) DepsDti = mi*(DxDti.dot(DvDti) + DepsDti);

      
      
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

      const auto Di = (damageRelieveRubble ? 
                       max(0.0, min(1.0, damage(nodeListi, i).Trace() - 1.0)) :
                       0.0);

      if (std::abs(localMi.Determinant()) > 1.0e-10 and
        numNeighborsi > Dimension::pownu(2)) {
        localMi = localMi.Inverse();
        localDvDxi = localDvDxi*localMi;
      }

      // Determine the deviatoric stress evolution.
      const auto deformation = localDvDxi.Symmetric();
      const auto spin = localDvDxi.SkewSymmetric();
      const auto deviatoricDeformation = deformation - (deformation.Trace()/3.0)*SymTensor::one;
      const auto spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;

      // In the presence of damage, add a term to reduce the stress on this point.
      DSDti = (1.0 - Di)*DSDti - 0.25/dt*Di*Si;
    }
  }
}


//===================================================================================
// Allows density sum to only be applied to some nodeLists
//===================================================================================
template<typename Dimension>
void
SolidFSISPHHydroBase<Dimension>::
computeFSISPHSumMassDensity(const ConnectivityMap<Dimension>& connectivityMap,
                            const TableKernel<Dimension>& W,
                            const FieldList<Dimension, typename Dimension::Vector>& position,
                            const FieldList<Dimension, typename Dimension::Scalar>& mass,
                            const FieldList<Dimension, typename Dimension::Scalar>& pressure,
                            const FieldList<Dimension, typename Dimension::Scalar>& bulkModulus,
                            const FieldList<Dimension, typename Dimension::SymTensor>& H,
                            FieldList<Dimension, typename Dimension::Scalar>& massDensity) {

  // Pre-conditions.
  const auto numNodeLists = massDensity.size();
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(mass.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  // Some useful variables.
  const auto W0 = W.kernelValue(0.0, 1.0);

  typedef typename Dimension::Scalar Scalar;
  FieldList<Dimension, Scalar> storeDensity(FieldStorageType::CopyFields);
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    storeDensity.appendNewField("storeDensity", massDensity[nodeListi]->nodeList(), 0.0);
  }


  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Now the pair contributions.
#pragma omp parallel
  {
    int i, j, nodeListi, nodeListj;
    auto storeMassDensity_thread = storeDensity.threadCopy();
    //auto massDensity_thread = massDensity.threadCopy();

#pragma omp for
    for (auto k = 0u; k < npairs; ++k) {
      i = pairs[k].i_node;
      j = pairs[k].j_node;
      nodeListi = pairs[k].i_list;
      nodeListj = pairs[k].j_list;

      if(mSumDensityNodeLists[nodeListi]==1 || mSumDensityNodeLists[nodeListj]==1){
        // State for node i
        const auto& ri = position(nodeListi, i);
        const auto  mi = mass(nodeListi, i);
        //const auto  Pi = pressure(nodeListi, i);
        //const auto  Ki = bulkModulus(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        //const auto rhoi = massDensity(nodeListi, i);
      
        // State for node j
        const auto& rj = position(nodeListj, j);
        const auto  mj = mass(nodeListj, j);
        //const auto  Pj = pressure(nodeListj, j);
        //const auto  Kj = bulkModulus(nodeListj, j);
        const auto& Hj = H(nodeListj, j);
        const auto  Hdetj = Hj.Determinant();
        //const auto rhoj = massDensity(nodeListj, j);
      
        // Kernel weighting and gradient.
        const auto rij = ri - rj;
        const auto etai = (Hi*rij).magnitude();
        const auto etaj = (Hj*rij).magnitude();
        const auto Wi = W.kernelValue(etai, Hdeti);
        const auto Wj = W.kernelValue(etaj, Hdetj);

        storeMassDensity_thread(nodeListi, i) += (mSumDensityNodeLists[nodeListi]==1 ? 
                                            (nodeListi == nodeListj ? mj : mi)*Wi : 
                                             0.0);
        storeMassDensity_thread(nodeListj, j) += (mSumDensityNodeLists[nodeListj]==1 ? 
                                            (nodeListi == nodeListj ? mi : mj)*Wj : 
                                             0.0);
      }
    }

#pragma omp critical
    {
      storeMassDensity_thread.threadReduce();
    }
  }

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = massDensity[nodeListi]->numInternalElements();
    if (mSumDensityNodeLists[nodeListi]==1){
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        const auto  mi = mass(nodeListi, i);
        const auto& Hi = H(nodeListi, i);
        const auto  Hdeti = Hi.Determinant();
        massDensity(nodeListi,i) = storeDensity(nodeListi, i) + mi*Hdeti*W0;
      }
    }
  }

}



}
