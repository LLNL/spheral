//---------------------------------Spheral++----------------------------------//
// SolidCRKSPHHydroBase -- The CRKSPH/ACRKSPH solid material hydrodynamic
// package for Spheral++.
//
// This is the area-weighted RZ specialization.
//
// Created by JMO, Fri May 13 10:50:36 PDT 2016
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
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicyRZ.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "Strength/RZDeviatoricStressPolicy.hh"
#include "Strength/RZPlasticStrainPolicy.hh"
#include "Strength/BulkModulusPolicy.hh"
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

#include "SolidCRKSPHHydroBaseRZ.hh"

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
SolidCRKSPHHydroBaseRZ::
SolidCRKSPHHydroBaseRZ(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
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
  SolidCRKSPHHydroBase<Dimension>(smoothingScaleMethod, 
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
  mDeviatoricStressTT(FieldSpace::FieldStorageType::CopyFields),
  mDdeviatoricStressTTDt(FieldSpace::FieldStorageType::CopyFields) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SolidCRKSPHHydroBaseRZ::
~SolidCRKSPHHydroBaseRZ() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
initializeProblemStartup(DataBase<Dim<2> >& dataBase) {

  // Call the ancestor.
  SolidCRKSPHHydroBase<Dimension>::initializeProblemStartup(dataBase);

  dataBase.isRZ = true;

  // Create storage for the state we're holding.
  mDeviatoricStressTT = dataBase.newSolidFieldList(0.0, SolidFieldNames::deviatoricStressTT);
  mDdeviatoricStressTTDt = dataBase.newSolidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + SolidFieldNames::deviatoricStressTT);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
registerState(DataBase<Dim<2> >& dataBase,
              State<Dim<2> >& state) {

  typedef State<Dimension>::PolicyPointer PolicyPointer;

  // Call the ancestor.
  SolidCRKSPHHydroBase<Dimension>::registerState(dataBase, state);

  // Create the local storage.
  dataBase.resizeFluidFieldList(mDeviatoricStressTT, 0.0, SolidFieldNames::deviatoricStressTT, false);

  // Reregister the deviatoric stress and plastic strain policies to the RZ specialized versions
  // that account for the theta-theta component of the stress.
  FieldList<Dimension, SymTensor> S = state.fields(SolidFieldNames::deviatoricStress, SymTensor::zero);
  FieldList<Dimension, Scalar> ps = state.fields(SolidFieldNames::plasticStrain, 0.0);
  PolicyPointer deviatoricStressPolicy(new RZDeviatoricStressPolicy());
  PolicyPointer plasticStrainPolicy(new RZPlasticStrainPolicy());
  state.enroll(S, deviatoricStressPolicy);
  state.enroll(ps, plasticStrainPolicy);
  state.enroll(mDeviatoricStressTT);

  // Are we using the compatible energy evolution scheme?
  // If so we need to override the ordinary energy registration with a specialized version.
  if (mCompatibleEnergyEvolution) {
    FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
    PolicyPointer thermalEnergyPolicy(new NonSymmetricSpecificThermalEnergyPolicyRZ(dataBase));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);

    // Get the policy for the position, and add the specific energy as a dependency.
    PolicyPointer positionPolicy = state.policy(state.buildFieldKey(HydroFieldNames::position, UpdatePolicyBase<Dimension>::wildcard()));
    positionPolicy->addDependency(HydroFieldNames::specificThermalEnergy);
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
registerDerivatives(DataBase<Dim<2> >& dataBase,
                    StateDerivatives<Dim<2> >& derivs) {

  // Call the ancestor method.
  SolidCRKSPHHydroBase<Dimension>::registerDerivatives(dataBase, derivs);

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  const string DSTTDtName = IncrementFieldList<Dimension, Vector>::prefix() + SolidFieldNames::deviatoricStressTT;
  dataBase.resizeFluidFieldList(mDdeviatoricStressTTDt, 0.0, DSTTDtName, false);

  derivs.enroll(mDdeviatoricStressTTDt);
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
evaluateDerivatives(const Dim<2>::Scalar time,
                    const Dim<2>::Scalar dt,
                    const DataBase<Dim<2> >& dataBase,
                    const State<Dim<2> >& state,
                    StateDerivatives<Dim<2> >& derivatives) const {

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  typedef Timing::Time Time;
  const double tiny = 1.0e-30;
  const bool compatibleEnergy = this->compatibleEnergyEvolution();
  const bool XSPH = this->XSPH();
  const Scalar epsTensile = this->epsilonTensile();
  const CRKOrder order = this->correctionOrder();
  const double correctionMin = this->correctionMin();
  const double correctionMax = this->correctionMax();

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
  const FieldList<Dimension, Scalar> STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
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
  CHECK(STT.size() == numNodeLists);
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
  CHECK(Adamage.size() == numNodeLists);
  CHECK(Bdamage.size() == numNodeLists);
  CHECK(Cdamage.size() == numNodeLists or order != CRKOrder::QuadraticOrder);
  CHECK(gradAdamage.size() == numNodeLists);
  CHECK(gradBdamage.size() == numNodeLists);
  CHECK(gradCdamage.size() == numNodeLists or order != CRKOrder::QuadraticOrder);

  const FieldList<Dimension, SymTensor>& Hfield0 = this->Hfield0();

  // Derivative FieldLists.
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
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
  FieldList<Dimension, Scalar> DSTTDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStressTT, 0.0);
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
  CHECK(DSTTDt.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy) {
    size_t nodeListi = 0;
    for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
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
  for (DataBase<Dimension>::ConstSolidNodeListIterator itr = dataBase.solidNodeListBegin();
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
    for (ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& posi = position(nodeListi, i);
      const Scalar ri = abs(posi.y());
      const Scalar circi = 2.0*M_PI*ri;
      const Scalar mi = mass(nodeListi, i);
      const Scalar mRZi = mi/circi;
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Scalar& epsi = specificThermalEnergy(nodeListi, i);
      const Scalar& Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar& ci = soundSpeed(nodeListi, i);
      const SymTensor& Si = S(nodeListi, i);
      const Scalar STTi = STT(nodeListi, i);
      const Scalar& mui = mu(nodeListi, i);
      const Scalar Ai = A(nodeListi, i);
      const Vector& gradAi = gradA(nodeListi, i);
      const Scalar Adami = Adamage(nodeListi, i);
      const Vector& gradAdami = gradAdamage(nodeListi, i);
      if (order != CRKOrder::ZerothOrder) {
        Bi = B(nodeListi, i);
        gradBi = gradB(nodeListi, i);
        Bdami = Bdamage(nodeListi, i);
        gradBdami = gradBdamage(nodeListi, i);
      }
      if (order == CRKOrder::QuadraticOrder) {
        Ci = C(nodeListi, i);
        gradCi = gradC(nodeListi, i);
        Cdami = Cdamage(nodeListi, i);
        gradCdami = gradCdamage(nodeListi, i);
      }
      const Scalar Hdeti = Hi.Determinant();
      const Scalar weighti = volume(nodeListi, i);  // Change CRKSPH weights here if need be!
      const Scalar zetai = abs((Hi*posi).y());
      const Scalar hri = ri*safeInv(zetai);
      const Scalar riInv = safeInv(ri, 0.25*hri);
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
      Scalar& DSTTDti = DSTTDt(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
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
              const Vector& posj = position(nodeListj, j);
              const Scalar rj = abs(posj.y());
              const Scalar circj = 2.0*M_PI*rj;
              const Scalar mj = mass(nodeListj, j);
              const Scalar mRZj = mj/circj;
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
              if (order != CRKOrder::ZerothOrder) {
                Bj = B(nodeListj, j);
                gradBj = gradB(nodeListj, j);
                Bdamj = Bdamage(nodeListj, j);
                gradBdamj = gradBdamage(nodeListj, j);
              }
              if (order == CRKOrder::QuadraticOrder) {
                Cj = C(nodeListj, j);
                gradCj = gradC(nodeListj, j);
                Cdamj = Cdamage(nodeListj, j);
                gradCdamj = gradCdamage(nodeListj, j);
              }
              const SymTensor& Sj = S(nodeListj, j);
              const Scalar STTj = STT(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar weightj = volume(nodeListj, j);     // Change CRKSPH weights here if need be!
              const Scalar zetaj = abs((Hj*posj).y());
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
              const Vector xij = posi - posj;
              const Vector etai = Hi*xij;
              const Vector etaj = Hj*xij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);
              const Vector vij = vi - vj;

              // Symmetrized kernel weight and gradient.
              Scalar gWi, gWj, Wi, Wj, gWdami, gWdamj, Wdami, Wdamj;
              Vector gradWi, gradWj, gradWdami, gradWdamj;
              CRKSPHKernelAndGradient(Wj, gWj, gradWj, W, CRKSPHHydroBase<Dimension>::correctionOrder(),  xij,  etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, Ci, gradAi, gradBi, gradCi, correctionMin, correctionMax);
              CRKSPHKernelAndGradient(Wi, gWi, gradWi, W, CRKSPHHydroBase<Dimension>::correctionOrder(), -xij, -etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, Cj, gradAj, gradBj, gradCj, correctionMin, correctionMax);
              CRKSPHKernelAndGradient(Wdamj, gWdamj, gradWdamj, W, CRKSPHHydroBase<Dimension>::correctionOrder(),  xij,  etai, Hi, Hdeti,  etaj, Hj, Hdetj, Adami, Bdami, Cdami, gradAdami, gradBdami, gradCdami, correctionMin, correctionMax);
              CRKSPHKernelAndGradient(Wdami, gWdami, gradWdami, W, CRKSPHHydroBase<Dimension>::correctionOrder(), -xij, -etaj, Hj, Hdetj, -etai, Hi, Hdeti, Adamj, Bdamj, Cdamj, gradAdamj, gradBdamj, gradCdamj, correctionMin, correctionMax);
              const Vector deltagrad = gradWj - gradWi;
              const Vector deltagraddam = gradWdamj - gradWdami;
              const Vector gradWSPHi = (Hi*etai.unitVector())*W.gradValue(etai.magnitude(), Hdeti);
              const Vector gradWSPHj = (Hj*etaj.unitVector())*W.gradValue(etaj.magnitude(), Hdetj);

              // Find the damaged pair weighting scaling.
              const double fij = coupling(nodeListi, i, nodeListj, j);
              CHECK(fij >= 0.0 and fij <= 1.0);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double fweightij = nodeListi == nodeListj ? 1.0 : mRZj*rhoi/(mRZi*rhoj);
              const double xij2 = xij.magnitude2();
              const auto thpt = xij.selfdyad()*safeInvVar(xij2*xij2*xij2);
              weightedNeighborSumi +=     fweightij*std::abs(gWi);
              weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
              massSecondMomenti +=     fweightij*gradWSPHi.magnitude2()*thpt;
              massSecondMomentj += 1.0/fweightij*gradWSPHj.magnitude2()*thpt;

              // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        posi, etai, vi, rhoi, ci, Hi,
                                                        posj, etaj, vj, rhoj, cj, Hj);
              const Vector Qaccij = (rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second).dot(deltagrad);
              // const Scalar workQij = 0.5*(vij.dot(Qaccij));
              const Scalar workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);                // CRK
              const Scalar workQj = rhoi*rhoi*QPiij.first .dot(vij).dot(deltagrad);                // CRK
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, 4.0*Qi);                                 // We need tighter timestep controls on the Q with CRK
              maxViscousPressurej = max(maxViscousPressurej, 4.0*Qj);
              effViscousPressurei += weightj * Qi * Wj;
              effViscousPressurej += weighti * Qj * Wi;
              viscousWorki += 0.5*weighti*weightj/mRZi*workQi;
              viscousWorkj += 0.5*weighti*weightj/mRZj*workQj;

              // Velocity gradient.
              DvDxi -= weightj*vij.dyad(gradWj);
              DvDxj += weighti*vij.dyad(gradWi);
              localDvDxi -= fij*weightj*vij.dyad(gradWdamj);
              localDvDxj += fij*weighti*vij.dyad(gradWdami);

              // We treat positive and negative pressures distinctly, so split 'em up.
              const Scalar Pposi = max(0.0, Pi),
                           Pnegi = min(0.0, Pi),
                           Pposj = max(0.0, Pj),
                           Pnegj = min(0.0, Pj);

              // Compute the stress tensors.
              SymTensor sigmai, sigmaj;
              if (nodeListi == nodeListj) {
                sigmai = -Pnegi*SymTensor::one + Si;
                sigmaj = -Pnegj*SymTensor::one + Sj;
              }

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
              const Vector forceij  = 0.5*weighti*weightj*((Pposi + Pposj)*deltagrad - fij*(sigmai + sigmaj)*deltagraddam + Qaccij);
              DvDti -= forceij/mRZi;
              DvDtj += forceij/mRZj;
              if (compatibleEnergy) {
                pairAccelerationsi.push_back(-forceij/mRZi);
                pairAccelerationsj.push_back( forceij/mRZj);
              }

              // Specific thermal energy evolution.
              DepsDti += 0.5*weighti*weightj*(Pposj*vij.dot(deltagrad) + fij*sigmaj.dot(vij).dot(deltagraddam) + workQi)/mRZi;
              DepsDtj += 0.5*weighti*weightj*(Pposi*vij.dot(deltagrad) + fij*sigmai.dot(vij).dot(deltagraddam) + workQj)/mRZj;

              // Estimate of delta v (for XSPH).
              XSPHDeltaVi -= fij*weightj*Wj*vij;
              XSPHDeltaVj += fij*weighti*Wi*vij;
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
                                                            posi,
                                                            DvDxi,
                                                            hmin,
                                                            hmax,
                                                            hminratio,
                                                            nPerh);
      Hideali = smoothingScaleMethod.newSmoothingScale(Hi,
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

      // If this node is damaged we begin to force it back to it's original H.
      const Scalar Di = max(0.0, min(1.0, damage(nodeListi, i).eigenValues().maxElement()));
      Hideali = (1.0 - Di)*Hideali + Di*Hfield0(nodeListi, i);

      // Finish the acceleration.
      const Vector deltaDvDti(Si(1,0)/rhoi*riInv,
                              (Si(1,1) - STTi)/rhoi*riInv);
      DvDti += deltaDvDti;
      pairAccelerationsi.push_back(deltaDvDti);

      // Determine the deviatoric stress evolution.
      const SymTensor deformation = localDvDxi.Symmetric();
      const Scalar deformationTT = vi.y()*riInv;
      const Tensor spin = localDvDxi.SkewSymmetric();
      const SymTensor deviatoricDeformation = deformation - ((deformation.Trace() + deformationTT)/3.0)*SymTensor::one;
      const SymTensor spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;
      DSTTDti = 2.0*mui*(deformationTT - (deformation.Trace() + deformationTT)/3.0);

      // In the presence of damage, add a term to reduce the stress on this point.
      // const Scalar Di = max(0.0, min(1.0, damage(nodeListi, i).eigenValues().maxElement()));
      DSDti = (1.0 - Di)*DSDti - 0.25/dt*Di*Si;
      DSTTDti = (1.0 - Di)*DSTTDti - 0.25/dt*Di*STTi;

      // Time evolution of the mass density.
      const Scalar vri = vi.y(); // + XSPHDeltaVi.y();
      DrhoDti = -rhoi*(localDvDxi.Trace() + vri*riInv);

      // We also adjust the density evolution in the presence of damage.
      if (rho0 > 0.0) DrhoDti = (1.0 - Di)*DrhoDti - 0.25/dt*Di*(rhoi - rho0);

      // Finish the specific thermal energy evolution.
      DepsDti += (STTi - Pi)/rhoi*vri*riInv;

      // If needed finish the total energy derivative.
      if (this->evolveTotalEnergy()) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

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
// Finalize the hydro.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
finalize(const Dim<2>::Scalar time,
         const Dim<2>::Scalar dt,
         DataBase<Dim<2> >& dataBase,
         State<Dim<2> >& state,
         StateDerivatives<Dim<2> >& derivs) {

  // If we're going to do the summation density, we need to convert the mass
  // to mass per unit length first.
  if (densityUpdate() == PhysicsSpace::MassDensityType::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::MassDensityType::CorrectedSumDensity) {
    FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> pos = state.fields(HydroFieldNames::position, Vector::zero);
    const unsigned numNodeLists = mass.numFields();
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = mass[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const Scalar circi = 2.0*M_PI*abs(pos(nodeListi, i).y());
        mass(nodeListi, i) /= circi;
      }
    }
  }

  // Base class finalization does most of the work.
  SolidCRKSPHHydroBase<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Now convert back to true masses and mass densities.  We also apply the RZ
  // correction factor to the mass density.
  if (densityUpdate() == PhysicsSpace::MassDensityType::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::MassDensityType::CorrectedSumDensity) {
    const TableKernel<Dimension>& W = this->kernel();
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const unsigned numNodeLists = massDensity.numFields();
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = massDensity[nodeListi]->numElements();
      for (unsigned i = 0; i != n; ++i) {
        const Vector& xi = position(nodeListi, i);
        const Scalar circi = 2.0*M_PI*abs(xi.y());
        mass(nodeListi, i) *= circi;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
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

  // Apply ordinary BCs.
  SolidCRKSPHHydroBase<Dim<2> >::applyGhostBoundaries(state, derivs);
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

  // Apply boundary conditions to our extra strength variables.
  FieldList<Dimension, Scalar> STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(STT);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
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

  // Apply ordinary BCs.
  SolidCRKSPHHydroBase<Dim<2> >::enforceBoundaries(state, derivs);

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

  // Enforce boundary conditions on the extra strength variables.
  FieldList<Dimension, Scalar> STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(STT);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
dumpState(FileIO& file, const string& pathName) const {

  // Ancestor method.
  SolidCRKSPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.write(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
void
SolidCRKSPHHydroBaseRZ::
restoreState(const FileIO& file, const string& pathName) {
 
  // Ancestor method.
  SolidCRKSPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.read(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

}
}

