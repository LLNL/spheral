//---------------------------------Spheral++----------------------------------//
// SolidSPHHydroBaseRZ -- The axisymmetric (RZ) SPH/ASPH solid material
//                        hydrodynamic package for Spheral++.
//
// This RZ version is a naive area-weighting implementation, nothing as
// highfalutin as the Garcia-Senz approach.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Mon May  9 11:01:51 PDT 2016
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "SolidSPHHydroBaseRZ.hh"
#include "DamagedNodeCoupling.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicyRZ.hh"
#include "Strength/SolidFieldNames.hh"
#include "Strength/SolidNodeList.hh"
#include "Strength/RZDeviatoricStressPolicy.hh"
#include "Strength/RZPlasticStrainPolicy.hh"
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
#include "FileIO/FileIO.hh"

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
// Compute the artificial tensile stress correction tensor for the given 
// stress tensor
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SolidSPHHydroBaseRZ::
SolidSPHHydroBaseRZ(const SmoothingScaleBase<Dim<2> >& smoothingScaleMethod,
                    ArtificialViscosity<Dim<2> >& Q,
                    const TableKernel<Dim<2> >& W,
                    const TableKernel<Dim<2> >& WPi,
                    const TableKernel<Dim<2> >& WGrad,
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
  SolidSPHHydroBase<Dim<2> >(smoothingScaleMethod, 
                             Q,
                             W,
                             WPi,
                             WGrad,
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
  mDeviatoricStressTT(FieldSpace::Copy),
  mDdeviatoricStressTTDt(FieldSpace::Copy) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SolidSPHHydroBaseRZ::
~SolidSPHHydroBaseRZ() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
initializeProblemStartup(DataBase<Dim<2> >& dataBase) {

  // Call the ancestor.
  SolidSPHHydroBase<Dim<2> >::initializeProblemStartup(dataBase);

  // Create storage for the state we're holding.
  mDeviatoricStressTT = dataBase.newFluidFieldList(0.0, SolidFieldNames::deviatoricStressTT);
  mDdeviatoricStressTTDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + SolidFieldNames::deviatoricStressTT);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
registerState(DataBase<Dim<2> >& dataBase,
              State<Dim<2> >& state) {

  typedef State<Dimension>::PolicyPointer PolicyPointer;

  // Call the ancestor.
  SolidSPHHydroBase<Dim<2> >::registerState(dataBase, state);

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
SolidSPHHydroBaseRZ::
registerDerivatives(DataBase<Dim<2> >& dataBase,
                    StateDerivatives<Dim<2> >& derivs) {

  // Call the ancestor method.
  SolidSPHHydroBase<Dim<2> >::registerDerivatives(dataBase, derivs);

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
SolidSPHHydroBaseRZ::
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
  const TableKernel<Dimension>& WG = this->GradKernel();
  const SmoothingScaleBase<Dimension>& smoothingScaleMethod = this->smoothingScaleMethod();

  // A few useful constants we'll use in the following loop.
  typedef Timing::Time Time;
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
  const FieldList<Dimension, Scalar> STT = state.fields(SolidFieldNames::deviatoricStressTT, 0.0);
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
  CHECK(STT.size() == numNodeLists);
  CHECK(mu.size() == numNodeLists);
  CHECK(damage.size() == numNodeLists);
  CHECK(gradDamage.size() == numNodeLists);
  CHECK(fragIDs.size() == numNodeLists);
  CHECK(pTypes.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
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
  FieldList<Dimension, Scalar> DSTTDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + SolidFieldNames::deviatoricStressTT, 0.0);
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

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
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

    // Build the functor we use to compute the effective coupling between nodes.
    DamagedNodeCoupling<Dimension> coupling(damage, gradDamage, H);

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
      const Scalar rhoi = massDensity(nodeListi, i);
      const Scalar epsi = specificThermalEnergy(nodeListi, i);
      const Scalar Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar ci = soundSpeed(nodeListi, i);
      const Scalar& omegai = omega(nodeListi, i);
      const SymTensor& Si = S(nodeListi, i);
      const Scalar STTi = STT(nodeListi, i);
      const Scalar& mui = mu(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar safeOmegai = safeInv(omegai, tiny);
      const int fragIDi = fragIDs(nodeListi, i);
      const int pTypei = pTypes(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      Scalar& rhoSumi = rhoSum(nodeListi, i);
      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      Tensor& Mi = M(nodeListi, i);
      Tensor& localMi = localM(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      Scalar& effViscousPressurei = effViscousPressure(nodeListi, i);
      Scalar& rhoSumCorrectioni = rhoSumCorrection(nodeListi, i);
      Scalar& viscousWorki = viscousWork(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Scalar& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
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
          const double fweightij = 1.0; // (nodeListi == nodeListj ? 1.0 : 0.2);
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

              // Get the state for node j.
              const Vector& posj = position(nodeListj, j);
              const Scalar rj = abs(posj.y());
              const Scalar circj = 2.0*M_PI*rj;
              const Scalar mj = mass(nodeListj, j);
              const Scalar mRZj = mj/circj;
              const Vector& vj = velocity(nodeListj, j);
              const Scalar rhoj = massDensity(nodeListj, j);
              const Scalar epsj = specificThermalEnergy(nodeListj, j);
              const Scalar Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar cj = soundSpeed(nodeListj, j);
              const Scalar& omegaj = omega(nodeListj, j);
              const SymTensor& Sj = S(nodeListj, j);
              const Scalar STTj = STT(nodeListj, j);
              const Scalar& muj = mu(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar safeOmegaj = safeInv(omegaj, tiny);
              const int fragIDj = fragIDs(nodeListj, j);
              const int pTypej = pTypes(nodeListj, j);
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Hdetj > 0.0);

              Scalar& rhoSumj = rhoSum(nodeListj, j);
              Vector& DxDtj = DxDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Tensor& Mj = M(nodeListj, j);
              Tensor& localMj = localM(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              Scalar& effViscousPressurej = effViscousPressure(nodeListj, j);
              Scalar& rhoSumCorrectionj = rhoSumCorrection(nodeListj, j);
              Scalar& viscousWorkj = viscousWork(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              Scalar& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
              Vector& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              Scalar& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              SymTensor& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Flag if this is a contiguous material pair or not.
              const bool sameMatij = (nodeListi == nodeListj and fragIDi == fragIDj);

              // Flag if at least one particle is free (0).
              const bool freeParticle = (pTypei == 0 or pTypej == 0);

              // Node displacement.
              const Vector xij = posi - posj;
              const Vector etai = Hi*xij;
              const Vector etaj = Hj*xij;
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
              const std::pair<double, double> WQWi = WQ.kernelAndGradValue(etaMagi, Hdeti);
              const Scalar WQi = WQWi.first;
              const Scalar gWQi = WQWi.second;
              const Vector gradWQi = gWQi*Hetai;
              const Vector gradWGi = WG.gradValue(etaMagi, Hdeti) * Hetai;

              const Vector Hetaj = Hj*etaj.unitVector();
              const std::pair<double, double> WWj = W.kernelAndGradValue(etaMagj, Hdetj);
              const Scalar Wj = WWj.first;
              const Scalar gWj = WWj.second;
              const Vector gradWj = gWj*Hetaj;
              const std::pair<double, double> WQWj = WQ.kernelAndGradValue(etaMagj, Hdetj);
              const Scalar WQj = WQWj.first;
              const Scalar gWQj = WQWj.second;
              const Vector gradWQj = gWQj*Hetaj;
              const Vector gradWGj = WG.gradValue(etaMagj, Hdetj) * Hetaj;

              // Determine how we're applying damage.
              const Scalar fDeffij = coupling(nodeListi, i, nodeListj, j);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double xij2 = xij.magnitude2();
              const SymTensor thpt = xij.selfdyad()/(xij2 + 1.0e-10) / FastMath::square(Dimension::pownu12(xij2 + 1.0e-10));
              weightedNeighborSumi += fweightij*abs(gWi);
              weightedNeighborSumj += fweightij*abs(gWj);
              massSecondMomenti += fweightij*gradWi.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWj.magnitude2()*thpt;

              // Contribution to the sum density (only if the same material).
              if (nodeListi == nodeListj) {
                rhoSumi += mRZj*Wi;
                rhoSumj += mRZi*Wj;
              }

              // Contribution to the sum density correction
              rhoSumCorrectioni += mRZj * WQi / rhoj ;
              rhoSumCorrectionj += mRZi * WQj / rhoi ;

              // Compute the pair-wise artificial viscosity.
              const Vector vij = vi - vj;
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Vector Qacci = 0.5*(QPiij.first *gradWQi);
              const Vector Qaccj = 0.5*(QPiij.second*gradWQj);
              const Scalar workQi = 0.5*(QPiij.first *vij).dot(gradWQi);
              const Scalar workQj = 0.5*(QPiij.second*vij).dot(gradWQj);
              // const Scalar workQi = vij.dot(Qacci);
              // const Scalar workQj = vij.dot(Qaccj);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += mRZj*Qi*WQi/rhoj;
              effViscousPressurej += mRZi*Qj*WQj/rhoi;
              viscousWorki += mRZj*workQi;
              viscousWorkj += mRZi*workQj;

              // Damage scaling of negative pressures.
              const Scalar Peffi = (Pi > 0.0 ? Pi : fDeffij*Pi);
              const Scalar Peffj = (Pj > 0.0 ? Pj : fDeffij*Pj);

              // Compute the stress tensors.
              SymTensor sigmai = -Peffi*SymTensor::one;
              SymTensor sigmaj = -Peffj*SymTensor::one;
              if (sameMatij) {
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
              const SymTensor sigmarhoi = safeOmegai*sigmai/(rhoi*rhoi);
              const SymTensor sigmarhoj = safeOmegaj*sigmaj/(rhoj*rhoj);
              const Vector deltaDvDt = fDeffij*(sigmarhoi*gradWi + sigmarhoj*gradWj) - Qacci - Qaccj;
              if (freeParticle) {
                DvDti += mRZj*deltaDvDt;
                DvDtj -= mRZi*deltaDvDt;
              }

              // Pair-wise portion of grad velocity.
              const Tensor deltaDvDxi = fDeffij*vij.dyad(gradWGi);
              const Tensor deltaDvDxj = fDeffij*vij.dyad(gradWGj);

              // Specific thermal energy evolution.
              DepsDti -= mRZj*(fDeffij*sigmarhoi.doubledot(deltaDvDxi.Symmetric()) - workQi);
              DepsDtj -= mRZi*(fDeffij*sigmarhoj.doubledot(deltaDvDxj.Symmetric()) - workQj);
              if (compatibleEnergy) {
                pairAccelerationsi.push_back( mRZj*deltaDvDt);
                pairAccelerationsj.push_back(-mRZi*deltaDvDt);
              }

              // Velocity gradient.
              DvDxi -= mRZj*deltaDvDxi;
              DvDxj -= mRZi*deltaDvDxj;
              if (sameMatij) {
                localDvDxi -= mRZj*deltaDvDxi;
                localDvDxj -= mRZi*deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (sameMatij) {
                const double wXSPHij = 0.5*(mRZi/rhoi*Wi + mRZj/rhoj*Wj);
                XSPHWeightSumi += wXSPHij;
                XSPHWeightSumj += wXSPHij;
                XSPHDeltaVi -= wXSPHij*vij;
                XSPHDeltaVj += wXSPHij*vij;
              }

              // Linear gradient correction term.
              Mi -= fDeffij*mRZj*xij.dyad(gradWGi);
              Mj -= fDeffij*mRZi*xij.dyad(gradWGj);
              if (sameMatij) {
                localMi -= fDeffij*mRZj*xij.dyad(gradWGi);
                localMj -= fDeffij*mRZi*xij.dyad(gradWGj);
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
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Add the self-contribution to density sum.
      rhoSumi += mRZi*W0*Hdeti;
      rhoSumi /= circi;

      // Add the self-contribution to density sum correction.
      rhoSumCorrectioni += mRZi*WQ0*Hdeti/rhoi ;

      // Correct the effective viscous pressure.
      effViscousPressurei /= rhoSumCorrectioni ;

      // Finish the acceleration.
      const Scalar zetai = abs((Hi*posi).y());
      const Scalar hri = ri*safeInvVar(zetai);
      const Scalar riInv = safeInv(ri, 0.01*hri);
      const Vector deltaDvDti(Si(1,0)/rhoi*riInv,
                              (Si(1,1) - STTi)/rhoi*riInv);
      DvDti += deltaDvDti;
      pairAccelerationsi.push_back(deltaDvDti);

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
      XSPHWeightSumi += Hdeti*mRZi/rhoi*W0;
      CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      XSPHDeltaVi /= XSPHWeightSumi;
      const Scalar vri = vi.y() + XSPHDeltaVi.y();
      DrhoDti = -rhoi*(DvDxi.Trace() + vri*riInv);

      // Finish the specific thermal energy evolution.
      DepsDti += (STTi - Pi)/rhoi*vri*riInv;

      // If needed finish the total energy derivative.
      if (this->mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

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

      // Determine the deviatoric stress evolution.
      const Scalar deformationTT = vi.y()*riInv;
      const SymTensor deformation = localDvDxi.Symmetric();
      const Tensor spin = localDvDxi.SkewSymmetric();
      const SymTensor deviatoricDeformation = deformation - ((deformation.Trace() + deformationTT)/3.0)*SymTensor::one;
      const SymTensor spinCorrection = (spin*Si + Si*spin).Symmetric();
      DSDti = spinCorrection + (2.0*mui)*deviatoricDeformation;
      DSTTDti = 2.0*mui*(deformationTT - (deformation.Trace() + deformationTT)/3.0);

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
SolidSPHHydroBaseRZ::
finalize(const Dim<2>::Scalar time,
         const Dim<2>::Scalar dt,
         DataBase<Dim<2> >& dataBase,
         State<Dim<2> >& state,
         StateDerivatives<Dim<2> >& derivs) {

  // If we're going to do the SPH summation density, we need to convert the mass
  // to mass per unit length first.
  if (densityUpdate() == PhysicsSpace::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::CorrectedSumDensity) {
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
  SPHHydroBase<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Now convert back to true masses and mass densities.  We also apply the RZ
  // correction factor to the mass density.
  if (densityUpdate() == PhysicsSpace::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::CorrectedSumDensity) {
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
SolidSPHHydroBaseRZ::
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

  // Apply ordinary SPH BCs.
  SolidSPHHydroBase<Dim<2> >::applyGhostBoundaries(state, derivs);
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
SolidSPHHydroBaseRZ::
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

  // Apply ordinary SPH BCs.
  SolidSPHHydroBase<Dim<2> >::enforceBoundaries(state, derivs);

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
SolidSPHHydroBaseRZ::
dumpState(FileIO& file, const string& pathName) const {

  // Ancestor method.
  SolidSPHHydroBase<Dimension>::dumpState(file, pathName);

  file.write(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.write(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
void
SolidSPHHydroBaseRZ::
restoreState(const FileIO& file, const string& pathName) {
 
  // Ancestor method.
  SolidSPHHydroBase<Dimension>::restoreState(file, pathName);

  file.read(mDeviatoricStressTT, pathName + "/deviatoricStressTT");
  file.read(mDdeviatoricStressTTDt, pathName + "/DdeviatoricStressTTDt");
}

}
}

