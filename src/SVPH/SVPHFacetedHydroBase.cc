//---------------------------------Spheral++----------------------------------//
// SVPHFacetedHydroBase -- The SVPH hydrodynamic package for Spheral++.
//
// Created by JMO, Sun Jul 28 20:57:01 PDT 2013
//----------------------------------------------------------------------------//
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "SVPHFacetedHydroBase.hh"
#include "computeSVPHCorrectionsOnFaces.hh"
#include "SVPHCorrectionsPolicy.hh"
#include "SPH/computeSumVoronoiCellMassDensity.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "SVPH/SpecificThermalEnergyVolumePolicy.hh"
#include "SVPH/CompatibleFaceSpecificThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/PositionPolicy.hh"
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
#include "FileIO/FileIO.hh"
#include "Mesh/Mesh.hh"

namespace Spheral {
namespace SVPHSpace {

using namespace std;
using PhysicsSpace::GenericHydro;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using NeighborSpace::Neighbor;
using MeshSpace::Mesh;

using PhysicsSpace::MassDensityType;
using PhysicsSpace::HEvolutionType;

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
SVPHFacetedHydroBase<Dimension>::
SVPHFacetedHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                     const TableKernel<Dimension>& W,
                     ArtificialViscosity<Dimension>& Q,
                     const double cfl,
                     const bool useVelocityMagnitudeForDt,
                     const bool compatibleEnergyEvolution,
                     const bool XSVPH,
                     const bool linearConsistent,
                     const MassDensityType densityUpdate,
                     const HEvolutionType HUpdate,
                     const Vector& xmin,
                     const Vector& xmax):
  GenericHydro<Dimension>(W, W, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mXSVPH(XSVPH),
  mLinearConsistent(linearConsistent),
  mXmin(xmin),
  mXmax(xmax),
  mMeshPtr(MeshPtr(new Mesh<Dimension>())),
  // mA(),
  // mB(),
  // mGradB(),
  mTimeStepMask(FieldList<Dimension, int>::Copy),
  mPressure(FieldList<Dimension, Scalar>::Copy),
  mSoundSpeed(FieldList<Dimension, Scalar>::Copy),
  mVolume(FieldList<Dimension, Scalar>::Copy),
  mSpecificThermalEnergy0(FieldList<Dimension, Scalar>::Copy),
  mVolume0(FieldList<Dimension, Scalar>::Copy),
  mHideal(FieldList<Dimension, SymTensor>::Copy),
  mMaxViscousPressure(FieldList<Dimension, Scalar>::Copy),
  mMassDensitySum(FieldList<Dimension, Scalar>::Copy),
  mWeightedNeighborSum(FieldList<Dimension, Scalar>::Copy),
  mMassSecondMoment(FieldList<Dimension, SymTensor>::Copy),
  mXSVPHDeltaV(FieldList<Dimension, Vector>::Copy),
  mDxDt(FieldList<Dimension, Vector>::Copy),
  mDvDt(FieldList<Dimension, Vector>::Copy),
  mDmassDensityDt(FieldList<Dimension, Scalar>::Copy),
  mDspecificThermalEnergyDt(FieldList<Dimension, Scalar>::Copy),
  mDHDt(FieldList<Dimension, SymTensor>::Copy),
  mDvDx(FieldList<Dimension, Tensor>::Copy),
  mInternalDvDx(FieldList<Dimension, Tensor>::Copy),
  mFaceVelocity(FieldList<Dimension, vector<Vector> >::Copy),
  mFaceVelocity0(FieldList<Dimension, vector<Vector> >::Copy),
  mFaceForce(FieldList<Dimension, vector<Vector> >::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SVPHFacetedHydroBase<Dimension>::
~SVPHFacetedHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Create storage for the pressure and sound speed.
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);

  // Initialize the pressure and sound speed.
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);

  // Construct the mesh and volumes.
  NodeList<Dimension> voidNodes("internal void", 0, 0);
  vector<const NodeList<Dimension>*> nodeLists(dataBase.nodeListBegin(), dataBase.nodeListEnd());
  nodeLists.push_back(&voidNodes);
  // std::sort(nodeLists.begin(), nodeLists.end(), typename NodeListRegistrar<Dimension>::NodeListComparator());
  MeshSpace::generateMesh<Dimension,
                          typename vector<const NodeList<Dimension>*>::iterator,
                          ConstBoundaryIterator>
    (nodeLists.begin(),
     nodeLists.end(),
     this->boundaryBegin(),
     this->boundaryEnd(),
     mXmin,
     mXmax,
     true,              // mesh ghost nodes
     false,             // generateVoid
     true,              // generateParallelConnectivity
     true,              // removeBoundaryZones
     2.0,               // voidThreshold
     *mMeshPtr,
     voidNodes);
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  for (unsigned nodeListi = 0; nodeListi != dataBase.numFluidNodeLists(); ++nodeListi) {
    const unsigned n = mVolume[nodeListi]->numInternalElements();
    const unsigned offset = mMeshPtr->offset(nodeListi);
    for (unsigned i = 0; i != n; ++i) {
      mVolume(nodeListi, i) = mMeshPtr->zone(offset + i).volume();
    }
  }

  // // Compute the SVPH normalization and corrections.
  // computeSVPHCorrectionsOnFaces<Dimension>(dataBase.connectivityMap(),
  //                                          this->kernel(),
  //                                          mVolume, dataBase.fluidPosition(), dataBase.fluidHfield(),
  //                                          mA, mB, mGradB);

}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  // dataBase.resizeFluidFieldList(mA, vector<Scalar>(), HydroFieldNames::A_CSPH);
  // dataBase.resizeFluidFieldList(mB, vector<Vector>(), HydroFieldNames::B_CSPH);
  // dataBase.resizeFluidFieldList(mGradB, vector<Tensor>(), HydroFieldNames::gradB_CSPH);
  dataBase.resizeFluidFieldList(mVolume, 0.0, HydroFieldNames::volume, false);
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);

  // Take a snapshot of the starting volume.
  mVolume0 = mVolume;
  mVolume0.copyFields();
  for (unsigned i = 0; i != mVolume0.size(); ++i) {
    mVolume0[i]->name(HydroFieldNames::volume + "0");
  }

  // Now register away.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {

    // Mass.
    state.enroll((*itr)->mass());

    // Mass density.
    if (densityUpdate() == PhysicsSpace::IntegrateDensity) {
      PolicyPointer rhoPolicy(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                           (*itr)->rhoMax()));
      state.enroll((*itr)->massDensity(), rhoPolicy);
    } else {
      PolicyPointer rhoPolicy(new ReplaceBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                         (*itr)->rhoMax()));
      state.enroll((*itr)->massDensity(), rhoPolicy);
    }

    // Mesh and volume.
    PolicyPointer meshPolicy(new MeshPolicy<Dimension>(*this, mXmin, mXmax, 2.0, true, false, true));
    PolicyPointer volumePolicy(new VolumePolicy<Dimension>());
    state.enrollMesh(mMeshPtr);
    state.enroll(HydroFieldNames::mesh, meshPolicy);
    state.enroll(*mVolume[nodeListi], volumePolicy);

    // // SVPH corrections.
    // // All of these corrections are computed in the same method/policy, so we register
    // // the A field with the update policy and the others just come along for the ride.
    // PolicyPointer Apolicy(new SVPHFaceCorrectionsPolicy<Dimension>(dataBase, this->kernel()));
    // state.enroll(*mA[nodeListi], Apolicy);
    // state.enroll(*mB[nodeListi]);
    // state.enroll(*mGradB[nodeListi]);

    // Register the position update.
    // PolicyPointer positionPolicy(new PositionPolicy<Dimension>());
    PolicyPointer positionPolicy(new IncrementState<Dimension, Vector>());
    state.enroll((*itr)->positions(), positionPolicy);

    // Are we using the compatible energy evolution scheme?
    // Register the H tensor.
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == PhysicsSpace::IntegrateH) {
      PolicyPointer Hpolicy(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
      state.enroll((*itr)->Hfield(), Hpolicy);
    } else {
      CHECK(HEvolution() == PhysicsSpace::IdealH);
      PolicyPointer Hpolicy(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
      state.enroll((*itr)->Hfield(), Hpolicy);
    }

    // Register the time step mask, initialized to 1 so that everything defaults to being
    // checked.
    state.enroll(*mTimeStepMask[nodeListi]);

    // Compute and register the pressure and sound speed.
    PolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
    PolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
    state.enroll(*mPressure[nodeListi], pressurePolicy);
    state.enroll(*mSoundSpeed[nodeListi], csPolicy);

    // Specific thermal energy.
    if (compatibleEnergyEvolution()) {
      PolicyPointer thermalEnergyPolicy(new CompatibleFaceSpecificThermalEnergyPolicy<Dimension>(this->kernel(), 
                                                                                                 dataBase,
                                                                                                 this->artificialViscosity(),
                                                                                                 mLinearConsistent));
      PolicyPointer velocityPolicy(new IncrementState<Dimension, Vector>(HydroFieldNames::position));
      state.enroll((*itr)->specificThermalEnergy(), thermalEnergyPolicy);
      state.enroll((*itr)->velocity(), velocityPolicy);
      // state.enroll(*mSpecificThermalEnergy0[nodeListi]);
    } else {
      // PolicyPointer thermalEnergyPolicy(new SpecificThermalEnergyVolumePolicy<Dimension>());
      PolicyPointer thermalEnergyPolicy(new IncrementState<Dimension, Scalar>());
      PolicyPointer velocityPolicy(new IncrementState<Dimension, Vector>());
      state.enroll((*itr)->specificThermalEnergy(), thermalEnergyPolicy);
      state.enroll((*itr)->velocity(), velocityPolicy);
      state.enroll(*mVolume0[nodeListi]);
    }

  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  typedef typename StateDerivatives<Dimension>::KeyType Key;
  const string DxDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const string DvDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  dataBase.resizeFluidFieldList(mMassDensitySum, 0.0, ReplaceState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSVPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  // dataBase.resizeFluidFieldList(mFaceVelocity, vector<Vector>(), HydroFieldNames::faceVelocity, false);
  // dataBase.resizeFluidFieldList(mFaceForce, vector<Vector>(), HydroFieldNames::faceForce, false);

  size_t i = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++i) {
    derivs.enroll(*mHideal[i]);
    derivs.enroll(*mMaxViscousPressure[i]);
    derivs.enroll(*mMassDensitySum[i]);
    derivs.enroll(*mWeightedNeighborSum[i]);
    derivs.enroll(*mMassSecondMoment[i]);
    derivs.enroll(*mXSVPHDeltaV[i]);

    // These two (the position and velocity updates) may be registered
    // by other physics packages as well, so we need to be careful
    // not to duplicate if so.
    const Key DxDtKey = State<Dimension>::buildFieldKey(DxDtName, (*itr)->name());
    const Key DvDtKey = State<Dimension>::buildFieldKey(DvDtName, (*itr)->name());
    if (not derivs.registered(DxDtKey)) derivs.enroll(*mDxDt[i]);
    if (not derivs.registered(DvDtKey)) derivs.enroll(*mDvDt[i]);

    derivs.enroll(*mDmassDensityDt[i]);
    derivs.enroll(*mDspecificThermalEnergyDt[i]);
    derivs.enroll(*mDHDt[i]);
    derivs.enroll(*mDvDx[i]);
    derivs.enroll(*mInternalDvDx[i]);
    // derivs.enroll(*mFaceVelocity[i]);
    // derivs.enroll(*mFaceForce[i]);
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro before trying to evaluateDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // Get the artificial viscosity and initialize it.
  // We depend on the caller knowing to finalize the ghost boundaries for the Q.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();
  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               this->kernel());
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  typedef typename Mesh<Dimension>::Zone Zone;
  typedef typename Mesh<Dimension>::Face Face;

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const Scalar W0 = W(0.0, 1.0);

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // The mesh.
  const Mesh<Dimension>& mesh = state.mesh();

  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);

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
  FieldList<Dimension, Vector> XSVPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  // FieldList<Dimension, vector<Vector> > faceVelocity = derivatives.fields(HydroFieldNames::faceVelocity, vector<Vector>());
  // FieldList<Dimension, vector<Vector> > faceForce = derivatives.fields(HydroFieldNames::faceForce, vector<Vector>());
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
  CHECK(XSVPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  // CHECK(faceVelocity.size() == numNodeLists);
  // CHECK(faceForce.size() == numNodeLists);

  // Prepare working arrays of face properties.
  const unsigned numFaces = mesh.numFaces();
  vector<Scalar> Pface(numFaces, 0.0);
  vector<Vector> dAface(numFaces, Vector::zero), posFace(numFaces, Vector::zero), velFace(numFaces, Vector::zero);
  vector<Tensor> Qface(numFaces, Tensor::zero);

  // Compute the SVPH corrections.
  vector<Scalar> A;
  vector<Vector> B;
  computeSVPHCorrectionsOnFaces(mesh, W,
                                volume, position, H,
                                A, B);
  if (not mLinearConsistent) B = vector<Vector>(numFaces, Vector::zero);

  // Walk the faces and sample their fluid properties.
  unsigned i, j, k, nodeListi, nodeListj, z1id, z2id;
  const SymTensor Hface = 1.0e100*SymTensor::one;
  for (k = 0; k != numFaces; ++k) {
    const Face& face = mesh.face(k);
    posFace[k] = face.position();
    dAface[k] = face.area() * face.unitNormal();

    const Scalar& Ai = A[k];
    const Vector& Bi = B[k];

    // Set the neighbors for this face.
    Neighbor<Dimension>::setMasterNeighborGroup(posFace[k], Hface,
                                                nodeLists.begin(), nodeLists.end(),
                                                W.kernelExtent());

    // Iterate over the NodeLists.
    for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      const NodeList<Dimension>& nodeList = *nodeLists[nodeListj];
      Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
      neighbor.setRefineNeighborList(posFace[k], Hface);
      for (typename Neighbor<Dimension>::const_iterator neighborItr = neighbor.refineNeighborBegin();
           neighborItr != neighbor.refineNeighborEnd();
           ++neighborItr) {
        j = *neighborItr;
      
        // Get the state for node j
        const Vector& rj = position(nodeListj, j);
        const Vector& vj = velocity(nodeListj, j);
        const Scalar& Pj = pressure(nodeListj, j);
        const SymTensor& Hj = H(nodeListj, j);
        const Scalar& Vj = volume(nodeListj, j);
        const Scalar Hdetj = Hj.Determinant();
        CHECK(Vj > 0.0);
        CHECK(Hdetj > 0.0);

        // Pair-wise kernel type stuff.
        const Vector rij = posFace[k] - rj;
        const Vector etaj = Hj*rij;
        const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);

        // Increment the face fluid properties.
        const Scalar VWRj = Vj*(1.0 + Bi.dot(rij))*Wj;
        velFace[k] += VWRj*vj;
        Pface[k] += VWRj*Pj;
      }
    }

    // Now for the Q...
    // Find the SVPH nodes on either side of this face.
    z1id = Mesh<Dimension>::positiveID(face.zone1ID());
    z2id = Mesh<Dimension>::positiveID(face.zone2ID());
    if (z1id != Mesh<Dimension>::UNSETID and
        z2id != Mesh<Dimension>::UNSETID) {
      mesh.lookupNodeListID(z1id, nodeListi, i);
      mesh.lookupNodeListID(z2id, nodeListj, j);

      // Net properties on the face.
      const Vector ri = 0.5*(position(nodeListi, i) + position(nodeListj, j));
      const Vector vi = 0.5*(velocity(nodeListi, i) + velocity(nodeListj, j));
      const Scalar rhoi = 0.5*(massDensity(nodeListi, i) + massDensity(nodeListj, j));
      const Scalar ci = 0.5*(soundSpeed(nodeListi, i) + soundSpeed(nodeListj, j));
      const SymTensor Hi = 0.5*(H(nodeListi, i) + H(nodeListj, j));

      for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const NodeList<Dimension>& nodeList = *nodeLists[nodeListj];
        Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
        neighbor.setRefineNeighborList(posFace[k], Hface);
        for (typename Neighbor<Dimension>::const_iterator neighborItr = neighbor.refineNeighborBegin();
             neighborItr != neighbor.refineNeighborEnd();
             ++neighborItr) {
          j = *neighborItr;
      
          // Get the state for node j
          const Vector& rj = position(nodeListj, j);
          const Vector& vj = velocity(nodeListj, j);
          const Scalar& rhoj = massDensity(nodeListj, j);
          const Scalar& cj = soundSpeed(nodeListj, j);
          const SymTensor& Hj = H(nodeListj, j);
          const Scalar& Vj = volume(nodeListj, j);
          const Scalar Hdetj = Hj.Determinant();
          CHECK(Vj > 0.0);
          CHECK(Hdetj > 0.0);

          // Pair-wise kernel type stuff.
          const Vector rij = posFace[k] - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);
          const Scalar VWRj = Vj*(1.0 + Bi.dot(rij))*Wj;

          // Get the face Q values (in this case P/rho^2).
          const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListj, j, nodeListj, j,
                                                    ri, etai, vi, rhoi, ci, Hi,
                                                    rj, etaj, vj, rhoj, cj, Hj);
          Qface[k] += 0.5*VWRj*(rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second);
          const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
          maxViscousPressure(nodeListj, j) = max(maxViscousPressure(nodeListj, j), Qj);
        }
      }
    }

    // Finish the face state.
    CHECK2(Ai >= 0.0, i << " " << Ai);
    velFace[k] *= Ai;
    Pface[k] *= Ai;
    Qface[k] *= Ai;

    // // Find the SVPH nodes on either side of this face.
    // z1id = Mesh<Dimension>::positiveID(face.zone1ID());
    // z2id = Mesh<Dimension>::positiveID(face.zone2ID());
    // if (z1id != Mesh<Dimension>::UNSETID and
    //     z2id != Mesh<Dimension>::UNSETID) {
    //   mesh.lookupNodeListID(z1id, nodeListi, i);
    //   mesh.lookupNodeListID(z2id, nodeListj, j);

    //   // Get the node properties.
    //   const Vector& ri = position(nodeListi, i);
    //   const Vector& vi = velocity(nodeListi, i);
    //   const Scalar& rhoi = massDensity(nodeListi, i);
    //   const Scalar& ci = soundSpeed(nodeListi, i);
    //   const SymTensor& Hi = H(nodeListi, i);

    //   const Vector& rj = position(nodeListj, j);
    //   const Vector& vj = velocity(nodeListj, j);
    //   const Scalar& rhoj = massDensity(nodeListj, j);
    //   const Scalar& cj = soundSpeed(nodeListj, j);
    //   const SymTensor& Hj = H(nodeListj, j);

    //   const Vector rij = ri - rj;
    //   const Vector etai = Hi*rij;
    //   const Vector etaj = Hj*rij;

    //   // Get the face Q values (in this case P/rho^2).
    //   const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
    //                                             ri, etai, vi, rhoi, ci, Hi,
    //                                             rj, etaj, vj, rhoj, cj, Hj);
    //   Qface[k] = 0.5*(rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second);
    //   const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
    //   const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
    //   maxViscousPressure(nodeListi, i) = max(maxViscousPressure(nodeListi, i), Qi);
    //   maxViscousPressure(nodeListj, j) = max(maxViscousPressure(nodeListj, j), Qj);
    // }
  }

  // Start our big loop over all FluidNodeLists.
  nodeListi = 0;
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

    // Iterate over the internal nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;
      const Zone& zonei = mesh.zone(nodeListi, i);
      const vector<int>& faceIDs = zonei.faceIDs();
      const unsigned nfaces = faceIDs.size();
      
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
      const Scalar& Vi = volume(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Vi > 0.0);

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
      Vector& XSVPHDeltaVi = XSVPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      // vector<Vector>& faceVeli = faceVelocity(nodeListi, i);
      // vector<Vector>& faceForcei = faceForce(nodeListi, i);

      Scalar Vsumi = 0.0;

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = (nodeListi == nodeListj ? 1.0 : 0.2);
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors and sum our properties across the faces.
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Get the state for node j
            const Vector& rj = position(nodeListj, j);
            const Scalar& mj = mass(nodeListj, j);
            const Vector& vj = velocity(nodeListj, j);
            const Scalar& rhoj = massDensity(nodeListj, j);
            const Scalar& epsj = specificThermalEnergy(nodeListj, j);
            const Scalar& Pj = pressure(nodeListj, j);
            const SymTensor& Hj = H(nodeListj, j);
            const Scalar& cj = soundSpeed(nodeListj, j);
            const Scalar& Vj = volume(nodeListj, j);
            const Scalar Hdetj = Hj.Determinant();
            CHECK(mj > 0.0);
            CHECK(rhoj > 0.0);
            CHECK(Vj > 0.0);

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

            // Zero'th and second moment of the node distribution -- used for the
            // ideal H calculation.
            const double rij2 = rij.magnitude2();
            const SymTensor thpt = rij.selfdyad()/(rij2 + 1.0e-10) / FastMath::square(Dimension::pownu12(rij2 + 1.0e-10));
            weightedNeighborSumi += fweightij*std::abs(gWi);
            massSecondMomenti += fweightij*gradWi.magnitude2()*thpt;

            // Contribution to the sum density (only if the same material).
            if (nodeListi == nodeListj) {
              rhoSumi += mj*Wi;
              Vsumi += Vj*Wi;
            }

          }
        }
      }

      // Walk the faces and measure the primary fluid derivatives.
      Tensor Qavg;
      Scalar Asum = 0.0;
      for (unsigned k = 0; k != nfaces; ++k) {  
        const unsigned fid = Mesh<Dimension>::positiveID(faceIDs[k]);
        const Vector dA = dAface[fid] * sgn(faceIDs[k]);
        const Vector fforce = Pface[fid]*dA + Qface[fid]*dA;
        Qavg += Qface[fid] * dA.magnitude();
        Asum += dA.magnitude();
        DvDti += fforce;
        DepsDti += fforce.dot(velFace[fid]);
        DvDxi += velFace[fid]*dA;
        // faceVeli.push_back(velFace[fid]);
        // faceForcei.push_back(-fforce);
      }
      CHECK(Asum > 0.0);
      Qavg /= Asum;
      // CHECK(faceVeli.size() == nfaces);
      // CHECK(faceForcei.size() == nfaces);

      // Finish the time derivatives.
      DvDti /= mi;
      DvDxi /= -Vi;
      localDvDxi = DvDxi;
      DrhoDti = -rhoi*DvDxi.Trace();
      if (mCompatibleEnergyEvolution) {
        DepsDti /= mi;
      } else {
        DepsDti = -(Pi + Qavg.Trace()/Dimension::nDim)/rhoi*DvDxi.Trace();
      }

      // Finish the density sum.
      rhoSumi = (rhoSumi + mi*W0*Hdeti)/(Vsumi + Vi*W0*Hdeti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSVPH or not.
      DxDti = vi;

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh,
                                                             maxNumNeighbors);
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                        weightedNeighborSumi,
                                                        massSecondMomenti,
                                                        numNeighborsi,
                                                        W,
                                                        hmin,
                                                        hmax,
                                                        hminratio,
                                                        nPerh,
                                                        maxNumNeighbors);

    }
  }
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    FieldList<Dimension, Vector> accelerations = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Base class finalization.
  GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (densityUpdate() == PhysicsSpace::RigorousSumDensity or
      densityUpdate() == PhysicsSpace::SumVoronoiCellDensity) {
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    SPHSpace::computeSumVoronoiCellMassDensity(connectivityMap, this->kernel(), position, mass, volume, H, massDensity);
  } else if (densityUpdate() == PhysicsSpace::SumDensity) {
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    FieldList<Dimension, Scalar> massDensitySum = derivs.fields(ReplaceState<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                                HydroFieldNames::massDensity, 0.0);
    massDensity.assignFields(massDensitySum);
  } else if (densityUpdate() == PhysicsSpace::VoronoiCellDensity) {
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    massDensity = mass / volume;
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Apply boundary conditions to the basic fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) {
    CHECK(state.fieldNameRegistered(HydroFieldNames::specificThermalEnergy + "0"));
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
  }

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    if (compatibleEnergyEvolution()) {
      (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  FieldList<Dimension, Scalar> volume = state.fields(HydroFieldNames::volume, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) {
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
  }

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    if (compatibleEnergyEvolution()) {
      (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
    }
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mVolume, pathName + "/volume");
  file.write(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.write(mHideal, pathName + "/Hideal");
  file.write(mMassDensitySum, pathName + "/massDensitySum");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mXSVPHDeltaV, pathName + "/XSVPHDeltaV");

  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SVPHFacetedHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mVolume, pathName + "/volume");
  file.read(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.read(mHideal, pathName + "/Hideal");
  file.read(mMassDensitySum, pathName + "/massDensitySum");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mXSVPHDeltaV, pathName + "/XSVPHDeltaV");

  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
}

}
}

