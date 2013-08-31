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
#include "computeSumVoronoiCellMassDensityFromFaces.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CopyState.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
#include "SVPH/SpecificThermalEnergyVolumePolicy.hh"
#include "SVPH/CompatibleFaceSpecificThermalEnergyPolicy.hh"
#include "SVPH/CellPressurePolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "SVPH/MeshIdealHPolicy.hh"
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
#include "Material/EquationOfState.hh"

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
using Material::EquationOfState;

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
                     const bool generateVoid,
                     const MassDensityType densityUpdate,
                     const HEvolutionType HUpdate,
                     const Scalar fcentroidal,
                     const Scalar fcellPressure,
                     const Vector& xmin,
                     const Vector& xmax):
  GenericHydro<Dimension>(W, W, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mXSVPH(XSVPH),
  mLinearConsistent(linearConsistent),
  mGenerateVoid(generateVoid),
  mfcentroidal(fcentroidal),
  mfcellPressure(fcellPressure),
  mXmin(xmin),
  mXmax(xmax),
  mMeshPtr(MeshPtr(new Mesh<Dimension>())),
  // mA(),
  // mB(),
  // mGradB(),
  mTimeStepMask(FieldList<Dimension, int>::Copy),
  mPressure(FieldList<Dimension, Scalar>::Copy),
  mCellPressure(FieldList<Dimension, Scalar>::Copy),
  mSoundSpeed(FieldList<Dimension, Scalar>::Copy),
  mVolume(FieldList<Dimension, Scalar>::Copy),
  mSpecificThermalEnergy0(FieldList<Dimension, Scalar>::Copy),
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
  mFaceForce(FieldList<Dimension, vector<Vector> >::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
  // Delegate range checking to our assignment methods.
  this->fcentroidal(fcentroidal);
  this->fcellPressure(fcellPressure);
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

  typedef typename Mesh<Dimension>::Zone Zone;

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
     mGenerateVoid,     // generateVoid
     true,              // generateParallelConnectivity
     (not mGenerateVoid),              // removeBoundaryZones
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

  // Make a pass through the H tensors and initialize them to the "ideal" value.
  if (Process::getRank() == 0) cout << "SVPHFacetedHydro initializing H tensors..." << endl;
  FieldList<Dimension, SymTensor> H = dataBase.globalHfield();
  const unsigned numNodeLists = H.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = H[nodeListi]->nodeList();
    const unsigned n = nodeList.numInternalNodes();
    const Scalar hmin = nodeList.hmin();
    const Scalar hmax = nodeList.hmax();
    const Scalar hminratio = nodeList.hminratio();
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();
    for (unsigned i = 0; i != n; ++i) {
      const Zone& zonei = mMeshPtr->zone(nodeListi, i);
      H(nodeListi, i) = mSmoothingScaleMethod.idealSmoothingScale(H(nodeListi, i),
                                                                  *mMeshPtr,
                                                                  zonei,
                                                                  hmin,
                                                                  hmax,
                                                                  hminratio,
                                                                  nPerh);
    }
  }

  // // Compute the SVPH normalization and corrections.
  // computeSVPHCorrectionsOnFaces<Dimension>(dataBase.connectivityMap(),
  //                                          this->kernel(),
  //                                          mVolume, dataBase.fluidPosition(), dataBase.fluidHfield(),
  //                                          mA, mB, mGradB);

  // If needed, initialize the cell pressure.
  mCellPressure = dataBase.newFluidFieldList(0.0, "Cell" + HydroFieldNames::pressure);
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    Field<Dimension, Scalar>& P = **mCellPressure.fieldForNodeList(**itr);
    const Field<Dimension, Scalar>& vol = **mVolume.fieldForNodeList(**itr);
    const Field<Dimension, Scalar>& mass = (*itr)->mass();
    const Field<Dimension, Scalar>& eps = (*itr)->specificThermalEnergy();
    Field<Dimension, Scalar> rho = mass/vol;
    const EquationOfState<Dimension>& eos = (*itr)->equationOfState();
    eos.setPressure(P, rho, eps);
  }
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
  typedef typename Mesh<Dimension>::Zone Zone;
  typedef typename Mesh<Dimension>::Face Face;

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  // dataBase.resizeFluidFieldList(mA, vector<Scalar>(), HydroFieldNames::A_CSPH);
  // dataBase.resizeFluidFieldList(mB, vector<Vector>(), HydroFieldNames::B_CSPH);
  // dataBase.resizeFluidFieldList(mGradB, vector<Tensor>(), HydroFieldNames::gradB_CSPH);
  dataBase.resizeFluidFieldList(mVolume, 0.0, HydroFieldNames::volume, false);
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);

  // If we're not using the cell by cell regularizing pressure, just make it point
  // at the normal pressure.
  if (mfcellPressure == 0.0) {
    mCellPressure = mPressure;
  }

  // If we're using the compatible energy discretization we need to copy initial state and
  // fill in the opposite node properties across faces.
  if (mCompatibleEnergyEvolution) {
    const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
    const FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
    const FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
    mSpecificThermalEnergy0.assignFields(dataBase.fluidSpecificThermalEnergy());
    mSpecificThermalEnergy0.copyFields();
    dataBase.resizeFluidFieldList(mSpecificThermalEnergy0, 0.0, HydroFieldNames::specificThermalEnergy + "0", false);
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
    PolicyPointer meshPolicy(new MeshPolicy<Dimension>(*this, mXmin, mXmax, 2.0, true, mGenerateVoid, (not mGenerateVoid)));
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

    // Velocity.
    PolicyPointer velocityPolicy(new IncrementState<Dimension, Vector>());
    state.enroll((*itr)->velocity(), velocityPolicy);

    // Are we using the compatible energy evolution scheme?
    // Register the H tensor.
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == PhysicsSpace::IntegrateH) {
      PolicyPointer Hpolicy(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
      state.enroll((*itr)->Hfield(), Hpolicy);
    } else {
      CHECK(HEvolution() == PhysicsSpace::IdealH);
      PolicyPointer Hpolicy(new MeshIdealHPolicy<Dimension>(mSmoothingScaleMethod,
                                                            (*itr)->hmin(),
                                                            (*itr)->hmax(),
                                                            (*itr)->hminratio(),
                                                            (*itr)->nodesPerSmoothingScale()));
      // PolicyPointer Hpolicy(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
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

    // The cell pressure for regularizing.
    if (mfcellPressure > 0.0) {
      PolicyPointer cellPressurePolicy(new CellPressurePolicy<Dimension>());
      state.enroll(*mCellPressure[nodeListi], cellPressurePolicy);
    } else {
      mCellPressure[nodeListi]->name("Cell" + HydroFieldNames::pressure); // Have to fix from the copy above.
      PolicyPointer cellPressurePolicy(new CopyState<Dimension, Scalar>(HydroFieldNames::pressure,
                                                                        "Cell" + HydroFieldNames::pressure));
      state.enroll(*mCellPressure[nodeListi], cellPressurePolicy);
    }

    // Specific thermal energy.
    if (compatibleEnergyEvolution()) {
      meshPolicy->addDependency(HydroFieldNames::specificThermalEnergy);
      velocityPolicy->addDependency(HydroFieldNames::position);
      velocityPolicy->addDependency(HydroFieldNames::specificThermalEnergy);
      PolicyPointer thermalEnergyPolicy(new CompatibleFaceSpecificThermalEnergyPolicy<Dimension>(this->kernel(), 
                                                                                                 dataBase,
                                                                                                 this->artificialViscosity(),
                                                                                                 mLinearConsistent));
      state.enroll((*itr)->specificThermalEnergy(), thermalEnergyPolicy);
      state.enroll(*mSpecificThermalEnergy0[nodeListi]);
    } else {
      PolicyPointer thermalEnergyPolicy(new IncrementState<Dimension, Scalar>());
      state.enroll((*itr)->specificThermalEnergy(), thermalEnergyPolicy);
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
  dataBase.resizeFluidFieldList(mFaceForce, vector<Vector>(), HydroFieldNames::faceForce, false);
  // dataBase.resizeFluidFieldList(mFaceAcceleration, vector<Vector>(), IncrementState<Dimension, Vector>::prefix() + "Face " + HydroFieldNames::velocity, false);

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
    derivs.enroll(*mFaceForce[i]);
    // derivs.enroll(*mFaceAcceleration[i]);
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

  // // Copy the starting specific thermal energy for compatible mode.
  // if (mCompatibleEnergyEvolution) {
  //   const FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  //   mSpecificThermalEnergy0.assignFields(dataBase.fluidSpecificThermalEnergy());
  //   mSpecificThermalEnergy0.copyFields();
  //   dataBase.resizeFluidFieldList(mSpecificThermalEnergy0, 0.0, HydroFieldNames::specificThermalEnergy + "0", false);
  // }
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

  // The set of NodeLists.
  const vector<const NodeList<Dimension>*> nodeLists(dataBase.fluidNodeListBegin(),
                                                     dataBase.fluidNodeListEnd());
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
  const FieldList<Dimension, Scalar> cellPressure = state.fields("Cell" + HydroFieldNames::pressure, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(volume.size() == numNodeLists);
  CHECK(cellPressure.size() == numNodeLists);

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
  FieldList<Dimension, vector<Vector> > faceForce = derivatives.fields(HydroFieldNames::faceForce, vector<Vector>());
  // FieldList<Dimension, vector<Vector> > faceAcceleration = derivatives.fields(IncrementState<Dimension, Vector>::prefix() + "Face " + HydroFieldNames::velocity, vector<Vector>());
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
  CHECK(faceForce.size() == numNodeLists);
  // CHECK(faceAcceleration.size() == numNodeLists);

  // Prepare working arrays of face properties.
  const unsigned numFaces = mesh.numFaces();
  vector<Scalar> rhoFace(numFaces, 0.0), csFace(numFaces, 0.0), Pface(numFaces, 0.0), PcellFace(numFaces, 0.0);
  vector<Vector> dAface(numFaces, Vector::zero), posFace(numFaces, Vector::zero), velFace(numFaces, Vector::zero);
  vector<Tensor> Qface(numFaces, Tensor::zero);
  vector<SymTensor> Hface(numFaces, SymTensor::zero);

  // Compute the SVPH corrections.
  vector<Scalar> A;
  vector<Vector> B;
  computeSVPHCorrectionsOnFaces(mesh, W,
                                volume, position, H,
                                this->boundaryBegin(),
                                this->boundaryEnd(),
                                A, B);
  if (not mLinearConsistent) B = vector<Vector>(numFaces, Vector::zero);

  // Walk the faces and sample their fluid properties.
  unsigned i, j, k, nodeListi, nodeListj, z1id, z2id;
  Scalar Hdetj;
  const SymTensor H0 = 1.0e100*SymTensor::one;
  for (k = 0; k != numFaces; ++k) {
    const Face& face = mesh.face(k);
    const Vector& Bi = B[k];

    // Look up the SVPH nodes either side of the Face.
    z1id = Mesh<Dimension>::positiveID(face.zone1ID());
    z2id = Mesh<Dimension>::positiveID(face.zone2ID());
    if (z1id > z2id) std::swap(z1id, z2id);
    CHECK(z1id != Mesh<Dimension>::UNSETID);
    mesh.lookupNodeListID(z1id, nodeListi, i);
    if (z2id == Mesh<Dimension>::UNSETID) {
      nodeListj = nodeListi;
      j = i;
    } else {
      mesh.lookupNodeListID(z2id, nodeListj, j);
    }

    // Properties on the face.
    dAface[k] = face.area() * face.unitNormal();
    posFace[k] = face.position();
    velFace[k] = 0.5*(velocity(nodeListi, i) + velocity(nodeListj, j));
    rhoFace[k] = 0.5*(massDensity(nodeListi, i) + massDensity(nodeListj, j));
    csFace[k] = 0.5*(soundSpeed(nodeListi, i) + soundSpeed(nodeListj, j));
    Hface[k] = 0.5*(H(nodeListi, i) + H(nodeListj, j));
    PcellFace[k] = 0.5*(cellPressure(nodeListi, i) + cellPressure(nodeListj, j));

    // Set the neighbors for this face.
    Neighbor<Dimension>::setMasterNeighborGroup(posFace[k], H0,
                                                nodeLists.begin(), nodeLists.end(),
                                                W.kernelExtent());

    // Iterate over the NodeLists.
    for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      const NodeList<Dimension>& nodeList = *nodeLists[nodeListj];
      Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
      neighbor.setRefineNeighborList(posFace[k], H0);
      for (typename Neighbor<Dimension>::const_iterator neighborItr = neighbor.refineNeighborBegin();
           neighborItr != neighbor.refineNeighborEnd();
           ++neighborItr) {
        j = *neighborItr;
      
        // Get the state for node j
        const Vector& rj = position(nodeListj, j);
        const Vector& vj = velocity(nodeListj, j);
        const Scalar& rhoj = massDensity(nodeListj, j);
        const Scalar& Pj = pressure(nodeListj, j);
        const Scalar& cj = soundSpeed(nodeListj, j);
        const SymTensor& Hj = H(nodeListj, j);
        const Scalar& Vj = volume(nodeListj, j);
        Hdetj = Hj.Determinant();
        CHECK(Vj > 0.0);
        CHECK(Hdetj > 0.0);

        // Pair-wise kernel type stuff.
        const Vector rij = posFace[k] - rj;
        const Vector etaj = Hj*rij;
        const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);
        const Scalar VWRj = Vj*(1.0 + Bi.dot(rij))*Wj;

        // Increment the face fluid properties.
        Pface[k] += VWRj*Pj;
      }
    }
  }

  // // Boundaries!
  // for (ConstBoundaryIterator itr = this->boundaryBegin();
  //      itr != this->boundaryEnd();
  //      ++itr) {
  //   (*itr)->enforceBoundary(rhoFace, mesh);
  //   (*itr)->enforceBoundary(csFace, mesh);
  //   (*itr)->enforceBoundary(Pface, mesh);
  //   (*itr)->enforceBoundary(velFace, mesh);
  //   (*itr)->enforceBoundary(Hface, mesh);
  // }

  // Determine the Q (requires correct boundary enforced velocities).
  for (k = 0; k != numFaces; ++k) {
    const Vector& Bi = B[k];

    // Set the neighbors for this face.
    Neighbor<Dimension>::setMasterNeighborGroup(posFace[k], H0,
                                                nodeLists.begin(), nodeLists.end(),
                                                W.kernelExtent());

    // Iterate over the NodeLists.
    for (nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
      const NodeList<Dimension>& nodeList = *nodeLists[nodeListj];
      Neighbor<Dimension>& neighbor = const_cast<Neighbor<Dimension>&>(nodeList.neighbor());
      neighbor.setRefineNeighborList(posFace[k], H0);
      for (typename Neighbor<Dimension>::const_iterator neighborItr = neighbor.refineNeighborBegin();
           neighborItr != neighbor.refineNeighborEnd();
           ++neighborItr) {
        j = *neighborItr;
      
        // Get the state for node j
        const Vector& rj = position(nodeListj, j);
        const Vector& vj = velocity(nodeListj, j);
        const Scalar& rhoj = massDensity(nodeListj, j);
        const Scalar& Pj = pressure(nodeListj, j);
        const Scalar& cj = soundSpeed(nodeListj, j);
        const SymTensor& Hj = H(nodeListj, j);
        const Scalar& Vj = volume(nodeListj, j);
        Hdetj = Hj.Determinant();
        CHECK(Vj > 0.0);
        CHECK(Hdetj > 0.0);

        // Pair-wise kernel type stuff.
        const Vector rij = posFace[k] - rj;
        const Vector etai = Hface[k]*rij;
        const Vector etaj = Hj*rij;
        const Scalar Wj = W.kernelValue(etaj.magnitude(), Hdetj);
        const Scalar VWRj = Vj*(1.0 + Bi.dot(rij))*Wj;

        // Get the face Q values (in this case P/rho^2).
        const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListj, j, nodeListj, j,
                                                  posFace[k], etai, velFace[k], rhoFace[k], csFace[k], Hface[k],
                                                  rj, etaj, vj, rhoj, cj, Hj);
        //Qface[k] += 0.5*VWRj*(rhoFace[k]*rhoFace[k]*QPiij.first + rhoj*rhoj*QPiij.second);
        Qface[k] += VWRj*rhoj*rhoj*QPiij.second;
        const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
        maxViscousPressure(nodeListj, j) = max(maxViscousPressure(nodeListj, j), Qj);
      }
    }
  }

  // // Boundaries!
  // for (ConstBoundaryIterator itr = this->boundaryBegin();
  //      itr != this->boundaryEnd();
  //      ++itr) {
  //   (*itr)->enforceBoundary(Qface, mesh);
  // }

  // Finish the face state.
  for (k = 0; k != numFaces; ++k) {
    const Scalar& Ai = A[k];
    CHECK2(Ai >= 0.0, i << " " << Ai);
    Pface[k] *= Ai;
    Qface[k] *= Ai;

    Pface[k] = (1.0 - mfcellPressure)*Pface[k] + mfcellPressure*PcellFace[k];
  }

  // Start our big loop over all FluidNodeLists.
  for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const NodeList<Dimension>& nodeList = *nodeLists[nodeListi];
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar hmin = nodeList.hmin();
    const Scalar hmax = nodeList.hmax();
    const Scalar hminratio = nodeList.hminratio();
    const int maxNumNeighbors = nodeList.maxNumNeighbors();
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();

    // Iterate over the internal nodes in this NodeList.
    const unsigned n = nodeList.numInternalNodes();
    for (unsigned i = 0; i != n; ++i) {
      const Zone& zonei = mesh.zone(nodeListi, i);
      const vector<int>& faceIDs = zonei.faceIDs();
      const unsigned nfaces = faceIDs.size();
      
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
      vector<Vector>& faceForcei = faceForce(nodeListi, i);

      // Walk the faces and measure the primary fluid derivatives.
      Tensor Qavg;
      Scalar Asum = 0.0;
      for (unsigned k = 0; k != nfaces; ++k) {  
        const unsigned fid = Mesh<Dimension>::positiveID(faceIDs[k]);
        const Vector dA = -(dAface[fid] * sgn(faceIDs[k]));
        CHECK(dA.dot(posFace[fid] - ri) < 0.0);
        const Vector fforce = Pface[fid]*dA + Qface[fid]*dA;
        const Scalar area = dAface[fid].magnitude();
        DvDti += fforce;
        DvDxi += velFace[fid]*dA;
        faceForcei.push_back(fforce);
        Asum += area;
        Qavg += area*Qface[fid];
        DxDti += area*(velFace[fid] - vi);
      }
      CHECK(faceForcei.size() == nfaces);
      CHECK(Asum > 0.0);
      Qavg /= Asum;

      // Finish the time derivatives.
      DvDti /= mi;
      DvDxi /= -Vi;
      localDvDxi = DvDxi;
      DrhoDti = -rhoi*DvDxi.Trace();
      DepsDti = -(Pi + Qavg.Trace()/Dimension::nDim)/rhoi*DvDxi.Trace();

      // Position update.
      if (this->XSVPH()) {
        DxDti = vi + DxDti/Asum;
      } else {
        DxDti = vi;
      }

      // Apply any centroidal filtering.
      const Vector drcent = mfcentroidal*(zonei.position() - ri);
      const Scalar flimitcent = min(1.0, DxDti.magnitude()*dt*safeInv(drcent.magnitude()));
      CHECK(flimitcent >= 0.0 and flimitcent <= 1.0);
      DxDti = (1.0 - mfcentroidal)*DxDti + drcent/dt*flimitcent;

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh,
                                                             maxNumNeighbors);
      Hideali = mSmoothingScaleMethod.idealSmoothingScale(Hi,
                                                          mesh,
                                                          zonei,
                                                          hmin,
                                                          hmax,
                                                          hminratio,
                                                          nPerh);
    }
  }

  // // We don't like high-frequency noise in the H field, so do a local filtering.
  // FieldList<Dimension, SymTensor> Havg = dataBase.newFluidFieldList(SymTensor::zero, "Havg");
  // for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //   const unsigned n = Hideal[nodeListi]->numInternalElements();
  //   for (i = 0; i != n; ++i) {
  //     const Zone& zone = mesh.zone(nodeListi, i);

  //     // Build the unique set of neighboring zones by nodes.
  //     vector<unsigned> otherZones;
  //     {
  //       const vector<unsigned>& nodeIDs = zone.nodeIDs();
  //       for (vector<unsigned>::const_iterator itr = nodeIDs.begin();
  //            itr != nodeIDs.end();
  //            ++itr) {
  //         const vector<unsigned>& nodeZoneIDs = mesh.node(*itr).zoneIDs();
  //         copy(nodeZoneIDs.begin(), nodeZoneIDs.end(), back_inserter(otherZones));
  //       }
  //       sort(otherZones.begin(), otherZones.end());
  //       otherZones.erase(unique(otherZones.begin(), otherZones.end()), otherZones.end());
  //     }
  //     CHECK(otherZones.size() > 0);

  //     // Build our averaged H tensor.
  //     for (vector<unsigned>::const_iterator itr = otherZones.begin();
  //          itr != otherZones.end();
  //          ++itr) {
  //       if (*itr == Mesh<Dimension>::UNSETID) {
  //         nodeListj = nodeListi;
  //         j = i;
  //       } else {
  //         mesh.lookupNodeListID(*itr, nodeListj, j);
  //       }
  //       Havg(nodeListi, i) += Hideal(nodeListj, j).Inverse();
  //     }
  //     Havg(nodeListi, i) = (Havg(nodeListi, i)/otherZones.size()).Inverse();
  //   }
  // }
  // Hideal.assignFields(Havg);

  // Finally, if we're using the compatible energy discretization we need to
  // fill in the opposite properties across faces.
  if (mCompatibleEnergyEvolution) {
    // for (nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    //   const NodeList<Dimension>& nodeList = *nodeLists[nodeListi];

    //   // Iterate over the internal nodes in this NodeList.
    //   const unsigned n = nodeList.numInternalNodes();
    //   for (unsigned i = 0; i != n; ++i) {
    //     const Zone& zonei = mesh.zone(nodeListi, i);
    //     const vector<int>& faceIDs = zonei.faceIDs();
    //     const unsigned nfaces = faceIDs.size();

    //     // Get the state for node i.
    //     vector<Vector>& faceAccelerationi = faceAcceleration(nodeListi, i);

    //     // Walk the faces.
    //     for (unsigned k = 0; k != nfaces; ++k) {  
    //       const unsigned fid = Mesh<Dimension>::positiveID(faceIDs[k]);
    //       const Face& face = mesh.face(fid);

    //       // Find the opposite node.
    //       const unsigned oppZoneID = Mesh<Dimension>::positiveID(face.oppositeZoneID(zonei.ID()));
    //       unsigned nodeListj = nodeListi, j = i;
    //       if (oppZoneID != Mesh<Dimension>::UNSETID) {
    //         mesh.lookupNodeListID(oppZoneID, nodeListj, j);
    //       }

    //       // Record the opposite node properties.
    //       faceAccelerationi.push_back(DvDt(nodeListj, j));
    //     }
    //     CHECK(faceAccelerationi.size() == nfaces);
    //   }
    // }

    // // Boundaries!
    // for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    //   for (ConstBoundaryIterator itr = this->boundaryBegin();
    //        itr != this->boundaryEnd();
    //        ++itr) {
    //     (*itr)->swapFaceValues(*faceAcceleration[nodeListi], *mMeshPtr);
    //   }
    // }

    // Boundaries!
    for (ConstBoundaryIterator itr = this->boundaryBegin();
         itr != this->boundaryEnd();
         ++itr) {
      (*itr)->applyFieldListGhostBoundary(DvDt);
    }
    for (ConstBoundaryIterator itr = this->boundaryBegin();
         itr != this->boundaryEnd();
         ++itr) {
      (*itr)->finalizeGhostBoundary();
    }
  }
}

//------------------------------------------------------------------------------
// Determine the timestep requirements for a hydro step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename SVPHFacetedHydroBase<Dimension>::TimeStepType
SVPHFacetedHydroBase<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   typename Dimension::Scalar currentTime) const {

  // Get some useful fluid variables from the DataBase.
  const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, Scalar());
  const FieldList<Dimension, int> mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> eps = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> maxViscousPressure = derivs.fields(HydroFieldNames::maxViscousPressure, 0.0);
  const FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const FieldList<Dimension, Vector> DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Initialize the return value to some impossibly high value.
  Scalar minDt = FLT_MAX;

  // Set up some history variables to track what set's our minimum Dt.
  Scalar lastMinDt = minDt;
  int lastNodeID;
  string lastNodeListName, reason;
  Scalar lastNodeScale, lastCs, lastRho, lastEps, lastDivVelocity, lastShearVelocity;
  Vector lastVelocity, lastAcc;

  // Loop over every fluid node.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator nodeListItr = dataBase.fluidNodeListBegin();
       nodeListItr != dataBase.fluidNodeListEnd();
       ++nodeListItr, ++nodeListi) {
    const FluidNodeList<Dimension>& fluidNodeList = **nodeListItr;
    const Scalar kernelExtent = fluidNodeList.neighbor().kernelExtent();
    CHECK(kernelExtent > 0.0);

    const unsigned n = fluidNodeList.numInternalNodes();
    for (unsigned i = 0; i != n; ++i) {

      // If this node is masked, don't worry about it.
      if (mask(nodeListi, i) == 1) {

        // Get this nodes length scale.  This is the only bit we specialize here from the
        // generic base class method.
        CHECK(vol(nodeListi, i) > 0.0);
        const Scalar nodeScale = 0.5*Dimension::rootnu(vol(nodeListi, i));

        // Sound speed limit.
        const double csDt = nodeScale/(soundSpeed(nodeListi, i) + FLT_MIN);
        if (csDt < minDt) {
          minDt = csDt;
          reason = "sound speed limit";
        }

        // Artificial viscosity effective sound speed.
        CHECK(rho(nodeListi, i) > 0.0);
        const Scalar csq = sqrt(maxViscousPressure(nodeListi, i)/rho(nodeListi, i));
        const double csqDt = nodeScale/(csq + FLT_MIN);
        if (csqDt < minDt) {
          minDt = csqDt;
          reason = "artificial viscosity sound speed limit";
        }

        // Velocity divergence limit.
        const Scalar divVelocity = DvDx(nodeListi, i).Trace();
        const double divvDt = 1.0/(std::abs(divVelocity) + FLT_MIN);
        if (divvDt < minDt) {
          minDt = divvDt;
          reason = "velocity divergence";
        }

        //     // Eigenvalues of the stress-strain tensor.
        //     Vector eigenValues = DvDx(nodeListi, i).Symmetric().eigenValues();
        //     Scalar maxValue = -1.0;
        //     for (int i = 0; i < Dimension::nDim; ++i) {
        //       maxValue = max(maxValue, fabs(eigenValues(i)));
        //     }
        //     const double strainDt = 1.0/(maxValue + FLT_MIN);
        //     if (strainDt < minDt) {
        //       minDt = strainDt;
        //       reason = "strain limit";
        //     }

        //     // Limit by the velocity shear.
        //     const Scalar shearVelocity = computeShearMagnitude(DvDx(nodeListi, i));
        //     const double shearDt = 1.0/(shearVelocity + FLT_MIN);
        //     if (shearDt < minDt) {
        //       minDt = shearDt;
        //       reason = "velocity shear limit";
        //     }

        // Total acceleration limit.
        const double dtAcc = sqrt(nodeScale/(DvDt(nodeListi, i).magnitude() + FLT_MIN));
        if (dtAcc < minDt) {
          minDt = dtAcc;
          reason = "total acceleration";
        }

        // If requested, limit against the absolute velocity.
        if (this->useVelocityMagnitudeForDt()) {
          const double velDt = nodeScale/(velocity(nodeListi, i).magnitude() + 1.0e-10);
          if (velDt < minDt) {
            minDt = velDt;
            reason = "velocity magnitude";
          }
        }

        if (minDt < lastMinDt) {
          lastMinDt = minDt;
          lastNodeID = i;
          lastNodeListName = fluidNodeList.name();
          lastNodeScale = nodeScale;
          lastCs = soundSpeed(nodeListi, i);
          lastAcc = DvDt(nodeListi, i);
          //      lastCsq = csq;
          lastRho = rho(nodeListi, i);
          lastEps = eps(nodeListi, i);
          lastVelocity = velocity(nodeListi, i);
          lastDivVelocity = divVelocity;
          //       lastShearVelocity = shearVelocity;
        }
      }
    }
  }

  stringstream reasonStream;
  reasonStream << "mindt = " << minDt << endl
	       << reason << endl
               << "  (nodeList, node) = (" << lastNodeListName << ", " << lastNodeID << ") | "
               << "  h = " << lastNodeScale << endl
               << "  cs = " << lastCs << endl
               << "  acc = " << lastAcc << endl
//               << " csq = " << lastCsq << endl
               << "  rho = " << lastRho << endl
               << "  eps = " << lastEps << endl
               << "  velocity = " << lastVelocity << endl
               << "  dtcs = " << lastNodeScale/(lastCs + FLT_MIN) << endl
               << "  divVelocity = " << lastDivVelocity << endl
//               << "  dtcsq = " << lastNodeScale/(lastCsq + FLT_MIN) << endl
//                << "  dtdivV = " << 1.0/(fabs(lastDivVelocity) + FLT_MIN) << endl
//                << "  shearVelocity = " << lastShearVelocity << endl
               << ends;

  // Now build the result.
  TimeStepType result(this->cfl()*minDt, reasonStream.str());

  return result;
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
  // if (compatibleEnergyEvolution()) {
  //   FieldList<Dimension, Vector> accelerations = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  // }
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
    const Mesh<Dimension>& mesh = state.mesh();
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    computeSumVoronoiCellMassDensityFromFaces(mesh, this->kernel(), dataBase, massDensity);
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
  FieldList<Dimension, Scalar> cellPressure = state.fields("Cell" + HydroFieldNames::pressure, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  FieldList<Dimension, Vector> DvDt;
  if (compatibleEnergyEvolution()) {
    CHECK(state.fieldNameRegistered(HydroFieldNames::specificThermalEnergy + "0"));
    CHECK(derivs.fieldNameRegistered(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity));
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
    DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
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
    (*boundaryItr)->applyFieldListGhostBoundary(cellPressure);
    if (compatibleEnergyEvolution()) {
      (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
      (*boundaryItr)->applyFieldListGhostBoundary(DvDt);
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
  FieldList<Dimension, Scalar> cellPressure = state.fields("Cell" + HydroFieldNames::pressure, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  FieldList<Dimension, Vector> DvDt;
  if (compatibleEnergyEvolution()) {
    CHECK(state.fieldNameRegistered(HydroFieldNames::specificThermalEnergy + "0"));
    CHECK(derivs.fieldNameRegistered(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity));
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
    DvDt = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
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
    (*boundaryItr)->enforceFieldListBoundary(cellPressure);
    if (compatibleEnergyEvolution()) {
      (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
      (*boundaryItr)->enforceFieldListBoundary(DvDt);
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
  file.write(mCellPressure, pathName + "/cellPressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.write(mHideal, pathName + "/Hideal");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
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

  file.write(mVolume, pathName + "/volume");
  file.write(mFaceForce, pathName + "/faceForce");
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
  file.read(mCellPressure, pathName + "/cellPressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.read(mHideal, pathName + "/Hideal");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
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

  file.read(mVolume, pathName + "/volume");
  file.read(mFaceForce, pathName + "/faceForce");
}

}
}

