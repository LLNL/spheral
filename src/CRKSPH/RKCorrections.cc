//---------------------------------Spheral++----------------------------------//
// RKCorrections
//
// Computes RK corrections for other physics packages
//----------------------------------------------------------------------------//
#include "RKCorrections.hh"

#include <limits>
#include "computeRKCorrections.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "Kernel/TableKernel.hh"
#include "LLNLPhysics/LLNLFieldNames.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
RKCorrections<Dimension>::
RKCorrections(const DataBase<Dimension>& dataBase,
              const TableKernel<Dimension>& W,
              const CRKOrder correctionOrder,
              const CRKVolumeType volumeType) :
  mDataBase(dataBase),
  mW(W),
  mCorrectionOrder(correctionOrder),
  mVolumeType(volumeType),
  mVolume(FieldStorageType::CopyFields),
  mA(FieldStorageType::CopyFields),
  mB(FieldStorageType::CopyFields),
  mC(FieldStorageType::CopyFields),
  mD(FieldStorageType::CopyFields),
  mGradA(FieldStorageType::CopyFields),
  mGradB(FieldStorageType::CopyFields),
  mGradC(FieldStorageType::CopyFields),
  mGradD(FieldStorageType::CopyFields),
  mHessA(FieldStorageType::CopyFields),
  mHessB(FieldStorageType::CopyFields),
  mHessC(FieldStorageType::CopyFields),
  mHessD(FieldStorageType::CopyFields),
  mSurfacePoint(FieldStorageType::CopyFields),
  mEtaVoidPoints(FieldStorageType::CopyFields),
  mCells(FieldStorageType::CopyFields),
  mCellFaceFlags(FieldStorageType::CopyFields),
  mVoidBoundary(mSurfacePoint, mEtaVoidPoints),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
RKCorrections<Dimension>::
~RKCorrections() {
}

//------------------------------------------------------------------------------
// Optional hook to initialize once when the problem is starting up
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  // Initialize state
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  
  // Initialize correction terms based on order
  switch (mCorrectionOrder) {
  case CRKOrder::CubicOrder:
    mD = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::D_CRKSPH);
    mGradD = dataBase.newFluidFieldList(FourthRankTensor::zero, HydroFieldNames::gradD_CRKSPH);
    mHessD = dataBase.newFluidFieldList(FifthRankTensor::zero, HydroFieldNames::hessD_CRKSPH);
  case CRKOrder::QuadraticOrder:
    mC = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::C_CRKSPH);
    mGradC = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::gradC_CRKSPH);
    mHessC = dataBase.newFluidFieldList(FourthRankTensor::zero, HydroFieldNames::hessC_CRKSPH);
  case CRKOrder::LinearOrder:
    mB = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::B_CRKSPH);
    mGradB = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::gradB_CRKSPH);
    mHessB = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::hessB_CRKSPH);
  default:
    mA = dataBase.newFluidFieldList(0.0, HydroFieldNames::A_CRKSPH);
    mGradA = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::gradA_CRKSPH);
    mHessA = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::hessA_CRKSPH);
  }

  // Initialize the Voronoi stuff
  mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::etaVoidPoints);
  
  // Get some more data
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();

  // Compute the volumes
  if (mVolumeType == CRKVolumeType::CRKMassOverDensity) {
    mVolume.assignFields(mass/massDensity);
  }
  else if (mVolumeType == CRKVolumeType::CRKSumVolume) {
    computeCRKSPHSumVolume(connectivityMap, mW, position, mass, H, mVolume);
  }
  else if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    mCells = dataBase.newFluidFieldList(FacetedVolume(), HydroFieldNames::cells);
    mCellFaceFlags = dataBase.newFluidFieldList(vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags);
    mVolume.assignFields(mass/massDensity);
    const FieldList<Dimension, typename Dimension::SymTensor> damage = dataBase.solidEffectiveDamage();
    computeVoronoiVolume(position, H, connectivityMap, damage,
                         vector<typename Dimension::FacetedVolume>(),               // no boundaries
                         vector<vector<typename Dimension::FacetedVolume> >(),      // no holes
                         vector<Boundary<Dimension>*>(this->boundaryBegin(),        // boundaries
                                                      this->boundaryEnd()),
                         FieldList<Dimension, typename Dimension::Scalar>(),        // no weights
                         mSurfacePoint, mVolume, mDeltaCentroid, mEtaVoidPoints,    // return values
                         mCells,                                                    // return cells
                         mCellFaceFlags);                                           // node cell multimaterial faces
  }
  else if (mVolumeType == CRKVolumeType::CRKHullVolume) {
    computeHullVolumes(connectivityMap, mW.kernelExtent(), position, H, mVolume);
  }
  else if (mVolumeType == CRKVolumeType::HVolume) {
    const Scalar nPerh = mVolume.nodeListPtrs()[0]->nodesPerSmoothingScale();
    computeHVolumes(nPerh, H, mVolume);
  }
  else {
    VERIFY2(false, "Unknown CRK volume weighting.");
  }

  // Apply boundaries to newly computed terms
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(mVolume);
    if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(mSurfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  
  // Compute the corrections
  computeRKCorrections(connectivityMap, mW, mVolume, position, H, mCorrectionOrder,
                       mA, mGradA, mHessA, mB, mGradB, mHessB,
                       mC, mGradC, mHessC, mD, mGradD, mHessD);
}

//------------------------------------------------------------------------------
// Register the state
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  state.enroll(mVolume);
  
  state.enroll(mA);
  state.enroll(mB);
  state.enroll(mC);
  state.enroll(mD);
  state.enroll(mGradA);
  state.enroll(mGradB);
  state.enroll(mGradC);
  state.enroll(mGradD);
  state.enroll(mHessA);
  state.enroll(mHessB);
  state.enroll(mHessC);
  state.enroll(mHessD);
  
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    state.enroll(mCells);
    state.enroll(mCellFaceFlags);
  }
  state.enroll(mSurfacePoint);
}

//------------------------------------------------------------------------------
// No derivatives to register
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  // 
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  //
}

//------------------------------------------------------------------------------
// No time step vote
//------------------------------------------------------------------------------
template<typename Dimension>
typename RKCorrections<Dimension>::TimeStepType
RKCorrections<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  return std::make_pair(std::numeric_limits<double>::max(), std::string("RKCorrections: no vote"));
}

//------------------------------------------------------------------------------
// Compute new RK corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // 
}

//------------------------------------------------------------------------------
// No derivatives to evaluate
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
}

//------------------------------------------------------------------------------
// Nothing to finalize
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
finalize(const Scalar time, 
         const Scalar dt,
         DataBase<Dimension>& dataBase, 
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {
}

//------------------------------------------------------------------------------
// Dump the current state to the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
dumpState(FileIO& file, const std::string& pathName) const {
  file.write(mVolume, pathName + "/Volume");
  
  file.write(mA, pathName + "/A");
  file.write(mA, pathName + "/B");
  file.write(mA, pathName + "/C");
  file.write(mA, pathName + "/D");
  file.write(mGradA, pathName + "/gradA");
  file.write(mGradB, pathName + "/gradB");
  file.write(mGradC, pathName + "/gradC");
  file.write(mGradD, pathName + "/gradD");
  file.write(mHessA, pathName + "/hessA");
  file.write(mHessB, pathName + "/hessB");
  file.write(mHessC, pathName + "/hessC");
  file.write(mHessD, pathName + "/hessD");
}

//------------------------------------------------------------------------------
// Restore the state from the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  file.read(mVolume, pathName + "/Volume");
  
  file.read(mA, pathName + "/A");
  file.read(mA, pathName + "/B");
  file.read(mA, pathName + "/C");
  file.read(mA, pathName + "/D");
  file.read(mGradA, pathName + "/gradA");
  file.read(mGradB, pathName + "/gradB");
  file.read(mGradC, pathName + "/gradC");
  file.read(mGradD, pathName + "/gradD");
  file.read(mHessA, pathName + "/hessA");
  file.read(mHessB, pathName + "/hessB");
  file.read(mHessC, pathName + "/hessC");
  file.read(mHessD, pathName + "/hessD");
}

} // end namespace Spheral
