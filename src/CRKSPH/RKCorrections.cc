//---------------------------------Spheral++----------------------------------//
// RKCorrections
//
// Computes RK corrections for other physics packages
//----------------------------------------------------------------------------//
#include "RKCorrections.hh"

#include <limits>
#include "computeRKCorrections.hh"
#include "computeRKVolumes.hh"
#include "Boundary/Boundary.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"

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
    mD = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::D_RK);
    mGradD = dataBase.newFluidFieldList(FourthRankTensor::zero, HydroFieldNames::gradD_RK);
    mHessD = dataBase.newFluidFieldList(FifthRankTensor::zero, HydroFieldNames::hessD_RK);
  case CRKOrder::QuadraticOrder:
    mC = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::C_RK);
    mGradC = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::gradC_RK);
    mHessC = dataBase.newFluidFieldList(FourthRankTensor::zero, HydroFieldNames::hessC_RK);
  case CRKOrder::LinearOrder:
    mB = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::B_RK);
    mGradB = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::gradB_RK);
    mHessB = dataBase.newFluidFieldList(ThirdRankTensor::zero, HydroFieldNames::hessB_RK);
  default:
    mA = dataBase.newFluidFieldList(0.0, HydroFieldNames::A_RK);
    mGradA = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::gradA_RK);
    mHessA = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::hessA_RK);
  }
  
  // Initialize the Voronoi stuff
  mSurfacePoint = dataBase.newFluidFieldList(0, HydroFieldNames::surfacePoint);
  mEtaVoidPoints = dataBase.newFluidFieldList(std::vector<Vector>(), HydroFieldNames::etaVoidPoints);
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    mCells = dataBase.newFluidFieldList(FacetedVolume(), HydroFieldNames::cells);
    mCellFaceFlags = dataBase.newFluidFieldList(std::vector<CellFaceFlag>(), HydroFieldNames::cellFaceFlags);
  }
  mDeltaCentroid = dataBase.newFluidFieldList(Vector::zero, "delta centroid");
  
  // Get some more data
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  const FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> damage;
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    damage = dataBase.solidEffectiveDamage();
  }
  
  // Compute the volumes
  computeRKVolumes(connectivityMap, mW,
                   position, mass, massDensity, H, damage,
                   this->boundaryConditions(), mVolumeType,
                   mSurfacePoint, mDeltaCentroid, mEtaVoidPoints, mCells, mCellFaceFlags,
                   mVolume);
  
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
       boundItr != this->boundaryEnd(); ++boundItr) {
    (*boundItr)->finalizeGhostBoundary();
  }
  
  // Compute corrections
  computeRKCorrections(connectivityMap, mW, mVolume, position, H, mCorrectionOrder,
                       mA, mB, mC, mD,
                       mGradA, mGradB, mGradC, mGradD,
                       mHessA, mHessB, mHessC, mHessD);
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
  
  state.enroll(mSurfacePoint);
  state.enroll(mEtaVoidPoints);
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    state.enroll(mCells);
    state.enroll(mCellFaceFlags);
  }
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
  // Get state variables
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto A = state.fields(HydroFieldNames::A_RK, 0.0);
  auto B = state.fields(HydroFieldNames::B_RK, Vector::zero);
  auto C = state.fields(HydroFieldNames::C_RK, Tensor::zero);
  auto D = state.fields(HydroFieldNames::D_RK, ThirdRankTensor::zero);
  auto gradA = state.fields(HydroFieldNames::gradA_RK, Vector::zero);
  auto gradB = state.fields(HydroFieldNames::gradB_RK, Tensor::zero);
  auto gradC = state.fields(HydroFieldNames::gradC_RK, ThirdRankTensor::zero);
  auto gradD = state.fields(HydroFieldNames::gradD_RK, FourthRankTensor::zero);
  auto hessA = state.fields(HydroFieldNames::hessA_RK, Tensor::zero);
  auto hessB = state.fields(HydroFieldNames::hessB_RK, ThirdRankTensor::zero);
  auto hessC = state.fields(HydroFieldNames::hessC_RK, FourthRankTensor::zero);
  auto hessD = state.fields(HydroFieldNames::hessD_RK, FifthRankTensor::zero);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  auto etaVoidPoints = state.fields(HydroFieldNames::etaVoidPoints, std::vector<Vector>());

  // Apply ghost boundary conditions
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(vol);
    (*boundaryItr)->applyFieldListGhostBoundary(A);
    (*boundaryItr)->applyFieldListGhostBoundary(B);
    (*boundaryItr)->applyFieldListGhostBoundary(C);
    (*boundaryItr)->applyFieldListGhostBoundary(D);
    (*boundaryItr)->applyFieldListGhostBoundary(gradA);
    (*boundaryItr)->applyFieldListGhostBoundary(gradB);
    (*boundaryItr)->applyFieldListGhostBoundary(gradC);
    (*boundaryItr)->applyFieldListGhostBoundary(gradD);
    (*boundaryItr)->applyFieldListGhostBoundary(hessA);
    (*boundaryItr)->applyFieldListGhostBoundary(hessB);
    (*boundaryItr)->applyFieldListGhostBoundary(hessC);
    (*boundaryItr)->applyFieldListGhostBoundary(hessD);
    (*boundaryItr)->applyFieldListGhostBoundary(surfacePoint);
    (*boundaryItr)->applyFieldListGhostBoundary(etaVoidPoints);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // Get state variables
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  auto A = state.fields(HydroFieldNames::A_RK, 0.0);
  auto B = state.fields(HydroFieldNames::B_RK, Vector::zero);
  auto C = state.fields(HydroFieldNames::C_RK, Tensor::zero);
  auto D = state.fields(HydroFieldNames::D_RK, ThirdRankTensor::zero);
  auto gradA = state.fields(HydroFieldNames::gradA_RK, Vector::zero);
  auto gradB = state.fields(HydroFieldNames::gradB_RK, Tensor::zero);
  auto gradC = state.fields(HydroFieldNames::gradC_RK, ThirdRankTensor::zero);
  auto gradD = state.fields(HydroFieldNames::gradD_RK, FourthRankTensor::zero);
  auto hessA = state.fields(HydroFieldNames::hessA_RK, Tensor::zero);
  auto hessB = state.fields(HydroFieldNames::hessB_RK, ThirdRankTensor::zero);
  auto hessC = state.fields(HydroFieldNames::hessC_RK, FourthRankTensor::zero);
  auto hessD = state.fields(HydroFieldNames::hessD_RK, FifthRankTensor::zero);

  // Enforce boundary conditions
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(vol);
    (*boundaryItr)->enforceFieldListBoundary(A);
    (*boundaryItr)->enforceFieldListBoundary(B);
    (*boundaryItr)->enforceFieldListBoundary(C);
    (*boundaryItr)->enforceFieldListBoundary(D);
    (*boundaryItr)->enforceFieldListBoundary(gradA);
    (*boundaryItr)->enforceFieldListBoundary(gradB);
    (*boundaryItr)->enforceFieldListBoundary(gradC);
    (*boundaryItr)->enforceFieldListBoundary(gradD);
    (*boundaryItr)->enforceFieldListBoundary(hessA);
    (*boundaryItr)->enforceFieldListBoundary(hessB);
    (*boundaryItr)->enforceFieldListBoundary(hessC);
    (*boundaryItr)->enforceFieldListBoundary(hessD);
  }    
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
// Compute new volumes
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  // Get data
  const auto& W = mW;
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  damage = state.fields(SolidFieldNames::effectiveTensorDamage, SymTensor::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto surfacePoint = state.fields(HydroFieldNames::surfacePoint, 0);
  FieldList<Dimension, FacetedVolume> cells;
  FieldList<Dimension, std::vector<CellFaceFlag>> cellFaceFlags;
  if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
    cells = state.fields(HydroFieldNames::cells, FacetedVolume());
    cellFaceFlags = state.fields(HydroFieldNames::cellFaceFlags, std::vector<CellFaceFlag>());
  }
  
  // Compute volumes
  computeRKVolumes(connectivityMap, W,
                   position, mass, massDensity, H, damage,
                   this->boundaryConditions(), mVolumeType,
                   surfacePoint, mDeltaCentroid, mEtaVoidPoints, cells, cellFaceFlags,
                   volume);
  
  // Apply ghost boundaries to Voronoi stuff
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(volume);
    if (mVolumeType == CRKVolumeType::CRKVoronoiVolume) {
      (*boundItr)->applyFieldListGhostBoundary(cells);
      (*boundItr)->applyFieldListGhostBoundary(surfacePoint);
      (*boundItr)->applyFieldListGhostBoundary(mEtaVoidPoints);
    }
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();
  
}

//------------------------------------------------------------------------------
// Compute new RK corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  // Get data
  const auto& W = mW;
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto A = state.fields(HydroFieldNames::A_RK, 0.0);
  auto B = state.fields(HydroFieldNames::B_RK, Vector::zero);
  auto C = state.fields(HydroFieldNames::C_RK, Tensor::zero);
  auto D = state.fields(HydroFieldNames::D_RK, ThirdRankTensor::zero);
  auto gradA = state.fields(HydroFieldNames::gradA_RK, Vector::zero);
  auto gradB = state.fields(HydroFieldNames::gradB_RK, Tensor::zero);
  auto gradC = state.fields(HydroFieldNames::gradC_RK, ThirdRankTensor::zero);
  auto gradD = state.fields(HydroFieldNames::gradD_RK, FourthRankTensor::zero);
  auto hessA = state.fields(HydroFieldNames::hessA_RK, Tensor::zero);
  auto hessB = state.fields(HydroFieldNames::hessB_RK, ThirdRankTensor::zero);
  auto hessC = state.fields(HydroFieldNames::hessC_RK, FourthRankTensor::zero);
  auto hessD = state.fields(HydroFieldNames::hessD_RK, FifthRankTensor::zero);
  
  // Compute corrections
  computeRKCorrections(connectivityMap, W, volume, position, H, mCorrectionOrder,
                       A, B, C, D,
                       gradA, gradB, gradC, gradD,
                       hessA, hessB, hessC, hessD);
  
  // Apply ghost boundaries to corrections
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(A);
    (*boundaryItr)->applyFieldListGhostBoundary(B);
    (*boundaryItr)->applyFieldListGhostBoundary(C);
    (*boundaryItr)->applyFieldListGhostBoundary(D);
    (*boundaryItr)->applyFieldListGhostBoundary(gradA);
    (*boundaryItr)->applyFieldListGhostBoundary(gradB);
    (*boundaryItr)->applyFieldListGhostBoundary(gradC);
    (*boundaryItr)->applyFieldListGhostBoundary(gradD);
    (*boundaryItr)->applyFieldListGhostBoundary(hessA);
    (*boundaryItr)->applyFieldListGhostBoundary(hessB);
    (*boundaryItr)->applyFieldListGhostBoundary(hessC);
    (*boundaryItr)->applyFieldListGhostBoundary(hessD);
  }
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
  
  // file.write(mA, pathName + "/A");
  // file.write(mA, pathName + "/B");
  // file.write(mA, pathName + "/C");
  // file.write(mA, pathName + "/D");
  // file.write(mGradA, pathName + "/gradA");
  // file.write(mGradB, pathName + "/gradB");
  // file.write(mGradC, pathName + "/gradC");
  // file.write(mGradD, pathName + "/gradD");
  // file.write(mHessA, pathName + "/hessA");
  // file.write(mHessB, pathName + "/hessB");
  // file.write(mHessC, pathName + "/hessC");
  // file.write(mHessD, pathName + "/hessD");
}

//------------------------------------------------------------------------------
// Restore the state from the given file
//------------------------------------------------------------------------------
template<typename Dimension>
void
RKCorrections<Dimension>::
restoreState(const FileIO& file, const std::string& pathName) {
  file.read(mVolume, pathName + "/Volume");
  
  // file.read(mA, pathName + "/A");
  // file.read(mA, pathName + "/B");
  // file.read(mA, pathName + "/C");
  // file.read(mA, pathName + "/D");
  // file.read(mGradA, pathName + "/gradA");
  // file.read(mGradB, pathName + "/gradB");
  // file.read(mGradC, pathName + "/gradC");
  // file.read(mGradD, pathName + "/gradD");
  // file.read(mHessA, pathName + "/hessA");
  // file.read(mHessB, pathName + "/hessB");
  // file.read(mHessC, pathName + "/hessC");
  // file.read(mHessD, pathName + "/hessD");
}

} // end namespace Spheral
