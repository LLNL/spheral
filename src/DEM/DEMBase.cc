//---------------------------------Spheral++----------------------------------//
// DEMBase -- The DEM package for Spheral++.
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/Physics.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/Timer.hh"

#include "DEM/DEMBase.hh"
#include "DEM/computeParticleRadius.hh"

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

// Declare timers
extern Timer TIME_DEM;
extern Timer TIME_DEMinitializeStartup;
extern Timer TIME_DEMregister;
extern Timer TIME_DEMregisterDerivs;
extern Timer TIME_DEMpreStepInitialize;
extern Timer TIME_DEMinitialize;
extern Timer TIME_DEMfinalizeDerivs;
extern Timer TIME_DEMghostBounds;
extern Timer TIME_DEMupdateVol;
extern Timer TIME_DEMenforceBounds;
extern Timer TIME_DEMevalDerivs;
extern Timer TIME_DEMevalDerivs_initial;
extern Timer TIME_DEMevalDerivs_pairs;
extern Timer TIME_DEMevalDerivs_final;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
DEMBase<Dimension>::
DEMBase(DataBase<Dimension>& dataBase,
        const TableKernel<Dimension>& W,
        const double cfl,
        const Vector& xmin,
        const Vector& xmax):
  Physics<Dimension>(),
  mKernel(W),
  mCfl(cfl),
  mxmin(xmin),
  mxmax(xmax),
  mTimeStepMask(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDomegaDt(FieldStorageType::CopyFields),
  mParticleRadius(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)){
    mTimeStepMask = dataBase.newFluidFieldList(int(0), "timeStepMask");
    mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::position);
    mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
    mDomegaDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Scalar>::prefix() + "angularVelocity");
    mParticleRadius = dataBase.newFluidFieldList( 0.0 , "particleRadius");
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DEMBase<Dimension>::
~DEMBase() {
}


//------------------------------------------------------------------------------
// Determine the timestep requirements for a hydro step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename DEMBase<Dimension>::TimeStepType
DEMBase<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   typename Dimension::Scalar /*currentTime*/) const {

    auto minDt = make_pair(1.0,("catPoop"));
    minDt.first*=this->mCfl;
    return minDt;
//   const auto& mask = state.fields(HydroFieldNames::timeStepMask, 1);
//   const auto& mass = state.fields(HydroFieldNames::mass, 0.0); 
//   const auto& position = state.fields(HydroFieldNames::position, Vector::zero);
//   const auto& velocity = state.fields(HydroFieldNames::velocity, Vector::zero);  
//   const auto& angularVelocity = state.fields("angularVelocity", Vector::zero);  

//   const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity(),
//                                                          this->requireOverlapConnectivity());

//   const float pi = 3.14159;
//   const float c1 = 1.0;
//   const float c2 = 1.0;

//   #pragma omp parallel
//   {

//     int i, j, nodeListi, nodeListj;

// #pragma omp for
//     for (auto kk = 0u; kk < npairs; ++kk) {

//       const auto mi = mass(nodeListi,i);
//       const auto mj = mass(nodeListj,j);
//       const auto mij = (mi*mj)/(mi+mj);

//       const auto stiffness = c1/mij;
//       const auto dissipation = c2/(2.0*mij);
//       const auto contactFrequency = std::sqrt(c1/mij - );
//       const auto contactTime = pi/contactFrequency;


//     }
//   }
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  TIME_DEMinitializeStartup.start();

  const auto H = dataBase.DEMHfield();
  auto particleRadius = this->particleRadius();
  computeParticleRadius(H,particleRadius);

  TIME_DEMinitializeStartup.stop();
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_DEMregister.start();
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  
  FieldList<Dimension, Vector> position = dataBase.DEMPosition();
  FieldList<Dimension, Vector> velocity = dataBase.DEMVelocity();
  FieldList<Dimension, Vector> angularVelocity = dataBase.DEMAngularVelocity();
  FieldList<Dimension, Scalar> mass = dataBase.DEMMass();
  FieldList<Dimension, SymTensor> Hfield = dataBase.DEMHfield();
  FieldList<Dimension, Scalar> particleRadius = this->particleRadius();

  PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
  PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,true));
  PolicyPointer angularVelocityPolicy(new IncrementFieldList<Dimension, Vector>());

  state.enroll(mass);
  state.enroll(Hfield);
  state.enroll(position, positionPolicy);
  state.enroll(velocity, velocityPolicy);
  state.enroll(angularVelocity, angularVelocityPolicy);
  state.enroll(mTimeStepMask);

  TIME_DEMregister.stop();
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_DEMregisterDerivs.start();

  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, false);
  dataBase.resizeFluidFieldList(mDomegaDt, Vector::zero, IncrementFieldList<Dimension, Scalar>::prefix() + "angularVelocity" , false);
  
  derivs.enroll(mDxDt);
  derivs.enroll(mDvDt);
  derivs.enroll(mDomegaDt);

  TIME_DEMregisterDerivs.stop();
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_DEMpreStepInitialize.start();


  TIME_DEMpreStepInitialize.stop();
}

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  TIME_DEMinitialize.start();


  // We depend on the caller knowing to finalize the ghost boundaries!
  TIME_DEMinitialize.stop();
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
  TIME_DEMevalDerivs.start();
//   TIME_DEMevalDerivs_initial.start();

//   // A few useful constants we'll use in the following loop.
//   const double tiny = 1.0e-30;
//   const auto c1 = 1.0;
//   const auto c2 = 1.0;

//   // The connectivity.
//   const auto& connectivityMap = dataBase.connectivityMap();
//   const auto& nodeLists = connectivityMap.nodeLists();
//   const auto numNodeLists = nodeLists.size();

//   // Get the state and derivative FieldLists.
//   // State FieldLists.
//   const auto mass = state.fields(HydroFieldNames::mass, 0.0);
//   const auto position = state.fields(HydroFieldNames::position, Vector::zero);
//   const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
//   const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
//   const auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
//   CHECK(mass.size() == numNodeLists);
//   CHECK(position.size() == numNodeLists);
//   CHECK(velocity.size() == numNodeLists);
//   CHECK(H.size() == numNodeLists);
//   CHECK(omega.size() == numNodeLists);

//   // Derivative FieldLists.
//   auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
//   auto  DvDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);
//   CHECK(DxDt.size() == numNodeLists);
//   CHECK(DvDt.size() == numNodeLists);

//   // The set of interacting node pairs.
//   const auto& pairs = connectivityMap.nodePairList();
//   const auto  npairs = pairs.size();

//   TIME_DEMevalDerivs_initial.stop();

//   // Walk all the interacting pairs.
//   TIME_DEMevalDerivs_pairs.start();
// #pragma omp parallel
//   {
//     // Thread private scratch variables
//     int i, j, nodeListi, nodeListj;

//     typename SpheralThreads<Dimension>::FieldListStack threadStack;
//     auto DvDt_thread = DvDt.threadCopy(threadStack);

// #pragma omp for
//     for (auto kk = 0u; kk < npairs; ++kk) {
//       i = pairs[kk].i_node;
//       j = pairs[kk].j_node;
//       nodeListi = pairs[kk].i_list;
//       nodeListj = pairs[kk].j_list;

//       // Get the state for node i.
//       const auto& ri = position(nodeListi, i);
//       const auto& mi = mass(nodeListi, i);
//       const auto& vi = velocity(nodeListi, i);
//       const auto& Hi = H(nodeListi, i);
//       const auto  Hdeti = Hi.Determinant();
      
//       auto& DvDti = DvDt_thread(nodeListi, i);

//       // Get the state for node j
//       const auto& rj = position(nodeListj, j);
//       const auto& mj = mass(nodeListj, j);
//       const auto& vj = velocity(nodeListj, j);
//       const auto& Hj = H(nodeListj, j);
//       const auto  Hdetj = Hj.Determinant();

//       auto& DvDtj = DvDt_thread(nodeListj, j);

//       CHECK(mi > 0.0);
//       CHECK(Hdeti > 0.0);
//       CHECK(mj > 0.0);
//       CHECK(Hdetj > 0.0);

//       const auto vij = vi-vj;
//       const auto rij = ri-rj;
//       const auto rhatij = rij.unitVector();
//       const auto rij2 = sqrt(rij.dot(rij));

//       const auto Ri = 1.0/Hdeti;
//       const auto Rj = 1.0/Hdetj;
//       const auto Rij2 = (Ri+Rj);

//       const auto delta = rij2-Rij2;  // negative will get ya a force

//       if (delta < 0.0){
//         const auto vn = vij.dot(rhatij);
//         const auto f = -(c1*delta - c2*vn);
//         DvDti += f/mi*rhatij;
//         DvDtj -= f/mj*rhatij;
//       }

//     } // loop over pairs

//     // Reduce the thread values to the master.
//     threadReduceFieldLists<Dimension>(threadStack);

//   }   // OpenMP parallel region
//   TIME_DEMevalDerivs_pairs.stop();

//   // Finish up the derivatives for each point.
//   TIME_DEMevalDerivs_final.start();
//   for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
//     const auto& nodeList = mass[nodeListi]->nodeList();
//     const auto  hmin = nodeList.hmin();
//     const auto  hmax = nodeList.hmax();
//     const auto  hminratio = nodeList.hminratio();
//     const auto  nPerh = nodeList.nodesPerSmoothingScale();

//     const auto ni = nodeList.numInternalNodes();
// #pragma omp parallel for
//     for (auto i = 0u; i < ni; ++i) {

//         //evaluate the node things

//     }
//   }
  //TIME_DEMevalDerivs_final.stop();
  TIME_DEMevalDerivs.stop();
}
//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_DEMfinalizeDerivs.start();

  TIME_DEMfinalizeDerivs.stop();
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {
  TIME_DEMghostBounds.start();
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Vector> angularVelocity = state.fields("angularVelocity", Vector::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(angularVelocity);
  }
  
  TIME_DEMghostBounds.stop();
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_DEMenforceBounds.start();
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Vector> angularVelocity = state.fields("angularVelocity", Vector::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(angularVelocity);
  }
 
  TIME_DEMenforceBounds.stop();
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  //file.write(mDxDt, pathName + "/DxDt");
  //file.write(mDvDt, pathName + "/DvDt");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  //file.read(mDxDt, pathName + "/DxDt");
  //file.read(mDvDt, pathName + "/DvDt");
}

}
