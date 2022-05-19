//---------------------------------Spheral++----------------------------------//
// DEMBase -- The DEM package for Spheral++.
//
//  NOT ADDING IF CONTACT IN ITITALIZED LIKE IN THE 2 PARTICLE TEST
//  WE're going to use allFields(dataType) to generalize our pairFielldList operations
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"

#include "Hydro/HydroFieldNames.hh"
#include "Hydro/PositionPolicy.hh"

#include "Physics/Physics.hh"

#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/DataBase.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"

#include "Boundary/Boundary.hh"

#include "Neighbor/ConnectivityMap.hh"

#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/Timer.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/registerWithRedistribution.hh"

#include "DEM/ContactStorageLocation.hh"
#include "DEM/IncrementPairFieldList.hh"
#include "DEM/ReplaceAndIncrementPairFieldList.hh"
#include "DEM/DEMBase.hh"
#include "DEM/DEMDimension.hh"
#include "DEM/DEMFieldNames.hh"

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
extern Timer TIME_DEMenforceBounds;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
DEMBase<Dimension>::
DEMBase(const DataBase<Dimension>& dataBase,
        const Scalar stepsPerCollision,
        const Vector& xmin,
        const Vector& xmax):
  Physics<Dimension>(),
  mDataBase(dataBase),
  mFirstCycle(true),
  mCyclesSinceLastKulling(0),
  mKullFrequency((int)stepsPerCollision),
  mStepsPerCollision(stepsPerCollision),
  mxmin(xmin),
  mxmax(xmax),
  mTimeStepMask(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mOmega(FieldStorageType::CopyFields),
  mDomegaDt(FieldStorageType::CopyFields),
  mUniqueIndices(FieldStorageType::CopyFields),
  mNeighborIndices(FieldStorageType::CopyFields),
  mEquilibriumOverlap(FieldStorageType::CopyFields),
  mShearDisplacement(FieldStorageType::CopyFields),
  mRollingDisplacement(FieldStorageType::CopyFields),
  mTorsionalDisplacement(FieldStorageType::CopyFields),
  mIsActiveContact(FieldStorageType::CopyFields),
  mDDtShearDisplacement(FieldStorageType::CopyFields),
  mNewShearDisplacement(FieldStorageType::CopyFields),
  mDDtRollingDisplacement(FieldStorageType::CopyFields),
  mNewRollingDisplacement(FieldStorageType::CopyFields),
  mDDtTorsionalDisplacement(FieldStorageType::CopyFields),
  mNewTorsionalDisplacement(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)),
  mRedistribute(registerWithRedistribution(*this,
                                           &DEMBase<Dimension>::initializeBeforeRedistribution,
                                           &DEMBase<Dimension>::finalizeAfterRedistribution)){
    
    mTimeStepMask = dataBase.newDEMFieldList(int(0), "timeStepMask");
    mDxDt = dataBase.newDEMFieldList(Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position);
    mDvDt = dataBase.newDEMFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
    mOmega = dataBase.newDEMFieldList(DEMDimension<Dimension>::zero, DEMFieldNames::angularVelocity);
    mDomegaDt = dataBase.newDEMFieldList(DEMDimension<Dimension>::zero, IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity);
    
    mUniqueIndices = dataBase.newDEMFieldList(int(0),  DEMFieldNames::uniqueIndices);
    mUniqueIndices += globalNodeIDs<Dimension>(dataBase.nodeListBegin(),dataBase.nodeListEnd()); 
    mIsActiveContact = dataBase.newDEMFieldList(std::vector<int>(), DEMFieldNames::isActiveContact);
    mNeighborIndices = dataBase.newDEMFieldList(std::vector<int>(), DEMFieldNames::neighborIndices);
    mShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), DEMFieldNames::shearDisplacement);
    mRollingDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), DEMFieldNames::rollingDisplacement);
    mTorsionalDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), DEMFieldNames::torsionalDisplacement);
    mDDtShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::shearDisplacement);
    mNewShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::shearDisplacement);
    mDDtRollingDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::rollingDisplacement);
    mNewRollingDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::rollingDisplacement);
    mDDtTorsionalDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::torsionalDisplacement);
    mNewTorsionalDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::torsionalDisplacement);
    mEquilibriumOverlap = dataBase.newDEMFieldList(std::vector<Scalar>(), DEMFieldNames::equilibriumOverlap);

}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
DEMBase<Dimension>::
~DEMBase() {
}

//------------------------------------------------------------------------------
// update numContacts for pair derivatives field lists
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
resizeDerivativePairFieldLists(StateDerivatives<Dimension>& derivs) const {
  
  auto DsDt = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() + DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto newShearDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() + DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto DDtRollingDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() + DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto newRollingDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() + DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto DDtTorsionalDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() + DEMFieldNames::torsionalDisplacement, std::vector<Vector>());
  auto newTorsionalDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() + DEMFieldNames::torsionalDisplacement, std::vector<Vector>());
  
  this->addContactsToPairFieldList(DsDt,Vector::zero);
  this->addContactsToPairFieldList(newShearDisp,Vector::zero);
  this->addContactsToPairFieldList(DDtRollingDisp,Vector::zero);
  this->addContactsToPairFieldList(newRollingDisp,Vector::zero);
  this->addContactsToPairFieldList(DDtTorsionalDisp,Vector::zero);
  this->addContactsToPairFieldList(newTorsionalDisp,Vector::zero);
}

//------------------------------------------------------------------------------
// update numContacts for pair state field lists
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
resizeStatePairFieldLists(State<Dimension>& state) const{

  auto eqOverlap = state.fields(DEMFieldNames::equilibriumOverlap, vector<Scalar>());
  auto shearDisp = state.fields(DEMFieldNames::shearDisplacement, vector<Vector>());
  auto rollingDisplacement = state.fields(DEMFieldNames::rollingDisplacement, vector<Vector>());
  auto torsionalDisplacement = state.fields(DEMFieldNames::torsionalDisplacement, vector<Vector>());

  this->addContactsToPairFieldList(eqOverlap,0.0);
  this->addContactsToPairFieldList(shearDisp,Vector::zero);
  this->addContactsToPairFieldList(rollingDisplacement,Vector::zero);
  this->addContactsToPairFieldList(torsionalDisplacement,Vector::zero);
}

//------------------------------------------------------------------------------
// remove old contacts
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
kullInactiveContactsFromStatePairFieldLists(State<Dimension>& state) const{

  auto eqOverlap = state.fields(DEMFieldNames::equilibriumOverlap, vector<Scalar>());
  auto shearDisp = state.fields(DEMFieldNames::shearDisplacement, vector<Vector>());
  auto rollingDisplacement = state.fields(DEMFieldNames::rollingDisplacement, vector<Vector>());
  auto torsionalDisplacement = state.fields(DEMFieldNames::torsionalDisplacement, vector<Vector>());

  this->kullInactiveContactsFromPairFieldList(eqOverlap);
  this->kullInactiveContactsFromPairFieldList(shearDisp);
  this->kullInactiveContactsFromPairFieldList(rollingDisplacement);
  this->kullInactiveContactsFromPairFieldList(torsionalDisplacement);

}


//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  TIME_DEMinitializeStartup.start();

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

  dataBase.resizeDEMFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeDEMFieldList(mOmega, DEMDimension<Dimension>::zero, DEMFieldNames::angularVelocity, false);
  dataBase.resizeDEMFieldList(mUniqueIndices, 0, DEMFieldNames::uniqueIndices, false);
  dataBase.resizeDEMFieldList(mIsActiveContact, vector<int>(), DEMFieldNames::isActiveContact, false);
  dataBase.resizeDEMFieldList(mNeighborIndices, vector<int>(), DEMFieldNames::neighborIndices, false);
  dataBase.resizeDEMFieldList(mShearDisplacement, vector<Vector>(), DEMFieldNames::shearDisplacement, false);
  dataBase.resizeDEMFieldList(mRollingDisplacement, vector<Vector>(), DEMFieldNames::rollingDisplacement, false);
  dataBase.resizeDEMFieldList(mTorsionalDisplacement, vector<Vector>(), DEMFieldNames::torsionalDisplacement, false);
  dataBase.resizeDEMFieldList(mEquilibriumOverlap, vector<Scalar>(), DEMFieldNames::equilibriumOverlap, false);

  auto position = dataBase.DEMPosition();
  auto velocity = dataBase.DEMVelocity();
  auto mass = dataBase.DEMMass();
  auto Hfield = dataBase.DEMHfield();
  auto radius = dataBase.DEMParticleRadius();
  auto compositeParticleIndex = dataBase.DEMCompositeParticleIndex();

  PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
  PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,true));
  PolicyPointer angularVelocityPolicy(new IncrementFieldList<Dimension, RotationType>());
  PolicyPointer shearDisplacementPolicy(new ReplaceAndIncrementPairFieldList<Dimension,std::vector<Vector>>());
  PolicyPointer rollingDisplacementPolicy(new ReplaceAndIncrementPairFieldList<Dimension,std::vector<Vector>>());
  PolicyPointer torsionalDisplacementPolicy(new ReplaceAndIncrementPairFieldList<Dimension,std::vector<Vector>>());

  state.enroll(mTimeStepMask);
  state.enroll(mass);
  state.enroll(radius);
  state.enroll(Hfield);
  state.enroll(compositeParticleIndex);
  state.enroll(position, positionPolicy);
  state.enroll(velocity, velocityPolicy);
  state.enroll(mOmega, angularVelocityPolicy);
  

  state.enroll(mIsActiveContact);
  state.enroll(mUniqueIndices);
  state.enroll(mNeighborIndices);
  state.enroll(mShearDisplacement, shearDisplacementPolicy);
  state.enroll(mRollingDisplacement, rollingDisplacementPolicy);
  state.enroll(mTorsionalDisplacement, torsionalDisplacementPolicy);
  state.enroll(mEquilibriumOverlap);

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

  dataBase.resizeDEMFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::position, false);
  dataBase.resizeDEMFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeDEMFieldList(mDomegaDt, DEMDimension<Dimension>::zero,  IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity , false);
  dataBase.resizeDEMFieldList(mDDtShearDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::shearDisplacement , false);
  dataBase.resizeDEMFieldList(mNewShearDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::shearDisplacement , false);
  dataBase.resizeDEMFieldList(mDDtRollingDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::rollingDisplacement , false);
  dataBase.resizeDEMFieldList(mNewRollingDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::rollingDisplacement , false);
  dataBase.resizeDEMFieldList(mDDtTorsionalDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::torsionalDisplacement , false);
  dataBase.resizeDEMFieldList(mNewTorsionalDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::torsionalDisplacement , false);
  
  derivs.enroll(mDxDt);
  derivs.enroll(mDvDt);
  derivs.enroll(mDomegaDt);
  derivs.enroll(mDDtShearDisplacement);
  derivs.enroll(mNewShearDisplacement);
  derivs.enroll(mDDtRollingDisplacement);
  derivs.enroll(mNewRollingDisplacement);
  derivs.enroll(mDDtTorsionalDisplacement);
  derivs.enroll(mNewTorsionalDisplacement);

  TIME_DEMregisterDerivs.stop();
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& /*dataBase*/, 
                  State<Dimension>& /*state*/,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_DEMpreStepInitialize.start();
  
  TIME_DEMpreStepInitialize.stop();
}

//------------------------------------------------------------------------------
// Call before deriv evaluation
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initialize(const Scalar  time,
           const Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs){

  // update our state pair fields for the current connectivity
  this->initializeContacts(dataBase);
  this->resizeStatePairFieldLists(state);

  // on the first pass we initialize our equilibrium overlap
  if(mFirstCycle){
    this->initializeOverlap(dataBase);
    mFirstCycle=false;
  }

  // every so often remove old contacts 
  mCyclesSinceLastKulling++;
  if (mCyclesSinceLastKulling % mKullFrequency == 0){ 
    this->kullInactiveContacts(dataBase);
    this->kullInactiveContactsFromStatePairFieldLists(state);
  }

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
                    StateDerivatives<Dimension>& /*derivs*/) const {
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
                     StateDerivatives<Dimension>& derivs) {
  TIME_DEMghostBounds.start();

  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto angularVelocity = state.fields(DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  auto radius = state.fields(DEMFieldNames::particleRadius, 0.0);
  auto compositeParticleIndex = state.fields(DEMFieldNames::compositeParticleIndex,int(0));
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mUniqueIndices);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(angularVelocity);
    (*boundaryItr)->applyFieldListGhostBoundary(radius);
    (*boundaryItr)->applyFieldListGhostBoundary(compositeParticleIndex);
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
                  StateDerivatives<Dimension>& derivs) {
  TIME_DEMenforceBounds.start();

  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto angularVelocity = state.fields(DEMFieldNames::angularVelocity, DEMDimension<Dimension>::zero);
  auto radius = state.fields(DEMFieldNames::particleRadius, 0.0);
  auto compositeParticleIndex = state.fields(DEMFieldNames::compositeParticleIndex,int(0));
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mUniqueIndices);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(angularVelocity);
    (*boundaryItr)->enforceFieldListBoundary(radius);
    (*boundaryItr)->enforceFieldListBoundary(compositeParticleIndex);
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

  file.write(mCyclesSinceLastKulling, pathName + "/cyclesSinceLastKulling");

  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mOmega, pathName + "/omega");
  file.write(mDomegaDt, pathName + "/DomegaDt");
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mUniqueIndices, pathName + "/uniqueIndices");

  file.write(mIsActiveContact, pathName + "/iActiveContact");
  file.write(mNeighborIndices, pathName + "/neighborIndices");
  file.write(mShearDisplacement, pathName + "/shearDisplacement");
  file.write(mRollingDisplacement, pathName + "/rollingDisplacement");
  file.write(mTorsionalDisplacement, pathName + "/torsionalDisplacement");
  
  file.write(mDDtShearDisplacement, pathName + "/DDtShearDisplacement");
  file.write(mNewShearDisplacement, pathName + "/newShearDisplacement");
  file.write(mDDtRollingDisplacement, pathName + "/DDtRollingDisplacement");
  file.write(mNewRollingDisplacement, pathName + "/newRollingDisplacement");
  file.write(mDDtTorsionalDisplacement, pathName + "/DDtTorsionalDisplacement");
  file.write(mNewTorsionalDisplacement, pathName + "/newTorsionalDisplacement");
  file.write(mEquilibriumOverlap, pathName + "/equilibriumOverlap");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  file.read(mCyclesSinceLastKulling, pathName + "/cyclesSinceLastKulling");

  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mOmega, pathName + "/omega");
  file.read(mDomegaDt, pathName + "/DomegaDt");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mUniqueIndices, pathName + "/uniqueIndices");

  file.read(mIsActiveContact, pathName + "/iActiveContact");
  file.read(mNeighborIndices, pathName + "/neighborIndices");
  file.read(mShearDisplacement, pathName + "/shearDisplacement");
  file.read(mRollingDisplacement, pathName + "/rollingDisplacement");
  file.read(mTorsionalDisplacement, pathName + "/torsionalDisplacement");

  file.read(mDDtShearDisplacement, pathName + "/DDtShearDisplacement");
  file.read(mNewShearDisplacement, pathName + "/newShearDisplacement");
  file.read(mDDtRollingDisplacement, pathName + "/DDtRollingDisplacement");
  file.read(mNewRollingDisplacement, pathName + "/newRollingDisplacement");
  file.read(mDDtTorsionalDisplacement, pathName + "/DDtTorsionalDisplacement");
  file.read(mNewTorsionalDisplacement, pathName + "/newTorsionalDisplacement");
  file.read(mEquilibriumOverlap, pathName + "/equilibriumOverlap");
}



//------------------------------------------------------------------------------
// when problem starts set our equilibrium overlap
//------------------------------------------------------------------------------

template<typename Dimension>
void
DEMBase<Dimension>::
initializeOverlap(const DataBase<Dimension>& dataBase){
  
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  numPairs = pairs.size();
  
  const auto& positions = dataBase.DEMPosition();
  const auto& particleRadius = dataBase.DEMParticleRadius();
  const auto& particleIndex = dataBase.DEMCompositeParticleIndex();

#pragma omp for
  for (auto kk = 0u; kk < numPairs; ++kk) {

    const auto i = pairs[kk].i_node;
    const auto j = pairs[kk].j_node;
    const auto nodeListi = pairs[kk].i_list;
    const auto nodeListj = pairs[kk].j_list;

    const auto storeNodeList = mContactStorageIndices[kk].storeNodeList;
    const auto storeNode = mContactStorageIndices[kk].storeNode;
    const auto storeContact = mContactStorageIndices[kk].storeContact;

    // node fields
    const auto ri = positions(nodeListi,i);
    const auto Ri = particleRadius(nodeListi,i);
    const auto pIDi = particleIndex(nodeListi,i);

    const auto rj = positions(nodeListj,j);
    const auto Rj = particleRadius(nodeListj,j);
    const auto pIDj = particleIndex(nodeListj,j);

    // composite particle 
    if(pIDi==pIDj){
      const auto rij = ri-rj;
      const auto rijMag = rij.magnitude();
      const auto delta0 = (Ri+Rj) - rijMag;
      mEquilibriumOverlap(storeNodeList,storeNode)[storeContact] = delta0;
    }
  }
}

//------------------------------------------------------------------------------
// REdistribution methods
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initializeBeforeRedistribution(){

  const auto& connectivityMap = mDataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // if the pair is internal copy stored pairwise values over
#pragma omp for
  for (auto kk = 0u; kk < npairs; ++kk) {
    const auto& contactkk = mContactStorageIndices[kk];
    const auto storageUniqueNode = mUniqueIndices(contactkk.storeNodeList,contactkk.storeNode);
    if (contactkk.pairNode < (int)mUniqueIndices[contactkk.pairNodeList]->numInternalElements()){
      mNeighborIndices(contactkk.pairNodeList,contactkk.pairNode).push_back(storageUniqueNode);
      mShearDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mShearDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mRollingDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mRollingDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mTorsionalDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mTorsionalDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mEquilibriumOverlap(contactkk.pairNodeList,contactkk.pairNode).push_back(mEquilibriumOverlap(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mDDtShearDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mDDtShearDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mNewShearDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mNewShearDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mDDtRollingDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mDDtRollingDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mNewRollingDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mNewRollingDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mDDtTorsionalDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mDDtTorsionalDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
      mNewTorsionalDisplacement(contactkk.pairNodeList,contactkk.pairNode).push_back(mNewTorsionalDisplacement(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
    } // if 
  }   // for
}     // method

template<typename Dimension>
void
DEMBase<Dimension>::
finalizeAfterRedistribution() {
  //std::cout << "SHABINGO!"<< std::endl;
}

//------------------------------------------------------------------------------
// set the size of our pairFieldLists based on mNeighborIndices 
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value1, typename Value2>
void
DEMBase<Dimension>::
addContactsToPairFieldList(      Value1& pairFieldList,
                           const Value2& newValue) const {

  const auto  numNodeLists = pairFieldList.numFields();
  const auto  nodeListPtrs = pairFieldList.nodeListPtrs();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = nodeListPtrs[nodeListi]->numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
      const auto numContacts = mNeighborIndices(nodeListi,i).size();
      pairFieldList(nodeListi,i).resize(numContacts,newValue);
    }
  }
}


//------------------------------------------------------------------------------
// cull dead contacts once isActive is initialized
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
void
DEMBase<Dimension>::
kullInactiveContactsFromPairFieldList(Value& pairFieldList) const {

  const auto  numNodeLists = pairFieldList.numFields();
  const auto  nodeListPtrs = pairFieldList.nodeListPtrs();                                     

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = nodeListPtrs[nodeListi]->numInternalNodes();
#pragma omp parallel for
    for (auto nodei = 0u; nodei < numNodes; ++nodei) {
      
      const auto& isActivei = mIsActiveContact(nodeListi,nodei);
      auto& pairFieldi = pairFieldList(nodeListi,nodei);

      const auto numContacts = pairFieldi.size();
      unsigned int activeContactCount = 0;

      // shift all active entries to begining of vector
      for (auto contacti = 0u; contacti < numContacts; ++contacti){
        pairFieldi[activeContactCount] = pairFieldi[contacti];
        if (isActivei[contacti]==1) activeContactCount++;
      }

      // resize to clip old entries
      pairFieldi.resize(activeContactCount);
    }
  }
}

//------------------------------------------------------------------------------
// get all of our state pair-fields up to date with the number of contacts
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initializeContacts(const DataBase<Dimension>& dataBase){

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  numPairs = pairs.size();
  
  mContactStorageIndices.resize(numPairs);

#pragma omp for
  for (auto kk = 0u; kk < numPairs; ++kk) {

    const auto i = pairs[kk].i_node;
    const auto j = pairs[kk].j_node;
    const auto nodeListi = pairs[kk].i_list;
    const auto nodeListj = pairs[kk].j_list;

    // get our unique global ID numbers
    const auto uIDi = mUniqueIndices(nodeListi,i);
    const auto uIDj = mUniqueIndices(nodeListj,j);
  
    // get our number of internal nodes
    const int numInternalNodesi = mUniqueIndices[nodeListi]->numInternalElements();
    const int numInternalNodesj = mUniqueIndices[nodeListj]->numInternalElements();

    // boolean operations to decide which pair-node maintains pair fields
    const auto nodeiIsInternal = i < numInternalNodesi;
    const auto nodejIsGhost = j >= numInternalNodesj;
    const auto nodeiIsMinUnique = uIDi <= uIDj;
    const auto selectNodei = nodeiIsInternal and (nodeiIsMinUnique or nodejIsGhost);

    // store info on storage/dummy locations
    auto& contactkk = mContactStorageIndices[kk];

    contactkk.storeNodeList = (selectNodei ? nodeListi : nodeListj);
    contactkk.storeNode     = (selectNodei ? i         : j);
    contactkk.pairNodeList  = (selectNodei ? nodeListj : nodeListi);
    contactkk.pairNode      = (selectNodei ? j         : i);

    // find contact if it already exists 
    const auto uniqueSearchIndex = (selectNodei ? uIDj : uIDi);
    const auto neighborContacts = mNeighborIndices(contactkk.storeNodeList,
                                                   contactkk.storeNode);
    const auto contactIndexPtr = std::find(neighborContacts.begin(),
                                           neighborContacts.end(),
                                           uniqueSearchIndex);
    const int  storageContactIndex = std::distance(neighborContacts.begin(),contactIndexPtr);
    contactkk.storeContact = storageContactIndex;

    // add the contact if it is new
    if (contactIndexPtr == neighborContacts.end()){
      mNeighborIndices(contactkk.storeNodeList,contactkk.storeNode).push_back(uniqueSearchIndex);
    }
  }
}


//------------------------------------------------------------------------------
// remove old contacts from our pair storage every so often
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
kullInactiveContacts(const DataBase<Dimension>& dataBase){

  const auto  numNodeLists = dataBase.numNodeLists();
  const auto  nodeListPtrs = dataBase.DEMNodeListPtrs();
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // initialize our tracker for contact activity as a vector of zeros
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = nodeListPtrs[nodeListi]->numInternalNodes();
#pragma omp parallel for
    for (auto nodei = 0u; nodei < numNodes; ++nodei) {  
      const auto numContacts = mNeighborIndices(nodeListi,nodei).size();
      mIsActiveContact(nodeListi,nodei).assign(numContacts,0);
    }
  }

  // set all our active contacts to 1
#pragma omp for
  for (auto kk = 0u; kk < npairs; ++kk) {
    const auto& contactkk = mContactStorageIndices[kk];
    mIsActiveContact(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact] = 1;
  }

  this->kullInactiveContactsFromPairFieldList(mNeighborIndices);
}


} //spheral namespace
