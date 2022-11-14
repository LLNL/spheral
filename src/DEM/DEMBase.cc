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
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/registerWithRedistribution.hh"

#include "DEM/ContactStorageLocation.hh"
#include "DEM/IncrementPairFieldList.hh"
#include "DEM/ReplaceAndIncrementPairFieldList.hh"
#include "DEM/DEMBase.hh"
#include "DEM/DEMDimension.hh"
#include "DEM/DEMFieldNames.hh"

#include "Utilities/Timer.hh"

#ifdef _OPENMP
#include "omp.h"
#endif


#include <iostream>
#include <stdexcept>
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
DEMBase<Dimension>::
DEMBase(const DataBase<Dimension>& dataBase,
        const Scalar stepsPerCollision,
        const Vector& xmin,
        const Vector& xmax):
  Physics<Dimension>(),
  mDataBase(dataBase),
  mFirstCycle(true),
  mCycle(0),
  mContactRemovalFrequency((int)stepsPerCollision),
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
    mTorsionalDisplacement = dataBase.newDEMFieldList(std::vector<Scalar>(), DEMFieldNames::torsionalDisplacement);
    mDDtShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::shearDisplacement);
    mNewShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::shearDisplacement);
    mDDtRollingDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::rollingDisplacement);
    mNewRollingDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::rollingDisplacement);
    mDDtTorsionalDisplacement = dataBase.newDEMFieldList(std::vector<Scalar>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::torsionalDisplacement);
    mNewTorsionalDisplacement = dataBase.newDEMFieldList(std::vector<Scalar>(), ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::torsionalDisplacement);
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
  auto DDtTorsionalDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::incrementPrefix() + DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  auto newTorsionalDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::replacePrefix() + DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  
  this->addContactsToPairFieldList(DsDt,Vector::zero);
  this->addContactsToPairFieldList(newShearDisp,Vector::zero);
  this->addContactsToPairFieldList(DDtRollingDisp,Vector::zero);
  this->addContactsToPairFieldList(newRollingDisp,Vector::zero);
  this->addContactsToPairFieldList(DDtTorsionalDisp,0.0);
  this->addContactsToPairFieldList(newTorsionalDisp,0.0);
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
  auto torsionalDisplacement = state.fields(DEMFieldNames::torsionalDisplacement, vector<Scalar>());

  this->addContactsToPairFieldList(eqOverlap,0.0);
  this->addContactsToPairFieldList(shearDisp,Vector::zero);
  this->addContactsToPairFieldList(rollingDisplacement,Vector::zero);
  this->addContactsToPairFieldList(torsionalDisplacement,0.0);
}

//------------------------------------------------------------------------------
// remove old contacts
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
removeInactiveContactsFromStatePairFieldLists(State<Dimension>& state) const{

  auto neighborIndices = state.fields(DEMFieldNames::neighborIndices,vector<int>());
  auto eqOverlap = state.fields(DEMFieldNames::equilibriumOverlap, vector<Scalar>());
  auto shearDisp = state.fields(DEMFieldNames::shearDisplacement, vector<Vector>());
  auto rollingDisplacement = state.fields(DEMFieldNames::rollingDisplacement, vector<Vector>());
  auto torsionalDisplacement = state.fields(DEMFieldNames::torsionalDisplacement, vector<Scalar>());

  this->removeInactiveContactsFromPairFieldList(neighborIndices);
  this->removeInactiveContactsFromPairFieldList(eqOverlap);
  this->removeInactiveContactsFromPairFieldList(shearDisp);
  this->removeInactiveContactsFromPairFieldList(rollingDisplacement);
  this->removeInactiveContactsFromPairFieldList(torsionalDisplacement);

}

//------------------------------------------------------------------------------
// remove old contacts
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
removeInactiveContactsFromDerivativePairFieldLists(StateDerivatives<Dimension>& derivs) const {
  
  auto DsDt = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() + DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto newShearDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() + DEMFieldNames::shearDisplacement, std::vector<Vector>());
  auto DDtRollingDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix() + DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto newRollingDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix() + DEMFieldNames::rollingDisplacement, std::vector<Vector>());
  auto DDtTorsionalDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::incrementPrefix() + DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());
  auto newTorsionalDisp = derivs.fields(ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::replacePrefix() + DEMFieldNames::torsionalDisplacement, std::vector<Scalar>());

  this->removeInactiveContactsFromPairFieldList(DsDt);
  this->removeInactiveContactsFromPairFieldList(newShearDisp);
  this->removeInactiveContactsFromPairFieldList(DDtRollingDisp);
  this->removeInactiveContactsFromPairFieldList(newRollingDisp);
  this->removeInactiveContactsFromPairFieldList(DDtTorsionalDisp);
  this->removeInactiveContactsFromPairFieldList(newTorsionalDisp);

}


//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  TIME_BEGIN("DEMinitializeProblemStartup");
  
  TIME_END("DEMinitializeProblemStartup");
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_BEGIN("DEMregister");
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  dataBase.resizeDEMFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeDEMFieldList(mOmega, DEMDimension<Dimension>::zero, DEMFieldNames::angularVelocity, false);
  dataBase.resizeDEMFieldList(mUniqueIndices, 0, DEMFieldNames::uniqueIndices, false);
  dataBase.resizeDEMFieldList(mIsActiveContact, vector<int>(), DEMFieldNames::isActiveContact, false);
  dataBase.resizeDEMFieldList(mNeighborIndices, vector<int>(), DEMFieldNames::neighborIndices, false);
  dataBase.resizeDEMFieldList(mShearDisplacement, vector<Vector>(), DEMFieldNames::shearDisplacement, false);
  dataBase.resizeDEMFieldList(mRollingDisplacement, vector<Vector>(), DEMFieldNames::rollingDisplacement, false);
  dataBase.resizeDEMFieldList(mTorsionalDisplacement, vector<Scalar>(), DEMFieldNames::torsionalDisplacement, false);
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
  PolicyPointer torsionalDisplacementPolicy(new ReplaceAndIncrementPairFieldList<Dimension,std::vector<Scalar>>());

  state.enroll(mTimeStepMask);
  state.enroll(mass);
  state.enroll(radius);
  state.enroll(Hfield);
  state.enroll(compositeParticleIndex);
  state.enroll(mUniqueIndices);
  state.enroll(position, positionPolicy);
  state.enroll(velocity, velocityPolicy);
  state.enroll(mOmega, angularVelocityPolicy);

  state.enroll(mIsActiveContact);
  state.enroll(mNeighborIndices);
  state.enroll(mEquilibriumOverlap);
  state.enroll(mShearDisplacement,     shearDisplacementPolicy);
  state.enroll(mRollingDisplacement,   rollingDisplacementPolicy);
  state.enroll(mTorsionalDisplacement, torsionalDisplacementPolicy);

  TIME_END("DEMregister");
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("DEMregisterDerivs");

  dataBase.resizeDEMFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::position, false);
  dataBase.resizeDEMFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeDEMFieldList(mDomegaDt, DEMDimension<Dimension>::zero,  IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::angularVelocity , false);
  dataBase.resizeDEMFieldList(mDDtShearDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::shearDisplacement , false);
  dataBase.resizeDEMFieldList(mNewShearDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::shearDisplacement , false);
  dataBase.resizeDEMFieldList(mDDtRollingDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::incrementPrefix()  + DEMFieldNames::rollingDisplacement , false);
  dataBase.resizeDEMFieldList(mNewRollingDisplacement, vector<Vector>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Vector>>::replacePrefix()  + DEMFieldNames::rollingDisplacement , false);
  dataBase.resizeDEMFieldList(mDDtTorsionalDisplacement, vector<Scalar>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::incrementPrefix()  + DEMFieldNames::torsionalDisplacement , false);
  dataBase.resizeDEMFieldList(mNewTorsionalDisplacement, vector<Scalar>(),  ReplaceAndIncrementPairFieldList<Dimension, std::vector<Scalar>>::replacePrefix()  + DEMFieldNames::torsionalDisplacement , false);
  
  derivs.enroll(mDxDt);
  derivs.enroll(mDvDt);
  derivs.enroll(mDomegaDt);
  derivs.enroll(mDDtShearDisplacement);
  derivs.enroll(mNewShearDisplacement);
  derivs.enroll(mDDtRollingDisplacement);
  derivs.enroll(mNewRollingDisplacement);
  derivs.enroll(mDDtTorsionalDisplacement);
  derivs.enroll(mNewTorsionalDisplacement);

  TIME_END("DEMregisterDerivs");
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivatives) {
  TIME_BEGIN("DEMpreStepInitialize");

  // update our state pair fields for the current connectivity
  this->updateContactMapAndNeighborIndices(dataBase);  // set our contactMap and neighborIndices pairFieldList
  this->resizeStatePairFieldLists(state);              // add entries for new contacts
  this->resizeDerivativePairFieldLists(derivatives);   // do same for derivs in case we're storing from last cycle
  
  if (mCycle % mContactRemovalFrequency == 0){ 

    // remove old contacts
    this->identifyInactiveContacts(dataBase);                              // create pairFieldList tracking active/inactive contacts
    this->removeInactiveContactsFromStatePairFieldLists(state);            // use it to remove old contacts from state fields
    this->removeInactiveContactsFromDerivativePairFieldLists(derivatives); // use it to remove old contacts from the derivatives

    this->updateContactMap(dataBase);

  }

  mCycle++;

  TIME_END("DEMpreStepInitialize");
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

  // This is a little ugly but we need the connectivity to be constructed
  // and the ghost nodes to be set to initialize our equilibrium overlap.
  if(mFirstCycle){
    this->initializeOverlap(dataBase);
    mFirstCycle=false;
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
  TIME_BEGIN("DEMfinalizeDerivs");

  TIME_END("DEMfinalizeDerivs");
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("DEMghostBounds");

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
  TIME_END("DEMghostBounds");
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("DEMenforceBounds");

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
  TIME_END("DEMenforceBounds");
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  file.write(mCycle, pathName + "/cycle");
  file.write(mFirstCycle, pathName + "/firstCycle");

  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mOmega, pathName + "/omega");
  file.write(mDomegaDt, pathName + "/DomegaDt");
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mUniqueIndices, pathName + "/uniqueIndices");

  file.write(mIsActiveContact, pathName + "/isActiveContact");
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

  file.read(mCycle, pathName + "/cycle");
  file.read(mFirstCycle, pathName + "/firstCycle");
  
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mOmega, pathName + "/omega");
  file.read(mDomegaDt, pathName + "/DomegaDt");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mUniqueIndices, pathName + "/uniqueIndices");

  file.read(mIsActiveContact, pathName + "/isActiveContact");
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
// Redistribution methods -- before we redistribute, we are going to make sure
// that each node in the pairwise interactions agrees regarding the stored
// pairwise variables. Note, if the kullFrequency -> infty then you will get 
// duplicate entries. It shouldn't break things but will lead to bloat pretty
// quickly if you're doing a lot of redistributing.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initializeBeforeRedistribution(){
  this->prepNeighborIndicesForRedistribution();

  this->prepPairFieldListForRedistribution(mShearDisplacement);
  this->prepPairFieldListForRedistribution(mRollingDisplacement);
  this->prepPairFieldListForRedistribution(mTorsionalDisplacement);
  this->prepPairFieldListForRedistribution(mEquilibriumOverlap);

  this->prepPairFieldListForRedistribution(mDDtShearDisplacement);
  this->prepPairFieldListForRedistribution(mNewShearDisplacement);
  this->prepPairFieldListForRedistribution(mDDtRollingDisplacement);
  this->prepPairFieldListForRedistribution(mNewRollingDisplacement);
  this->prepPairFieldListForRedistribution(mDDtTorsionalDisplacement);
  this->prepPairFieldListForRedistribution(mNewTorsionalDisplacement);
}

template<typename Dimension>
void
DEMBase<Dimension>::
finalizeAfterRedistribution() {
}

//------------------------------------------------------------------------------
// redistribution sub function
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
void
DEMBase<Dimension>::
prepPairFieldListForRedistribution(Value& pairFieldList) {
  const auto  nContacts = mContactStorageIndices.size();
#pragma omp for
  for (auto kk = 0u; kk < nContacts; ++kk) {
    const auto& contactkk = mContactStorageIndices[kk];
    const int numInternalNodes = pairFieldList[contactkk.pairNodeList]->numInternalElements();
    if (contactkk.pairNode < numInternalNodes){
      pairFieldList(contactkk.pairNodeList,contactkk.pairNode).push_back(pairFieldList(contactkk.storeNodeList,contactkk.storeNode)[contactkk.storeContact]);
    } 
  }   
}

template<typename Dimension>
void
DEMBase<Dimension>::
prepNeighborIndicesForRedistribution() {

  const auto  nContacts = mContactStorageIndices.size();
#pragma omp for
  for (auto kk = 0u; kk < nContacts; ++kk) {
    const auto& contactkk = mContactStorageIndices[kk];
    const int numInternalNodes = mUniqueIndices[contactkk.pairNodeList]->numInternalElements();
    const auto storageUniqueNode = mUniqueIndices(contactkk.storeNodeList,contactkk.storeNode);
    if (contactkk.pairNode < numInternalNodes){
      mNeighborIndices(contactkk.pairNodeList,contactkk.pairNode).push_back(storageUniqueNode);
    } 
  }   
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
removeInactiveContactsFromPairFieldList(Value& pairFieldList) const {

  const auto  numNodeLists = pairFieldList.numFields();
  const auto  nodeListPtrs = pairFieldList.nodeListPtrs();                                     

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = nodeListPtrs[nodeListi]->numInternalNodes();
#pragma omp parallel for
    for (auto nodei = 0u; nodei < numNodes; ++nodei) {
      
      const auto& isActivei = mIsActiveContact(nodeListi,nodei);
      auto& pairFieldi = pairFieldList(nodeListi,nodei);

      if(not(isActivei.size()==pairFieldi.size())) throw std::invalid_argument("wrong sizes");

      const auto numContacts = isActivei.size();
      unsigned int activeContactCount = 0;

      // shift all active entries to begining of vector
      for (auto contacti = 0u; contacti < numContacts; ++contacti){
        pairFieldi[activeContactCount] = pairFieldi[contacti];
        if (isActivei[contacti]==1) ++activeContactCount;
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
updateContactMapAndNeighborIndices(const DataBase<Dimension>& dataBase){

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
// reset contacts
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
updateContactMap(const DataBase<Dimension>& dataBase){

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

  }
}

//------------------------------------------------------------------------------
// remove old contacts from our pair storage every so often
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
identifyInactiveContacts(const DataBase<Dimension>& dataBase){

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
}


} //spheral namespace
