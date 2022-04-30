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

#include "DEM/IncrementPairFieldList.hh"
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
  mIsActiveContact(FieldStorageType::CopyFields),
  mDDtShearDisplacement(FieldStorageType::CopyFields),
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
    mDDtShearDisplacement = dataBase.newDEMFieldList(std::vector<Vector>(), IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::shearDisplacement);
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
  
  auto DsDt = derivs.fields(IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::shearDisplacement, std::vector<Vector>());
  
  this->addContactsToPairFieldList(DsDt,Vector::zero);

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

  this->addContactsToPairFieldList(eqOverlap,0.0);
  this->addContactsToPairFieldList(shearDisp,Vector::zero);

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

  this->kullInactiveContactsFromPairFieldList(eqOverlap);
  this->kullInactiveContactsFromPairFieldList(shearDisp);

}


//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& /*dataBase*/) {

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
  dataBase.resizeDEMFieldList(mEquilibriumOverlap, vector<Scalar>(), DEMFieldNames::uniqueIndices, false);

  auto position = dataBase.DEMPosition();
  auto velocity = dataBase.DEMVelocity();
  auto mass = dataBase.DEMMass();
  auto Hfield = dataBase.DEMHfield();
  auto radius = dataBase.DEMParticleRadius();

  PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
  PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,true));
  PolicyPointer angularVelocityPolicy(new IncrementFieldList<Dimension, RotationType>());
  PolicyPointer shearDisplacementPolicy(new IncrementPairFieldList<Dimension, std::vector<Vector>>());

  state.enroll(mTimeStepMask);
  state.enroll(mass);
  state.enroll(radius);
  state.enroll(Hfield);
  state.enroll(position, positionPolicy);
  state.enroll(velocity, velocityPolicy);
  state.enroll(mOmega, angularVelocityPolicy);
  
  state.enroll(mIsActiveContact);
  state.enroll(mUniqueIndices);
  state.enroll(mNeighborIndices);
  state.enroll(mShearDisplacement, shearDisplacementPolicy);
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
  dataBase.resizeDEMFieldList(mDDtShearDisplacement, vector<Vector>(),  IncrementFieldList<Dimension, Scalar>::prefix() + DEMFieldNames::shearDisplacement , false);
  
  derivs.enroll(mDxDt);
  derivs.enroll(mDvDt);
  derivs.enroll(mDomegaDt);
  derivs.enroll(mDDtShearDisplacement);

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
  this->addContacts(dataBase,state,derivs);

  mCyclesSinceLastKulling++;
  if (mCyclesSinceLastKulling % mKullFrequency == 0){ 
    this->kullInactiveContacts(dataBase,state,derivs);
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

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mUniqueIndices);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(angularVelocity);
    (*boundaryItr)->applyFieldListGhostBoundary(radius);
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

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mUniqueIndices);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(angularVelocity);
    (*boundaryItr)->enforceFieldListBoundary(radius);
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
  file.write(mDDtShearDisplacement, pathName + "/DDtShearDisplacement");
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
  file.read(mDDtShearDisplacement, pathName + "/DDtShearDisplacement");
  file.read(mEquilibriumOverlap, pathName + "/equilibriumOverlap");
}


//------------------------------------------------------------------------------
// REdistribution methods
//------------------------------------------------------------------------------

template<typename Dimension>
void
DEMBase<Dimension>::
initializeBeforeRedistribution(){

  //const auto  numNodeLists = mDataBase.numNodeLists();
  //const auto  nodeListPtrs = mDataBase.DEMNodeListPtrs();
  const auto& connectivityMap = mDataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // set all our active contacts to 1
#pragma omp for
  for (auto kk = 0u; kk < npairs; ++kk) {

    const auto i = pairs[kk].i_node;
    const auto j = pairs[kk].j_node;
    const auto nodeListi = pairs[kk].i_list;
    const auto nodeListj = pairs[kk].j_list;

    const auto storageLocation = this->findContactIndex(nodeListi,i,nodeListj,j);
    const auto storageNodeList = storageLocation[0];
    const auto storageNode = storageLocation[1];
    const auto storageContact = storageLocation[2];

    const auto storageUniqueNode = mUniqueIndices(storageNodeList,storageNode);

    int pairNodeList, pairNode;
    if (nodeListi == storageNodeList){
      pairNodeList = nodeListj;
      pairNode = j;
    } else{
      pairNodeList = nodeListi;
      pairNode = i;
    }
    mNeighborIndices(pairNodeList,pairNode).push_back(storageUniqueNode);
    mShearDisplacement(pairNodeList,pairNode).push_back(mShearDisplacement(storageNodeList,storageNode)[storageContact]);
    mEquilibriumOverlap(pairNodeList,pairNode).push_back(mEquilibriumOverlap(storageNodeList,storageNode)[storageContact]);
    mDDtShearDisplacement(pairNodeList,pairNode).push_back(mDDtShearDisplacement(storageNodeList,storageNode)[storageContact]);
  }
}

template<typename Dimension>
void
DEMBase<Dimension>::
finalizeAfterRedistribution() {
  std::cout << "SHABINGO!"<< std::endl;
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
addContacts(const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/){

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();
  
  auto neighborIds = state.fields(DEMFieldNames::neighborIndices, vector<int>());

#pragma omp for
  for (auto kk = 0u; kk < npairs; ++kk) {

    const auto i = pairs[kk].i_node;
    const auto j = pairs[kk].j_node;
    const auto nodeListi = pairs[kk].i_list;
    const auto nodeListj = pairs[kk].j_list;

    const auto dataLocation = this->storageNodeSelection(nodeListi,i,nodeListj,j);
    const auto storageNodeListIndex = dataLocation[0];
    const auto storageNodeIndex = dataLocation[1];
    const auto uniqueIndex = dataLocation[2];
 
    // get our contact indices 
    const auto neighborContacts = neighborIds(storageNodeListIndex,storageNodeIndex);
    const auto contactIndexPtr = std::find(neighborContacts.begin(),neighborContacts.end(),uniqueIndex);

    // add the contact if it is new
    if (contactIndexPtr == neighborContacts.end()) neighborIds(storageNodeListIndex,storageNodeIndex).push_back(uniqueIndex);

  }

  this->resizeStatePairFieldLists(state);

}


//------------------------------------------------------------------------------
// remove old contacts from our pair storage every so often
//------------------------------------------------------------------------------
template<typename Dimension>
void
DEMBase<Dimension>::
kullInactiveContacts(const DataBase<Dimension>& dataBase,
                           State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs){

  auto neighborIds = state.fields(DEMFieldNames::neighborIndices, vector<int>());
  auto isActive = state.fields(DEMFieldNames::isActiveContact, vector<int>());

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
      const auto numContacts = neighborIds(nodeListi,nodei).size();
      isActive(nodeListi,nodei).assign(numContacts,0);
    }
  }

  // set all our active contacts to 1
#pragma omp for
  for (auto kk = 0u; kk < npairs; ++kk) {

    const auto i = pairs[kk].i_node;
    const auto j = pairs[kk].j_node;
    const auto nodeListi = pairs[kk].i_list;
    const auto nodeListj = pairs[kk].j_list;

    const auto storageLocation = this->findContactIndex(nodeListi,i,nodeListj,j);
    const auto nodeListStore = storageLocation[0];
    const auto nodeStore = storageLocation[1];
    const auto contactStore = storageLocation[2];

    // add the contact if it is new
    if (contactStore == (int)isActive(nodeListStore,nodeStore).size()){
      isActive(nodeListStore,nodeStore).push_back(1);
    }
    isActive(nodeListStore,nodeStore)[contactStore] = 1;
  }

  this->kullInactiveContactsFromPairFieldList(neighborIds);
  this->kullInactiveContactsFromStatePairFieldLists(state);

}


//-------;-----------------------------------------------------------------------
// method to find where we are storing the pair contact
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<int>
DEMBase<Dimension>::
findContactIndex(int nodeListi,
                 int i,
                 int nodeListj,
                 int j) const {

  const auto dataLocation = this->storageNodeSelection(nodeListi,i,nodeListj,j);
  const auto storageNodeListIndex = dataLocation[0];
  const auto storageNodeIndex = dataLocation[1];
  const auto uniqueIndex = dataLocation[2];

  // get our contact indices 
  auto& neighborContacts = mNeighborIndices(storageNodeListIndex,storageNodeIndex);
  const auto contactIndexPtr = std::find(neighborContacts.begin(),neighborContacts.end(),uniqueIndex);
  const int  storageContactIndex = std::distance(neighborContacts.begin(),contactIndexPtr);

  // these contacts should all exist by the time this method is called
  if (contactIndexPtr == neighborContacts.end()) std::cout<<"Oopsie doodle" << std::endl;

  std::vector<int> storageIndices = {storageNodeListIndex,storageNodeIndex,storageContactIndex};
  return storageIndices; 
}

//------------------------------------------------------------------------------
// where are we storing this?
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<int>
DEMBase<Dimension>::
storageNodeSelection(int nodeListi,
                     int i,
                     int nodeListj,
                     int j) const {

  // get our unique global ID numbers
  const auto uIDi = mUniqueIndices(nodeListi,i);
  const auto uIDj = mUniqueIndices(nodeListj,j);

  const auto numInternalNodesi = mUniqueIndices[nodeListi]->numInternalElements();
  const auto numInternalNodesj = mUniqueIndices[nodeListj]->numInternalElements();

  // first we'll pick the internal node over a ghost. 
  // if both are internal we'll use the one with the lowerst global id
  const auto selectNodei = (i < (int)numInternalNodesi) and (uIDi < uIDj or j > (int)numInternalNodesj);

  int storageNodeListIndex, storageNodeIndex, uniqueIndex;
  if (selectNodei) {
    storageNodeListIndex = nodeListi;
    storageNodeIndex = i;
    uniqueIndex = uIDj;
  } else{
    storageNodeListIndex = nodeListj;
    storageNodeIndex = j;
    uniqueIndex = uIDi;
  }

  // return the nodelist,node ids plus the unique index we'll search for
  std::vector<int> storageIndices = {storageNodeListIndex,storageNodeIndex,uniqueIndex};
  return storageIndices; 
}

}
