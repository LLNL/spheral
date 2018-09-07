#include "Utilities/DBC.hh"
#include "NodeList/NodeListRegistrar.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Access the current time.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Integrator<Dimension>::currentTime() const {
  return mCurrentTime;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
currentTime(typename Dimension::Scalar time) {
  mCurrentTime = time;
}

//------------------------------------------------------------------------------
// Access the current cycle.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
Integrator<Dimension>::currentCycle() const {
  return mCurrentCycle;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
currentCycle(int cycle) {
  mCurrentCycle = cycle;
}

//------------------------------------------------------------------------------
// Access the minimum allowed time step.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Integrator<Dimension>::dtMin() const {
  return mDtMin;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
dtMin(typename Dimension::Scalar dt) {
  mDtMin = dt;
}

//------------------------------------------------------------------------------
// Access the maximum allowed time step.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Integrator<Dimension>::dtMax() const {
  return mDtMax;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
dtMax(typename Dimension::Scalar dt) {
  mDtMax = dt;
}

//------------------------------------------------------------------------------
// Access the last time step.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Integrator<Dimension>::lastDt() const {
  return mLastDt;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
lastDt(typename Dimension::Scalar dt) {
  mLastDt = dt;
}

//------------------------------------------------------------------------------
// Access the time step growth limit.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
Integrator<Dimension>::dtGrowth() const {
  return mDtGrowth;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
dtGrowth(typename Dimension::Scalar fraction) {
  mDtGrowth = fraction;
}

//------------------------------------------------------------------------------
// Access the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const DataBase<Dimension>& 
Integrator<Dimension>::dataBase() const {
  CHECK(mDataBasePtr);
  return *mDataBasePtr;
}

//------------------------------------------------------------------------------
// Access the physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const std::vector<Physics<Dimension>*>&
Integrator<Dimension>::physicsPackages() const {
  return mPhysicsPackages;
}

//------------------------------------------------------------------------------
// Provide iterators over over the physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<Physics<Dimension>*>::iterator
Integrator<Dimension>::physicsPackagesBegin() {
  return mPhysicsPackages.begin();
}

template<typename Dimension>
inline
typename std::vector<Physics<Dimension>*>::iterator
Integrator<Dimension>::physicsPackagesEnd() {
  return mPhysicsPackages.end();
}

//------------------------------------------------------------------------------
// Provide const iterators over over the physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename std::vector<Physics<Dimension>*>::const_iterator
Integrator<Dimension>::physicsPackagesBegin() const {
  return mPhysicsPackages.begin();
}

template<typename Dimension>
inline
typename std::vector<Physics<Dimension>*>::const_iterator
Integrator<Dimension>::physicsPackagesEnd() const {
  return mPhysicsPackages.end();
}

//------------------------------------------------------------------------------
// Determine whether or not to treat the boundaries rigorously.  This really
// means whether or not we call setGhostNodes for every call to 
// enforceBoundaries.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Integrator<Dimension>::rigorousBoundaries() const {
  return mRigorousBoundaries;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
rigorousBoundaries(bool value) {
  mRigorousBoundaries = value;
}

//------------------------------------------------------------------------------
// In the case of non-rigorous boundaries, we specify a frequency for updating
// the boundaries and connectivity.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
Integrator<Dimension>::
updateBoundaryFrequency() const {
  return mUpdateBoundaryFrequency;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
updateBoundaryFrequency(const int value) {
  mUpdateBoundaryFrequency = value;
}

//------------------------------------------------------------------------------
// Select whether the integrator is verbose or not.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Integrator<Dimension>::verbose() const {
  return mVerbose;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
verbose(bool value) {
  mVerbose = value;
}

//------------------------------------------------------------------------------
// Select whether we're enforcing culling of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
Integrator<Dimension>::
cullGhostNodes() const {
  return mCullGhostNodes;
}

template<typename Dimension>
inline
void
Integrator<Dimension>::
cullGhostNodes(const bool x) {
  mCullGhostNodes = x;
}

//------------------------------------------------------------------------------
// Descendent classes can get write access to the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
DataBase<Dimension>& 
Integrator<Dimension>::accessDataBase() {
  CHECK(mDataBasePtr);
  return *mDataBasePtr;
}

}
