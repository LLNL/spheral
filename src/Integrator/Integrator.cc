//---------------------------------Spheral++----------------------------------//
// Integrator -- The topmost abstract base class for all integrator classes
// Spheral++.  Integrator classes take a list of Physics packages and advance
// them in time.
//
// Created by JMO, Wed May 31 21:58:08 PDT 2000
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>

#include "Integrator.hh"
#include "DataOutput/Restart.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"
#include "FileIO/FileIO.hh"
#include "Hydro/HydroFieldNames.hh"
// #include "Utilities/timingUtilities.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/allReduce.hh"
#include "Distributed/Communicator.hh"

#include "Utilities/DBC.hh"

namespace Spheral {
namespace IntegratorSpace {

using namespace std;

using NodeSpace::NodeList;
using DataBaseSpace::DataBase;
using PhysicsSpace::Physics;
using FileIOSpace::FileIO;
using FieldSpace::FieldList;
using BoundarySpace::Boundary;
using NeighborSpace::ConnectivityMap;

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Integrator<Dimension>::Integrator():
  mDtMin(0.0),
  mDtMax(FLT_MAX),
  mDtGrowth(2.0),
  mLastDt(1e-5),
  mCurrentTime(0.0),
  mCurrentCycle(0),
  mVerbose(false),
  mRequireConnectivity(true),
  mDtThreshold(1.0e-10),
  mDataBasePtr(0),
  mPhysicsPackages(0),
  mRigorousBoundaries(false),
  mUpdateBoundaryFrequency(1),
  mCullGhostNodes(true),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
Integrator<Dimension>::
Integrator(DataBase<Dimension>& dataBase):
  mDtMin(0.0),
  mDtMax(FLT_MAX),
  mDtGrowth(2.0),
  mLastDt(1e-5),
  mCurrentTime(0.0),
  mCurrentCycle(0),
  mVerbose(false),
  mRequireConnectivity(true),
  mDtThreshold(1.0e-10),
  mDataBasePtr(&dataBase),
  mPhysicsPackages(0),
  mRigorousBoundaries(false),
  mUpdateBoundaryFrequency(1),
  mCullGhostNodes(true),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
Integrator<Dimension>::
Integrator(DataBase<Dimension>& dataBase,
           const vector<Physics<Dimension>*>& physicsPackages):
  mDtMin(0.0),
  mDtMax(FLT_MAX),
  mDtGrowth(2.0),
  mLastDt(1e-5),
  mCurrentTime(0.0),
  mCurrentCycle(0),
  mVerbose(false),
  mRequireConnectivity(true),
  mDtThreshold(1.0e-10),
  mDataBasePtr(&dataBase),
  mPhysicsPackages(physicsPackages),
  mRigorousBoundaries(false),
  mUpdateBoundaryFrequency(1),
  mCullGhostNodes(true),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
Integrator<Dimension>::~Integrator() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
Integrator<Dimension>&
Integrator<Dimension>::
operator=(const Integrator<Dimension>& rhs) {
  if (this != &rhs) {
    mDtMin = rhs.mDtMin;
    mDtMax = rhs.mDtMax;
    mDtGrowth = rhs.mDtGrowth;
    mLastDt = rhs.mLastDt;
    mCurrentTime = rhs.mCurrentTime;
    mCurrentCycle = rhs.mCurrentCycle;
    mDataBasePtr = rhs.mDataBasePtr;
    mPhysicsPackages = rhs.mPhysicsPackages;
    mRigorousBoundaries = rhs.mRigorousBoundaries;
    mUpdateBoundaryFrequency = rhs.mUpdateBoundaryFrequency;
    mCullGhostNodes = rhs.mCullGhostNodes;
    mVerbose = rhs.mVerbose;
    mRequireConnectivity = rhs.mRequireConnectivity;
    mDtThreshold = rhs.mDtThreshold;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Loop over the stored physics packages and pick the minimum timestep.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
Integrator<Dimension>::
selectDt(const typename Dimension::Scalar dtMin,
         const typename Dimension::Scalar dtMax,
         const State<Dimension>& state,
         const StateDerivatives<Dimension>& derivs) const {

  REQUIRE(dtMin >= 0 and dtMax > 0);
  REQUIRE(dtMin <= dtMax);

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Get the current time and data base.
  Scalar t = currentTime();
  const DataBase<Dimension>& db = dataBase();

  // Loop over each package, and pick their timesteps.
  TimeStepType dt(dtMax, "");
  for (ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd(); 
       ++physicsItr) {
    TimeStepType dtVote = (*physicsItr)->dt(db, state, derivs, t);
    if (dtVote.first > 0.0 and dtVote.first < dt.first) dt = dtVote;
  }

  // In the parallel case we need to find the minimum timestep across all
  // processors.
  const Scalar globalDt = allReduce(dt.first, MPI_MIN, Communicator::communicator());
  if (dt.first == globalDt and 
      (verbose() or globalDt < mDtThreshold)) {
    cout << "Selected timestep of "
	 << dt.first << endl
         << dt.second << endl;
  }
  cout.flush();
  dt.first = globalDt;

  // We also require that the timestep is not allowed to grow faster than a
  // prescribed fraction.
  dt.first = min(dt.first, dtGrowth()*lastDt());

  // Enforce the timestep boundaries.
  dt.first = min(dtMax, max(dtMin, dt.first));

  CHECK(dt.first >= 0.0 and
        dt.first >= dtMin and dt.first <= dtMax);

  return dt.first;
}

//------------------------------------------------------------------------------
// Perform basic initializations for the integrators.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::initialize(State<Dimension>& state,
                                  StateDerivatives<Dimension>& derivs) {

  // Check if we need to construct connectivity.
  mRequireConnectivity = false;
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) mRequireConnectivity = (mRequireConnectivity or
                                             (*physicsItr)->requireConnectivity());

  // Set the boundary conditions.
  if ((not mRigorousBoundaries) and (mCurrentCycle % mUpdateBoundaryFrequency == 0)) {
    setGhostNodes();
  }
  applyGhostBoundaries(state, derivs);

  // Register the now updated connectivity with the state.
  if (mRequireConnectivity) {
    DataBase<Dimension>& db = accessDataBase();
    state.enrollConnectivityMap(db.connectivityMapPtr());
  }

  // Prepare individual packages and node lists.
  preStepInitialize(currentTime(), 0.0, state, derivs);
}

//------------------------------------------------------------------------------
// Initializations that need to be called before a time advance step begins.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::preStepInitialize(const double t,
                                         const double dt,
                                         State<Dimension>& state,
                                         StateDerivatives<Dimension>& derivs) {

  // Loop over the NodeLists in the DataBase and initialize them.  Also make sure
  // they know about the smoothing scale constraints.
  DataBase<Dimension>& db = accessDataBase();
  for (typename DataBase<Dimension>::NodeListIterator nodeListItr = db.nodeListBegin();
       nodeListItr < db.nodeListEnd();
       ++nodeListItr) {
    (*nodeListItr)->work() = 0.0;
  }

  // Loop over the physics packages and perform any necessary initializations.
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) {
    (*physicsItr)->initialize(t, dt, db, state, derivs);
  }

  // Physics packages may have called boundary conditions as well, so finalize any
  // outstanding boundary conditions here.
  this->finalizeGhostBoundaries();
}

//------------------------------------------------------------------------------
// Finalize at the completion of a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::finalize(const double t,
                                const double dt,
                                State<Dimension>& state,
                                StateDerivatives<Dimension>& derivs) {

  // Loop over the physics packages and perform any necessary finalizations.
  DataBase<Dimension>& db = accessDataBase();
//   for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = db.fluidNodeListBegin();
//        nodeListItr != db.fluidNodeListEnd(); 
//        ++nodeListItr) {
//     (*nodeListItr)->neighbor().updateNodes();
//   }
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) {
    (*physicsItr)->finalize(t, dt, db, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Iterate over all physics packages and call evaluateDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::evaluateDerivatives(const Scalar t,
                                           const Scalar dt,
                                           const DataBase<Dimension>& dataBase,
                                           const State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivs) const {

  // Loop over the physics packages and have them evaluate their derivatives.
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) {
    (*physicsItr)->evaluateDerivatives(t, dt, dataBase, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Iterate over all physics packages and call finalizeDerivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::finalizeDerivatives(const Scalar t,
                                           const Scalar dt,
                                           const DataBase<Dimension>& dataBase,
                                           const State<Dimension>& state,
                                           StateDerivatives<Dimension>& derivs) const {

  // Loop over the physics packages and have them finalize their derivatives.
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) {
    (*physicsItr)->finalizeDerivatives(t, dt, dataBase, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Iterate over all physics packages and let them do any post state update 
// stuff.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::postStateUpdate(const DataBase<Dimension>& dataBase, 
                                       State<Dimension>& state,
                                       const StateDerivatives<Dimension>& derivs) const {

  // Loop over the physics packages.
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) {
    (*physicsItr)->postStateUpdate(dataBase, state, derivs);
  }
}

//------------------------------------------------------------------------------
// Add a physics package.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
appendPhysicsPackage(Physics<Dimension>& package) {
  if (!havePhysicsPackage(package)) {
    mPhysicsPackages.push_back(&package);
  } else {
    cerr << "Warning: attempt to append Physics package " << &package
         << "to Integrator " << this << " which already has it." << endl;
  }
}

//------------------------------------------------------------------------------
// Test if the given physics package is listed in the integrator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Integrator<Dimension>::
havePhysicsPackage(const Physics<Dimension>& package) const {
  return count(mPhysicsPackages.begin(), mPhysicsPackages.end(), &package) > 0;
}

//------------------------------------------------------------------------------
// Get the unique set of boundary conditions across all physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<Boundary<Dimension>*>
Integrator<Dimension>::
uniqueBoundaryConditions() const {

  // The result we're going to build up.
  vector<Boundary<Dimension>*> result;

  // Iterate over each physics package.
  for (ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) {

    // Iterate over the boundary conditions associated with this package.
    for (ConstBoundaryIterator boundaryItr = (*physicsItr)->boundaryBegin();
         boundaryItr != (*physicsItr)->boundaryEnd();
         ++boundaryItr) {
      if (find(result.begin(), result.end(), *boundaryItr) == result.end())
        result.push_back(*boundaryItr);
    }

  }

  BEGIN_CONTRACT_SCOPE;
  // Ensure that all boundary conditions are included in the result
  for (ConstPackageIterator physicsItr = physicsPackagesBegin();
       physicsItr != physicsPackagesEnd();
       ++physicsItr) {
    for (ConstBoundaryIterator boundaryItr = (*physicsItr)->boundaryBegin();
         boundaryItr != (*physicsItr)->boundaryEnd();
         ++boundaryItr) {
      ENSURE(find(result.begin(), result.end(), *boundaryItr) != result.end());
    }
  }

  // Also ensure that there are no duplicates.
  for (ConstBoundaryIterator boundaryItr = result.begin();
       boundaryItr != result.end();
       ++boundaryItr) {
    ENSURE(count(result.begin(), result.end(), *boundaryItr) == 1);
  }
  END_CONTRACT_SCOPE;

  // That's it.
  return result;
}

//------------------------------------------------------------------------------
// Set up the ghost nodes for boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::setGhostNodes() {

//   // Start our work timer.
//   typedef Timing::Time Time;
//   const Time start = Timing::currentTime();

  if (mRequireConnectivity) {

    // Get that DataBase.
    DataBase<Dimension>& db = accessDataBase();

    // Remove any old ghost node information from the NodeLists.
    for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = db.fluidNodeListBegin();
         nodeListItr != db.fluidNodeListEnd(); 
         ++nodeListItr) {
      (*nodeListItr)->numGhostNodes(0);
      (*nodeListItr)->neighbor().updateNodes();
    }

    // Get the complete set of unique boundary conditions.
    const vector<Boundary<Dimension>*> boundaries = uniqueBoundaryConditions();

    // Iterate over the boundaries and set their ghost node info.
    for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
         boundaryItr != boundaries.end();
         ++boundaryItr) {
      (*boundaryItr)->setAllGhostNodes(db);
      (*boundaryItr)->finalizeGhostBoundary();
      for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = db.fluidNodeListBegin();
           nodeListItr != db.fluidNodeListEnd(); 
           ++nodeListItr) {
        (*nodeListItr)->neighbor().updateNodes();
      }
    }

    //   // Update the work per node fields.
    //   double deltaLocal = Timing::difference(start, Timing::currentTime());
    //   int numLocalNodes = 0;
    //   for (typename DataBase<Dimension>::NodeListIterator nodeListItr = db.nodeListBegin();
    //        nodeListItr != db.nodeListEnd(); 
    //        ++nodeListItr) numLocalNodes += (*nodeListItr)->numInternalNodes();
    //   const double workPerNode = deltaLocal/(numLocalNodes + 1.0e-30);
    //   FieldList<Dimension, Scalar> work = db.globalWork();
    //   work += workPerNode;

    // Update the connectivity.
    db.updateConnectivityMap();

    // If we're culling ghost nodes, do it now.
    if (mCullGhostNodes and not this->domainDecompositionIndependent()) {

      // First build the set of flags indicating which nodes are used.
      FieldList<Dimension, int> flags = db.newGlobalFieldList(int(0), "active nodes");
      const ConnectivityMap<Dimension>& cm = db.connectivityMap();
      int nodeListi = 0;
      for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = db.nodeListBegin();
           nodeListItr != db.nodeListEnd(); 
           ++nodeListItr, ++nodeListi) {
        const NodeList<Dimension>& nodeList = **nodeListItr;
        for (int i = 0; i != nodeList.numInternalNodes(); ++i) {
          flags(nodeListi, i) = 1;
          const vector<vector<int> >& fullConnectivity = cm.connectivityForNode(&nodeList, i);
          for (int nodeListj = 0; nodeListj != fullConnectivity.size(); ++nodeListj) {
            const vector<int>& connectivity = fullConnectivity[nodeListj];
            for (vector<int>::const_iterator jItr = connectivity.begin();
                 jItr != connectivity.end();
                 ++jItr) flags(nodeListj, *jItr) = 1;
          }
        }

        // Ghost nodes that are control nodes for other ghost nodes we're keeping must
        // be kept as well.
        const size_t firstGhostNode = nodeList.firstGhostNode();
        for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
             boundaryItr != boundaries.end();
             ++boundaryItr) {
          const Boundary<Dimension>& boundary = **boundaryItr;
          const vector<int>& controlNodes = boundary.controlNodes(nodeList);
          const vector<int>& ghostNodes = boundary.ghostNodes(nodeList);
          // CHECK(controlNodes.size() == ghostNodes.size());  // Not true if this is a DistributedBoundary!
          for (size_t k = 0; k != controlNodes.size(); ++k) {
            if (controlNodes[k] >= firstGhostNode) { //  and flags(nodeListi, ghostNodes[k]) == 1)
              flags(nodeListi, controlNodes[k]) = 1;
            }
          }
        }
      }

      // Create the index mapping from old to new node orderings.
      FieldList<Dimension, int> old2newIndexMap = db.newGlobalFieldList(int(0), "index map");
      nodeListi = 0;
      for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = db.nodeListBegin();
           nodeListItr != db.nodeListEnd(); 
           ++nodeListItr, ++nodeListi) {
        const int numNodes = (**nodeListItr).numNodes();
        for (int i = 0; i != numNodes; ++i) old2newIndexMap(nodeListi, i) = i;
      }

      // Now use these flags to cull the boundary conditions.
      const size_t numNodeLists = db.numNodeLists();
      vector<int> numNodesRemoved(numNodeLists, 0);
      for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
           boundaryItr != boundaries.end();
           ++boundaryItr) {
        (*boundaryItr)->cullGhostNodes(flags, old2newIndexMap, numNodesRemoved);
      }

      // Patch up the connectivity map.
      db.patchConnectivityMap(flags, old2newIndexMap);

      // Now the boundary conditions have been updated, so we can go ahead and remove
      // the ghost nodes themselves from the NodeLists.
      nodeListi = 0;
      for (typename DataBase<Dimension>::ConstNodeListIterator nodeListItr = db.nodeListBegin();
           nodeListItr != db.nodeListEnd(); 
           ++nodeListItr, ++nodeListi) {
        NodeList<Dimension>& nodeList = **nodeListItr;
        vector<int> nodesToRemove;
        for (size_t i = nodeList.firstGhostNode(); i != nodeList.numNodes(); ++i) {
          if (flags(nodeListi, i) == 0) nodesToRemove.push_back(i);
        }
        nodeList.deleteNodes(nodesToRemove);
        nodeList.neighbor().updateNodes();
      }

      // All nodes should now be labeled as keepers.
      BEGIN_CONTRACT_SCOPE;
      {
        for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
          ENSURE(flags[nodeListi]->numElements() == 0 or
                 *min_element(flags[nodeListi]->begin(), flags[nodeListi]->end()) == 1);
        }
      }
      END_CONTRACT_SCOPE;

      // The ConnectivityMap should be valid too.
      ENSURE(db.connectivityMap().valid());
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::applyGhostBoundaries(State<Dimension>& state,
                                            StateDerivatives<Dimension>& derivs) {

//   // Start our work timer.
//   typedef Timing::Time Time;
//   const Time start = Timing::currentTime();

  if (mRequireConnectivity) {

    // Get that DataBase.
    DataBase<Dimension>& db = accessDataBase();

    // If we're being rigorous about boundaries, we have to reset the ghost nodes.
    const vector<Boundary<Dimension>*> boundaries = uniqueBoundaryConditions();
    if (mRigorousBoundaries) {
      setGhostNodes();

    } else {

      // If we didn't call setGhostNodes, then make each boundary update it's 
      // ghost node info (position and H).
      for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
           boundaryItr != boundaries.end();
           ++boundaryItr) {
        for (typename DataBase<Dimension>::NodeListIterator nodeListItr = db.nodeListBegin();
             nodeListItr != db.nodeListEnd(); 
             ++nodeListItr) {
          (*boundaryItr)->updateGhostNodes(**nodeListItr);
        }
        (*boundaryItr)->finalizeGhostBoundary();
      }
      for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = db.fluidNodeListBegin();
           nodeListItr != db.fluidNodeListEnd(); 
           ++nodeListItr) {
        (*nodeListItr)->neighbor().updateNodes();
      }
    }

    // Iterate over the physics packages, and have them apply ghost boundaries
    // for their state.
    for (ConstPackageIterator physicsItr = physicsPackagesBegin();
         physicsItr != physicsPackagesEnd();
         ++physicsItr) {
      (*physicsItr)->applyGhostBoundaries(state, derivs);
    }

    //   // Update the work per node fields.
    //   double deltaLocal = Timing::difference(start, Timing::currentTime());
    //   int numLocalNodes = 0;
    //   for (typename DataBase<Dimension>::NodeListIterator nodeListItr = db.nodeListBegin();
    //        nodeListItr != db.nodeListEnd(); 
    //        ++nodeListItr) numLocalNodes += (*nodeListItr)->numInternalNodes();
    //   const double workPerNode = deltaLocal/(numLocalNodes + 1.0e-30);
    //   FieldList<Dimension, Scalar> work = db.globalWork();
    //   work += workPerNode;

    // Finalize the boundaries.
    this->finalizeGhostBoundaries();
  }
}

//------------------------------------------------------------------------------
// Finalize the ghost boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::finalizeGhostBoundaries() {

//   // Start our work timer.
//   typedef Timing::Time Time;
//   const Time start = Timing::currentTime();

  if (mRequireConnectivity) {

    // Get that DataBase.
    DataBase<Dimension>& db = accessDataBase();

    // If we're being rigorous about boundaries, we have to reset the ghost nodes.
    const vector<Boundary<Dimension>*> boundaries = uniqueBoundaryConditions();
    for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
         boundaryItr != boundaries.end();
         ++boundaryItr) {
      (*boundaryItr)->finalizeGhostBoundary();
    }

    //   // Update the work per node fields.
    //   double deltaLocal = Timing::difference(start, Timing::currentTime());
    //   int numLocalNodes = 0;
    //   for (typename DataBase<Dimension>::NodeListIterator nodeListItr = db.nodeListBegin();
    //        nodeListItr != db.nodeListEnd(); 
    //        ++nodeListItr) numLocalNodes += (*nodeListItr)->numInternalNodes();
    //   const double workPerNode = deltaLocal/(numLocalNodes + 1.0e-30);
    //   FieldList<Dimension, Scalar> work = db.globalWork();
    //   work += workPerNode;
  }
}

//------------------------------------------------------------------------------
// Identify the internal nodes in violation of the boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::setViolationNodes() {

  if (mRequireConnectivity) {

    // Get that DataBase.
    DataBase<Dimension>& db = accessDataBase();

    // Get the complete set of unique boundary conditions.
    const vector<Boundary<Dimension>*> boundaries = uniqueBoundaryConditions();

    // Have each boundary identify the set of nodes that violate it.
    for (ConstBoundaryIterator boundaryItr = boundaries.begin(); 
         boundaryItr != boundaries.end();
         ++boundaryItr) {
      (*boundaryItr)->setAllViolationNodes(db);
    }

    // Fix neighbor information.
    for (typename DataBase<Dimension>::FluidNodeListIterator nodeListItr = db.fluidNodeListBegin();
         nodeListItr != db.fluidNodeListEnd(); 
         ++nodeListItr) {
      (*nodeListItr)->neighbor().updateNodes();
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::enforceBoundaries(State<Dimension>& state,
                                         StateDerivatives<Dimension>& derivs) {

  if (mRequireConnectivity) {

    // Have each boundary identify the set of nodes in violation.  This also resets
    // the positions and H's of the nodes to be in compliance.
    setViolationNodes();

    // Iterate over the physics packages, and have them apply ghost boundaries
    // for their state.
    for (ConstPackageIterator physicsItr = physicsPackagesBegin();
         physicsItr != physicsPackagesEnd();
         ++physicsItr) (*physicsItr)->enforceBoundaries(state, derivs);
  }
}

//------------------------------------------------------------------------------
// Copy the ghost positions and H tensor from one state to another.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::copyGhostState(const State<Dimension>& state0,
                                      State<Dimension>& state1) const {
  const FieldList<Dimension, Vector> x0 = state0.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H0 = state0.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Vector> x1 = state1.fields(HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, SymTensor> H1 = state1.fields(HydroFieldNames::H, SymTensor::zero);
  for (GhostNodeIterator<Dimension> itr = x0.ghostNodeBegin();
       itr != x0.ghostNodeEnd();
       ++itr) {
    x1(itr) = x0(itr);
    H1(itr) = H0(itr);
  }
}

//------------------------------------------------------------------------------
// Advance the set of physics packages to the given time.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
advance(typename Dimension::Scalar goalTime) {

  // Loop and advance the system until the goal time is achieved.
  while (currentTime() < goalTime) {
    step(goalTime);
  }
}

//------------------------------------------------------------------------------
// Flag for whether we should try to run in a domain decomposition independent/
// reproducing mode.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
Integrator<Dimension>::
domainDecompositionIndependent() const {
  return NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent();
}

template<typename Dimension>
void
Integrator<Dimension>::
domainDecompositionIndependent(const bool x) {
  NodeListRegistrar<Dimension>::instance().domainDecompositionIndependent(x);
  // if (x) mCullGhostNodes = false;
}

//------------------------------------------------------------------------------
// Dump the current state of the Integrator to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  // file.write(dtMin(), pathName + "/dtMin");
  // file.write(dtMax(), pathName + "/dtMax");
  // file.write(dtGrowth(), pathName + "/dtGrowth");
  file.write(lastDt(), pathName + "/lastDt");
  file.write(currentTime(), pathName + "/currentTime");
  file.write(currentCycle(), pathName + "/currentCycle");
}  

//------------------------------------------------------------------------------
// Restore the state of the NodeList from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Integrator<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  // file.read(mDtMin, pathName + "/dtMin");
  // file.read(mDtMax, pathName + "/dtMax");
  // file.read(mDtGrowth, pathName + "/dtGrowth");
  file.read(mLastDt, pathName + "/lastDt");
  file.read(mCurrentTime, pathName + "/currentTime");
  file.read(mCurrentCycle, pathName + "/currentCycle");
}

}
}

