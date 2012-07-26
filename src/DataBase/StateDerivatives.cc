//---------------------------------Spheral++----------------------------------//
// State -- Accumulate and cart the state for a set of physics packages around.
//
// Created by JMO, Fri Aug 27 10:56:40 2004
//----------------------------------------------------------------------------//

#include "StateDerivatives.hh"
#include "StateBase.hh"
#include "DataBase.hh"
#include "Physics/Physics.hh"
#include "Field/Field.hh"

#include "TAU.h"

namespace Spheral {

using namespace std;
using FieldSpace::FieldBase;
using FieldSpace::Field;
using PhysicsSpace::Physics;

//------------------------------------------------------------------------------
// Construct with the derivatives for the given set of Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
StateDerivatives<Dimension>::
StateDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                 typename StateDerivatives<Dimension>::PackageList& physicsPackages):
  StateBase<Dimension>(),
  mCalculatedNodePairs(),
  mNumSignificantNeighbors() {
  TAU_PROFILE("StateDerivatives", "::StateDerivatives(db, packages)", TAU_USER);

  // Iterate over the physics packages, and have them register their derivatives.
  for (PackageIterator itr = physicsPackages.begin();
       itr != physicsPackages.end();
       ++itr) (*itr)->registerDerivatives(dataBase, *this);
}

//------------------------------------------------------------------------------
// Construct with the derivatives for the given set of Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
StateDerivatives<Dimension>::
StateDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                 typename StateDerivatives<Dimension>::PackageIterator physicsPackageBegin,
                 typename StateDerivatives<Dimension>::PackageIterator physicsPackageEnd):
  StateBase<Dimension>(),
  mCalculatedNodePairs(),
  mNumSignificantNeighbors() {
  TAU_PROFILE("StateDerivatives", "::StateDerivatives(db, packageBegin, packageEnd)", TAU_USER);

  // Iterate over the physics packages, and have them register their derivatives.
  for (PackageIterator itr = physicsPackageBegin;
       itr != physicsPackageEnd;
       ++itr) (*itr)->registerDerivatives(dataBase, *this);
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateDerivatives<Dimension>::
StateDerivatives(const StateDerivatives<Dimension>& rhs):
  StateBase<Dimension>(rhs),
  mCalculatedNodePairs(rhs.mCalculatedNodePairs),
  mNumSignificantNeighbors(rhs.mNumSignificantNeighbors) {
  TAU_PROFILE("StateDerivatives", "::StateDerivatives(StateDerivatives&)", TAU_USER);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateDerivatives<Dimension>::
~StateDerivatives() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
StateDerivatives<Dimension>&
StateDerivatives<Dimension>::
operator=(const StateDerivatives<Dimension>& rhs) {
  TAU_PROFILE("StateDerivatives", "::operator=(StateDerivatives&)", TAU_USER);
  if (this != &rhs) {
    StateBase<Dimension>::operator=(rhs);
    mCalculatedNodePairs = rhs.mCalculatedNodePairs;
    mNumSignificantNeighbors = rhs.mNumSignificantNeighbors;
  }
  return *this;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateDerivatives<Dimension>::
operator==(const StateDerivatives<Dimension>& rhs) const {
  return StateBase<Dimension>::operator==(rhs);
}

//------------------------------------------------------------------------------
// (Re)initialize the internal data structure for tracking calculated node 
// pairs.  This also initializes the number of significant neighbor tracking.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateDerivatives<Dimension>::
initializeNodePairInformation() {
  TAU_PROFILE("StateDerivatives", "::initializeNodePairInformation()", TAU_USER);

  // Clear out any existing info.
  mCalculatedNodePairs = CalculatedPairType();
  mNumSignificantNeighbors = SignificantNeighborMapType();

}

//------------------------------------------------------------------------------
// Check to see if the node interaction map is symmetric.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateDerivatives<Dimension>::
calculatedNodePairsSymmetric() const {
  TAU_PROFILE("StateDerivatives", "::calculatedNodePairsSymmetric()", TAU_USER);
  bool result = true;
  typename CalculatedPairType::const_iterator itr = mCalculatedNodePairs.begin();
  while (result && itr != mCalculatedNodePairs.end()) {
    const NodeIteratorBase<Dimension> nodeI = itr->first;
    const vector<NodeIteratorBase<Dimension> > neighbors = itr->second;
    for (typename vector<NodeIteratorBase<Dimension> >::const_iterator nodeJItr = neighbors.begin();
         (nodeJItr != neighbors.end()) && result;
         ++nodeJItr) {
      typename CalculatedPairType::const_iterator itr2 = mCalculatedNodePairs.find(*nodeJItr);
      CHECK(itr2 != mCalculatedNodePairs.end());
      const vector<NodeIteratorBase<Dimension> > neighborsJ = itr->second;
      result = result && (find(neighborsJ.begin(), neighborsJ.end(), nodeI) != neighborsJ.end());
    }
    ++itr;
  }
  return result;
}

//------------------------------------------------------------------------------
// Zero out all the stored derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateDerivatives<Dimension>::
Zero() {

  // Walk the state fields and zero them.
  for (typename StateBase<Dimension>::StorageType::iterator itr = this->mStorage.begin();
       itr != this->mStorage.end();
       ++itr) itr->second->Zero();

  // Reinitialize the node pair interaction information.
  initializeNodePairInformation();
}

}

//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class StateDerivatives<Dim<1> >;
  template class StateDerivatives<Dim<2> >;
  template class StateDerivatives<Dim<3> >;
}
