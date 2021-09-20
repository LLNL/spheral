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

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateDerivatives<Dimension>::
StateDerivatives():
  StateBase<Dimension>(),
  mCalculatedNodePairs(),
  mNumSignificantNeighbors() {
}

//------------------------------------------------------------------------------
// Construct with the derivatives for the given set of Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
StateDerivatives<Dimension>::
StateDerivatives(DataBase<Dimension>& dataBase,
                 typename StateDerivatives<Dimension>::PackageList& physicsPackages):
  StateBase<Dimension>(),
  mCalculatedNodePairs(),
  mNumSignificantNeighbors() {

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
StateDerivatives(DataBase<Dimension>& dataBase,
                 typename StateDerivatives<Dimension>::PackageIterator physicsPackageBegin,
                 typename StateDerivatives<Dimension>::PackageIterator physicsPackageEnd):
  StateBase<Dimension>(),
  mCalculatedNodePairs(),
  mNumSignificantNeighbors() {

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
operator==(const StateBase<Dimension>& rhs) const {
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
  bool result = true;
  typename CalculatedPairType::const_iterator itr = mCalculatedNodePairs.begin();
  while (result && itr != mCalculatedNodePairs.end()) {
    const NodeIteratorBase<Dimension> nodeI = itr->first;
    const vector<NodeIteratorBase<Dimension> > neighbors = itr->second;
    for (typename vector<NodeIteratorBase<Dimension> >::const_iterator nodeJItr = neighbors.begin();
         (nodeJItr != neighbors.end()) && result;
         ++nodeJItr) {
      typename CalculatedPairType::const_iterator itr2 = mCalculatedNodePairs.find(*nodeJItr);
      CONTRACT_VAR(itr2);
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
       ++itr) {

    try {
      auto ptr = boost::any_cast<FieldBase<Dimension>*>(itr->second);
      ptr->Zero();

    } catch (const boost::bad_any_cast&) {
      try {
        auto ptr = boost::any_cast<vector<Vector>*>(itr->second);
        ptr->clear();

      } catch (const boost::bad_any_cast&) {
        try {
          auto ptr = boost::any_cast<vector<Scalar>*>(itr->second);
          ptr->clear();

        } catch (const boost::bad_any_cast&) {
          VERIFY2(false, "StateDerivatives::Zero ERROR: unknown type for key " << itr->first << "\n");
        }
      }
    }
  }

  // Reinitialize the node pair interaction information.
  initializeNodePairInformation();
}

}

