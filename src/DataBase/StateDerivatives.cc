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
#include "Utilities/AnyVisitor.hh"

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
  for (auto pkg: physicsPackages) pkg->registerDerivatives(dataBase, *this);
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
  for (auto pkg: range(physicsPackageBegin, physicsPackageEnd)) pkg->registerDerivatives(dataBase, *this);
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

  // Build a visitor to zero each data type
  AnyVisitor<void, std::any&> ZERO;
  ZERO.addVisitor<std::reference_wrapper<FieldBase<Dimension>>>         ([](const std::any& x) { std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(x).get().Zero(); });
  ZERO.addVisitor<std::reference_wrapper<Scalar>>                       ([](const std::any& x) { std::any_cast<std::reference_wrapper<Scalar>>(x).get() = 0.0; });
  ZERO.addVisitor<std::reference_wrapper<Vector>>                       ([](const std::any& x) { std::any_cast<std::reference_wrapper<Vector>>(x).get() = Vector::zero; });
  ZERO.addVisitor<std::reference_wrapper<Tensor>>                       ([](const std::any& x) { std::any_cast<std::reference_wrapper<Tensor>>(x).get() = Tensor::zero; });
  ZERO.addVisitor<std::reference_wrapper<SymTensor>>                    ([](const std::any& x) { std::any_cast<std::reference_wrapper<SymTensor>>(x).get() = SymTensor::zero; });
  ZERO.addVisitor<std::reference_wrapper<vector<Scalar>>>               ([](const std::any& x) { std::any_cast<std::reference_wrapper<vector<Scalar>>>(x).get().clear(); });
  ZERO.addVisitor<std::reference_wrapper<vector<Vector>>>               ([](const std::any& x) { std::any_cast<std::reference_wrapper<vector<Vector>>>(x).get().clear(); });
  ZERO.addVisitor<std::reference_wrapper<vector<Tensor>>>               ([](const std::any& x) { std::any_cast<std::reference_wrapper<vector<Tensor>>>(x).get().clear(); });
  ZERO.addVisitor<std::reference_wrapper<vector<SymTensor>>>            ([](const std::any& x) { std::any_cast<std::reference_wrapper<vector<SymTensor>>>(x).get().clear(); });
  ZERO.addVisitor<std::reference_wrapper<set<int>>>                     ([](const std::any& x) { std::any_cast<std::reference_wrapper<set<int>>>(x).get().clear(); });
  ZERO.addVisitor<std::reference_wrapper<set<RKOrder>>>                 ([](const std::any& x) { std::any_cast<std::reference_wrapper<set<int>>>(x).get().clear(); });
  ZERO.addVisitor<std::reference_wrapper<ReproducingKernel<Dimension>>> ([](const std::any& x) { } );

  // Walk the state values and zero them
  for (auto itr: mStorage) {
    ZERO.visit(itr.second);
  }

  // Reinitialize the node pair interaction information.
  initializeNodePairInformation();
}

}

