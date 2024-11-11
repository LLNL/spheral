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
  AnyVisitor<void, std::any&, std::any&> ZERO;
  ZERO.addVisitor<FieldBase<Dimension>*>         ([](const std::any& x, const std::any& y) { std::any_cast<FieldBase<Dimension>*>(x)->Zero(); });
  ZERO.addVisitor<Scalar*>                       ([](const std::any& x, const std::any& y) { *std::any_cast<Scalar*>(x) = 0.0; });
  ZERO.addVisitor<Vector*>                       ([](const std::any& x, const std::any& y) { *std::any_cast<Vector*>(x) = Vector::zero; });
  ZERO.addVisitor<Tensor*>                       ([](const std::any& x, const std::any& y) { *std::any_cast<Tensor*>(x) = Tensor::zero; });
  ZERO.addVisitor<SymTensor*>                    ([](const std::any& x, const std::any& y) { *std::any_cast<SymTensor*>(x) = SymTensor::zero; });
  ZERO.addVisitor<vector<Scalar>*>               ([](const std::any& x, const std::any& y) { std::any_cast<vector<Scalar>*>(x)->clear(); });
  ZERO.addVisitor<vector<Vector>*>               ([](const std::any& x, const std::any& y) { std::any_cast<vector<Vector>*>(x)->clear(); });
  ZERO.addVisitor<vector<Tensor>*>               ([](const std::any& x, const std::any& y) { std::any_cast<vector<Tensor>*>(x)->clear(); });
  ZERO.addVisitor<vector<SymTensor>*>            ([](const std::any& x, const std::any& y) { std::any_cast<vector<SymTensor>*>(x)->clear(); });
  ZERO.addVisitor<set<int>*>                     ([](const std::any& x, const std::any& y) { std::any_cast<set<int>*>(x)->clear(); });
  ZERO.addVisitor<set<RKOrder>*>                 ([](const std::any& x, const std::any& y) { std::any_cast<set<int>*>(x)->clear(); });
  ZERO.addVisitor<ReproducingKernel<Dimension>*> ([](const std::any& x, const std::any& y) { });

  // Walk the state fields and zero them.
  for (auto [key, anyvalptr]: mStorage) {
    ZERO.visit(anyvalptr, anyvalptr);
  }

  // Reinitialize the node pair interaction information.
  initializeNodePairInformation();
}

}

