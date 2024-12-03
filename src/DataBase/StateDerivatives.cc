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
#include "Neighbor/PairwiseField.hh"
#include "Utilities/AnyVisitor.hh"
#include "Utilities/DataTypeTraits.hh"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;
using std::any;
using std::any_cast;
using std::reference_wrapper;

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Add zero methods to a Visitor
//------------------------------------------------------------------------------
template<typename VisitorType, typename T>
void addZero(VisitorType& visitor) {
  visitor.template addVisitor<reference_wrapper<T>>([](std::any& x) {
                                                      any_cast<reference_wrapper<T>>(x).get() = DataTypeTraits<T>::zero();
                                                    });
}

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
  auto cmp = dataBase.connectivityMapPtr();
  if (cmp) this->enrollConnectivityMap(cmp);
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
  auto cmp = dataBase.connectivityMapPtr();
  if (cmp) this->enrollConnectivityMap(cmp);
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
  using VisitorType = AnyVisitor<void, std::any&>;
  VisitorType ZERO;
  ZERO.addVisitor<std::reference_wrapper<FieldBase<Dimension>>>                          ([](const any& x) { any_cast<reference_wrapper<FieldBase<Dimension>>>(x).get().Zero(); });
  addZero<VisitorType, Scalar>           (ZERO);
  addZero<VisitorType, Vector>           (ZERO);
  addZero<VisitorType, Tensor>           (ZERO);
  addZero<VisitorType, SymTensor>        (ZERO);
  addZero<VisitorType, vector<Scalar>>   (ZERO);
  addZero<VisitorType, vector<Vector>>   (ZERO);
  addZero<VisitorType, vector<Tensor>>   (ZERO);
  addZero<VisitorType, vector<SymTensor>>(ZERO);
  addZero<VisitorType, set<int>>         (ZERO);
  addZero<VisitorType, set<RKOrder>>     (ZERO);
  ZERO.addVisitor<std::reference_wrapper<ReproducingKernel<Dimension>>>                  ([](const any& x) { } );
  ZERO.addVisitor<std::reference_wrapper<PairwiseField<Dimension, Vector>>>              ([](const any& x) { any_cast<reference_wrapper<PairwiseField<Dimension, Vector>>>(x).get().Zero(); });
  ZERO.addVisitor<std::reference_wrapper<PairwiseField<Dimension, pair<Vector, Vector>>>>([](const any& x) { any_cast<reference_wrapper<PairwiseField<Dimension, pair<Vector, Vector>>>>(x).get().Zero(); });
  ZERO.addVisitor<std::reference_wrapper<PairwiseField<Dimension, Scalar>>>              ([](const any& x) { any_cast<reference_wrapper<PairwiseField<Dimension, Scalar>>>(x).get().Zero(); });
  ZERO.addVisitor<std::reference_wrapper<PairwiseField<Dimension, pair<Scalar, Scalar>>>>([](const any& x) { any_cast<reference_wrapper<PairwiseField<Dimension, pair<Scalar, Scalar>>>>(x).get().Zero(); });

  // Walk the state values and zero them
  for (auto itr: mStorage) {
    ZERO.visit(itr.second);
  }

  // Reinitialize the node pair interaction information.
  initializeNodePairInformation();
}

}

