//---------------------------------Spheral++----------------------------------//
// StateDerivatives -- Accumulate and cart around the derivatives/changes to
// the state for a set of physics packages.
//
// Created by JMO, Fri Aug 27 15:50:37 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StateDerivatives_hh__
#define __Spheral_StateDerivatives_hh__

#include "StateBase.hh"
#include "Field/Field.hh"
#include "Field/NodeIteratorBase.hh"

#include <vector>
#include <map>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class DataBase;
template<typename Dimension> class Physics;

template<typename Dimension>
class StateDerivatives: public StateBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Vector3d = typename Dimension::Vector3d;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PackageList = std::vector<Physics<Dimension>*>;
  using PackageIterator = typename PackageList::iterator;

  using KeyType = typename StateBase<Dimension>::KeyType;

  // Constructors, destructor.
  StateDerivatives();
  StateDerivatives(DataBase<Dimension>& dataBase, PackageList& physicsPackage);
  StateDerivatives(DataBase<Dimension>& dataBase,
                   PackageIterator physicsPackageBegin,
                   PackageIterator physicsPackageEnd);
  StateDerivatives(const StateDerivatives& rhs);
  virtual ~StateDerivatives();

  // Assignment.
  StateDerivatives& operator=(const StateDerivatives& rhs);

  // Test if two StateDerivatives have equivalent fields.
  virtual bool operator==(const StateBase<Dimension>& rhs) const override;

  // Methods for setting/interrogating if a given pair of nodes has been 
  // calculated.
  bool nodePairCalculated(const NodeIteratorBase<Dimension>& node1,
                          const NodeIteratorBase<Dimension>& node2) const;
  void flagNodePairCalculated(const NodeIteratorBase<Dimension>& node1,
                              const NodeIteratorBase<Dimension>& node2);
  void initializeNodePairInformation();
  bool calculatedNodePairsSymmetric() const;

  // Convenience bookkeeping methods for maintaining a record of how many significant
  // neighbors a node interacts with.
  int numSignificantNeighbors(const NodeIteratorBase<Dimension>& node) const;
  void incrementSignificantNeighbors(const NodeIteratorBase<Dimension>& node);

  // Force all derivative FieldLists to zero.
  void Zero();

private:
  //--------------------------- Private Interface ---------------------------//
  // Map for storing information about pairs of nodes that have already been
  // calculated.
  using CalculatedPairType = std::map<NodeIteratorBase<Dimension>,
                                      std::vector<NodeIteratorBase<Dimension>>>;
  CalculatedPairType mCalculatedNodePairs;

  // Map for maintaining the number of significant neighbors per node.
  using SignificantNeighborMapType = std::map<NodeIteratorBase<Dimension>, int>;
  SignificantNeighborMapType mNumSignificantNeighbors;

  using StateBase<Dimension>::mFieldStorage;
  using StateBase<Dimension>::mMiscStorage;
};

}

#include "StateDerivativesInline.hh"

#endif

