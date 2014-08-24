//---------------------------------Spheral++----------------------------------//
// StateDerivatives -- Accumulate and cart around the derivatives/changes to
// the state for a set of physics packages.
//
// Created by JMO, Fri Aug 27 15:50:37 2004
//----------------------------------------------------------------------------//
#ifndef __Spheral_StateDerivatives_hh__
#define __Spheral_StateDerivatives_hh__

#include <vector>
#include <map>

#include "StateBase.hh"
#include "Field/Field.hh"
#include "Field/NodeIteratorBase.hh"

namespace Spheral {

// Forward declarations.
namespace DataBaseSpace {
  template<typename Dimension> class DataBase;
}
namespace PhysicsSpace {
  template<typename Dimension> class Physics;
}

template<typename Dimension>
class StateDerivatives: public StateBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  // Useful typedefs
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Vector3d Vector3d;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef std::vector<PhysicsSpace::Physics<Dimension>*> PackageList;
  typedef typename PackageList::iterator PackageIterator;

  typedef typename StateBase<Dimension>::KeyType KeyType;

  // Constructors, destructor.
  StateDerivatives();
  StateDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase, PackageList& physicsPackage);
  StateDerivatives(DataBaseSpace::DataBase<Dimension>& dataBase,
                   PackageIterator physicsPackageBegin,
                   PackageIterator physicsPackageEnd);
  StateDerivatives(const StateDerivatives& rhs);
  virtual ~StateDerivatives();

  // Assignment.
  StateDerivatives& operator=(const StateDerivatives& rhs);

  // Test if two StateDerivatives have equivalent fields.
  virtual bool operator==(const StateBase<Dimension>& rhs) const;

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
  typedef std::map<NodeIteratorBase<Dimension>,
                   std::vector<NodeIteratorBase<Dimension> > > CalculatedPairType;
  CalculatedPairType mCalculatedNodePairs;

  // Map for maintaining the number of significant neighbors per node.
  typedef std::map<NodeIteratorBase<Dimension>, int> SignificantNeighborMapType;

  SignificantNeighborMapType mNumSignificantNeighbors;

  using typename StateBase<Dimension>::StorageType;
  using StateBase<Dimension>::mStorage;
};

}

#ifndef __GCCXML__
#include "StateDerivativesInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class StateDerivatives;
}

#endif

