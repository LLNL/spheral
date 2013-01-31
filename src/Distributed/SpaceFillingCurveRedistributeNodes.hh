//---------------------------------Spheral++----------------------------------//
// SpaceFillingCurveRedistributeNodes
//
// This is an abstract base for the space filling curve family of 
// repartitioners.  The assumption is that the descendent classes will provide
// the computeHashedIndicies method to assign unique keys to each point in the
// order that that algorithm wants the points distributed.
//
// Created by JMO, Wed Apr  9 13:13:46 PDT 2008
//----------------------------------------------------------------------------//
#ifndef SpaceFillingCurveRedistributeNodes_HH
#define SpaceFillingCurveRedistributeNodes_HH

#ifndef __GCCXML__

#include <vector>
#include <map>
#ifdef USE_MPI
#include "mpi.h"
#endif

#else

#include "fakestl.hh"

#endif

#include "RedistributeNodes.hh"
#include "Utilities/KeyTraits.hh"

namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace BoundarySpace {
    template<typename Dimension> class Boundary;
  }
}

namespace Spheral {
namespace PartitionSpace {

template<typename Dimension>
class SpaceFillingCurveRedistributeNodes: public RedistributeNodes<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef KeyTraits::Key Key;

  // Constructors
  SpaceFillingCurveRedistributeNodes(const double dummy,
                                     const double minNodesPerDomainFraction,
                                     const double maxNodesPerDomainFraction,
                                     const bool workBalance = true,
                                     const bool localReorderOnly = false);

  // Destructor
  virtual ~SpaceFillingCurveRedistributeNodes();

  // This is the required method for all descendant classes.
  virtual
  FieldSpace::FieldList<Dimension, Key> 
  computeHashedIndicies(const DataBaseSpace::DataBase<Dimension>& dataBase) const = 0;

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  // This is the method required of all descendent classes.
  virtual void redistributeNodes(DataBaseSpace::DataBase<Dimension>& dataBase,
                                 std::vector<BoundarySpace::Boundary<Dimension>*> boundaries =
                                 std::vector<BoundarySpace::Boundary<Dimension>*>());

  // Compute the cell size in each dimension.
  Vector computeStepSize(const std::pair<Vector, Vector>& box) const;
  
  // Stitch together the given indicies and DomainNode list.
  // This returns the set sorted by the index.
  std::vector<std::pair<Key, DomainNode<Dimension> > >
  buildIndex2IDPairs(const FieldSpace::FieldList<Dimension, Key>& indicies,
                     const std::vector<DomainNode<Dimension> >& domainNodes) const;

  // Find the hashed index the given amount of work above the specified lower bound.
  Key findUpperKey(const std::vector<Key>& indicies,
                   const std::vector<int>& count,
                   const std::vector<Scalar>& work,
                   const Key lowerBound,
                   const Key maxUpperBound,
                   const Scalar workTarget,
                   const int minNodes,
                   const int maxNodes,
                   Key& upperKey,
                   int& numNodes) const;

  // Compute the global number of nodes in the given index range.
  int numIndiciesInRange(const std::vector<Key>& indicies,
                         const std::vector<int>& count,
                         const Key lowerBound,
                         const Key upperBound) const;

  // Compute the global work for the nodes in the given index range.
  Scalar workInRange(const std::vector<Key>& indicies,
                     const std::vector<Scalar>& work,
                     const Key lowerBound,
                     const Key upperBound) const;

  // Combines the above.
  void workAndNodesInRange(const std::vector<Key>& indicies,
                           const std::vector<int>& count,
                           const std::vector<Scalar>& work,
                           const Key lowerBound,
                           const Key upperBound,
                           int& countInRange,
                           Scalar& workInRange) const;

  // Compute the number of nodes we want per process.
  int targetNumNodes(const int numGlobal,
                     const int numProcs,
                     const int targetProc) const;

  // Get the next (global) index following the given value.
  Key findNextIndex(const std::vector<Key>& indicies,
                    const Key index,
                    const Key maxIndex) const;

  // Find the domain for the given index given the set of index ranges for 
  // each processors.
  int domainForIndex(const Key index,
                     const std::vector<std::pair<Key, Key> >& indexRanges) const;

  // The allowed range of nodes per domain as a fraction of the average.
  double minNodesPerDomainFraction() const;
  double maxNodesPerDomainFraction() const;

  void minNodesPerDomainFraction(const double x);
  void maxNodesPerDomainFraction(const double x);

  // Flag for whether we should compute the work per node or strictly balance by
  // node count.
  bool workBalance() const;
  void workBalance(const bool val);

  // Flag that will cause us not to repartition between domains, but only sort locally on 
  // each domain.
  bool localReorderOnly() const;
  void localReorderOnly(const bool val);

private:
  //--------------------------- Private Interface ---------------------------//
  double mMinNodesPerDomainFraction, mMaxNodesPerDomainFraction;
  bool mWorkBalance;
  bool mLocalReorderOnly;

  // No copy or assignment operations.
  SpaceFillingCurveRedistributeNodes(const SpaceFillingCurveRedistributeNodes& nodes);
  SpaceFillingCurveRedistributeNodes& operator=(const SpaceFillingCurveRedistributeNodes& rhs);

};

}
}

#else
// Forward declare the SpaceFillingCurveRedistributeNodes class.
namespace Spheral {
  namespace PartitionSpace {
    template<typename Dimension> class SpaceFillingCurveRedistributeNodes;
  }
}

#endif
