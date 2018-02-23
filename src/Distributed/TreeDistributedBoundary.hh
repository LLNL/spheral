//---------------------------------Spheral++----------------------------------//
// TreeDistributedBoundary -- Implementation of the Distributed Boundary
// condition for use with TreeNeighbor based NodeLists.
//
// Created by JMO, Mon Aug 27 21:57:51 PDT 2001
//----------------------------------------------------------------------------//
#ifndef __Spheral_TreeDistributedBoundary__
#define __Spheral_TreeDistributedBoundary__

#include "DistributedBoundary.hh"
#include "Neighbor/TreeNeighbor.hh"

#include <vector>

namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class TreeDistributedBoundary: public DistributedBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename NeighborSpace::TreeNeighbor<Dimension>::CellKey CellKey;
  typedef typename DistributedBoundary<Dimension>::DomainBoundaryNodes DomainBoundaryNodes;

  // This method returns the singleton instance of the object.
  static TreeDistributedBoundary& instance();
  static TreeDistributedBoundary* instancePtr();

  // Destructor.
  virtual ~TreeDistributedBoundary();

  //**********************************************************************
  // Set the ghost nodes based on the NodeLists in the given DataBase.
  virtual void setAllGhostNodes(DataBaseSpace::DataBase<Dimension>& dataBase) override;
  //**********************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Singleton instance pointer.
  static TreeDistributedBoundary* mInstance;

  // Disabled methods.
  TreeDistributedBoundary();
  TreeDistributedBoundary(const TreeDistributedBoundary&);
  TreeDistributedBoundary& operator=(const TreeDistributedBoundary&);

  // Private methods.
  const NeighborSpace::TreeNeighbor<Dimension>* getTreeNeighborPtr(const NodeSpace::NodeList<Dimension>* nodeListPtr) const;
  std::vector<std::vector<CellKey>> flattenTrees(const DataBaseSpace::DataBase<Dimension>& dataBase) const;
  void broadcastTree(const std::vector<std::vector<CellKey>>& localTree);
  void buildSendNodes(const DataBaseSpace::DataBase<Dimension>& dataBase,
                      const std::vector<std::vector<CellKey>>& localTree);
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    template<typename Dimension> class TreeDistributedBoundary;
  }
}

#endif
