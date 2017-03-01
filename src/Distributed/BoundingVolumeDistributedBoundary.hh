//---------------------------------Spheral++----------------------------------//
// BoundingVolumeDistributedBoundary
//
// Build a distributed boundary based on testing for intersecting bounding 
// volumes of domains.
//
// Created by JMO, Tue Jan 19 09:22:37 PST 2010
//----------------------------------------------------------------------------//

#ifndef BoundingVolumeDistributedBoundary_HH
#define BoundingVolumeDistributedBoundary_HH

#include <string>

#include "DistributedBoundary.hh"

namespace Spheral {
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
  namespace NeighborSpace {
    template<typename Dimension> class BoundingVolumeNeighbor;
    template<typename Dimension> class GridCellIndex;
  }
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class BoundingVolumeDistributedBoundary: public DistributedBoundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename DistributedBoundary<Dimension>::DomainBoundaryNodes DomainBoundaryNodes;

  // This method returns the singleton instance of the object.
  static BoundingVolumeDistributedBoundary& instance();
  static BoundingVolumeDistributedBoundary* instancePtr();

  // Destructor.
  virtual ~BoundingVolumeDistributedBoundary();

  //**********************************************************************
  // Apply the boundary condition to the given Field.
  // Set the ghost nodes based on the NodeLists in the given DataBase.
  virtual void setAllGhostNodes(DataBaseSpace::DataBase<Dimension>& dataBase);
  //**********************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // Singleton instance pointer.
  static BoundingVolumeDistributedBoundary* mInstance;

  // Disabled methods.
  BoundingVolumeDistributedBoundary();
  BoundingVolumeDistributedBoundary(const BoundingVolumeDistributedBoundary&);
  BoundingVolumeDistributedBoundary& operator=(const BoundingVolumeDistributedBoundary&);

  // Private methods.
  void buildSendNodes(const DataBaseSpace::DataBase<Dimension>& dataBase);
  void packNodeListBuffers(const DataBaseSpace::DataBase<Dimension>& dataBase,
                           std::vector<int>& numNodesPerNodes,
                           std::vector<std::string>& positionBuffers,
                           std::vector<std::string>& Hbuffers) const;
};

}
}

#else

// Forward declaration.
namespace Spheral {
  namespace BoundarySpace {
    template<typename Dimension> class BoundingVolumeDistributedBoundary;
  }
}

#endif
