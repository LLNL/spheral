//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_Neighbor_hh__
#define __Spheral_Neighbor_hh__

#ifndef __GCCXML__
#include <vector>
#include "Field/Field.hh"
#else
#include "fakestl.hh"
#endif

#include "Geometry/Dimension.hh"

namespace Spheral {
  template<typename Dimension> class GeomPlane;
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
  }
}

namespace Spheral {
namespace NeighborSpace {

enum NeighborSearchType {
  None = 0,
  Gather = 1,
  Scatter = 2,
  GatherScatter = 3
};

template<typename Dimension>
class Neighbor {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef std::vector<int>::iterator iterator;
  typedef std::vector<int>::const_iterator const_iterator;

  // Constructors and destructors
  Neighbor(NodeSpace::NodeList<Dimension>& nodeList, 
           const NeighborSearchType searchType,
           const double kernelExtent);
  virtual ~Neighbor();

  // Choose the type of neighbor search we wish to use.
  NeighborSearchType neighborSearchType() const;
  void neighborSearchType(NeighborSearchType searchType);

  // All neighboring classes need to now how far to sample.
  double kernelExtent() const;
  void kernelExtent(double kernelExtent);

  // Determine the maximum extent of a given H smoothing scale along the
  // Cartesian axes.
  static Vector HExtent(const Scalar& H, const double kernelExtent);
  static Vector HExtent(const SymTensor& H, const double kernelExtent);

  // Allow access to the field of node extents.
  const FieldSpace::Field<Dimension, Vector>& nodeExtentField() const;

  // Return the current number of (master, coarse, fine) nodes.
  int numMaster() const;
  int numCoarse() const;
  int numRefine() const;

  // Return the current list of potential neighbor indicies.
  const std::vector<int>& masterList() const;
  const std::vector<int>& coarseNeighborList() const;
  const std::vector<int>& refineNeighborList() const;

  // Provide iterators to go over the neighbor indicies.
  const_iterator masterBegin() const;
  const_iterator masterEnd() const;
  const_iterator coarseNeighborBegin() const;
  const_iterator coarseNeighborEnd() const;
  const_iterator refineNeighborBegin() const;
  const_iterator refineNeighborEnd() const;

  // Access the node list.
  const NodeSpace::NodeList<Dimension>& nodeList() const;
  void nodeList(NodeSpace::NodeList<Dimension>& nodeList);

  const NodeSpace::NodeList<Dimension>* nodeListPtr() const;
  void nodeListPtr(NodeSpace::NodeList<Dimension>* nodeListPtr);

  void unregisterNodeList();

  // Determine the node extent for an individual node.
  Vector nodeExtent(int nodeID) const;

  // Force the node extent field to be computed.
  void setNodeExtents();
  void setNodeExtents(const std::vector<int>& nodeIDs);
  void setInternalNodeExtents();
  void setGhostNodeExtents();

  // Set or refine the neighbor lists for a given node ID.
  virtual void setMasterList(int nodeID);
  virtual void setRefineNeighborList(int nodeID);

  // Helper method to cull lists of neighbors based on min/max positions and 
  // min/max extents.
  std::vector<int> 
  precullList(const Vector& minMasterPosition, const Vector& maxMasterPostion,
              const Vector& minMasterExtent, const Vector& maxMasterExtent,
              const std::vector<int>& coarseList) const;

  // Cull the local (to this NodeList) neighbor info based on the current master
  // state.
  // *NOTE* -- this is not safe to do when you want to use this neighbor info 
  // with different NodeLists!
  void precullForLocalNodeList();

  // ********** Descendent Neighbor types must provide these methods. **********
  // Set or refine the neighbor lists for the given position and smoothing 
  // scale.
  virtual void setMasterList(const Vector& position,
                             const Scalar& H) = 0;
  virtual void setMasterList(const Vector& position,
                             const SymTensor& H) = 0;

  virtual void setRefineNeighborList(const Vector& position,
                                     const Scalar& H) = 0;
  virtual void setRefineNeighborList(const Vector& position,
                                     const SymTensor& H) = 0;

  // Set Neighbors for the given position.
  virtual void setMasterList(const Vector& position) = 0;
  virtual void setRefineNeighborList(const Vector& position) = 0;

  // Set the neighbor lists based on proximity to planes.
  virtual void setMasterList(const GeomPlane<Dimension>& enterPlane,
                             const GeomPlane<Dimension>& exitPlane) = 0;

  // Force the update of internal data for the NodeList.
  virtual void updateNodes() = 0;
  virtual void updateNodes(const std::vector<int>& nodeIDs) = 0;
  //****************************************************************************

  // Determine if the Neighbor is in a valid, ready to use state.
  virtual bool valid() const;

  // A static method for set master/coarse neighbor sets on multiple NodeLists
  // with culling of the coarse set.  This method is static because it is 
  // intended to deal with multiple Neighbor objects (not just this particular
  // one), and must access private data of the Neighbors.
  template<typename NodeListIteratorType>
  static void setMasterNeighborGroup(const Vector& position,
                                     const SymTensor& H,
                                     const NodeListIteratorType& nodeListBegin,
                                     const NodeListIteratorType& nodeListEnd,
                                     const double kernelExtent);

protected:
  //-------------------------- Protected Interface --------------------------//
  // Provide read/write access to the node index vectors for descendent classes.
  std::vector<int>& accessMasterList();
  std::vector<int>& accessCoarseNeighborList();
  std::vector<int>& accessRefineNeighborList();
  FieldSpace::Field<Dimension, Vector>& accessNodeExtentField();

private:
  //--------------------------- Private Interface ---------------------------//
  NeighborSearchType mSearchType;
  double mKernelExtent;

  std::vector<int>* mMasterListPtr;
  std::vector<int>* mCoarseNeighborListPtr;
  std::vector<int>* mRefineNeighborListPtr;

#ifndef __GCCXML__
  NodeSpace::NodeList<Dimension>* mNodeListPtr;
  FieldSpace::Field<Dimension, Vector> mNodeExtent;
#endif
};

// We explicitly specialize the HExtent method for 1, 2, & 3 dimensions.
template<> Dim<1>::Vector Neighbor< Dim<1> >::HExtent(const Dim<1>::SymTensor&, const double kernelExtent);
template<> Dim<2>::Vector Neighbor< Dim<2> >::HExtent(const Dim<2>::SymTensor&, const double kernelExtent);
template<> Dim<3>::Vector Neighbor< Dim<3> >::HExtent(const Dim<3>::SymTensor&, const double kernelExtent);

}
}

#ifndef __GCCXML__
#include "NeighborInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace NeighborSpace {
    template<typename Dimension> class Neighbor;
  }
}

#endif
