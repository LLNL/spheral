//---------------------------------Spheral++----------------------------------//
// Boundary -- Abstract base class for the boundary condition classes.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//
// Modified by:
// JMO, Mon Aug 20 23:09:01 PDT 2001
//    Switching to using DataBase and FieldLists.
//----------------------------------------------------------------------------//
#ifndef Boundary_HH
#define Boundary_HH

#include <vector>
#include <map>
#include "boost/shared_ptr.hpp"

namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class Field;
    template<typename Dimension, typename DataType> class FieldList;
  }
  namespace DataBaseSpace {
    template<typename Dimension> class DataBase;
  }
}

namespace Spheral {
namespace BoundarySpace {

template<typename Dimension>
class Boundary {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Vector3d Vector3d;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // An internal type to hold the paired control/ghost node indicies.
  // Also maintains a list of any internal nodes that are in violation
  // of the boundary condition.
  struct BoundaryNodes {
    std::vector<int> controlNodes;
    std::vector<int> ghostNodes;
    std::vector<int> violationNodes;
  };

  // Constructors and destructors.
  Boundary();
  virtual ~Boundary();

  // Boundary conditions define three types of nodes:
  //   Control nodes -- the real nodes which are being shadowed as ghosts.
  //   Ghost nodes -- the indicies of the ghost nodes corresponding to the 
  //                  controls.
  //   Violation nodes -- any internal nodes in the NodeList that our out
  //                      of bounds for this boundary condition.

#ifndef __GCCXML__
  // Allow read access to the map of NodeList->BoundaryNodes.
  const std::map<NodeSpace::NodeList<Dimension>*, BoundaryNodes>& boundaryNodeMap() const;
#endif

  // Check if we have entries for the given NodeList.
  bool haveNodeList(const NodeSpace::NodeList<Dimension>& nodeList) const;

  // Control, Ghost, and Violation nodes for a given NodeList.
  const std::vector<int>& controlNodes(const NodeSpace::NodeList<Dimension>& nodeList) const;
  const std::vector<int>& ghostNodes(const NodeSpace::NodeList<Dimension>& nodeList) const;
  const std::vector<int>& violationNodes(const NodeSpace::NodeList<Dimension>& nodeList) const;

  // Provide iterators over the control, ghost, and violation nodes for a
  // given NodeList.
  std::vector<int>::const_iterator controlBegin(const NodeSpace::NodeList<Dimension>& nodeList) const;
  std::vector<int>::const_iterator controlEnd(const NodeSpace::NodeList<Dimension>& nodeList) const;

  std::vector<int>::const_iterator ghostBegin(const NodeSpace::NodeList<Dimension>& nodeList) const;
  std::vector<int>::const_iterator ghostEnd(const NodeSpace::NodeList<Dimension>& nodeList) const;

  std::vector<int>::const_iterator violationBegin(const NodeSpace::NodeList<Dimension>& nodeList) const;
  std::vector<int>::const_iterator violationEnd(const NodeSpace::NodeList<Dimension>& nodeList) const;

  // Set the ghost nodes based on the FluidNodeLists in the given DataBase.
  virtual void setAllGhostNodes(DataBaseSpace::DataBase<Dimension>& dataBase);

  // Apply the boundary condition to the ghost nodes in the given FieldList.
  template<typename DataType>
  void applyFieldListGhostBoundary(FieldSpace::FieldList<Dimension, DataType>& fieldList) const;

  // Select any nodes based in the FluidNodeLists in the given DataBase that
  // are in violation of boundary condition.
  virtual void setAllViolationNodes(DataBaseSpace::DataBase<Dimension>& dataBase);

  // Enforce the boundary condition on the values of the FieldList for the nodes in
  // violation in the given FieldList.
  template<typename DataType>
  void enforceFieldListBoundary(FieldSpace::FieldList<Dimension, DataType>& fieldList) const;

  // Use a set of flags to cull out inactive ghost nodes.
  virtual void cullGhostNodes(const FieldSpace::FieldList<Dimension, int>& flagSet,
                              FieldSpace::FieldList<Dimension, int>& old2newIndexMap,
                              std::vector<int>& numNodesRemoved);

  //**********************************************************************
  // All Boundary conditions must provide the following methods:
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList) = 0;

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList) = 0;

  // Apply the boundary condition to the ghost node values in the given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const = 0;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const = 0;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const = 0;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector3d>& field) const = 0;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const = 0;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const = 0;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const = 0;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Scalar> >& field) const = 0;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList) = 0;

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList) = 0;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const = 0;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const = 0;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const = 0;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector3d>& field) const = 0;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const = 0;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const = 0;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const = 0;
  //**********************************************************************

  // Provide an optional hook that is to be called when all ghost boundaries are
  // to have been set.
  virtual void finalizeGhostBoundary() const {};

  // Overridable hook for clearing out the boundary condition.
  virtual void reset(const DataBaseSpace::DataBase<Dimension>& dataBase);

  // Report the number of ghost nodes in this boundary.
  virtual int numGhostNodes() const;

#ifndef __GCCXML__
  // protected:
  //--------------------------- Protected Interface ---------------------------//
  // Descendent classes are allowed to access the BoundaryNodes for the
  // given NodeList.
  std::map<NodeSpace::NodeList<Dimension>*, BoundaryNodes>& accessBoundaryNodes();
  BoundaryNodes& accessBoundaryNodes(NodeSpace::NodeList<Dimension>& nodeList);

  virtual void addNodeList(NodeSpace::NodeList<Dimension>& nodeList);

private:
  //--------------------------- Private Interface ---------------------------//
  std::map<NodeSpace::NodeList<Dimension>*, BoundaryNodes> mBoundaryNodes;
#endif
};

}
}

#ifndef __GCCXML__
#include "BoundaryInline.hh"
#endif

#else

namespace Spheral {
  namespace BoundarySpace {
    // Forward declaration.
    template<typename Dimension> class Boundary;
  }
}

#endif
