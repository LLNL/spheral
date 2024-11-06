//---------------------------------Spheral++----------------------------------//
// Boundary -- Abstract base class for the boundary condition classes.
//
// Boundary conditions define three types of nodes:
//   Control nodes -- the real nodes which are being shadowed as ghosts.
//   Ghost nodes -- the indices of the ghost nodes corresponding to the 
//                  controls.
//   Violation nodes -- any internal nodes in the NodeList that our out
//                      of bounds for this boundary condition.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//
// Modified by:
// JMO, Mon Aug 20 23:09:01 PDT 2001
//    Switching to using DataBase and FieldLists.
//----------------------------------------------------------------------------//
#ifndef Boundary_HH
#define Boundary_HH

#include "Field/Field.hh"

#include <vector>
#include <map>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class NodeList;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;
template<typename Dimension> class Mesh;
template<typename Dimension> struct RKCoefficients;

template<typename Dimension>
class Boundary {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // An internal type to hold the paired control/ghost node indices.
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

  //****************************************************************************
  // All Boundary conditions must provide the following methods:
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) = 0;

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) = 0;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) = 0;

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) = 0;

  //****************************************************************************
  // Apply the boundary condition to the ghost node values in the given Field.
  // Catch all method for Fields, assumes copy control->ghost
  virtual void applyGhostBoundary(FieldBaseView<Dimension>& fieldBase) const;                

  // Specific Field ValueTypes -- default to just calling base method
  virtual void applyGhostBoundary(FieldView<Dimension, int>& field) const                       { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, Scalar>& field) const                    { this->applyGhostBoundary(field); }
  virtual void applyGhostBoundary(FieldView<Dimension, Vector>& field) const                    { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, Tensor>& field) const                    { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, SymTensor>& field) const                 { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, ThirdRankTensor>& field) const           { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, FourthRankTensor>& field) const          { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, FifthRankTensor>& field) const           { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, FacetedVolume>& field) const             { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, std::vector<Scalar>>& field) const       { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, std::vector<Vector>>& field) const       { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void applyGhostBoundary(FieldView<Dimension, RKCoefficients<Dimension>>& field) const { this->applyGhostBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  
  //****************************************************************************
  // Apply the boundary condition to the violation node values in the given Field.
  // Catch all method for Fields, assumes noop
  virtual void enforceBoundary(FieldBaseView<Dimension>& ) const {};

  // Specific Field ValueTypes -- default to just calling base method
  virtual void enforceBoundary(FieldView<Dimension, int>& field) const                       { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, Scalar>& field) const                    { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, Vector>& field) const                    { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, Tensor>& field) const                    { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, SymTensor>& field) const                 { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, ThirdRankTensor>& field) const           { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, FourthRankTensor>& field) const          { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, FifthRankTensor>& field) const           { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, FacetedVolume>& field) const             { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, std::vector<Scalar>>& field) const       { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, std::vector<Vector>>& field) const       { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }
  virtual void enforceBoundary(FieldView<Dimension, RKCoefficients<Dimension>>& field) const { this->enforceBoundary(dynamic_cast<FieldBaseView<Dimension>&>(field)); }

  //****************************************************************************
  // Set the ghost and violation nodes based on the FluidNodeLists in the given DataBase.
  virtual void setAllGhostNodes(DataBase<Dimension>& dataBase);
  virtual void setAllViolationNodes(DataBase<Dimension>& dataBase);

  // Apply the boundary condition to the ghost nodes in the given FieldList.
  template<typename DataType>
  void applyFieldListGhostBoundary(FieldList<Dimension, DataType>& fieldList) const;

  // Enforce the boundary condition on the values of the FieldList for the nodes in
  // violation in the given FieldList.
  template<typename DataType>
  void enforceFieldListBoundary(FieldList<Dimension, DataType>& fieldList) const;

  //****************************************************************************
  // Apply the boundary condition to face centered fields on a tessellation.
  virtual void enforceBoundary(std::vector<int>&              /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void enforceBoundary(std::vector<Scalar>&           /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void enforceBoundary(std::vector<Vector>&           /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void enforceBoundary(std::vector<Tensor>&           /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void enforceBoundary(std::vector<SymTensor>&        /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void enforceBoundary(std::vector<ThirdRankTensor>&  /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void enforceBoundary(std::vector<FourthRankTensor>& /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void enforceBoundary(std::vector<FifthRankTensor>&  /*faceField*/, const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }

  // Fill in faces on this boundary with effective opposite face values.
  virtual void swapFaceValues(Field<Dimension, std::vector<Scalar> >& /*field*/,
                              const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }
  virtual void swapFaceValues(Field<Dimension, std::vector<Vector> >& /*field*/,
                              const Mesh<Dimension>& /*mesh*/) const { VERIFY2(false, "Not implemented"); }

  //****************************************************************************
  // Use a set of flags to cull out inactive ghost nodes.
  virtual void cullGhostNodes(const FieldList<Dimension, int>& flagSet,
                              FieldList<Dimension, int>& old2newIndexMap,
                              std::vector<int>& numNodesRemoved);

  // Some boundaries need to know when a problem is starting up and all the physics
  // packages have been initialized.
  virtual void initializeProblemStartup(const bool)                              {};

  // Provide an optional hook that is to be called when all ghost boundaries are
  // to have been set.
  virtual void finalizeGhostBoundary() const                                     {};

  // Overridable hook for clearing out the boundary condition.
  virtual void reset(const DataBase<Dimension>& dataBase);

  // Report the number of ghost nodes in this boundary.
  virtual int numGhostNodes() const;

  // Optionally the boundary can clip an input box range.
  // Defaults to no-op.
  virtual void clip(Vector& /*xmin*/, Vector& /*xmax*/) const                            {};

  // Optionally opt-out of ghost node culling.
  virtual bool allowGhostCulling() const                                         { return true; }

  // Some boundaries have ghosts we should exclude from tessellations.
  // Provide a hook to note such cases.
  virtual bool meshGhostNodes() const                                            { return true; };

  // Allow read access to the map of NodeList->BoundaryNodes.
  const std::map<NodeList<Dimension>*, BoundaryNodes>& boundaryNodeMap() const;

  // Check if we have entries for the given NodeList.
  bool haveNodeList(const NodeList<Dimension>& nodeList) const;

  // Control, Ghost, and Violation nodes for a given NodeList.
  const std::vector<int>& controlNodes(const NodeList<Dimension>& nodeList) const;
  const std::vector<int>& ghostNodes(const NodeList<Dimension>& nodeList) const;
  const std::vector<int>& violationNodes(const NodeList<Dimension>& nodeList) const;

  // Provide iterators over the control, ghost, and violation nodes for a
  // given NodeList.
  std::vector<int>::const_iterator controlBegin(const NodeList<Dimension>& nodeList) const;
  std::vector<int>::const_iterator controlEnd(const NodeList<Dimension>& nodeList) const;

  std::vector<int>::const_iterator ghostBegin(const NodeList<Dimension>& nodeList) const;
  std::vector<int>::const_iterator ghostEnd(const NodeList<Dimension>& nodeList) const;

  std::vector<int>::const_iterator violationBegin(const NodeList<Dimension>& nodeList) const;
  std::vector<int>::const_iterator violationEnd(const NodeList<Dimension>& nodeList) const;

  // Descendent classes are allowed to access the BoundaryNodes for the
  // given NodeList.
  std::map<NodeList<Dimension>*, BoundaryNodes>&       accessBoundaryNodes();
  const std::map<NodeList<Dimension>*, BoundaryNodes>& accessBoundaryNodes() const;

  BoundaryNodes&       accessBoundaryNodes(NodeList<Dimension>& nodeList);
  const BoundaryNodes& accessBoundaryNodes(NodeList<Dimension>& nodeList) const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Descendent classes need to add NodeLists as they build ghosts.
  virtual void addNodeList(NodeList<Dimension>& nodeList);

private:
  //--------------------------- Private Interface ---------------------------//
  std::map<NodeList<Dimension>*, BoundaryNodes> mBoundaryNodes;
};

}

#include "BoundaryInline.hh"

#endif
