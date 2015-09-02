//---------------------------------Spheral++----------------------------------//
// ConstantBoundary -- Take a snapshot of the state of a set of nodes on 
// NodeList, and create ghost nodes with that state from that time on.
//
// Created by JMO, Fri Nov 21 13:57:34 PST 2003
//
// Modified by:
//----------------------------------------------------------------------------//
#ifndef ConstantBoundary_HH
#define ConstantBoundary_HH

#include "Boundary.hh"
#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeList.hh"
#include "DataBase/StateBase.hh" // For constructing Field keys.

namespace Spheral {
  namespace NodeSpace {
    template<typename Dimension> class NodeList;
  }
  namespace FieldSpace {
    template<typename Dimension> class FieldBase;
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
class ConstantBoundary: public Boundary<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename StateBase<Dimension>::KeyType KeyType;

  // Constructors and destructors.
  ConstantBoundary(const NodeSpace::NodeList<Dimension>& nodeList,
                   const std::vector<int>& nodeIDs,
                   const GeomPlane<Dimension>& denialPlane);
  virtual ~ConstantBoundary();

  //**********************************************************************
  // All Boundary conditions must provide the following methods:
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the ghost node values in the given Field.
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, std::vector<Scalar> >& field) const;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeSpace::NodeList<Dimension>& nodeList);

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(FieldSpace::Field<Dimension, int>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Scalar>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;
  //**********************************************************************

  // After physics have been initialized we take a snapshot of the node state.
  virtual void initializeProblemStartup();

  // Minimal valid test.
  virtual bool valid() const;

  // Accessor methods.
  std::vector<int> nodeIndices() const;
  int numConstantNodes() const;
  const NodeSpace::NodeList<Dimension>& nodeList() const;
  Tensor reflectOperator() const;

private:
  //--------------------------- Private Interface ---------------------------//
  const NodeSpace::NodeList<Dimension>* mNodeListPtr;
  std::vector<int> mNodeFlags;
  int mNumConstantNodes;
  GeomPlane<Dimension> mDenialPlane;
  Tensor mReflectOperator;

  std::map<KeyType, std::vector<int> > mIntValues;
  std::map<KeyType, std::vector<Scalar> > mScalarValues;
  std::map<KeyType, std::vector<Vector> > mVectorValues;
  std::map<KeyType, std::vector<Tensor> > mTensorValues;
  std::map<KeyType, std::vector<SymTensor> > mSymTensorValues;
  std::map<KeyType, std::vector<ThirdRankTensor> > mThirdRankTensorValues;
  std::map<KeyType, std::vector<std::vector<Scalar> > > mVectorScalarValues;

  // No default or copy constructors.
  ConstantBoundary();
  ConstantBoundary(ConstantBoundary&);
};

}
}

#include "ConstantBoundaryInline.hh"

#else

// Forward declaration.
namespace Spheral {
  namespace ConstantBoundarySpace {
    template<typename Dimension> class ConstantBoundary;
  }
}

#endif
