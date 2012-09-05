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
#include "NodeList/NodeList.hh"

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
  typedef typename Dimension::Vector3d Vector3d;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Constructors and destructors.
  ConstantBoundary(const NodeSpace::NodeList<Dimension>& nodeList,
                   const std::vector<int>& nodeIDs);
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
  virtual void applyGhostBoundary(FieldSpace::Field<Dimension, Vector3d>& field) const;
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
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Vector3d>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, Tensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, SymTensor>& field) const;
  virtual void enforceBoundary(FieldSpace::Field<Dimension, ThirdRankTensor>& field) const;
  //**********************************************************************

  // Minimal valid test.
  virtual bool valid() const;

  // Accessor methods.
  int numConstantNodes() const;
  const NodeSpace::NodeList<Dimension>& nodeList() const;

#ifndef __GCCXML__
  // Helper method for inserting Field values in the internal map data 
  // structures.
  template<typename DataType>
  void storeFieldValues(const NodeSpace::NodeList<Dimension>& nodeList,
                        const std::vector<int>& nodeIDs,
                        std::map<const FieldSpace::FieldBase<Dimension>*, 
                                 std::vector<DataType> >& values) const;

  // Set the ghost values in the given field using the given map.
  template<typename DataType>
  void setGhostValues(FieldSpace::Field<Dimension, DataType>& field,
                      const std::map<const FieldSpace::FieldBase<Dimension>*, 
                                     std::vector<DataType> >& values) const;
#endif

private:
  //--------------------------- Private Interface ---------------------------//
#ifndef __GCCXML__
  int mNumConstantNodes;
  const NodeSpace::NodeList<Dimension>* mNodeListPtr;

  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<int> > mIntValues;
  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<Scalar> > mScalarValues;
  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<Vector> > mVectorValues;
  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<Vector3d> > mVector3dValues;
  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<Tensor> > mTensorValues;
  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<SymTensor> > mSymTensorValues;
  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<ThirdRankTensor> > mThirdRankTensorValues;
  std::map<const FieldSpace::FieldBase<Dimension>*, std::vector<std::vector<Scalar> > > mVectorScalarValues;
#endif

  // No default or copy constructors.
  ConstantBoundary();
  ConstantBoundary(ConstantBoundary&);
};

}
}

#ifndef __GCCXML__
#include "ConstantBoundaryInline.hh"
#endif

#else

// Forward declaration.
namespace Spheral {
  namespace ConstantBoundarySpace {
    template<typename Dimension> class ConstantBoundary;
  }
}

#endif
