//------------------------------------------------------------------------------
// Trampoline classes for Boundary interface.
// This one inherits from the abstract interface and just overrides the pure
// virtual methods.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyBoundary__
#define __Spheral_PyBoundary__

#include "Geometry/GeomPlane.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
using Spheral::MeshSpace::Mesh;

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// PyBoundary
//------------------------------------------------------------------------------
template<typename Dimension, class BoundaryBase>
class PyBoundary: public BoundaryBase {
public:
  using BoundaryBase::BoundaryBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Spheral::GeomPlane<Dimension> Plane;

  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           setGhostNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           updateGhostNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, int>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, std::vector<Scalar>>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           setViolationNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           updateViolationNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, int>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, Vector>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

};

}
}

#endif
