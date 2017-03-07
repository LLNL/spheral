//------------------------------------------------------------------------------
// Trampoline classes for abstract Boundary interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyAbstractBoundary__
#define __Spheral_PyAbstractBoundary__

#include "Geometry/GeomPlane.hh"

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// PyAbstractBoundary
//------------------------------------------------------------------------------
template<typename Dimension, class BoundaryBase>
class PyAbstractBoundary: public BoundaryBase {
public:
  using BoundaryBase::BoundaryBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef Spheral::GeomPlane<Dimension> Plane;
  using Spheral::NodeSpace::NodeList;
  using Spheral::FieldSpace::Field;
  using Spheral::FieldSpace::FieldList;

  virtual void setAllGhostNodes(DataBase<Dimension>& dataBase) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      setAllGhostNodes,// name of method
                      dataBase      // arguments
      );
  }

  virtual void setAllViolationNodes(DataBase<Dimension>& dataBase) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      setAllViolationNodes,// name of method
                      dataBase      // arguments
      );
  }

  virtual void cullGhostNodes(const FieldList<Dimension, int>& flagSet,
                              FieldList<Dimension, int>& old2newIndexMap,
                              std::vector<int>& numNodesRemoved) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      cullGhostNodes,// name of method
                      flagSet, old2newIndexMap, numNodesRemoved      // arguments
      );
  }

  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           setGhostNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           updateGhostNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, int>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void applyGhostBoundary(Field<Dimension, std::vector<Scalar>>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           applyGhostBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           setViolationNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           updateViolationNodes,// name of method
                           nodeList      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, int>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, Vector>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const override {
    PYBIND11_OVERLOAD_PURE(void,         // Return type
                           BoundaryBase, // Parent class
                           enforceBoundary,// name of method
                           field      // arguments
      );
  }

  virtual void enforceBoundary(std::vector<int>& faceField, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,            // Return type
                      BoundaryBase,    // Parent class
                      enforceBoundary, // name of method
                      faceField, mesh  // arguments
      );
  }

  virtual void enforceBoundary(std::vector<Scalar>& faceField, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,            // Return type
                      BoundaryBase,    // Parent class
                      enforceBoundary, // name of method
                      faceField, mesh  // arguments
      );
  }

  virtual void enforceBoundary(std::vector<Vector>& faceField, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,            // Return type
                      BoundaryBase,    // Parent class
                      enforceBoundary, // name of method
                      faceField, mesh  // arguments
      );
  }

  virtual void enforceBoundary(std::vector<Tensor>& faceField, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,            // Return type
                      BoundaryBase,    // Parent class
                      enforceBoundary, // name of method
                      faceField, mesh  // arguments
      );
  }

  virtual void enforceBoundary(std::vector<SymTensor>& faceField, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,            // Return type
                      BoundaryBase,    // Parent class
                      enforceBoundary, // name of method
                      faceField, mesh  // arguments
      );
  }

  virtual void enforceBoundary(std::vector<ThirdRankTensor>& faceField, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,            // Return type
                      BoundaryBase,    // Parent class
                      enforceBoundary, // name of method
                      faceField, mesh  // arguments
      );
  }

  virtual void swapFaceValues(Field<Dimension, std::vector<Scalar>>& field, cons Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,           // Return type
                      BoundaryBase,   // Parent class
                      swapFaceValues, // name of method
                      field, mesh     // arguments
      );
  }

  virtual void swapFaceValues(Field<Dimension, std::vector<Vector>>& field, cons Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,           // Return type
                      BoundaryBase,   // Parent class
                      swapFaceValues, // name of method
                      field, mesh     // arguments
                      );
  }

  virtual void initializeProblemStartup() override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      initializeProblemStartup, // name of method
      );
  }

  virtual void finalizeGhostBoundary() const override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      finalizeGhostBoundary,// name of method
      );
  }

  virtual void reset(const DataBase<Dimension>& dataBase) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      reset,        // name of method
                      dataBase      // arguments
      );
  }

  virtual int numGhostNodes() const override {
    PYBIND11_OVERLOAD(int,         // Return type
                      BoundaryBase, // Parent class
                      numGhostNodes,// name of method
      );
  }

  virtual void clip(Vector& xmin, Vector& xmax) const override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      clip,         // name of method
                      xmin, xmax    // arguments
      );
  }

  virtual bool meshGhostNodes() const override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      meshGhostNodes// name of method
      );
  }

};

}
}

#endif
