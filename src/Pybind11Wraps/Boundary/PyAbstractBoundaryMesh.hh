//------------------------------------------------------------------------------
// Trampoline classes for abstract Boundary interface.
// This one supports the Mesh methods.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyAbstractBoundaryMesh__
#define __Spheral_PyAbstractBoundaryMesh__

#include "Geometry/GeomPlane.hh"
#include "NodeList/NodeList.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "Mesh/Mesh.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
using Spheral::MeshSpace::Mesh;

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// PyAbstractBoundaryMesh
//------------------------------------------------------------------------------
template<typename Dimension, class BoundaryBase>
class PyAbstractBoundaryMesh: public BoundaryBase {
public:
  using BoundaryBase::BoundaryBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef Spheral::GeomPlane<Dimension> Plane;

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

  virtual void swapFaceValues(Field<Dimension, std::vector<Scalar>>& field, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,           // Return type
                      BoundaryBase,   // Parent class
                      swapFaceValues, // name of method
                      field, mesh     // arguments
      );
  }

  virtual void swapFaceValues(Field<Dimension, std::vector<Vector>>& field, const Mesh<Dimension>& mesh) const override {
    PYBIND11_OVERLOAD(void,           // Return type
                      BoundaryBase,   // Parent class
                      swapFaceValues, // name of method
                      field, mesh     // arguments
                      );
  }

  virtual bool meshGhostNodes() const override {
    PYBIND11_OVERLOAD(bool,         // Return type
                      BoundaryBase, // Parent class
                      meshGhostNodes// name of method
      );
  }

};

}
}

#endif
