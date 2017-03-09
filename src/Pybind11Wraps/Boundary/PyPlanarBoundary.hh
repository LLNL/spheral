//------------------------------------------------------------------------------
// Trampoline classes for PlanarBoundary interface.
// This one inherits from the abstract interface and just overrides the pure
// virtual methods.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyPlanarBoundary__
#define __Spheral_PyPlanarBoundary__

#include "Geometry/GeomPlane.hh"

using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
using Spheral::MeshSpace::Mesh;

namespace Spheral {
namespace BoundarySpace {

//------------------------------------------------------------------------------
// PyPlanarBoundary
//------------------------------------------------------------------------------
template<typename Dimension, class BoundaryBase>
class PyPlanarBoundary: public BoundaryBase {
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

};

}
}

#endif
