//------------------------------------------------------------------------------
// Trampoline classes for the Boundary types.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyBoundary__
#define __Spheral_PyBoundary__

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
  typedef typename Dimension::SymTensor SymTensor;
  typedef Spheral::GeomPlane<Dimension> Plane;

  virtual void setMasterList(int nodeID) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      BoundaryBase, // Parent class
                      setMasterList,// name of method
                      nodeID        // arguments
      );
  }

  virtual void setRefineBoundaryList(int nodeID) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      BoundaryBase,         // Parent class
                      setRefineBoundaryList,// name of method
                      nodeID                // arguments
      );
  }

  virtual void setMasterList(const Vector& position, const Scalar& H) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           setMasterList,        // name of method
                           position, H           // arguments
      );
  }

  virtual void setMasterList(const Vector& position, const SymTensor& H) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           setMasterList,        // name of method
                           position, H           // arguments
      );
  }

  virtual void setRefineBoundaryList(const Vector& position, const Scalar& H) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           setRefineBoundaryList,// name of method
                           position, H           // arguments
      );
  }

  virtual void setRefineBoundaryList(const Vector& position, const SymTensor& H) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           setRefineBoundaryList,// name of method
                           position, H           // arguments
      );
  }

  virtual void setMasterList(const Vector& position) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           setMasterList,        // name of method
                           position              // arguments
      );
  }

  virtual void setRefineBoundaryList(const Vector& position) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           setRefineBoundaryList,// name of method
                           position              // arguments
      );
  }

  virtual void setMasterList(const Plane& enterPlane, const Plane& exitPlane) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           setMasterList,        // name of method
                           enterPlane, exitPlane // arguments
      );
  }

  virtual void updateNodes() override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           updateNodes,        // name of method
      );
  }

  virtual void updateNodes(const std::vector<int>& nodeIDs) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           BoundaryBase,         // Parent class
                           updateNodes,          // name of method
                           nodeIDs               // arguments
      );
  }

  virtual bool valid() const override {
    PYBIND11_OVERLOAD(bool,                 // Return type
                      BoundaryBase,         // Parent class
                      valid                 // name of method
      );
  }
};

}
}

#endif
