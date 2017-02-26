#include <vector>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"
#include "Neighbor/GridCellIndex.hh"
#include "Neighbor/GridCellPlane.hh"
#include "Neighbor/Neighbor.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::NeighborSpace;

namespace Spheral {
namespace NeighborSpace {

//------------------------------------------------------------------------------
// PyNeighbor
//------------------------------------------------------------------------------
template<typename Dimension, class NeighborBase>
class PyNeighbor: public NeighborBase {
public:
  using NeighborBase::NeighborBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef Spheral::GeomPlane<Dimension> Plane;

  virtual void setMasterList(int nodeID) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      NeighborBase, // Parent class
                      setMasterList,// name of method
                      nodeID        // arguments
      );
  }

  virtual void setRefineNeighborList(int nodeID) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setRefineNeighborList,// name of method
                      nodeID                // arguments
      );
  }

  virtual void setMasterList(const Vector& position, const Scalar& H) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setMasterList,        // name of method
                      position, H           // arguments
      );
  }

  virtual void setMasterList(const Vector& position, const SymTensor& H) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setMasterList,        // name of method
                      position, H           // arguments
      );
  }

  virtual void setRefineNeighborList(const Vector& position, const Scalar& H) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setRefineNeighborList,// name of method
                      position, H           // arguments
      );
  }

  virtual void setRefineNeighborList(const Vector& position, const SymTensor& H) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setRefineNeighborList,// name of method
                      position, H           // arguments
      );
  }

  virtual void setMasterList(const Vector& position) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setMasterList,        // name of method
                      position              // arguments
      );
  }

  virtual void setRefineNeighborList(const Vector& position) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setRefineNeighborList,// name of method
                      position              // arguments
      );
  }

  virtual void setMasterList(const Plane& enterPlane, const Plane& exitPlane) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      setMasterList,        // name of method
                      enterPlane, exitPlane // arguments
      );
  }

  virtual void updateNodes() override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      updateNodes,        // name of method
      );
  }

  virtual void updateNodes(const std::vector<int>& nodeIDs) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NeighborBase,         // Parent class
                      updateNodes,          // name of method
                      nodeIDs               // arguments
      );
  }

  virtual bool valid() const override {
    PYBIND11_OVERLOAD(bool,                 // Return type
                      NeighborBase,         // Parent class
                      valid                 // name of method
      );
  }
};

}
}

namespace {  // anonymous

//------------------------------------------------------------------------------
// Common virtual methods of Neighbor objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void virtualNeighborBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef Spheral::GeomPlane<Dimension> Plane;

  obj

    // Methods
    .def("setMasterList", (void (Obj::*)(int)) &Obj::setMasterList)
    .def("setMasterList", (void (Obj::*)(const Vector&, const Scalar&)) &Obj::setMasterList, "position"_a, "h"_a)
    .def("setMasterList", (void (Obj::*)(const Vector&, const SymTensor&)) &Obj::setMasterList, "position"_a, "H"_a)
    .def("setMasterList", (void (Obj::*)(const Vector&)) &Obj::setMasterList, "position"_a)
    
    .def("setRefineNeighborList", (void (Obj::*)(int)) &Obj::setRefineNeighborList)
    .def("setRefineNeighborList", (void (Obj::*)(const Vector&, const Scalar&)) &Obj::setRefineNeighborList, "position"_a, "h"_a)
    .def("setRefineNeighborList", (void (Obj::*)(const Vector&, const SymTensor&)) &Obj::setRefineNeighborList, "position"_a, "H"_a)
    .def("setRefineNeighborList", (void (Obj::*)(const Vector&)) &Obj::setRefineNeighborList, "position"_a)

    .def("setMasterList", (void (Obj::*)(const Plane&, const Plane&)) &Obj::setMasterList, "enterPlane"_a, "exitPlane"_a)

    .def("updateNodes", (void (Obj::*)()) &Obj::updateNodes)
    .def("updateNodes", (void (Obj::*)(const std::vector<int>&)) &Obj::updateNodes, "nodeIDs"_a)

    .def("valid", &Obj::valid)
    ;
}

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  using Spheral::NodeSpace::NodeList;

  //............................................................................
  // GridCellIndex
  typedef GridCellIndex<Dimension> GCI;
  py::class_<GCI>(m, ("GridCellIndex" + suffix).c_str())

    // Constructors
    .def(py::init<>())
    .def(py::init<int>(), "xIndex"_a)
    .def(py::init<int, int>(), "xIndex"_a, "yIndex"_a)
    .def(py::init<int, int, int>(), "xIndex"_a, "yIndex"_a, "zIndex"_a)
    .def(py::init<const GCI&>(), "rhs"_a)

    // Methods
    .def("setIndices", (void (GCI::*)(int)) &GCI::setIndices, "xIndex"_a)
    .def("setIndices", (void (GCI::*)(int, int)) &GCI::setIndices, "xIndex"_a, "yIndex"_a)
    .def("setIndices", (void (GCI::*)(int, int, int)) &GCI::setIndices, "xIndex"_a, "yIndex"_a, "zIndex"_a)
    .def("dot", &GCI::dot)
    .def("compare", &GCI::compare)
    .def("inRange", &GCI::inRange)
    .def("magnitude", &GCI::magnitude)
    .def("minElement", &GCI::minElement)
    .def("maxElement", &GCI::maxElement)
    .def("sumElements", &GCI::sumElements)
    .def("productElements", &GCI::productElements)
    .def("indexMin", &GCI::indexMin)
    .def("indexMax", &GCI::indexMax)

    // Operators
    .def("__call__", (int (GCI::*)(int) const) &GCI::operator(), "index"_a, py::is_operator())
    .def(-py::self)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self + int())
    .def(py::self - int())
    .def(py::self * int())
    .def(py::self / int())

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self)
    .def(py::self >  py::self)
    .def(py::self <= py::self)
    .def(py::self >= py::self)

    // Attributes
    .def_property("xIndex", (int (GCI::*)() const) &GCI::xIndex, (void (GCI::*)(int)) &GCI::xIndex)
    .def_property("yIndex", (int (GCI::*)() const) &GCI::yIndex, (void (GCI::*)(int)) &GCI::yIndex)
    .def_property("zIndex", (int (GCI::*)() const) &GCI::zIndex, (void (GCI::*)(int)) &GCI::zIndex)
    ;

  //............................................................................
  // GridCellPlane
  typedef GridCellPlane<Dimension> GCP;
  py::class_<GCP>(m, ("GridCellPlane" + suffix).c_str())

    // Constructors
    .def(py::init<>())
    .def(py::init<const GCP&>(), "rhs"_a)
    .def(py::init<const GCI&, const GCI&>(), "point"_a, "normal"_a)

    // Methods
    .def("minimumDistance", &GCP::minimumDistance)
    .def("coplanar", &GCP::coplanar)
    .def("parallel", &GCP::parallel)
    .def("valid", &GCP::valid)

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self >  GCI())
    .def(py::self <  GCI())
    .def(py::self >= GCI())
    .def(py::self <= GCI())

    // Attributes
    .def_property("point", &GCP::point, &GCP::setPoint)
    .def_property("normal", &GCP::normal, &GCP::setNormal)
    ;

  //............................................................................
  // Neighbor
  typedef Neighbor<Dimension> NT;
  py::class_<NT, PyNeighbor<Dimension, NT>> neighborPB11(m, ("Neighbor" + suffix).c_str());
  virtualNeighborBindings<Dimension, NT>(m, suffix, neighborPB11);
  neighborPB11
    
    // Constructors
    .def(py::init<NodeList<Dimension>&, const NeighborSearchType, const double>(), "nodeList"_a, "searchType"_a, "kernelExtent"_a)

    // Methods
    .def("nodeExtentField", &NT::nodeExtentField)

    ;

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralNeighbor) {
  py::module m("SpheralNeighbor", "Spheral Neighbor module.");

  //............................................................................
  // NeighborSearchType
  py::enum_<Spheral::NeighborSpace::NeighborSearchType>(m, "NeighborSearchType")
    .value("None",          Spheral::NeighborSpace::NeighborSearchType::None)
    .value("Gather",        Spheral::NeighborSpace::NeighborSearchType::Gather)
    .value("Scatter",       Spheral::NeighborSpace::NeighborSearchType::Scatter)
    .value("GatherScatter", Spheral::NeighborSpace::NeighborSearchType::GatherScatter)
    .export_values();

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif
// #ifdef SPHERAL2D
//   dimensionBindings<Spheral::Dim<2>>(m, "2d");
// #endif
// #ifdef SPHERAL3D
//   dimensionBindings<Spheral::Dim<3>>(m, "3d");
// #endif

  return m.ptr();
}
