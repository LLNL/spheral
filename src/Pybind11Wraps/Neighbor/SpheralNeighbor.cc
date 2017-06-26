// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

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

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::NeighborSpace::GridCellIndex<Spheral::Dim<1>>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::NeighborSpace::GridCellIndex<Spheral::Dim<1>>>>);

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::NeighborSpace::GridCellIndex<Spheral::Dim<2>>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::NeighborSpace::GridCellIndex<Spheral::Dim<2>>>>);

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::NeighborSpace::GridCellIndex<Spheral::Dim<3>>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::NeighborSpace::GridCellIndex<Spheral::Dim<3>>>>);

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
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           setMasterList,        // name of method
                           position, H           // arguments
      );
  }

  virtual void setMasterList(const Vector& position, const SymTensor& H) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           setMasterList,        // name of method
                           position, H           // arguments
      );
  }

  virtual void setRefineNeighborList(const Vector& position, const Scalar& H) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           setRefineNeighborList,// name of method
                           position, H           // arguments
      );
  }

  virtual void setRefineNeighborList(const Vector& position, const SymTensor& H) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           setRefineNeighborList,// name of method
                           position, H           // arguments
      );
  }

  virtual void setMasterList(const Vector& position) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           setMasterList,        // name of method
                           position              // arguments
      );
  }

  virtual void setRefineNeighborList(const Vector& position) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           setRefineNeighborList,// name of method
                           position              // arguments
      );
  }

  virtual void setMasterList(const Plane& enterPlane, const Plane& exitPlane) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           setMasterList,        // name of method
                           enterPlane, exitPlane // arguments
      );
  }

  virtual void updateNodes() override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           NeighborBase,         // Parent class
                           updateNodes,        // name of method
      );
  }

  virtual void updateNodes(const std::vector<int>& nodeIDs) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
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
    .def("setMasterList", (void (Obj::*)(const Vector&, const Scalar&)) &Obj::setMasterList, "position"_a, "h"_a)
    .def("setMasterList", (void (Obj::*)(const Vector&, const SymTensor&)) &Obj::setMasterList, "position"_a, "H"_a)
    .def("setMasterList", (void (Obj::*)(const Vector&)) &Obj::setMasterList, "position"_a)
    
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
  typedef Spheral::GeomPlane<Dimension> Plane;
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
    .def("nodeList", (const NodeList<Dimension>& (NT::*)() const) &NT::nodeList)
    .def("nodeList", (void (NT::*)(NodeList<Dimension>&)) &NT::nodeList)
    .def("unregisterNodeList", &NT::unregisterNodeList)
    .def("nodeExtent", &NT::nodeExtent)
    .def("setNodeExtents", (void (NT::*)()) &NT::setNodeExtents)
    .def("setNodeExtents", (void (NT::*)(const std::vector<int>&)) &NT::setNodeExtents)
    .def("setInternalNodeExtents", &NT::setInternalNodeExtents)
    .def("setGhostNodeExtents", &NT::setGhostNodeExtents)
    .def("precullList", &NT::precullList)
    .def("precullForLocalNodeList", &NT::precullForLocalNodeList)
    .def_static("HExtent", (Vector (*)(const Scalar&, const double)) &NT::HExtent) //, "H"_a, "kernelExtent"_a)
    .def_static("HExtent", (Vector (*)(const SymTensor&, const double)) &NT::HExtent, "H"_a, "kernelExtent"_a)

    // Virtual methods
    .def("setMasterList", (void (NT::*)(int)) &NT::setMasterList, "nodeID"_a)
    .def("setRefineNeighborList", (void (NT::*)(int)) &NT::setRefineNeighborList, "nodeID"_a)

    // Attributes
    .def_property("neighborSearchType",
                  (NeighborSearchType (NT::*)() const) &NT::neighborSearchType,
                  (void (NT::*)(NeighborSearchType)) &NT::neighborSearchType)
    .def_property("kernelExtent",
                  (double (NT::*)() const) &NT::kernelExtent,
                  (void (NT::*)(double)) &NT::kernelExtent)
    .def_property_readonly("numMaster", &NT::numMaster)
    .def_property_readonly("numCoarse", &NT::numCoarse)
    .def_property_readonly("numRefine", &NT::numRefine)
    .def_property_readonly("masterList", &NT::masterList)
    .def_property_readonly("coarseNeighborList", &NT::coarseNeighborList)
    .def_property_readonly("refineNeighborList", &NT::refineNeighborList)
    ;

  //............................................................................
  // NestedGridNeighbor
  typedef NestedGridNeighbor<Dimension> NGT;
  py::class_<NGT, NT, PyNeighbor<Dimension, NGT>> nestedgridneighborPB11(m, ("NestedGridNeighbor" + suffix).c_str());
  virtualNeighborBindings<Dimension, NGT>(m, suffix, nestedgridneighborPB11);
  nestedgridneighborPB11
    
    // Constructors
    .def(py::init<NodeList<Dimension>&, const NeighborSearchType, int, double, Vector, double, int>(),
         "nodeList"_a,
         "searchType"_a = NeighborSearchType::GatherScatter,
         "numGridLevels"_a = 31,
         "topGridCellSize"_a = 100.0,
         "origin"_a = Vector::zero,
         "kernelExtent"_a = 2.0,
         "gridCellInfluenceRadius"_a = 1)

    // Methods
    .def("gridLevel", (int (NGT::*)(int) const) &NGT::gridLevel)
    .def("gridCellIndex", (GCI (NGT::*)(int, int) const) &NGT::gridCellIndex, "nodeID"_a, "gridLevel"_a)
    .def("gridCellIndex", (GCI (NGT::*)(const Vector&, int) const) &NGT::gridCellIndex, "position"_a, "gridLevel"_a)
    .def("translateGridCellRange", &NGT::translateGridCellRange)
    .def("cellOccupied", &NGT::cellOccupied)
    .def("occupiedGridCells", (const std::vector<std::vector<GCI>>& (NGT::*)() const) &NGT::occupiedGridCells)
    .def("occupiedGridCells", (const std::vector<GCI>& (NGT::*)(const int) const) &NGT::occupiedGridCells, "gridLevelID")
    .def("headOfGridCell", &NGT::headOfGridCell)
    .def("nextNodeInCell", &NGT::nextNodeInCell, "nodeID"_a)
    .def("internalNodesInCell", &NGT::internalNodesInCell)
    .def("nodesInCell", &NGT::nodesInCell)
    .def("appendNodesInCell", &NGT::appendNodesInCell)
    .def("occupiedGridCellsInRange", &NGT::occupiedGridCellsInRange)
    .def("gridNormal", &NGT::gridNormal)
    .def("mapGridCell", &NGT::mapGridCell)
    .def("setNestedMasterList", (void (NGT::*)(const GCI&, const int)) &NGT::setNestedMasterList, "gridCell"_a, "gridLevel"_a)
    .def("findNestedNeighbors", &NGT::findNestedNeighbors)

    // Virtual methods
    .def("setMasterList", (void (NGT::*)(int)) &NT::setMasterList, "nodeID"_a)
    .def("setRefineNeighborList", (void (NGT::*)(int)) &NT::setRefineNeighborList, "nodeID"_a)

    // Attributes
    .def_property("numGridLevels",
                  (int (NGT::*)() const) &NGT::numGridLevels,
                  (void (NGT::*)(int)) &NGT::numGridLevels)

    ;

  //............................................................................
  // TreeNeighbor
  typedef TreeNeighbor<Dimension> TN;
  py::class_<TN, NT, PyNeighbor<Dimension, TN>> treeneighborPB11(m, ("TreeNeighbor" + suffix).c_str());
  virtualNeighborBindings<Dimension, TN>(m, suffix, treeneighborPB11);
  treeneighborPB11
    
    // Constructors
    .def(py::init<NodeList<Dimension>&, const NeighborSearchType, double, Vector, Vector>(),
         "nodeList"_a,
         "searchType"_a,
         "kernelExtent"_a,
         "xmin"_a,
         "xmax"_a)

    // Methods
    .def("gridLevel", (unsigned (TN::*)(const double&) const) &TN::gridLevel, "h"_a)
    .def("gridLevel", (unsigned (TN::*)(const SymTensor&) const) &TN::gridLevel, "H"_a)
    .def("dumpTree", &TN::dumpTree)
    .def("dumpTreeStatistics", &TN::dumpTreeStatistics)

    // Attributes
    .def_property_readonly("xmin", &TN::xmin)
    .def_property_readonly("xmax", &TN::xmax)
    .def_property_readonly("boxLength", &TN::boxLength)
    ;

  //............................................................................
  // ConnectivityMap
  typedef ConnectivityMap<Dimension> CM;
  py::class_<CM>(m, ("ConnectivityMap" + suffix).c_str())

    // Constructors
    .def(py::init<>())

    // Methods
    .def("patchConnectivity", &CM::patchConnectivity, "flags"_a, "old2new"_a)
    .def("connectivityForNode", (const std::vector<std::vector<int>>& (CM::*)(const NodeList<Dimension>*, const int) const) &CM::connectivityForNode, "nodeList"_a, "nodeID"_a)
    .def("connectivityForNode", (const std::vector<std::vector<int>>& (CM::*)(const int, const int) const) &CM::connectivityForNode, "nodeListID"_a, "nodeID"_a)
    .def("connectivityIntersectionForNodes", &CM::connectivityIntersectionForNodes)
    .def("connectivityUnionForNodes", &CM::connectivityUnionForNodes)
    .def("numNeighborsForNode", (int (CM::*)(const NodeList<Dimension>*, const int) const) &CM::numNeighborsForNode, "nodeList"_a, "nodeID"_a)
    .def("numNeighborsForNode", (int (CM::*)(const int, const int) const) &CM::numNeighborsForNode, "nodeListID"_a, "nodeID"_a)
    .def("calculatePairInteraction", &CM::calculatePairInteraction)
    .def("nodeList", &CM::nodeList)
    .def("nodeListIndex", &CM::nodeListIndex)
    .def("valid", &CM::valid)

    // Attributes
    .def_property_readonly("buildGhostConnectivity", &CM::buildGhostConnectivity)
    ;

  //............................................................................
  // The STL containers of Neighbor objects.
  py::bind_vector<std::vector<GCI>>(m, "vector_of_GridCellIndex" + suffix);
  py::bind_vector<std::vector<std::vector<GCI>>>(m, "vector_of_vector_of_GridCellIndex" + suffix);
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
#ifdef SPHERAL2D
  dimensionBindings<Spheral::Dim<2>>(m, "2d");
#endif
#ifdef SPHERAL3D
  dimensionBindings<Spheral::Dim<3>>(m, "3d");
#endif

  return m.ptr();
}
