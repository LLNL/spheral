// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Boundary/Boundary.hh"
#include "Boundary/PlanarBoundary.hh"
#include "Boundary/ReflectingBoundary.hh"
#include "Boundary/RigidBoundary.hh"
#include "Boundary/PeriodicBoundary.hh"
#include "Boundary/ConstantVelocityBoundary.hh"
#include "Boundary/ConstantXVelocityBoundary.hh"
#include "Boundary/ConstantYVelocityBoundary.hh"
#include "Boundary/ConstantZVelocityBoundary.hh"
#include "Boundary/ConstantRVelocityBoundary.hh"
#include "Boundary/ConstantBoundary.hh"
#include "Boundary/SphericalBoundary.hh"
#include "Boundary/CylindricalBoundary.hh"
#include "Boundary/AxialSymmetryBoundary.hh"
#include "Boundary/AxisBoundaryRZ.hh"

#include "PyAbstractBoundary.hh"
#include "PyBoundary.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::BoundarySpace;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Boundary<Spheral::Dim<1>>*>);

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Boundary<Spheral::Dim<2>>*>);

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Boundary<Spheral::Dim<3>>*>);

namespace {  // anonymous

//------------------------------------------------------------------------------
// Common virtual methods of Boundary objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void virtualBoundaryBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef Spheral::GeomPlane<Dimension> Plane;

  obj

    // Methods
    .def("setAllGhostNodes", &Obj::setAllGhostNodes)
    .def("setAllViolationNodes", &Obj::setAllViolationNodes)
    .def("cullGhostNodes", &Obj::cullGhostNodes)
    .def("setGhostNodes", &Obj::setGhostNodes)
    .def("updateGhostNodes", &Obj::updateGhostNodes)

    .def("applyGhostBoundary", (void (Obj::*)(Field<Dimension, int>&) const) &Obj::applyGhostBoundary, "field"_a)
    .def("applyGhostBoundary", (void (Obj::*)(Field<Dimension, Scalar>&) const) &Obj::applyGhostBoundary, "field"_a)
    .def("applyGhostBoundary", (void (Obj::*)(Field<Dimension, Vector>&) const) &Obj::applyGhostBoundary, "field"_a)
    .def("applyGhostBoundary", (void (Obj::*)(Field<Dimension, Tensor>&) const) &Obj::applyGhostBoundary, "field"_a)
    .def("applyGhostBoundary", (void (Obj::*)(Field<Dimension, SymTensor>&) const) &Obj::applyGhostBoundary, "field"_a)
    .def("applyGhostBoundary", (void (Obj::*)(Field<Dimension, ThirdRankTensor>&) const) &Obj::applyGhostBoundary, "field"_a)
    .def("applyGhostBoundary", (void (Obj::*)(Field<Dimension, std::vector<Scalar>>&) const) &Obj::applyGhostBoundary, "field"_a)

    .def("setViolationNodes", &Obj::setViolationNodes)
    .def("updateViolationNodes", &Obj::updateViolationNodes)

    .def("enforceBoundary", (void (Obj::*)(Field<Dimension, int>&) const) &Obj::enforceBoundary, "field"_a)
    .def("enforceBoundary", (void (Obj::*)(Field<Dimension, Scalar>&) const) &Obj::enforceBoundary, "field"_a)
    .def("enforceBoundary", (void (Obj::*)(Field<Dimension, Vector>&) const) &Obj::enforceBoundary, "field"_a)
    .def("enforceBoundary", (void (Obj::*)(Field<Dimension, Tensor>&) const) &Obj::enforceBoundary, "field"_a)
    .def("enforceBoundary", (void (Obj::*)(Field<Dimension, SymTensor>&) const) &Obj::enforceBoundary, "field"_a)
    .def("enforceBoundary", (void (Obj::*)(Field<Dimension, ThirdRankTensor>&) const) &Obj::enforceBoundary, "field"_a)

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
  // Boundary
  typedef Boundary<Dimension> NT;
  py::class_<NT, PyAbstractBoundary<Dimension, NT>> boundaryPB11(m, ("Boundary" + suffix).c_str());
  virtualBoundaryBindings<Dimension, NT>(m, suffix, boundaryPB11);
  boundaryPB11
    
    // Constructors
    .def(py::init<NodeList<Dimension>&, const BoundarySearchType, const double>(), "nodeList"_a, "searchType"_a, "kernelExtent"_a)

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
    .def("setRefineBoundaryList", (void (NT::*)(int)) &NT::setRefineBoundaryList, "nodeID"_a)

    // Attributes
    .def_property("neighborSearchType",
                  (BoundarySearchType (NT::*)() const) &NT::neighborSearchType,
                  (void (NT::*)(BoundarySearchType)) &NT::neighborSearchType)
    .def_property("kernelExtent",
                  (double (NT::*)() const) &NT::kernelExtent,
                  (void (NT::*)(double)) &NT::kernelExtent)
    .def_property_readonly("numMaster", &NT::numMaster)
    .def_property_readonly("numCoarse", &NT::numCoarse)
    .def_property_readonly("numRefine", &NT::numRefine)
    .def_property_readonly("masterList", &NT::masterList)
    .def_property_readonly("coarseBoundaryList", &NT::coarseBoundaryList)
    .def_property_readonly("refineBoundaryList", &NT::refineBoundaryList)
    ;

  //............................................................................
  // The STL containers of Boundary objects.
  py::bind_vector<std::vector<GCI>>(m, "vector_of_GridCellIndex" + suffix);
  py::bind_vector<std::vector<std::vector<GCI>>>(m, "vector_of_vector_of_GridCellIndex" + suffix);
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralBoundary) {
  py::module m("SpheralBoundary", "Spheral Boundary module.");

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
