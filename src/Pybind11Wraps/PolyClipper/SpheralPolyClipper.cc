// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "Geometry/polyclipper.hh"

namespace {

//------------------------------------------------------------------------------
// Plane methods
//------------------------------------------------------------------------------
template<typename Dimension,
         typename Plane,
         typename PlanePB11>
void planeBindings(PlanePB11& x) {

  typedef typename Plane::Vector Vector;
  namespace py = pybind11;
  using namespace pybind11::literals;

  // Constructors
  x.def(py::init<>());
  x.def(py::init<const Plane&>());
  x.def(py::init<const double, const Vector&>(), "dist"_a, "normal"_a);
  x.def(py::init<const Vector&, const Vector&>(), "point"_a, "normal"_a);

  // Attributes
  x.def_readwrite("dist", &Plane::dist);
  x.def_readwrite("normal", &Plane::normal);
}

//------------------------------------------------------------------------------
// Vertex methods
//------------------------------------------------------------------------------
template<typename Dimension,
         typename Vertex,
         typename VertexPB11>
void vertexBindings(VertexPB11& x) {

  typedef typename Vertex::Vector Vector;
  namespace py = pybind11;
  using namespace pybind11::literals;

  // Constructors
  x.def(py::init<>());
  x.def(py::init<const Vertex&>());
  x.def(py::init<const Vector&>(), "position"_a);
  x.def(py::init<const Vector&, int>(), "position"_a, "c"_a);

  // Attributes
  x.def_readwrite("position", &Vertex::position);
  x.def_readwrite("neighbors", &Vertex::neighbors);
  x.def_readwrite("comp", &Vertex::comp);
  x.def_readwrite("ID", &Vertex::ID);

  // Comparisons
  x.def(py::self == py::self);
}

}

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralOpenMP, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;

  m.doc() = "Spheral OpenMP wrappers.";

  //............................................................................
  // 2D plane
  py::class_<PolyClipper::Plane2d> plane2dPB11(m, "PolyClipperPlane2d");
  planeBindings<Spheral::Dim<2>, PolyClipper::Plane2d>(plane2dPB11);
  py::bind_vector<std::vector<PolyClipper::Plane2d>>(m, "vector_of_PolyClipperPlane2d");

  //............................................................................
  // 3D plane
  py::class_<PolyClipper::Plane3d> plane3dPB11(m, "PolyClipperPlane3d");
  planeBindings<Spheral::Dim<3>, PolyClipper::Plane3d>(plane3dPB11);
  py::bind_vector<std::vector<PolyClipper::Plane3d>>(m, "vector_of_PolyClipperPlane3d");

  //............................................................................
  // Vertex2d
  py::class_<PolyClipper::Vertex2d> vertex2dPB11(m, "Vertex2d");
  vertexBindings<Spheral::Dim<2>, PolyClipper::Vertex2d>(vertex2dPB11);
  py::bind_vector<std::vector<PolyClipper::Vertex2d>>(m, "Polygon");

  //............................................................................
  // Vertex3d
  py::class_<PolyClipper::Vertex3d> vertex3dPB11(m, "Vertex3d");
  vertexBindings<Spheral::Dim<3>, PolyClipper::Vertex3d>(vertex3dPB11);
  py::bind_vector<std::vector<PolyClipper::Vertex3d>>(m, "Polyhedron");

  //............................................................................
  // Polygon methods
  m.def("initializePolygon", &PolyClipper::initializePolygon,
        "poly"_a, "positions"_a, "neighbors"_a,
        "Initialize a PolyClipper::Polygon from vertex positions and vertex neighbors.");
  m.def("polygon2string", &PolyClipper::polygon2string, "poly"_a,
        "Return a formatted string representation for a PolyClipper::Polygon.");
}
