// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>
#include <utility>

#include "Geometry/Dimension.hh"
#include "Utilities/packElement.hh"
#include "boost/math/special_functions/legendre.hpp"
#include "Utilities/Functors.hh"
#include "Utilities/erff.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/iterateIdealH.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/peanoHilbertOrderIndices.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/segmentIntersectEdges.hh"
#include "Utilities/integrateThroughMeshAlongSegment.hh"
#include "Utilities/numberDensity.hh"
//#include "Utilities/writeRectilinearMesh.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/KeyTraits.hh"
#include "Utilities/pointInPolygon.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/pointInPolyhedron.hh"
#include "Utilities/pointOnPolyhedron.hh"
#include "Utilities/refinePolyhedron.hh"
#include "Utilities/overlayRemapFields.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "Utilities/Timer.hh"

#include "Pybind11Wraps/Utilities/PyAbstractSpheralFunctor.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

#ifndef NOR3D
#include "Utilities/r3d_utils.hh"
#else
//------------------------------------------------------------------------------
// Stub these methods out when we're not allowing R3D.
//------------------------------------------------------------------------------
namespace Spheral {
  inline Dim<2>::FacetedVolume clipFacetedVolume(const Dim<2>::FacetedVolume& poly,
                                                 const std::vector<GeomPlane<Dim<2> > >& planes) {
    VERIFY2(false, "ERROR: clipFacetedVolume unavailable without R3D.");
  }

  inline Dim<3>::FacetedVolume clipFacetedVolume(const Dim<3>::FacetedVolume& poly,
                                                 const std::vector<GeomPlane<Dim<3> > >& planes) {
    VERIFY2(false, "ERROR: clipFacetedVolume unavailable without R3D.");
  }
}
#endif

namespace py = pybind11;
using namespace pybind11::literals;
using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::FieldList;
using Spheral::PythonBoundFunctors::SpheralFunctor;

using namespace Spheral;

namespace {  // anonymous

//------------------------------------------------------------------------------
// toString
//------------------------------------------------------------------------------
template<typename DataType>
std::string toString(const DataType& x) {
  vector<char> buffer;
  packElement(x, buffer);
  return string(buffer.begin(), buffer.end());
}

//------------------------------------------------------------------------------
// fromString
//------------------------------------------------------------------------------
template<typename DataType>
DataType fromString(const string& x) {
  DataType result;
  vector<char> buffer(x.begin(), x.end());
  vector<char>::const_iterator itr = buffer.begin();
  unpackElement(result, itr, buffer.end());
  return result;
}

//------------------------------------------------------------------------------
// Bind the objects templated on Dimension.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  //............................................................................
  // VectorScalarFunctor
  {
    typedef SpheralFunctor<Vector, double> VSF;
    py::class_<VSF, PyAbstractSpheralFunctor<Vector, double, VSF>> x(m, ("VectorScalarFunctor" + suffix).c_str());

    // Constructors
    x.def(py::init<>());

    // Methods
    x.def("__call__", &VSF::__call__, "x"_a);
  }

  //............................................................................
  // VectorVectorFunctor
  {
    typedef SpheralFunctor<Vector, Vector> VVF;
    py::class_<VVF, PyAbstractSpheralFunctor<Vector, Vector, VVF>> x(m, ("VectorVectorFunctor" + suffix).c_str());

    // Constructors
    x.def(py::init<>());

    // Methods
    x.def("__call__", &VVF::__call__, "x"_a);
  }

  //............................................................................
  // VectorPairScalarFunctor
  {
    typedef SpheralFunctor<Vector, std::pair<double, double>> VPSF;
    py::class_<VPSF, PyAbstractSpheralFunctor<Vector, std::pair<double, double>, VPSF>> x(m, ("VectorPairScalarFunctor" + suffix).c_str());

    // Constructors
    x.def(py::init<>());

    // Methods
    x.def("__call__", &VPSF::__call__, "x"_a);
  }

  //............................................................................
  m.def("boundingBox", &boundingBox<Vector>, "positions"_a, "xmin"_a, "xmax"_a);
  m.def("globalBoundingBox", (void (*)(const FieldList<Dimension, Vector>&, Vector&, Vector&, bool)) &globalBoundingBox<Dimension>, "positions"_a, "xmin"_a, "xmax"_a, "ghost"_a=false);
  m.def("globalBoundingVolumes", &globalBoundingVolumes<Dimension>, "dataBase"_a, "nodeVolume"_a, "sampleVolume"_a);

  m.def("collinear", &collinear<Vector>, "a"_a, "b"_a, "c"_a, "tol"_a=1.0e-10, "Test if three points are collinear.");
  m.def("between", &between<Vector>, "a"_a, "b"_a, "c"_a, "tol"_a=1.0e-10, "Test if point c is between (a,b).");
  m.def("segmentSegmentDistance", (double (*)(const Vector&, const Vector&, const Vector&, const Vector&)) &segmentSegmentDistance,
        "a0"_a, "a1"_a, "b0"_a, "b1"_a, "Find the distance between line segements (a0,a1) -> (b0,b1)");
  m.def("segmentSegmentIntersection", (char (*)(const Vector&, const Vector&, const Vector&, const Vector&, Vector&, Vector&, const double)) &segmentSegmentIntersection,
        "a0"_a, "a1"_a, "b0"_a, "b1"_a, "result1"_a, "result2"_a, "tol"_a=1e-8, "Compute the intersection of two line segments (a0,a1) (b0,b1), returned as last two args.");
  m.def("segmentSegmentIntersection", &segmentSegmentIntersection<Vector>,
        "a0"_a, "a1"_a, "b0"_a, "b1"_a, "tol"_a=1e-8, "Test if two line segments (a0,a1) (b0,b1) intersect.");

  //............................................................................
  m.def("numGlobalNodes", (int (*)(const NodeList<Dimension>&)) &NodeSpace::numGlobalNodes<Dimension>, "nodes"_a, 
        "Global number of nodes in the NodeList.");
  m.def("numGlobalNodes", (int (*)(const DataBase<Dimension>&)) &NodeSpace::numGlobalNodes<Dimension>, "dataBase"_a, 
        "Global number of nodes in the DataBase.");
  m.def("globalNodeIDs", (Field<Dimension, int> (*)(const NodeList<Dimension>&)) &NodeSpace::globalNodeIDs<Dimension>, "nodes"_a, 
        "Determine unique global node IDs for the nodes in a NodeList.");
  m.def("globalNodeIDs", (FieldList<Dimension, int> (*)(const DataBase<Dimension>&)) &NodeSpace::globalNodeIDs<Dimension>, "dataBase"_a, 
        "Determine unique global node IDs for the nodes in a DataBase.");

  m.def("iterateIdealH", &Spheral::iterateIdealH<Dimension>,
        "dataBase"_a, "boundaries"_a, "W"_a, "smoothingScaleMethod"_a,
        "maxIterations"_a = 100,
        "tolerance"_a = 1.0e-10,
        "nPerhForIteration"_a = 0.0,
        "sphericalStart"_a = false,
        "fixDeterminant"_a = false,
        "Iterate the Hfield for NodeLists in the DataBase using the ideal H algorithm.");

  m.def("mortonOrderIndices",
        (FieldList<Dimension, KeyTraits::Key> (*)(const FieldList<Dimension, Vector>&)) &Spheral::mortonOrderIndices<Dimension>,
        "positions"_a,
        "Compute indices for nodes obeying Morton ordering given the positions.");
  m.def("mortonOrderIndices",
        (FieldList<Dimension, KeyTraits::Key> (*)(const DataBase<Dimension>&)) &Spheral::mortonOrderIndices<Dimension>,
        "dataBase"_a,
        "Compute indices for nodes obeying Morton ordering in the given DataBase.");
  m.def("mortonOrderIndices",
        (FieldList<Dimension, KeyTraits::Key> (*)(const DataBase<Dimension>&, const FieldList<Dimension, int>&)) &Spheral::mortonOrderIndices<Dimension>,
        "dataBase"_a, "mask"_a,
        "Compute indices for nodes obeying Morton ordering in the given DataBase and mask.");
        
  m.def("peanoHilbertOrderIndices",
        (FieldList<Dimension, KeyTraits::Key> (*)(const FieldList<Dimension, Vector>&)) &Spheral::peanoHilbertOrderIndices<Dimension>,
        "positions"_a,
        "Compute indices for nodes obeying Peano-Hilbert ordering given the positions.");
  m.def("peanoHilbertOrderIndices",
        (FieldList<Dimension, KeyTraits::Key> (*)(const DataBase<Dimension>&)) &Spheral::peanoHilbertOrderIndices<Dimension>,
        "dataBase"_a,
        "Compute indices for nodes obeying Peano-Hilbert ordering in the given DataBase.");
        
  m.def("numberDensity", &Spheral::numberDensity<Dimension>, "dataBase"_a, "W"_a,
        "Compute the ASPH sum number density for each node in a DataBase.");

  m.def("integrateThroughMeshAlongSegment", &Spheral::integrateThroughMeshAlongSegment<Dimension, Scalar>,
        "values"_a, "xmin"_a, "xmax"_a, "ncells"_a, "s0"_a, "s1"_a,
        "Integrate through a lattice sampled field along a line segment.");

  //............................................................................
  // Shepards interpolation methods
  m.def("computeShepardsInterpolation", &Spheral::computeShepardsInterpolation<Dimension, Scalar>,
        "fieldList"_a, "connectivityMap"_a, "W"_a, "position"_a, "H"_a, "weight"_a, 
        "Interpolate a FieldList using a Shepards function approach.");
  m.def("computeShepardsInterpolation", &Spheral::computeShepardsInterpolation<Dimension, Vector>,
        "fieldList"_a, "connectivityMap"_a, "W"_a, "position"_a, "H"_a, "weight"_a, 
        "Interpolate a FieldList using a Shepards function approach.");
  m.def("computeShepardsInterpolation", &Spheral::computeShepardsInterpolation<Dimension, Tensor>,
        "fieldList"_a, "connectivityMap"_a, "W"_a, "position"_a, "H"_a, "weight"_a, 
        "Interpolate a FieldList using a Shepards function approach.");
  m.def("computeShepardsInterpolation", &Spheral::computeShepardsInterpolation<Dimension, SymTensor>,
        "fieldList"_a, "connectivityMap"_a, "W"_a, "position"_a, "H"_a, "weight"_a, 
        "Interpolate a FieldList using a Shepards function approach.");

  //............................................................................
  m.def("rotationMatrix", (Tensor (*)(const Vector&)) &Spheral::rotationMatrix, "runit"_a,
        "Rotational transformation to align with the given unit vector.");
  m.def("testBoxIntersection", (bool (*)(const Vector&, const Vector&, const Vector&, const Vector&, const double)) &Spheral::testBoxIntersection,
        "xmin1"_a, "xmax1"_a, "xmin2"_a, "xmax2"_a, "tol"_a = 1.0e-10,
        "Test if the two boxes intersect.");
  m.def("testPointInBox", (bool (*)(const Vector&, const Vector&, const Vector&, const double)) &Spheral::testPointInBox,
        "point"_a, "xmin"_a, "xmax"_a, "tol"_a = 1.0e-10,
        "Test if the point is in the box.");
  m.def("planarReflectingOperator", &Spheral::planarReflectingOperator<GeomPlane<Dimension>>, "plane"_a,
        "Generate the planar reflection transformation for th given plane.");
}

//------------------------------------------------------------------------------
// These methods are only valid in (2D,3D).
//------------------------------------------------------------------------------
template<typename Dimension>
void dimension23Bindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  m.def("pointPlaneDistance", &pointPlaneDistance<Vector>, "point"_a, "origin"_a, "unitNormal"_a,
        "Compute the distance from (point) to the plane defined by (origin, unitNormal).");
  m.def("closestPointOnSegment", &closestPointOnSegment<Vector>, "p"_a, "a0"_a, "a1"_a,
        "Find the closest point on a line segment (a0,a1) to point (p).");
  m.def("overlayRemapFields", &overlayRemapFields<Dimension>, 
        "boundaries"_a, "scalarDonorFields"_a, "vectorDonorFields"_a, "tensorDonorFields"_a, "symTensorDonorFields"_a, 
        "scalarAcceptorFields"_a, "vectorAcceptorFields"_a, "tensorAcceptorFields"_a, "symTensorAcceptorFields"_a,
        "Do a simple donor overlay using geometric intersection.");
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralUtilities, m) {
  using namespace Spheral;

  m.doc() = "Spheral Utilities module: a grab-bag of stuff.";

  //............................................................................
  // KeyTraits
  {
    py::class_<KeyTraits> x(m, "KeyTraits");

    // Attributes
    x.def_readonly_static("numbits", &KeyTraits::numbits);
    x.def_readonly_static("numbits1d", &KeyTraits::numbits1d);
    x.def_readonly_static("zero", &KeyTraits::zero);
    x.def_readonly_static("one", &KeyTraits::one);
    x.def_readonly_static("two", &KeyTraits::two);
    x.def_readonly_static("maxKey1d", &KeyTraits::maxKey1d);
    x.def_readonly_static("maxKey", &KeyTraits::maxKey);
  }

  //............................................................................
  // Timer
  {
    py::class_<Timer> x(m, "Timer");

    // Constructors
    x.def(py::init<>());
    x.def(py::init<std::string>(), "name"_a);
    x.def(py::init<std::string, Timer&>(), "name"_a, "base"_a);

    // Methods
    x.def("setup", &Timer::setup);
    x.def("start", &Timer::start);
    x.def("stop", &Timer::stop);
    x.def("clear", &Timer::clear);
    x.def("getTimeStampWC", &Timer::getTimeStampWC);
    x.def("wc_time", &Timer::wc_time);
    x.def("setup", &Timer::setup);
    x.def("Name", &Timer::Name);
    x.def("Count", &Timer::Count);
    x.def_static("TimerSummary", (void (*)()) &Timer::TimerSummary);
  }

  //............................................................................
  // ScalarScalarFunctor
  {
    typedef SpheralFunctor<double, double> SSF;
    py::class_<SSF, PyAbstractSpheralFunctor<double, double, SSF>> x(m, "ScalarScalarFunctor");

    // Constructors
    x.def(py::init<>());

    // Methods
    x.def("__call__", &SSF::__call__, "x"_a);
  }

  //............................................................................
  // ScalarPairScalarFunctor
  {
    typedef SpheralFunctor<double, std::pair<double, double>> SPSF;
    py::class_<SPSF, PyAbstractSpheralFunctor<double, std::pair<double, double>, SPSF>> x(m, "ScalarPairScalarFunctor");

    // Constructors
    x.def(py::init<>());

    // Methods
    x.def("__call__", &SPSF::__call__, "x"_a);
  }

  //............................................................................
  // Module functions
  m.def("erff", &Spheral::erff, "You know, the error function.");
  m.def("newtonRaphsonFindRoot", &newtonRaphson<SpheralFunctor<double, std::pair<double, double>>>,
        "function"_a, "x1"_a, "x2"_a, "xaccuracy"_a=1.0e-15, "yaccuracy"_a=1.0e-15, "maxIterations"_a=100,
        "Newton-Raphson root finder.");
  m.def("simpsonsIntegrationDouble", &simpsonsIntegration<SpheralFunctor<double, double>, double, double>,
        "function"_a, "x0"_a, "x1"_a, "numBins"_a,
        "Simpsons rule integration for double functions.");

  // packElement
  m.def("packElement", (void (*)(const int&, std::vector<char>&)) &packElement<int>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const double&, std::vector<char>&)) &packElement<double>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const uint32_t&, std::vector<char>&)) &packElement<uint32_t>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const uint64_t&, std::vector<char>&)) &packElement<uint64_t>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const std::vector<unsigned>&, std::vector<char>&)) &packElement<unsigned>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const std::vector<int>&, std::vector<char>&)) &packElement<int>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const std::vector<float>&, std::vector<char>&)) &packElement<float>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const std::vector<double>&, std::vector<char>&)) &packElement<double>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const std::vector<uint32_t>&, std::vector<char>&)) &packElement<uint32_t>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const std::vector<uint64_t>&, std::vector<char>&)) &packElement<uint64_t>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const std::vector<std::vector<unsigned>>&, std::vector<char>&)) &packElement<vector<unsigned>>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const Dim<1>::FacetedVolume&, std::vector<char>&)) &packElement<Dim<1>::FacetedVolume>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const Dim<2>::FacetedVolume&, std::vector<char>&)) &packElement<Dim<2>::FacetedVolume>, "x"_a, "buffer"_a);
  m.def("packElement", (void (*)(const Dim<3>::FacetedVolume&, std::vector<char>&)) &packElement<Dim<3>::FacetedVolume>, "x"_a, "buffer"_a);

  // toString
  m.def("toString", &toString<int>, "x"_a);
  m.def("toString", &toString<double>, "x"_a);
  m.def("toString", &toString<uint32_t>, "x"_a);
  m.def("toString", &toString<uint64_t>, "x"_a);
  m.def("toString", &toString<vector<unsigned>>, "x"_a);
  m.def("toString", &toString<vector<int>>, "x"_a);
  m.def("toString", &toString<vector<float>>, "x"_a);
  m.def("toString", &toString<vector<double>>, "x"_a);
  m.def("toString", &toString<vector<uint32_t>>, "x"_a);
  m.def("toString", &toString<vector<uint64_t>>, "x"_a);
  m.def("toString", &toString<vector<vector<unsigned>>>, "x"_a);
  m.def("toString", &toString<Dim<1>::FacetedVolume>, "x"_a);
  m.def("toString", &toString<Dim<2>::FacetedVolume>, "x"_a);
  m.def("toString", &toString<Dim<3>::FacetedVolume>, "x"_a);

  // fromString
  m.def("toInt", &fromString<int>, "x"_a);
  m.def("toDouble", &fromString<double>, "x"_a);
  m.def("toUint32_t", &fromString<uint32_t>, "x"_a);
  m.def("toUint64_t", &fromString<uint64_t>, "x"_a);
  m.def("toVectorUnsigned", &fromString<vector<unsigned>>, "x"_a);
  m.def("toVectorInt", &fromString<vector<int>>, "x"_a);
  m.def("toVectorFloat", &fromString<vector<float>>, "x"_a);
  m.def("toVectorDouble", &fromString<vector<double>>, "x"_a);
  m.def("toVectorUint32_t", &fromString<vector<uint32_t>>, "x"_a);
  m.def("toVectorUint64_t", &fromString<vector<uint64_t>>, "x"_a);
  m.def("toVectorVectorUnsigned", &fromString<vector<vector<unsigned>>>, "x"_a);
  m.def("toBox1d", &fromString<Dim<1>::FacetedVolume>, "x"_a);
  m.def("toPolygon", &fromString<Dim<2>::FacetedVolume>, "x"_a);
  m.def("toPolyhedron", &fromString<Dim<3>::FacetedVolume>, "x"_a);

  m.def("closestPointOnPlane", &closestPointOnPlane, "p"_a, "origin"_a, "unitNormal"_a,
        "Find the closest point in the plane (origin,normal) to point (p).");
  m.def("segmentPlaneIntersection", (char (*)(const Dim<3>::Vector&, const Dim<3>::Vector&, const Dim<3>::Vector&, const Dim<3>::Vector&, Dim<3>::Vector&, double)) &segmentPlaneIntersection,
        "s0"_a, "s1"_a, "point"_a, "normal"_a, "result"_a, "tol"_a=1.0e-8,
        "Compute the intesection of a line segment (s0,s1) with a plane (point,normal).");
  m.def("segmentPlaneIntersection", (char (*)(const Dim<3>::Vector&, const Dim<3>::Vector&, const Dim<3>::Vector&, const Dim<3>::Vector&, const Dim<3>::Vector&, Dim<3>::Vector&, const double)) &segmentPlaneIntersection,
        "s0"_a, "s1"_a, "p0"_a, "p1"_a, "p2"_a, "result"_a, "tol"_a=1.0e-8,
        "Compute the intesection of a line segment (s0,s1) with a plane (p0,p1,p2).");
  m.def("segmentPlanarSectionIntersection", &segmentPlanarSectionIntersection,
        "s0"_a, "s1"_a, "pverts"_a, "result"_a, "tol"_a=1.0e-8,
        "Compute the intesection of a line segment (s0,s1) with a polygonal planar section (pverts).");
  m.def("segmentPlaneIntersection", (bool (*)(const Dim<3>::Vector&, const Dim<3>::Vector&, const vector<Dim<3>::Vector>&, const Dim<3>::Vector&, const double)) &segmentPlaneIntersection,
        "a0"_a, "a1"_a, "vertices"_a, "normal"_a, "tol"_a=1.0e-10,
        "Test if a line segment intersects a planar section.");
  m.def("segmentPlaneIntersection", (bool (*)(const Dim<3>::Vector&, const Dim<3>::Vector&, const vector<Dim<3>::Vector>&, const vector<unsigned>&, const Dim<3>::Vector&, const double)) &segmentPlaneIntersection,
        "a0"_a, "a1"_a, "vertices"_a, "ipoints"_a, "normal"_a, "tol"_a=1.0e-10,
        "Test if a line segment intersects a planar section.");

  // Polygon utilities
  m.def("pointOnPolygon", (bool (*)(const Dim<2>::Vector&, const Dim<2>::FacetedVolume&, const double)) &pointOnPolygon, 
        "p"_a, "poly"_a, "tol"_a=1.0e-10,
        "Test if the given point is on the boundary of a polygon specified by it's vertices.");
  m.def("pointOnPolygon", (bool (*)(const Dim<2>::Vector&, const vector<Dim<2>::Vector>&, const vector<unsigned>&, const double)) &pointOnPolygon, 
        "p"_a, "vertices"_a, "ipoints"_a, "tol"_a=1.0e-10,
        "Test if the given point is on the boundary of a polygon specified by its vertices.");
  m.def("pointInPolygon", (bool (*)(const Dim<2>::Vector&, const vector<Dim<2>::Vector>&)) &pointInPolygon, 
        "p"_a, "vertices"_a, 
        "Test if the given point is contained within a polygon.");
  m.def("pointInPolygon", (bool (*)(const Dim<2>::Vector&, const Dim<2>::FacetedVolume&, const bool, const double)) &pointInPolygon, 
        "p"_a, "poly"_a, "countBoundary"_a=false, "tol"_a=1.0e-10,
        "Test if the given point is contained within a polygon.");
  m.def("segmentIntersectEdges", (bool (*)(const Dim<2>::Vector&, const Dim<2>::Vector&, const Dim<2>::FacetedVolume&, const double)) &segmentIntersectEdges, 
        "a0"_a, "a1"_a, "poly"_a, "tol"_a=1.0e-8, 
        "Test if the give line segment intersects any edges/vertices of the given polygon.");

  m.def("refinePolyhedron", &refinePolyhedron, "poly0"_a, "numLevels"_a,
        "Return a new Polyhedron based on refining an existing one a given number of levels.");

  // R2D/R3D
  m.def("clipFacetedVolume", (Dim<2>::FacetedVolume (*)(const Dim<2>::FacetedVolume&, const vector<GeomPlane<Dim<2>>>&)) &clipFacetedVolume,
        "poly"_a, "planes"_a, 
        "Clip a polygon with a set of planes.");
  m.def("clipFacetedVolume", (Dim<3>::FacetedVolume (*)(const Dim<3>::FacetedVolume&, const vector<GeomPlane<Dim<3>>>&)) &clipFacetedVolume,
        "poly"_a, "planes"_a, 
        "Clip a polyhedron with a set of planes.");
  m.def("clippedVolume", (double (*)(const Dim<2>::FacetedVolume&, const vector<GeomPlane<Dim<2>>>&)) &clippedVolume,
        "Return the area of a polygon clipped with a set of planes.");
  m.def("clippedVolume", (double (*)(const Dim<3>::FacetedVolume&, const vector<GeomPlane<Dim<3>>>&)) &clippedVolume,
        "Return the volume of a polyhedron clipped with a set of planes.");

  // Boost math functions
  m.def("legendre_p", (double (*)(int, int, double)) &boost::math::legendre_p, "l"_a, "m"_a, "x"_a,
        "Compute the associated Legendre polynomial.");

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif

#ifdef SPHERAL2D
  dimensionBindings<Spheral::Dim<2>>(m, "2d");
  dimension23Bindings<Spheral::Dim<2>>(m, "2d");
#endif

#ifdef SPHERAL3D
  dimensionBindings<Spheral::Dim<3>>(m, "3d");
  dimension23Bindings<Spheral::Dim<3>>(m, "3d");
#endif

}
