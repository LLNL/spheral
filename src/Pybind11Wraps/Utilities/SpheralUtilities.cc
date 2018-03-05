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

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralGravity, m) {
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

}
