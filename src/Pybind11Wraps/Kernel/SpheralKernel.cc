// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "Kernel/BSplineKernel.hh"
#include "Kernel/W4SplineKernel.hh"
#include "Kernel/GaussianKernel.hh"
#include "Kernel/SuperGaussianKernel.hh"
#include "Kernel/PiGaussianKernel.hh"
#include "Kernel/HatKernel.hh"
#include "Kernel/SincKernel.hh"
#include "Kernel/NSincPolynomialKernel.hh"
#include "Kernel/NBSplineKernel.hh"
#include "Kernel/QuarticSplineKernel.hh"
#include "Kernel/QuinticSplineKernel.hh"
#include "Kernel/TableKernel.hh"
#include "Kernel/WendlandC2Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"
#include "Kernel/ExpInvKernel.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::KernelSpace;

namespace {  // anonymous

//------------------------------------------------------------------------------
// Common methods to Kernel objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void genericKernelBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  obj

    // Methods
    .def("__call__", (double (Obj::*)(double, const double&) const) &Obj::operator(), py::is_operator())
    .def("__call__", (double (Obj::*)(const Vector&, const double&) const) &Obj::operator(), py::is_operator())
    .def("__call__", (double (Obj::*)(double, const SymTensor&) const) &Obj::operator(), py::is_operator())
    .def("__call__", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::operator(), py::is_operator())

    .def("grad", (double (Obj::*)(double, const double&) const) &Obj::grad)
    .def("grad", (double (Obj::*)(const Vector&, const double&) const) &Obj::grad)
    .def("grad", (double (Obj::*)(double, const SymTensor&) const) &Obj::grad)
    .def("grad", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::grad)

    .def("grad2", (double (Obj::*)(double, const double&) const) &Obj::grad2)
    .def("grad2", (double (Obj::*)(const Vector&, const double&) const) &Obj::grad2)
    .def("grad2", (double (Obj::*)(double, const SymTensor&) const) &Obj::grad2)
    .def("grad2", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::grad2)

    .def("gradh", (double (Obj::*)(double, const double&) const) &Obj::gradh)
    .def("gradh", (double (Obj::*)(const Vector&, const double&) const) &Obj::gradh)
    .def("gradh", (double (Obj::*)(double, const SymTensor&) const) &Obj::gradh)
    .def("gradh", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::gradh)

    .def("kernelValue", &Obj::kernelValue)
    .def("gradValue", &Obj::gradValue)
    .def("grad2Value", &Obj::grad2Value)
    .def("gradhValue", &Obj::gradhValue)
    .def("valid", &Obj::valid)

    // These methods are in protected scope
    // .def("setVolumeNormalization", &Obj::setVolumeNormalization)
    // .def("setKernelExtent", &Obj::setKernelExtent)
    // .def("setInflectionPoint", &Obj::setInflectionPoint)

    // Attributes
    .def_property_readonly("volumeNormalization", &Obj::volumeNormalization)
    .def_property_readonly("kernelExtent", &Obj::kernelExtent)
    .def_property_readonly("inflectionPoint", &Obj::inflectionPoint)
    ;
}

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  using namespace Spheral::KernelSpace;

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  //............................................................................
  // BSpline
  py::class_<BSplineKernel<Dimension>> bsplinePB11(m, ("BSplineKernel" + suffix).c_str());
  genericKernelBindings<Dimension, BSplineKernel<Dimension>>(m, suffix, bsplinePB11);
  bsplinePB11.def(py::init<>());

  //............................................................................
  // W4Spline
  py::class_<W4SplineKernel<Dimension>> w4splinePB11(m, ("W4SplineKernel" + suffix).c_str());
  genericKernelBindings<Dimension, W4SplineKernel<Dimension>>(m, suffix, w4splinePB11);
  w4splinePB11.def(py::init<>());

  //............................................................................
  // Gaussian
  py::class_<GaussianKernel<Dimension>> gaussianPB11(m, ("GaussianKernel" + suffix).c_str());
  genericKernelBindings<Dimension, GaussianKernel<Dimension>>(m, suffix, gaussianPB11);
  gaussianPB11.def(py::init<double>(), "extent"_a);

  //............................................................................
  // SuperGaussian
  py::class_<SuperGaussianKernel<Dimension>> supergaussianPB11(m, ("SuperGaussianKernel" + suffix).c_str());
  genericKernelBindings<Dimension, SuperGaussianKernel<Dimension>>(m, suffix, supergaussianPB11);
  supergaussianPB11.def(py::init<>());

  //............................................................................
  // PiGaussian
  py::class_<PiGaussianKernel<Dimension>> pigaussianPB11(m, ("PiGaussianKernel" + suffix).c_str());
  genericKernelBindings<Dimension, PiGaussianKernel<Dimension>>(m, suffix, pigaussianPB11);
  pigaussianPB11
    .def(py::init<>())
    .def(py::init<double>(), "K"_a)
    .def_property("K", &PiGaussianKernel<Dimension>::getK, &PiGaussianKernel<Dimension>::setK)
    ;

  //............................................................................
  // Hat
  py::class_<HatKernel<Dimension>> hatPB11(m, ("HatKernel" + suffix).c_str());
  genericKernelBindings<Dimension, HatKernel<Dimension>>(m, suffix, hatPB11);
  hatPB11
    .def(py::init<double, double>(), "eta0"_a, "W0"_a)
    .def_property_readonly("eta0", &HatKernel<Dimension>::eta0)
    .def_property_readonly("W0", &HatKernel<Dimension>::eta0)
    ;

  //............................................................................
  // Sinc
  py::class_<SincKernel<Dimension>> sincPB11(m, ("SincKernel" + suffix).c_str());
  genericKernelBindings<Dimension, SincKernel<Dimension>>(m, suffix, sincPB11);
  sincPB11.def(py::init<double>(), "extent"_a);

  //............................................................................
  // NSinc
  py::class_<NSincPolynomialKernel<Dimension>> nsincPB11(m, ("NSincPolynomialKernel" + suffix).c_str());
  genericKernelBindings<Dimension, NSincPolynomialKernel<Dimension>>(m, suffix, nsincPB11);
  nsincPB11.def(py::init<int>(), "order"_a);

  //............................................................................
  // NBSpline
  py::class_<NBSplineKernel<Dimension>> nbsplinePB11(m, ("NBSplineKernel" + suffix).c_str());
  genericKernelBindings<Dimension, NBSplineKernel<Dimension>>(m, suffix, nbsplinePB11);
  nbsplinePB11
    .def(py::init<int>(), "order"_a)
    .def("factorial", &NBSplineKernel<Dimension>::factorial)
    .def("binomialCoefficient", &NBSplineKernel<Dimension>::binomialCoefficient)
    .def("oneSidedPowerFunction", &NBSplineKernel<Dimension>::oneSidedPowerFunction)
    .def_property("order", &NBSplineKernel<Dimension>::order, &NBSplineKernel<Dimension>::setOrder)
    ;

  //............................................................................
  // QuarticSpline
  py::class_<QuarticSplineKernel<Dimension>> quarticsplinePB11(m, ("QuarticSplineKernel" + suffix).c_str());
  genericKernelBindings<Dimension, QuarticSplineKernel<Dimension>>(m, suffix, quarticsplinePB11);
  quarticsplinePB11.def(py::init<>());

  //............................................................................
  // QuinticSpline
  py::class_<QuinticSplineKernel<Dimension>> quinticsplinePB11(m, ("QuinticSplineKernel" + suffix).c_str());
  genericKernelBindings<Dimension, QuinticSplineKernel<Dimension>>(m, suffix, quinticsplinePB11);
  quinticsplinePB11.def(py::init<>());

  //............................................................................
  // WendlandC2
  py::class_<WendlandC2Kernel<Dimension>> wc2PB11(m, ("WendlandC2Kernel" + suffix).c_str());
  genericKernelBindings<Dimension, WendlandC2Kernel<Dimension>>(m, suffix, wc2PB11);
  wc2PB11.def(py::init<>());

  //............................................................................
  // WendlandC4
  py::class_<WendlandC4Kernel<Dimension>> wc4PB11(m, ("WendlandC4Kernel" + suffix).c_str());
  genericKernelBindings<Dimension, WendlandC4Kernel<Dimension>>(m, suffix, wc4PB11);
  wc4PB11.def(py::init<>());

  //............................................................................
  // WendlandC6
  py::class_<WendlandC6Kernel<Dimension>> wc6PB11(m, ("WendlandC6Kernel" + suffix).c_str());
  genericKernelBindings<Dimension, WendlandC6Kernel<Dimension>>(m, suffix, wc6PB11);
  wc6PB11.def(py::init<>());

  //............................................................................
  // ExpInv
  py::class_<ExpInvKernel<Dimension>> expinvPB11(m, ("ExpInvKernel" + suffix).c_str());
  genericKernelBindings<Dimension, ExpInvKernel<Dimension>>(m, suffix, expinvPB11);
  expinvPB11.def(py::init<>());

  //............................................................................
  // TableKernel
  py::class_<TableKernel<Dimension>> tablePB11(m, ("TableKernel" + suffix).c_str());
  genericKernelBindings<Dimension, TableKernel<Dimension>>(m, suffix, tablePB11);
  tablePB11

    // Constructors
    .def(py::init<BSplineKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<W4SplineKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<GaussianKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<SuperGaussianKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<PiGaussianKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<HatKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<SincKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<NSincPolynomialKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<QuarticSplineKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<QuinticSplineKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<NBSplineKernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<WendlandC2Kernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<WendlandC4Kernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)
    .def(py::init<WendlandC6Kernel<Dimension>, int, double>(), "kernel"_a, "numPoints"_a=1000, "hmult"_a=1.0)

    // Methods
    .def("kernelAndGradValue", &TableKernel<Dimension>::kernelAndGradValue)
    .def("kernelAndGradValues", &TableKernel<Dimension>::kernelAndGradValues)
    .def("equivalentNodesPerSmoothingScale", &TableKernel<Dimension>::equivalentNodesPerSmoothingScale)
    .def("equivalentWsum", &TableKernel<Dimension>::equivalentWsum)
    .def("f1", &TableKernel<Dimension>::f1)
    .def("f2", &TableKernel<Dimension>::f2)
    .def("gradf1", &TableKernel<Dimension>::gradf1)
    .def("gradf2", &TableKernel<Dimension>::gradf2)
    .def("f1Andf2", &TableKernel<Dimension>::f1Andf2)
    .def("lowerBound", &TableKernel<Dimension>::lowerBound)
    .def("valid", &TableKernel<Dimension>::valid)

    // Attributes
    .def_property_readonly("nperhValues", &TableKernel<Dimension>::nperhValues)
    .def_property_readonly("WsumValues", &TableKernel<Dimension>::WsumValues)
    .def_property_readonly("numPoints", &TableKernel<Dimension>::numPoints)
    .def_property_readonly("stepSize", &TableKernel<Dimension>::stepSize)
    .def_property_readonly("stepSizeInv", &TableKernel<Dimension>::stepSizeInv)
    ;
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralKernel) {
  py::module m("SpheralKernel", "Spheral Kernel module.");

  // //............................................................................
  // // imports
  // py::object vector_of_double = (py::object) py::module::import("SpheralCXXTypes").attr("vector_of_double");

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
