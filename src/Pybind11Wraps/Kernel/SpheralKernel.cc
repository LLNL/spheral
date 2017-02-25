#include <vector>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

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

  py::class_<SincKernel<Dimension>> sincPB11(m, ("SincKernel" + suffix).c_str());
  py::class_<NSincPolynomialKernel<Dimension>> nsincPB11(m, ("NSincPolynomialKernel" + suffix).c_str());
  py::class_<NBSplineKernel<Dimension>> nbsplinePB11(m, ("NBSplineKernel" + suffix).c_str());
  py::class_<QuarticSplineKernel<Dimension>> quarticsplinePB11(m, ("QuarticSplineKernel" + suffix).c_str());
  py::class_<QuinticSplineKernel<Dimension>> quinticsplinePB11(m, ("QuinticSplineKernel" + suffix).c_str());
  py::class_<WendlandC2Kernel<Dimension>> wc2PB11(m, ("WendlandC2Kernel" + suffix).c_str());
  py::class_<WendlandC4Kernel<Dimension>> wc4PB11(m, ("WendlandC4Kernel" + suffix).c_str());
  py::class_<WendlandC6Kernel<Dimension>> wc6PB11(m, ("WendlandC6Kernel" + suffix).c_str());
  py::class_<ExpInvKernel<Dimension>> expinvPB11(m, ("ExpInvKernel" + suffix).c_str());
  py::class_<TableKernel<Dimension>> tablePB11(m, ("TableKernel" + suffix).c_str());

  //

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralKernel) {
  py::module m("SpheralKernel", "Spheral Kernel module.");

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
