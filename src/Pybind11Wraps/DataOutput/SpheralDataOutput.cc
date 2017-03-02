// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "DataOutput/RestartRegistrar.hh"
#include "RestartableObject.hh"
#include "FileIO/FileIO.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::DataOutput;

namespace Spheral {
namespace DataOutput {

//------------------------------------------------------------------------------
// PyRestartableObject
//------------------------------------------------------------------------------
class PyRestartableObject: public RestartableObject {
public:
  using RestartableObject::RestartableObject;  // inherit constructors

  virtual std::string label() const override {
    PYBIND11_OVERLOAD_PURE(std::string,          // Return type
                           RestartableObject,    // Parent class
                           label,                // name of method
      );                        // arguments
  }

  virtual void dumpState(FileIOSpace::FileIO& file, const std::string pathName) const override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           RestartableObject,    // Parent class
                           dumpState,            // name of method
                           file,                 // arguments
                           pathName
      );
  }

  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string pathName) override {
    PYBIND11_OVERLOAD_PURE(void,                 // Return type
                           RestartableObject,    // Parent class
                           restoreState,         // name of method
                           file,                 // arguments
                           pathName
      );
  }
};

}
}

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralDataOutput) {
  py::module m("SpheralDataOutput", "Spheral DataOutput module.");

  //............................................................................
  // RestartRegistrar
  py::class_<RestartRegistrar, std::unique_ptr<RestartRegistrar, py::nodelete>>(m, "RestartRegistrar", py::metaclass())
    .def_property_readonly_static("instance", &RestartRegistrar::instance)
    .def("removeExpiredPointers", &RestartRegistrar::removeExpiredPointers)
    .def("uniqueLabels", &RestartRegistrar::uniqueLabels)
    .def("printLabels", &RestartRegistrar::printLabels)
    .def("dumpState", &RestartRegistrar::dumpState)
    .def("restoreState", &RestartRegistrar::restoreState)
    ;

  //............................................................................
  // RestartableObject
  py::class_<RestartableObject, PyRestartableObject>(m, "RestartableObject")
    .def(py::init<unsigned>(), "priority"_a)
    .def("label", &RestartableObject::label)
    .def("dumpState", &RestartableObject::dumpState)
    .def("restoreState", &RestartableObject::restoreState)
    ;

  return m.ptr();
}
