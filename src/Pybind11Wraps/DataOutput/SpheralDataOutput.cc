#include <vector>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "Geometry/Dimension.hh"
#include "DataOutput/RestartRegistrar.hh"
#include "RestartableObject.hh"
#include "FileIO/FileIO.hh"

namespace py = pybind11;
using namespace pybind11::literals;

namespace Spheral {
namespace DataOutput {

//------------------------------------------------------------------------------
// PyRestartableObject
//------------------------------------------------------------------------------
class PyRestartableObject: public RestartableObject {
public:
  using RestartableObject::RestartableObject;  // inherit constructors

  virtual std::string label() const override {
    PYBIND11_OVERLOAD(std::string,          // Return type
                      RestartableObject,    // Parent class
                      label,                // name of method
      );                        // arguments
  }

  virtual void dumpState(FileIOSpace::FileIO& file, const std::string pathName) const override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      RestartableObject,    // Parent class
                      dumpState,            // name of method
                      file,                 // arguments
                      pathName
      );
  }

  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string pathName) override {
    PYBIND11_OVERLOAD(void,                 // Return type
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
  py::class_<Spheral::DataOutput::RestartRegistrar, std::unique_ptr<Spheral::DataOutput::RestartRegistrar, py::nodelete>>(m, "RestartRegistrar", py::metaclass())
    .def_property_readonly_static("instance", &Spheral::DataOutput::RestartRegistrar::instance)
    .def("removeExpiredPointers", &Spheral::DataOutput::RestartRegistrar::removeExpiredPointers)
    .def("uniqueLabels", &Spheral::DataOutput::RestartRegistrar::uniqueLabels)
    .def("printLabels", &Spheral::DataOutput::RestartRegistrar::printLabels)
    .def("dumpState", &Spheral::DataOutput::RestartRegistrar::dumpState)
    .def("restoreState", &Spheral::DataOutput::RestartRegistrar::restoreState)
    ;

  //............................................................................
  // RestartableObject
  py::class_<Spheral::DataOutput::RestartableObject, Spheral::DataOutput::PyRestartableObject>(m, "RestartableObject")
    .def(py::init<unsigned>(), "priority"_a)
    .def("label", &Spheral::DataOutput::RestartableObject::label)
    .def("dumpState", &Spheral::DataOutput::RestartableObject::dumpState)
    .def("restoreState", &Spheral::DataOutput::RestartableObject::restoreState)
    ;

  return m.ptr();
}
