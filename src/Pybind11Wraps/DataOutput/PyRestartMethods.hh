//------------------------------------------------------------------------------
// Trampoline classes for restart interface.
//------------------------------------------------------------------------------
#ifndef __Spheral_PyRestartMethods__
#define __Spheral_PyRestartMethods__

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "FileIO/FileIO.hh"

namespace py = pybind11;
using namespace pybind11::literals;

namespace Spheral {

//------------------------------------------------------------------------------
// PyRestartMethods
//------------------------------------------------------------------------------
template<class Base>
class PyRestartMethods: public Base {
public:

  virtual std::string label() const override {
    PYBIND11_OVERLOAD(std::string,  // Return type
                      Base,         // Parent class
                      label         // name of method
      );
  }

  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const override {
    PYBIND11_OVERLOAD(void,          // Return type
                      Base,          // Parent class
                      dumpState,     // name of method
                      file, pathName // arguments
      );
  }

  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName) override {
    PYBIND11_OVERLOAD(void,          // Return type
                      Base,          // Parent class
                      restoreState,  // name of method
                      file, pathName // arguments
      );
  }
};

//------------------------------------------------------------------------------
// A function to add the restart bindings to a class.
//------------------------------------------------------------------------------
template<typename Obj, typename PB11Obj>
void restartMethodBindings(py::module& m, PB11Obj& obj) {
  obj
    .def("label", (std::string (Obj::*)() const) &Obj::label)
    .def("dumpState", (void (Obj::*)(FileIOSpace::FileIO&, const std::string&) const) &Obj::dumpState, "file"_a, "pathName"_a)
    .def("restoreState", (void (Obj::*)(const FileIOSpace::FileIO&, const std::string&)) &Obj::restoreState, "file"_a, "pathName"_a)
    ;
}

}

#endif

