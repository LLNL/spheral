#include <vector>
#include <map>
#include <set>
#include <string>

#include <limits.h>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"

namespace py = pybind11;

// vector
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);

// map
typedef std::map<std::string, double> map_string_double;
typedef std::map<int, std::string> map_int_string;
PYBIND11_MAKE_OPAQUE(map_string_double);
PYBIND11_MAKE_OPAQUE(map_int_string);

// Make the module
PYBIND11_PLUGIN(SpheralCXXTypes) {
  py::module m("SpheralCXXTypes", "Spheral C++ types module.");

  // vector
  py::bind_vector<std::vector<int>>(m, "vector_of_int");
  py::bind_vector<std::vector<double>>(m, "vector_of_double");
  py::bind_vector<std::vector<std::string>>(m, "vector_of_string");

  // map
  py::bind_map<std::map<std::string, double>>(m, "map_string_double");
  py::bind_map<std::map<int, std::string>>(m, "map_int_string");

  return m.ptr();
}
