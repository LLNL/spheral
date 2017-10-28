// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <map>
#include <set>
#include <string>

#include <limits.h>

namespace py = pybind11;

// vector
PYBIND11_MAKE_OPAQUE(std::vector<char>);
PYBIND11_MAKE_OPAQUE(std::vector<unsigned>);
PYBIND11_MAKE_OPAQUE(std::vector<uint64_t>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);

// vector<vector<>>
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<char>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<unsigned>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<uint64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::string>>);

typedef std::pair<double, double> pair_double_double;
typedef std::pair<double, std::string> pair_double_string;
typedef std::pair<unsigned, unsigned> pair_unsigned_unsigned;
typedef std::pair<uint64_t, uint64_t> pair_ULL_ULL;
typedef std::pair<std::string, std::string> pair_string_string;
PYBIND11_MAKE_OPAQUE(std::vector<pair_double_double>);
PYBIND11_MAKE_OPAQUE(std::vector<pair_double_string>);
PYBIND11_MAKE_OPAQUE(std::vector<pair_unsigned_unsigned>);
PYBIND11_MAKE_OPAQUE(std::vector<pair_ULL_ULL>);
PYBIND11_MAKE_OPAQUE(std::vector<pair_string_string>);

// map
typedef std::map<std::string, double> map_string_double;
typedef std::map<int, std::string> map_int_string;
PYBIND11_MAKE_OPAQUE(map_string_double);
PYBIND11_MAKE_OPAQUE(map_int_string);

// Make the module
PYBIND11_PLUGIN(SpheralCXXTypes) {
  py::module m("SpheralCXXTypes", "Spheral C++ types module.");

  // vector
  py::bind_vector<std::vector<char>>(m, "vector_of_char");
  py::bind_vector<std::vector<unsigned int>>(m, "vector_of_unsigned");
  py::bind_vector<std::vector<uint64_t>>(m, "vector_of_ULL");
  py::bind_vector<std::vector<int>>(m, "vector_of_int");
  py::bind_vector<std::vector<float>>(m, "vector_of_float");
  py::bind_vector<std::vector<double>>(m, "vector_of_double");
  py::bind_vector<std::vector<std::string>>(m, "vector_of_string");

  // vector<vector<>>
  py::bind_vector<std::vector<std::vector<char>>>(m, "vector_of_vector_of_char");
  py::bind_vector<std::vector<std::vector<unsigned int>>>(m, "vector_of_vector_of_unsigned");
  py::bind_vector<std::vector<std::vector<uint64_t>>>(m, "vector_of_vector_of_ULL");
  py::bind_vector<std::vector<std::vector<int>>>(m, "vector_of_vector_of_int");
  py::bind_vector<std::vector<std::vector<float>>>(m, "vector_of_vector_of_float");
  py::bind_vector<std::vector<std::vector<double>>>(m, "vector_of_vector_of_double");
  py::bind_vector<std::vector<std::vector<std::string>>>(m, "vector_of_vector_of_string");

  // vector<pair>
  py::bind_vector<std::vector<std::pair<double, double>>>(m, "vector_of_pair_double_double");
  py::bind_vector<std::vector<std::pair<double, std::string>>>(m, "vector_of_pair_double_string");
  py::bind_vector<std::vector<std::pair<unsigned, unsigned>>>(m, "vector_of_pair_unsigned_unsigned");
  py::bind_vector<std::vector<std::pair<uint64_t, uint64_t>>>(m, "vector_of_pair_ULL_ULL");
  py::bind_vector<std::vector<std::pair<std::string, std::string>>>(m, "vector_of_pair_string_string");

  // map
  py::bind_map<std::map<std::string, double>>(m, "map_string_double");
  py::bind_map<std::map<int, std::string>>(m, "map_int_string");

  return m.ptr();
}
