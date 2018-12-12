#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/operators.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"

namespace py = pybind11;
using namespace pybind11::literals;

#include <iostream>
#include <vector>

PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>)

void printStuff(std::vector<std::vector<int>>& stuff) {
  std::cout << "{\n";
  for (auto& vec: stuff) {
    std::cout << "    {";
    for (auto x: vec) std::cout << " " << x;
    std::cout << " }\n";
  }
  std::cout << "}\n";
}

void modifyStuff(std::vector<std::vector<int>>& stuff) {
  stuff.clear();
  stuff.push_back({1,2,3});
  stuff.push_back({10, 20, 30, 40, 50});
}

PYBIND11_MODULE(test_vector_vector, m) {
  m.def("printStuff", &printStuff);
  m.def("modifyStuff", &modifyStuff);
  py::bind_vector<std::vector<int>>(m, "vector_of_int", py::module_local(false));
  py::bind_vector<std::vector<std::vector<int>>>(m, "vector_of_vector_of_int", py::module_local(false));
}
