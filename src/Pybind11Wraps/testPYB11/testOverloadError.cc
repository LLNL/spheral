#include "pybind11/pybind11.h"
#include "pybind11/functional.h"
#include "pybind11/operators.h"

namespace py = pybind11;
using namespace pybind11::literals;

#include <iostream>

class A {
public:
  virtual void func()               { std::cout << "A::func()" << std::endl; }
  virtual void func(int x) = 0;
};

class B: public A{
public:
  virtual void func(int x) override { std::cout << "B::func(" << x << ")" << std::endl; }
};

class Atrampoline: public A {
  using A::A;
  virtual void func() override      { PYBIND11_OVERLOAD     (void, A, func,); }
  virtual void func(int x) override { PYBIND11_OVERLOAD_PURE(void, A, func, x); }
};

class Btrampoline: public B {
  using B::B;
  virtual void func() override      { PYBIND11_OVERLOAD     (void, A, func,); }            // <--- ugliness I'm not sure is right!
  virtual void func(int x) override { PYBIND11_OVERLOAD     (void, B, func, x); }
};

PYBIND11_MODULE(testOverloadError, m) {
  py::class_<A, Atrampoline>(m, "A")
    .def(py::init<>())
    .def("func", (void (A::*)()) &A::func)
    .def("func", (void (A::*)(int)) &A::func)
    ;

  py::class_<B, Btrampoline>(m, "B")
    .def(py::init<>())
    .def("func", (void (B::*)()) &B::func)                   // <-- this shouldn't be necessary, but overload doesn't show up in python otherwise
    .def("func", (void (B::*)(int)) &B::func)
    ;
}
