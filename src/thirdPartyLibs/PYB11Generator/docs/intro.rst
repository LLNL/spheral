Introduction to PYB11Generator
==============================

Using PYB11Generator to write the interface for a C++ binding is intended to emulate writing that same interface in Python, so if you're familiar with Python it should be easy to get started with PYB11Generator.  As an example, if we have a header for a C++ function that looks like::

   // A really important function!
   int func();

We can define the PYB11Generator prescription for binding this method by writing a python method as::

  def func():
      "A really important function!"
      return "int"

Wherever possible we try to use orindary python syntax to correspond to pybind11/C++ constructs: python functions correspond to and generate binding code for C++ functions as above; a python class generates binding code for pybind11 to bind a C++ class; arguments for functions and methods in python generate corresponding argument specifications in pybind11 function pointer syntax.  Because Python is not a strongly typed language, we specify C++ types using strings (if needed).  We also use Python decorators to annotate Python methods with uniquely C++ concepts such as ``const``, ``virtual``, etc.

-------------------------------
A first example start to finish
-------------------------------

To explicitly demonstrate the stages of creating bindings using PYB11Generator, here we recreate an early example in the pybind11 documentation: `Creating bindings for a simple function <https://pybind11.readthedocs.io/en/stable/basics.html>`_.  Say we want to wrap the following C++ method::

  int add(int i, int j) {
    return i + j;
  }

We can use PYB11Generator to create the same pybind11 code used to wrap this method in the pybind11 tutorial by writing a file ``simple_example.py`` containing::

  "pybind11 example plugin"
  
  PYB11preamble = """
  int add(int i, int j) {
    return i + j;
  }
  """
  
  def add():
      "A function which adds two numbers"
      return

Now executing the command::

  python -c 'from PYB11Generator import *; import simple_example; PYB11generateModule(simple_example, "example")'

creates a file ``example.cc`` with (ignoring the boilerplate preamble code with ``#include``'s)::

  int add(int i, int j) {
    return i + j;
  }

  //------------------------------------------------------------------------------
  // Make the module
  //------------------------------------------------------------------------------
  PYBIND11_MODULE(example, m) {

    m.doc() = "pybind11 example plugin"  ;

    //...........................................................................
    // Methods
    m.def("add", &add, "A function which adds two numbers");
  }

This is identical to the pybind11 code shown for this case in the `pybind11 tutorial <https://pybind11.readthedocs.io/en/stable/basics.html>`_, modulo some comments.  This code can now be compiled to the final Python shared module as described in the pybind11 tutorial::

  $ c++ -O3 -Wall -shared -std=c++11 -fPIC `python -m pybind11 --includes` example.cc -o example.so
