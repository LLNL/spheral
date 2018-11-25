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
