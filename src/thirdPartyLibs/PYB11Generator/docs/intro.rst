.. _intro:

Introduction to PYB11Generator
==============================

Using PYB11Generator to write the interface for a C++ binding is intended to emulate writing that same interface in Python, so if you're familiar with Python it should be easy to get started with PYB11Generator.  As an example, if we have a header for a C++ function that looks like::

   // A really important function!
   int func();

We can define the PYB11Generator prescription for binding this method by writing a python method as::

  def func():
      "A really important function!"
      return "int"

Wherever possible we try to use ordinary python syntax to correspond to pybind11/C++ constructs: python functions correspond to and generate binding code for C++ functions as above; a python class generates binding code for pybind11 to bind a C++ class; arguments for functions and methods in python generate corresponding argument specifications in C++ function pointer syntax.  Because Python is not a strongly typed language, we specify C++ types using strings (if needed) as above, where we specify the return ``int`` type by returning the string ``"int"`` from ``func``.  We also use Python decorators to annotate Python methods with uniquely C++ concepts such as ``const``, ``virtual``, etc., as will be discussed in succeeding sections.

.. _first-example:

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

creates a file ``example.cc``, which looks like (omitting the boilerplate preamble code with ``#include``'s)::

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

A few things worth noting in this example:

* This example uses the fact that if the function being wrapped is unambiguous, the C++ function pointer does not need the full explicit function prescription.  This is reflected in the PYB11Generator syntax when we write the ``def add()`` function in python without arguments or a return type.
* In order to directly insert the C++ function definition into the resulting C++ file, we have used the special variable ``PYB11preamble`` variable.  A more typical use case will require ``#include`` ing the necessary C++ header files in the generated code, which is accomplished through another special variable, ``PYB11includes``, described later.

  * In general special variables and commands to PYB11Generator use the prefix ``PYB11`` such as ``PYB11preamble`` in this example.

* Note also that ordinary Python doc strings (both for the module and function) are picked up from ``simple_example.py`` and propagated to the pybind11 bindings.

This example demonstrates the steps necessary to create a usable python module using PYB11Generator:

* Create a python file describing the desired interface using ordinary python syntax, based on the C++ methods and classes to be bound.
* Run a python line like above to generate the ``pybind11`` C++ code from this python input.
* Compile the resulting ``pybind11`` C++ code to create the python shared module.

In the following sections we describe the nuances of creating the PYB11 python input files in much more detail; we will not show the compilation examples beyond this point since it is no different than using ``pybind11`` directly, and the above example pretty much covers it.
