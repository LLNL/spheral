.. _functions:

Binding functions
=================

We have already introduced a quick example of binding a function in :ref:`first-example`; this section will go through much more detail on how to generate ``pybind11`` bindings for functions, including complications such as overloaded methods and C++ templates.

.. _functions-overload:

Ordinary and overloaded functions
---------------------------------

Suppose we have a header defining the following functions that we wish to bind for Python::

  double unique_function(int a, double b);
  double overloaded_function(int a, int b, int c);
  double overloaded_function(double a, double b, double c);

We can use PYB11Generator to bind these functions with a file containing the following code::

  def unique_function():
      "This is a unique function prescription, and so requires no details about arguments or return types"
      return

  def overloaded_function(a = "int",
                          b = "int",
                          c = "int"):
      "This is the version of overloaded_function that takes ints"
      return "double"

  @PYB11pycppname("overloaded_function")
  def overloaded_function1(a = "double",
                           b = "double",
                           c = "double"):
      "This is the version of overloaded_function that takes doubles"
      return "double"

The first function ``unique_function`` is trivial, since it is unambiguous and can be wrapped with an unadorned C++ function pointer as shown in :ref:`first-example`.  In this case PYB11Generator assumes the C++ function name is the same as the Python function name, and all is simple.

The overloaded functions take a bit more work.  The first challenge we have is that Python does not support the concept of function overloading: two Python functions cannot have the same name.  Therefore we need to use unique Python names for the C++ ``overloaded_function`` Python descriptions, and therefore we define ``overloaded_function`` and ``overloaded_function1`` in the python source.  In order to tell PYB11Generator that we really want to call ``overloaded_function1`` ``overloaded_function`` in both the C++ and python bindings, we have used our first PYB11 decorator: ``PYB11pycppname``.  This decorator tells PYB11Generator that that function in question is really called ``overloaded_function`` in C++, and we wish the python name in the resulting binding code to also call this function ``overloaded_function`` in Python as well.  This is actually two statements, and there are two PYB11 decorators that can do these individual tasks independently if needed (``PYB11cppname`` and ``PYB11pyname``): ``PYB11pycppname`` is simply a convenient shorthand combination to cover the common case of wanting to simulataneously rename the bound method for C++ and Python.  For a full listing of the PYB11 decorators see :ref:`decorators`.

Note we have also now specified the arguments and return types for both bindings of ``overloaded_function``.  This is required since the C++ functions are overloaded, and in order for the C++ compiler to distinguish which one we want it is necessary to fully prescribe the function signatures for the function pointers in the ``pybind11`` binding code.  PYB11Generator always checks the return value for a wrapped function: if a return value is present, it should be a string describing the C++ return type (as shown here, with both ''overloaded_function`` and ``overloaded_function1`` returning the string value ``"double"``).  If such a return value is specified PYB11Generator assumes a fully qualified C++ function pointer signature is required, and will also look for and generate the argument types as well.  The function arguments should be named what the argument name will be in the resulting Python code, and set equal to a string with the C++ type of the argument as shown above for the ``overloaded_function`` descriptions.  Note, a C++ ``void`` return value or argument should be set to the string ``"void"`` for PYB11Generator for such explicit specifications.
