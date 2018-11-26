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

The overloaded functions take a bit more work.  The first challenge is that Python does not support the concept of function overloading: two Python functions cannot have the same name.  Therefore we need to use unique Python names for the C++ ``overloaded_function`` Python descriptions, and therefore we define ``overloaded_function`` and ``overloaded_function1`` in the python source.  In order to tell PYB11Generator that we really want to call ``overloaded_function1`` ``overloaded_function`` in both the C++ and python bindings, we have used our first PYB11 decorator: ``PYB11pycppname``.  This decorator tells PYB11Generator that that function in question is really called ``overloaded_function`` in C++, and we wish the python name in the resulting binding code to also call this function ``overloaded_function`` in Python as well.  This is actually two statements, and there are two PYB11 decorators that can do these individual tasks independently if needed (``PYB11cppname`` and ``PYB11pyname``): ``PYB11pycppname`` is simply a convenient shorthand combination to cover the common case of wanting to simulataneously rename the bound method for C++ and Python.  For a full listing of the PYB11 decorators see :ref:`decorators`.

Note we have also now specified the arguments and return types for both bindings of ``overloaded_function``.  This is required since the C++ functions are overloaded, and in order for the C++ compiler to distinguish which one we want it is necessary to fully prescribe the function signatures for the function pointers in the ``pybind11`` binding code.  PYB11Generator always checks the return value for a wrapped function: if a return value is present, it should be a string describing the C++ return type (as shown here, with both ``overloaded_function`` and ``overloaded_function1`` returning the string value ``"double"``).  If such a return value is specified PYB11Generator assumes a fully qualified C++ function pointer signature is required, and will also look for and generate the argument types as well.  The function arguments should be named what the argument name will be in the resulting Python code, and set equal to a string with the C++ type of the argument as shown above for the ``overloaded_function`` descriptions.  Note, a C++ ``void`` return value or argument should be set to the string ``"void"`` for PYB11Generator for such explicit specifications.

.. _functions-default-args:

Default argument values
-----------------------

Another useful feature of ``pybind11`` is the ability to specify default values for arguments to functions/methods in Python, and naturally PYB11Generator supports this feature as well.  In order to specify a default value for an argument, we set the value of the argument in the python binding code as a tuple, where the first element is a string describing the C++ type, and the second a string with the C++ default value.  As an example suppose we wish to bind the following function that has two arguments (an ``int`` and a ``std::string``)::

  void howToDrawADragon(int numberOfBeefyArms, std::string label name);

and we want to use the default values ``1`` and ``"Trogdor"`` for these arguments.  The PYB11Generator code would then look like::

  def howToDrawADragon(numberOfBeefyArms = ("int", "1"),
                       name = ("std::string", "Trogdor")):
      return "void"

.. _functions-template:

C++ template functions
----------------------

C++ templates present another challenge, as this another concept not found in Python.  Suppose we wish to expose several instantiations of the following method::

  template<typename ValueA, typename ValueB, typename ValueC>
  ValueC
  transmogrify(const ValueA& x, const ValueB& y);

It is always possible to explicitly (and repetitively) define the function over and over again for each template instantiation combination of (``ValueA``, ``ValueB``, ``ValueC``), but we would rather write the prescription once and have the computer generate the necessary redundant code.  PYB11Generator has such a facility: a template method can be defined with the ``@PYB11template`` decorator, which takes the template arguments as a set of string arguments.  The function can then be instantiated as many times as needed using the function ``PYB11TemplateFunction``.  The complete PYB11Generator binding code then might look like::

  from PYB11Generator import *     # Necessary to get decorators and PYB11TemplateFunction

  @PYB11template("ValueA", "ValueB", "ValueC")
  def transmogrify(x = "const ValueA&",
                   y = "const ValueB&"):
      "I'm sure this does something useful..."
      return "ValueC"

  transmogrifyIntIntDouble = PYB11TemplateFunction(transmogrify, ("int", "int", "double"),             pyname="transmogrify")
  transmogrifyI32I32I64    = PYB11TemplateFunction(transmogrify, ("uint32_t", "uint32_t", "uint64_t"), pyname="transmogrify")

Note we have made these different instantiations overloaded in python by forcing them all to have the name ``transmogrify`` via the ``pyname="transmogrify"`` argument.  This is not necessarily required: we must give each instantiation of the template a unique name in Python (``transmogrifyIntIntDouble`` and ``transmogrifyI32I32I64`` in this case), and if we are happy with those being the Python names of the wrapped results we would not need to specify ``pyname``.  Such unique names in Python are safest, in that which instantiation the user wants to call down the line in the wrapped library call is unambiguous, but often it is nicer to force the Python names to match the C++ as we do in this case.  It is also notable that functions decorated by ``@PYB11template`` do not generate any binding code by themselves -- PYB11Generator ignores these functions until a template instantiation is created with ``PYB11TemplateFunction``.

The full list of allowed arguments to ``PYB11TemplateFunction`` is::

  PYB11TemplateFunction(func_template, template_parameters, cppname=None, pyname=None, pyext="")

``func_template``
  The function description decorated by ``@PYB11template``.

``template_parameters``
  A tuple of C++ strings, one for each of the template parameters specified in the template function spec of ``@PYB11template``.

``cppname``
  Optional -- override the C++ name of the function.  Defaults to the name of ``func_template``.

``pyname``
  Optional -- override the Python name of the wrapped function.  Defaults to the python name of the instantion (``transmogrifyIntIntDouble`` in the first case above had we not specified ``pyname``)

``docext``
  Optional -- a string to tack onto the documentation string specified in ``func_template``, if any.
