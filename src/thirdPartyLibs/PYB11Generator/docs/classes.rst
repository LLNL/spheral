.. _classes:

Binding classes
===============

Binding classes in PYB11Generator is based on writing the desired interface as a Python class, similar to the process for :ref:`functions`.  As a first example consider the example struct used as the first such example in the pybind11 class documentation :ref:`pybind11:classes`:

.. code-block:: cpp

  struct Pet {
    Pet(const std::string &name) : name(name) { }
    void setName(const std::string &name_) { name = name_; }
    const std::string &getName() const { return name; }

    std::string name;
  };

This struct can be wrapped in straighforward fashion in PYB11Generator as:

.. code-block:: python

  class Pet:

      def pyinit(self,
                 name = "const std::string&"):
         return

      def setName(self,
                 name = "const std::string&"):
          return "void"

      def getName(self):
          return "const std::string"

Processing this Python class definition through PYB11Generator results in the following (omitting generic preamble code):

.. code-block:: cpp

  // Class Pet
  {
    py::class_<Pet> obj(m, "Pet");

    // Constructors
    obj.def(py::init<const std::string&>(), "name"_a);

    // Methods
    obj.def("setName", (void (Pet::*)(const std::string&)) &Pet::setName, "name"_a);
    obj.def("getName", (const std::string (Pet::*)()) &Pet::getName);
  }

which is very similar to the native pybind11 code presented in :ref:`pybind11:classes`.  This example demonstrates a few important aspects of generating class bindings with PYB11Generator:

* A python class results in the generation of a pybind11 ``class_<>`` declaration.

* Binding class methods with PYB11Generator is directly analogous to binding free functions: we write the method signature in python syntax, with the arguments set equal to the C++ type as a string.

  * If the C++ class method is unambiguous (not overloaded), then just as with functions we can specify the method in python with no arguments and an empty return value.

  * If a default value for an argument is desirable, simply set the argument equal to a tuple of two strings: ``arg = ("C++ type", "C++ default value")``, identically to the treatment of functions in :ref:`functions-default-args`.

* Constructors are specified by any class method starting with the string ``pyinit``.

.. _class-constructors:

Class constructors
------------------

In general PYB11Generator interprets methods of classes as ordinary methods to exposed via pybind11 -- the one exception to this rule is class constructors.  Any method that begins with the name ``pyinit`` is interpreted as a class constructor, allowing the specification of an arbitrary number of constructors.  For instance, if we have a C++ class with the following constructors:

.. code-block:: cpp

  class A {
  public:

    A();                                             // Default constructor
    A(const std::string name);                       // Build with a name, default priority
    A(const std::string name, const int priority);   // Build with a name and priority
  };

We can bind these three different constructors using the following Python specification:

.. code-block:: py

  class A:

      def pyinit(self):
          "Default constructor"

      def pyinit1(self, name="const std::string"):
          "Build with a name, default priority"

      def pyinit2(self, name="const std::string", priority="const int"):
          "Build with a name and priority"

.. _class-methods:

Class methods
-------------

Class methods are wrapped much like free functions using PYB11Generator: we simply define a python class method with the desired name.

* If the method is unambiguous (not overloaded), we do not necessarily have to specify the return types and arguments (though full specifications are always allowed, and at times preferable to generate more explicit help in Python).

* Just as with :ref:`function-overloads`, Overloaded methods require full call specifications, as well as unique names in python.  We use the PYB11 decorators ``@PYB11pyname``, ``@PYB11cppname``, or ``@PYB11pycppname`` to link the proper C++/Python names as needed.
