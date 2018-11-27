.. _classes:

Binding classes
===============

Binding classes in PYB11Generator is based on writing the desired interface as a Python class, similar to the process for :ref:`functions`.  As a first example consider the example struct used as the first such example in the `pybind11 documentation <https://pybind11.readthedocs.io/en/stable/classes.html>`_::

  struct Pet {
    Pet(const std::string &name) : name(name) { }
    void setName(const std::string &name_) { name = name_; }
    const std::string &getName() const { return name; }

    std::string name;
  };

This struct can be wrapped in straighforward fashion in PYB11Generator as::

  class Pet:

      def pyinit(self,
                 name = "const std::string&"):
         return

      def setName(self,
                 name = "const std::string&"):
          return "void"

      def getName(self):
          return "const std::string"

Processing this Python class definition through PYB11Generator results in the following (omitting generic preamble code)::

  // Class Pet
  {
    py::class_<Pet> obj(m, "Pet");

    // Constructors
    obj.def(py::init<const std::string&>(), "name"_a);

    // Methods
    obj.def("setName", (void (Pet::*)(const std::string&)) &Pet::setName, "name"_a);
    obj.def("getName", (const std::string (Pet::*)()) &Pet::getName);
  }

which is very similar to the native pybind11 `pybind11 code <https://pybind11.readthedocs.io/en/stable/classes.html>`_.  This example demonstrates a few important aspects of generating class bindings with PYB11Generator:

* A python class results in the generation of a pybind11 ``class_<>`` declaration.

* Binding class methods with PYB11Generator is directly analogous to binding free functions: we write the method signature in python syntax, with the arguments set equal to the C++ type as a string.

  * If the C++ class method is unambiguous (not overloaded), then just as with functions we can specify the method in python with no arguments and an empty return value.

  * If a default value for an argument is desirable, simply set the argument equal to a tuple of two strings: ``arg = ("C++ type", "C++ default value")``, identically to the treatment of functions in :ref:`functions-default-args`.

* Constructors are specified by any class method starting with the string ``pyinit``.
