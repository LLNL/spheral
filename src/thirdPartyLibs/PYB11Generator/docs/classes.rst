.. _classes:

===============
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

------------------
Class constructors
------------------

In general PYB11Generator interprets methods of classes as ordinary methods to exposed via pybind11 -- the one exception to this rule is class constructors.  Any method that begins with the name ``pyinit`` is interpreted as a class constructor, allowing the specification of an arbitrary number of constructors.  For instance, if we have a C++ class with the following constructors::

  class A {
  public:
    A();                                             // Default constructor
    A(const std::string name);                       // Build with a name, default priority
    A(const std::string name, const int priority);   // Build with a name and priority
  };

We can bind these three different constructors using the following Python specification::

  class A:
      "A class that does something with a string and an int..."

      def pyinit(self):
          "Default constructor"

      def pyinit1(self, name="const std::string"):
          "Build with a name, default priority"

      def pyinit2(self, name="const std::string", priority="const int"):
          "Build with a name and priority"

.. _class-methods:

-------------
Class methods
-------------

Class methods are wrapped much like free functions using PYB11Generator: we simply define a python class method with the desired name.  If the method is unambiguous (not overloaded), we do not necessarily have to specify the return types and arguments (though full specifications are always allowed, and at times preferable to generate more explicit help in Python).  The syntax for specifying C++ return types and arguments for methods is identical to that used for for :ref:`functions`, as is evident in the examples below.

.. _overloaded-class-methods:

Overloaded class methods
------------------------

Just as with :ref:`function-overloads`, overloaded methods require full call specifications, as well as unique names in python.  We use the PYB11 decorators ``@PYB11pyname``, ``@PYB11cppname``, or ``@PYB11pycppname`` to link the proper C++/Python names as needed.  As an example, consider the following C++ class:

.. code-block:: cpp

    class A {
    public:
      int process(const int x);                     // Process the internal state somehow to answer this query
      std::string label();                          // Return a string label
      std::string label(const std::string suffix);  // Return a string label including a specified suffix
    }

In this case we have one unambiguous method (``process``), and two overloaded methods (``label``).  We can write PYB11Generator bindings for these methods as::

  class A:

      def process(self):
          "Process the internal state somehow to answer this query"
          return

      def label(self):
          "Return a string label"
          return "std::string"

      @PYB11pycppname("label")
      def label1(self, suffix="const std::string"):
          "Return a string label including a specified suffix"
          return "std::string"

We have chosen to bind the unambiguous ``A::process`` method using no method signature (i.e., no return type or arguments) for brevity.  The overloaded ``A::label`` methods however require the complete method prescriptions be specified in order for the compiler to know which C++ ``A::label`` we are referring to.  Because Python does not allow class methods with the same name however, we must use unique method names in our Python class binding (hence ``A.label`` and ``A.label1``).  We use the PYB11 decorator ``@PYB11pycppname`` on ``A.label1`` to indicate we want the bound Python and C++ names to be ``label``.   This is identical to how this overloading problem is handled for :ref:`function-overloads`.

.. Note::

   In this example we have made the typical choice to overload the ``label`` method in Python just as in C++.  We could, however, decide to leave the Python ``label`` and ``label1`` methods with unique names, removing the unpythonic overloading concept from the python interface.  If we want to leave the Python name of the second binding of ``A::label`` as ``A.label1``, we still need to tell PYB11Generator that the C++ name is ``A::label`` rather than ``A::label1``.  In this case we would simply change the decorator to specify the C++ name alone::

      @PYB11cppname("label")
      def label1(self, suffix="const std::string"):
          "Return a string label including a specified suffix"
          return "std::string"

.. _const-methods:

Const class methods
-------------------

Const'ness is a concept in C++ not shared by Python, so we use a decorator (``@PYB11const``) to denote a const method when needed.  For instance, the following C++ class definition:

.. code-block:: cpp

  class A {
  public:
    int square(const int x) const { return x*x; }  // Return the square of the argument
  };

can be specfied in PYB11 using::

  class A:

      @PYB11const
      def square(self, x="const int"):
          "Return the square of the argument"
          return "int"

.. _protected-methods:

Protected class methods
-----------------------

It is possible to bind protected class methods in pybind11 as described in `the pybind11 documentation <https://pybind11.readthedocs.io/en/stable/advanced/classes.html#binding-protected-member-functions>`_.  In the pybind11 code this requires writing an intermediate C++ class to publish the protected methods.  PYB11Generator automates the production of such publisher classes as needed, however, so all that is required to expose a protected class method is to decorate the PYB11 binding with ``@PYB11protected``.  In order to expose the protected method of the following example:

.. code-block:: cpp

  class A {
  protected:
    void some_protected_method(const int x);     // A protected method to apply x->A somehow
  }

we simply provide a decorated PYB11 binding as::

  class A:

      @PYB11protected
      def some_protected_method(self, x="int"):
          "A protected method to apply x->A somehow"
          return "void"

.. _static-methods:

Static class methods
--------------------

Static C++ methods are denoted to PYB11Generator using the ``@PYB11static`` decorator as in the following example.

C++ class with a static method:

.. code-block:: cpp

  class A {
  public:
    static int func(int x);    // This method does something with x
  };

PYB11 binding code::

  class A:

      @PYB11static
      def func(x = "int"):
          "This method does something with x"
          return "int"

.. _class-inheritance:

-----------------
Class inheritance
-----------------

Class inheritance hierarchies in C++ are simple to reflect in PYB11Generator, as this is an OO concept shared by both C++ and Python: all that is required is to reflect the inheritance hierarchy in the Python PYB11 code.  In order to expose the following C++ classes:

.. code-block:: cpp

  class A {
    A();                    // Default constructor
    int func(int x);        // Some useful function of A
  };

  class B: public A {
    B();                    // Default constructor
    double dfunc(double x); // Some useful function of B
  };

we can simply reflect this object hiearchy in the PYB11Generator code::

  class A:

      def pyinit(self):
          "Default constructor"

      def func(self, x="int"):
          "Some useful function of A"
          return "int"

  class B(A):

      def pyinit(self):
          "Default constructor"

      def dfunc(self, x="double"):
          "Some useful function of B"
          return "double"
