.. _classes:

=======
Classes
=======

Binding classes in PYB11Generator is based on writing the desired interface as a Python class, similar to the process for :ref:`functions`.  As a first example consider the example struct used as the first such example in the pybind11 class documentation :ref:`pybind11:classes`:

.. code-block:: cpp

  struct Pet {
    Pet(const std::string &name) : name(name) { }
    void setName(const std::string &name_) { name = name_; }
    const std::string &getName() const { return name; }

    std::string name;
  };

This struct can be wrapped in straightforward fashion in PYB11Generator as:

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

------------
Constructors
------------

In general PYB11Generator interprets methods of classes as ordinary methods to exposed via pybind11 -- the one exception to this rule is class constructors.  Any method that begins with the name ``pyinit`` is interpreted as a class constructor, allowing the specification of an arbitrary number of constructors.  For instance, if we have a C++ class with the following constructors:

.. code-block:: cpp

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

For constructors it does not matter what names are used past the ``pyinit`` string: any such name will be interpreted as a constructor.  All that is required is that any class ``pyinit*`` name be unique -- remember, python does not allow overloading, so defining successive methods with the same name simply causes the earlier method definitions to be lost.  Not that the author has made such mistakes in creating my own binding code...

.. _class-inheritance:

-----------
Inheritance
-----------

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

.. Note::

   Cross module inheritance (binding a class in one module that inherits from a class bound in another) is a slightly trickier case.  See the discussion in :ref:`cross-module-inheritance` for an example of how to do this.

.. Note::

   Another esoteric case is having a non-templated class inherit from a templated one.  A method of handling this situation is discussed in :ref:`non-template-to-template-inheritance`.

.. _class-methods:

-------
Methods
-------

Class methods are wrapped much like free functions using PYB11Generator: we simply define a python class method with the desired name.  If the method is unambiguous (not overloaded), we do not necessarily have to specify the return types and arguments (though full specifications are always allowed, and at times preferable to generate more explicit help in Python).  The syntax for specifying C++ return types and arguments for methods is identical to that used for for :ref:`functions`, as is evident in the examples below.

.. _overloaded-class-methods:

Overloaded methods
------------------

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

Const methods
-------------

Const'ness is a concept in C++ not shared by Python, so we use a decorator (``@PYB11const``) to denote a const method when needed.  For instance, the following C++ class definition:

.. code-block:: cpp

  class A {
  public:
    int square(const int x) const { return x*x; }  // Return the square of the argument
  };

can be specified in PYB11 using::

  class A:

      @PYB11const
      def square(self, x="const int"):
          "Return the square of the argument"
          return "int"

.. _virtual-methods:

Virtual methods
---------------

If we simply wish to expose C++ virtual methods as ordinary class methods in Python (i.e., not allowing overriding the implementation of such methods from Python), then nothing extra need be done in the method binding for PYB11.  However, in pybind11 it is also possible to expose C++ virtual methods such that they *can* be overridden from Python descendants, which is a very powerful capability.  Exposing such overridable virtual methods in pybind11 involves writing an intermediate "trampoline" class as described in the pybind11 documentation :ref:`pybind11:overriding_virtuals`.  PYB11Generator automates the generation of such intermediate redundant code (this was in fact the motivating factor in the creation of PYB11Generator), removing much of the bookkeeping necessary to maintain such coding in face of a changing interface.  In PYB11Generator all that is required for making a virtual method overridable from Python is decorating such virtual methods with ``@PYB11virtual``/``@PYB11pure_virtual`` as appropriate.  Consider binding the C++ example from the pybind11 documentation :ref:`pybind11:overriding_virtuals`:

.. code-block:: cpp

    class Animal {
    public:
        virtual ~Animal() { }
        virtual std::string go(int n_times) = 0;
    };

    class Dog : public Animal {
    public:
        virtual std::string go(int n_times) override {
            std::string result;
            for (int i=0; i<n_times; ++i)
                result += "woof! ";
            return result;
        }
    };

All that is necessary to bind this code using PYB11Generator is the following::

  class Animal:

      def pyinit(self):
          "Default constructor"

      @PYB11pure_virtual
      def go(self, n_times="int"):
          return "std::string"

  class Dog(Animal):

      def pyinit(self):
          "Default constructor"

      @PYB11virtual
      def go(self, n_times="int"):
          return "std::string"

Now both ``Animal`` and ``Dog`` are accessible from Python, and PYB11Generator automatically generates the necessary trampoline classes such that the ``go`` method can be overridden by descendant Python classes as desired.  Note the use of the PYB11 decorators: ``PYB11virtual`` and ``PYB11pure_virtual``.  The use of these two should evident from their names and uses in this example:

* ``PYB11virtual`` decorates C++ methods that are virtual (such as ``Dog::go``).

* ``PYB11pure_virtual`` decorates C++ methods are pure virtual (such as ``Animal::go``), marking such classes as abstract.

.. _protected-methods:

Protected methods
-----------------

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

Static methods
--------------

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

.. _class-operators:

-----------------------------------
Special class operators and methods
-----------------------------------

Python has a number of `special methods for classes <https://docs.python.org/2/reference/datamodel.html#special-method-names>`_, such as ``__len__``, ``__add__``, etc., which allow the object behavior to be controlled for operations such as +, +=, ``len()``, and so forth.  pybind11 supports `these operators <https://pybind11.readthedocs.io/en/stable/advanced/classes.html#operator-overloading>`_, so naturally PYB11Generator does as well.  In keeping with PYB11Generators interface, these are specified by providing these special method names in your Python class description.

Numeric operators
-----------------

The numeric operators supported by PYB11Generator are ``__add__``, ``__sub__``, ``__mul__``, ``__div__``, ``__mod__``, ``__and__``, ``__xor__``, ``__or__``, ``__radd__``, ``__rsub__``, ``__rmul__``, ``__rdiv__``, ``__rmod__``, ``__rand__``, ``__rxor__``, ``__ror__``, 
``__iadd__``, ``__isub__``, ``__imul__``, ``__idiv__``, ``__imod__``, ``__iand__``, ``__ixor__``, ``__ior__``, ``__neg__``, and ``__invert__``.  Python descriptions of these methods are available `here <https://docs.python.org/2/reference/datamodel.html#emulating-numeric-types>`_.

In the common case for binary operators where the argument is of the same type as the class we're binding, we can omit the the argument specification and return type.  However, in the case where the binary operator accepts a different C++ type, we need to specify this argument type in the usual PYB11 syntax for arguments and return types.

It is also important to remember that Python does not allow us to define a method name more than once in a class, so if we have overloaded C++ math operators (say ``operator+`` can accept more than one type), we must give each binding a unique name, but then use decorators such as ``@PYB11pyname`` to force the special operator name for the method.

As an example, consider the following C++ class which supports addition with itself or a ``double``, multiplication by a ``double``, and the unary negative operator:

.. code-block:: cpp

  class Vector3d {
  public:
     Vector3d  operator-() const;

     Vector3d& operator+=(const Vector3d& rhs);
     Vector3d  operator+ (const Vector3d& rhs) const;

     Vector3d& operator+=(const double rhs);
     Vector3d  operator+ (const double rhs) const;

     Vector3d& operator*=(const double rhs);
     Vector3d  operator* (const double rhs) const;
  };

We can bind these numeric operations for the Python version of ``Vector3d`` with PYB11Generator using normal Python operator syntax::

  class Vector3d:

      def __neg__(self):
          return

      def __iadd__(self):
          return

      def __add__(self):
          return

      @PYB11pyname("__iadd__")
      def __iadd__double(self, rhs="const double"):
          return

      @PYB11pyname("__add__")
      def __add__double(self, rhs="const double"):
          return

      def __imul__(self, rhs="const double"):
          return "Vector3d&"

      def __mul__(self, rhs="const double"):
          return "Vector3d"
    

Comparison operators
--------------------

The `comparison operators <https://docs.python.org/2/reference/datamodel.html#object.__lt__>`_ supported are ``__lt__``, ``__le__``, ``__eq__``, ``__ne__``, ``__gt__``, and ``__ge__``.  Usage of these methods (naturally all binary operators in this case) follow the same pattern as the numeric binary operators.  As an example, suppose our ``Vector3d`` class in the previous example also defined comparisons with with either ``Vector3d`` or ``double``:

.. code-block:: cpp

  class Vector3d {
  public:
    bool operator==(const Vector3d& rhs) const;
    bool operator!=(const Vector3d& rhs) const;
    bool operator< (const Vector3d& rhs) const;

    bool operator==(const double rhs) const;
    bool operator!=(const double rhs) const;
    bool operator< (const double rhs) const;
  };

We can expose these operations to Python similarly to the binary math operators::

  class Vector3d:

      def __eq__(self):
          return

      def __ne__(self):
          return

      def __lt__(self):
          return

      @PYB11pyname("__eq__")
      def __eq__double(self, rhs="const double"):
          return "bool"
          
      @PYB11pyname("__ne__")
      def __ne__double(self, rhs="const double"):
          return "bool"
          
      @PYB11pyname("__lt__")
      def __lt__double(self, rhs="const double"):
          return "bool"


Functor (call) operator
-----------------------

A special class operator in Python is the ``__call__`` operator (corresponding to the C++ ``operator()`` method), which allows a class to operate like a function.  If we have a C++ functor class, we can expose this functor behavior by binding the C++ ``operator()`` call as ``__call__``.  As an example, suppose we have C++ functor like the following:

.. code-block:: cpp

  class Transmute {
  public:
    double operator()(const double x);
  };

we can expose this functor nature of ``Transmute`` via this sort of PYB11 binding::

  class Transmute:

      def __call__(self, x="const double"):
          return "double"

PYB11Generator automatically associates ``__call__`` with the C++ method ``operator()``, unless overridden with something like ``@PYB11implementation``.

.. _class-misc-operators:

Miscellaneous operators
-----------------------

Another pair other useful operators supported are ``__repr__`` and ``__str__``.  These are used to create string representations of objects for slightly different purposes, as explained in the official `Python documentation <https://docs.python.org/2/reference/datamodel.html#object.__repr__>`_ -- essentially ``__repr__`` should return a string representation of the object such that it could be reconstructed, vs. ``__str__`` which should produce a human friendly string.

Any function or method that produces such strings is fine to bind to these names (often via renaming such as ``@PYB11pyname("__str__")``), but a very common pattern is to use lambda functions with the :func:`PYB11implementation` decorator to implement these methods directly in the binding code.  As one example, we might bind useful versions of these operators for the example C++ class ``Vector3d`` above as::

  class Vector3d:

      @PYB11implementation("[](const Vector3d& self) -> std::string { return "[" + self.x + ", " + self.y + ", " + self.z + "]" }")
      def __repr__(self):
          return "std::string"

      @PYB11implementation("[](const Vector3d& self) -> std::string { return "Vector3d(" + self.x + " " + self.y + " " + self.z + ")" }")
      def __str__(self):
          return "std::string"

Sequence methods
----------------

Probably the first thing to point out here is this section is *not* necessary for handling STL containers: pybind11 has built-in support for :ref:`pybind11:stl_bind`, which PYB11Generator provides convenient wrappers for.  In fact, so long as implicit copying of STL containers through the Python-C++ interface is acceptable, nothing need be done with STL containers at all -- they will automatically be handled by pybind11 transparently.

Binding the Python sequence methods for your own C++ types can at times be a complicated process, and there is not necessarily a single solution that fits all cases.  There are `several methods <https://docs.python.org/2/reference/datamodel.html#emulating-container-types>`_ in Python you can override to provide sequence information: ``__len__``, ``__getitem__``, ``__setitem__``, ``__getslice__``, ``__setslice__``, ``__iter__``, etc.  PYB11Generator allows all these methods to be used via pybind11, but it definitely behooves the interested user to thoroughly understand the `pybind11 <https://pybind11.readthedocs.io/en/stable/advanced/misc.html#binding-sequence-data-types-iterators-the-slicing-protocol-etc>`_ and `Python <https://docs.python.org/2/reference/datamodel.html#emulating-container-types>`_ documentation on this subject.  It will often require writing some lightweight interstitial code to translate C++ container information to Python and back, for which lambda functions and the :py:func:`PYB11implementation` decorator are handy.

As the bare beginning of an example, here is a version of one of the pybind11 test C++ sequence classes (stripped to just the interface) drawn from the ``pybind11/tests/test_sequences_and_iterators.cpp`` test code:

.. code-block:: cpp

  class Sequence {
    public:
        Sequence(size_t size);
        Sequence(const std::vector<float> &value);
        Sequence(const Sequence &s);

        bool operator==(const Sequence &s) const;
        bool operator!=(const Sequence &s) const;

        float operator[](size_t index) const;
        float &operator[](size_t index);

        bool contains(float v) const;

        Sequence reversed() const;

        size_t size() const;

        const float *begin() const;
        const float *end() const;
    };

and here is an example binding for these methods translated from the pybind11 test code in ``pybind11/tests/test_sequences_and_iterators.cpp`` to PYB11Generator Python syntax:

.. code-block:: py

  class Sequence:

     def pyinit0(self, size="size_t"):
         return
     def pyinit1(self, value="const std::vector<float>&"):
         return
     def pyinit2(self, s="const Sequence&"):
         return

     def __eq__(self):
         return
     def __ne__(self):
         return

     # Sequence methods
     @PYB11cppname("size")
     def __len__(self):
         return "size_t"

     @PYB11implementation("[](const Sequence &s, size_t i) { if (i >= s.size()) throw py::index_error(); return s[i]; }")
     def __getitem__(self, i="size_t"):
         return "float"

     @PYB11implementation("[](Sequence &s, size_t i, float v) { if (i >= s.size()) throw py::index_error(); s[i] = v; }")
     def __setitem__(self, i="size_t", v="float"):
         return "void"

     # Optional sequence methods
     @PYB11keepalive(0, 1)   # Essential: keep object alive while iterator exists
     @PYB11implementation("[](const Sequence &s) { return py::make_iterator(s.begin(), s.end()); }")
     def __iter__(self):
         return "py::iterator"

     @PYB11cppname("contains")
     @PYB11const
     def __contains__(self, v="float"):
         return "bool"

     @PYB11cppname("reversed")
     @PYB11const
     def __reversed__(self):
         return "Sequence"

     # Slicing protocol
     @PYB11pyname("__getitem__")
     @PYB11implementation("""[](const Sequence &s, py::slice slice) -> Sequence* {
                            size_t start, stop, step, slicelength;
                            if (!slice.compute(s.size(), &start, &stop, &step, &slicelength)) throw py::error_already_set();
                             Sequence *seq = new Sequence(slicelength);
                             for (size_t i = 0; i < slicelength; ++i) {
                               (*seq)[i] = s[start]; start += step;
                             }
                            return seq;
                          }""")
     def __getitem__slice(self, slice="py::slice"):
         return "Sequence*"

     @PYB11pyname("__setitem__")
     @PYB11implementation("""[](Sequence &s, py::slice slice, const Sequence &value) {
                            size_t start, stop, step, slicelength;
                            if (!slice.compute(s.size(), &start, &stop, &step, &slicelength))
                                throw py::error_already_set();
                            if (slicelength != value.size())
                                throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
                            for (size_t i = 0; i < slicelength; ++i) {
                                s[start] = value[i]; start += step;
                            }"""
                          }""")
     def __getitem__slice(self, slice="py::slice", value="const Sequence&"):
         return "void"

This rather in-depth example uses a few concepts not introduced yet (such as ``@PYB11keepalive``) which are discussed later, but hopefully gives a flavor of what is needed.  Mapping types are also supported through the same sort of overriding of built-in Python methods analogous to above.

.. _template_methods:

-----------------
Templated methods
-----------------

Templated methods are handled in a very similar manner to :ref:`function-templates`.  Suppose we want to bind the templated method in the following C++ class:

.. code-block:: cpp

  class A {
  public:

    template<typename ValueA, typename ValueB, typename ValueC>
    ValueC
    transmogrify(const ValueA& x, const ValueB& y);

  };

In order to bind this method we first create a python class method and decorate it with ``@PYB11template`` and the template types as strings.  We then create however many instantiations of this method as we like using :func:`PYB11TemplateMethod`::

  class A:

      @PYB11template("ValueA", "ValueB", "ValueC")
      def transmogrify(self, x="%(ValueA)s", y="%(ValueB)s"):
          "I'm sure this does something useful..."
          return "%(ValueC)s"

      transmogrifyIntIntDouble = PYB11TemplateMethod(transmogrify, ("int", "int", "double"),             pyname="transmogrify")
      transmogrifyI32I32I64    = PYB11TemplateMethod(transmogrify, ("uint32_t", "uint32_t", "uint64_t"), pyname="transmogrify")

Comparing this with the example in :ref:`function-templates`, we see that handling template class methods is nearly identical to template functions.  The only real difference is we instantiate the template class method using ``PYB11TemplateMethod`` (assigned to class attributes) instead of ``PYB11TemplateFunction``.

.. _class-attributes:

----------
Attributes
----------

C++ structs and classes can have attributes, such as:

.. code-block:: cpp

  struct A {
    double x;                // An ordinary attribute
    const double y;          // A readonly attribute
    static double xstatic;   // A static attribute
  };

Attributes in pybind11 are discussed in :ref:`pybind11:properties`; PYB11Generator exposes these kinds of attributes via the special PYB11 types ``PYB1readwrite`` and ``PYB11readonly``.  We can expose the attributes of ``A`` in this case via PYB11Generator using::

  class A:
      x = PYB11readwrite(doc="An ordinary attribute")
      y = PYB11readonly(doc="A readonly attribute")
      xstatic = PYB11readwrite(static=True, doc="A static attribute")

In this example we have used the optional arguments ``doc`` to add document strings to our attributes, and ``static`` to indicate a static attribute -- for the full set of options to these functions see :func:`PYB11readwrite` and :func:`PYB11readonly`.

.. _class-properties:

----------
Properties
----------

A related concept to attributes is class properties, where we use getter and setter methods for data of classes as though they were attributes.  Consider the following C++ class definition:

.. code-block:: cpp

  class A {
    public:
    double getx() const;     // Getter for a double named "x"
    void setx(double val);   // Setter for a double named "x"
  };

There are at least two ways we can go about creating ``A.x`` as a property.

Option 1: use ``PYB11property``
-------------------------------

The most convenient method (or at least most succinct) to treat ``A.x`` as a property is via the ``PYB11property`` helper type.  In this example we could simply write::

  class A:
      x = PYB11property(getter="getx", settter="setx",
                        doc="Some helpful description of x for this class")

This minimal example demonstrates that using ``PYB11property`` we can expose properties in a single line like this -- see full description of :func:`PYB11property`.

Option 2: use an ordinary python property definition
----------------------------------------------------

Python has native support for properties via the built-in :py:func:`property`; PYB11Generator is able to interpret use of this function to define pybind11 properties as well.  We can use this method to create ``A.x`` as follows::

  class A:

      def getx(self):
          return

      def setx(self):
          return

      x = property(getx, setx, doc="Some helpful description of x for this class")

This method has the advantage we are using all ordinary python constructs, which PYB11Generator is able to parse and create the property as desired.

.. Note::

   In this second example we have also exposed the ``getx`` and ``setx`` methods to be bound in pybind11.  If this is not desired, we can decorate these methods with ``@PYB11ignore``, allowing these methods to be used in the :py:func:`property` definition while preventing them from being directly exposed themselves.

.. _dynamic-attributes:

------------------
Dynamic attributes
------------------

By default pybind11 classes are immutable from Python, so it is an error to try and insert new attributes to an instance of a pybind11 bound C++ class.  This is different than default behavior in Python however, which allows instances of classes to be modified with new attributes.  For example, the following is legal Python code:

.. code-block:: pycon

   >>> class Strongbadia:
   ...     headOfState = "Strong Bad"
   ...
   >>> country = Strongbadia()
   >>> country.population = "Tire"    # Valid, we just dynamically added a new attribute

pybind11 allows us to specify if we want classes to be modifiable in this way (see `pybind11 docs <https://pybind11.readthedocs.io/en/stable/classes.html#dynamic-attributes>`_), which is reflected in PYB11Generator by using the decorator ``@PYB11dynamic_attr``.  So if we wanted to modify one of our class bindings for ``A`` above to allow dynamic attributes, we can simply decorate the class declaration like::

  @PYB11dynamic_attr
  class A:
  ...

.. _class-singletons:

----------
Singletons
----------

Suppose we have declared a C++ class to be a singleton object (i.e., declared all constructors and destructors private) like so:

.. code-block:: cpp

  class Asingleton {
  public:
    static A* instance() { return instanceptr; }

  private:
    static A* instanceptr;
    A();
    A(const A&);
    A& operator=(const A&);
    ~A();
  };

pybind11 (via its use of ``std::unique_ptr`` to hold Python instances) assumes bound objects are destructible, but for singletons such as ``Asingleton`` above the destructor is private.  We must notify pybind11 that singletons such as this are different (as discussed in pybind11 for :ref:`pybind11:classes_with_non_public_destructors`) -- PYB11Generator accomplishes this via the decorator ``@PYB11singleton`` like so::

  @PYB11singleton
  class Asingleton:

      @PYB11static
      @PYB11returnpolicy("reference")
      def instance(self):
          return "Asingleton*"

This example also involves setting a policy for handling the memory of the ``Asingleton*`` returned by ``A.instance``: these sorts of memory management details are discussed in :ref:`return-policies`.

.. _class-templates:

-----------------
Templated classes
-----------------

PYB11 handles C++ class templates similarly to :ref:`function-templates`: first, we decorate a class definition with ``@PYB11template``, which takes an arbitrary number of string arguments representing the template parameters; second, we use the :func:`PYB11TemplateClass` function to create instantiations of the template class.  Consider a C++ template class definition:

.. code-block:: cpp

  template<typename Scalar>
  class Vector {
  public:
    Scalar x, y, z;                        // Coordinate attributes

    Vector(Scalar x, Scalar y, Scalar z);  // Constructor
    Scalar magnitude() const;              // Compute the magnitude (norm)
  };

We can create PYB11Generator instantiations of this class for ``double`` and ``float`` types using::

  @PYB11template("Scalar")
  class Vector:
      "A simple three-dimensional Vector type using %(Scalar)s coordinates"

      x = PYB11readwrite()
      y = PYB11readwrite()
      z = PYB11readwrite()

      def pyinit(self, x="%(Scalar)s", y="%(Scalar)s", z="%(Scalar)s"):
          "Construct with specified coordinates"

      @PYB11const
      def magnitude(self):
          "Compute the magnitude (norm)"
          return "%(Scalar)s"

  FloatVector = PYB11TemplateClass(Vector, template_parameters="float")
  DoubleVector = PYB11TemplateClass(Vector, template_parameters="double")

Just as is the case with template functions, classes decorated with ``@PYB11template`` are implicitly ignored by PYB11Generator until an instantiation is created with :func:`PYB11TemplateClass`.  Additionally, template parameters specified in ``@PYB11template`` become named patterns which can be substituted with the types used to instantiate the templates.  So, in the ``Vector`` example above, ``%(Scalar)s`` becomes ``float`` in the first instantiation and ``double`` in the second.  See :func:`PYB11template` and :func:`PYB11TemplateClass` for further details.
