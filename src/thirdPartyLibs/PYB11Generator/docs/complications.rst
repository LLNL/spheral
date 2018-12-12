.. _complications:

==============================
Complications and corner cases
==============================

.. _cross-module-inheritance:

------------------------
Cross-module inheritance
------------------------

For the most part having C++ types exposed in different modules is transparent: so long as you import all the necessary modules once they are all compiled and bound, everything just works.  However, the exception to this rule is if you want to bind a class in one module that inherits from a class bound in another.   Suppose for instance we have two C++ classes (``A`` and ``B``), defined in two different headers ``A.hh`` and ``B.hh``, as follows.

A.hh:

.. code-block:: cpp

  struct A {
    A()                          { printf("A::A()\n"); }
    ~A()                          { printf("A::~A()\n"); }
    virtual int func(const int x) { printf("A::func(%d)\n", x); return x; }
  };

B.hh:

.. code-block:: cpp

  #include "A.hh"

  struct B: public A {
    B(): A()                               { printf("B::B()\n"); }
    ~B()                                   { printf("B::~B()\n"); }
    virtual int func(const int x) override { printf("B::func(%d)\n", x); return x*x; }
  };


We want to expose these two class in two different modules, ``Amodule`` and ``Bmodule``.  We will now need to annotate the bindings for the ``A`` class with one piece of new information -- the module it will be bound in.  This is accomplished with a new decorator: ``PYB11module``, and the bindings for ``Amodule`` might look like::

  from PYB11Generator import *

  PYB11includes = ['"A.hh"']

  @PYB11module("Amodule")       # <--- This is our new annotation
  class A:

      def pyinit(self):
          "Default constructor"

      @PYB11virtual
      def func(self, x="int"):
          "A::func"
          return "int"

Let's suppose the above binding source is stored in file ``Amodule_bindings.py``.  We can now write our binding source for ``Bmodule`` as normal, but we need to import ``Amodule_bindings`` so we can express the inheritance relation between ``B`` and ``A``::

  from PYB11Generator import *

  import Amodule_bindings

  PYB11includes = ['"B.hh"']

  class B(Amodule_bindings.A):

      def pyinit(self):
          "Default constructor"

      @PYB11virtual
      def func(self, x="int"):
          "B::func"
          return "int"

The ``@PYB11module`` decoration on ``A`` tells PYB11Generator how to generate the pybind11 code to correctly import ``A`` rather than generate ``A`` locally, as described in the `pybind11 documentation <https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules>`_.

.. Note::

   It is critical here in the bindings for ``Bmodule`` that we use ``import Amodule_bindings``, and do *not* import ``A`` into the local scope using ``from Amodule_bindings import A``!  If we put ``A`` in the top-level scope of our bindings for ``B``, the binding code for ``A`` will be generated redundantly in the new bindings, and cause a conflict when we try to import the two modules together.

.. _non-template-to-template-inheritance:

-----------------------------------------------------
Non-templated class inheriting from a templated class
-----------------------------------------------------

PYB11Generator needs to know template parameters for templated classes in order to create concrete instantiations, but since Python does not have the concept of templates we have adopted a two-stage process for creating template class instantiations in PYB11 as described in :ref:`class-templates`.  However, if we have a non-templated class which inherits from a templated base, there is no longer the second-stage of this procedure using :func:`PYB11TemplateClass` to instantiate the base with the proper template parameters.

It is possible to handle this situation, but it requires two decorations be applied to the non-templated descendant:

#. Because the descendant will inherit the template decoration of the base class, we must explicitly state that the descendant has no template parameters with ``@PYB11template()``.

#. We still need to specify what template parameters should be used for the base class.  Template parameters in PYB11Generator are specified using python dictionary matching, so we can directly insert the proper template parameter choices in the appropriate dictionary for our non-templated descendant using ``@PYB11template_dict``.

These two steps are best demonstrated by an example -- consider the following C++ class hierarchy:

.. code-block:: cpp

  template<typename Value1, typename Value2>
  class A {
  public:
    A();
    virtual ~A();
    virtual std::string func(const Value1& x, const Value2& y) const;
  };

  class B: public A<double, int> {
  public:
    B();
    virtual ~B();
    virtual std::string func(const double& x, const int& y) const;
  };

PYB11Generator can represent this hierarchy with::

  @PYB11template("Value1", "Value2")
  class A:

      def pyinit(self):
          "Default A()"

      @PYB11virtual
      @PYB11const
      def func(self, x="const %(Value1)s&", y="const %(Value2)s&"):
          "Default A::func"
          return "std::string"

  @PYB11template()                                             # <--- force not to inherit template parameters from A
  @PYB11template_dict({"Value1" : "double", "Value2" : "int"}) # <--- specify the template parameter substitutions
  class B(A):

      def pyinit(self):
          "Default B()"

      @PYB11virtual
      @PYB11const
      def func(self, x="const double&", y="const int&"):
          "B::func override"
          return "std::string"

  # We still need to instantiate any versions of A that we need/use.
  A_double_int = PYB11TemplateClass(A, template_parameters=("double", "int"))
