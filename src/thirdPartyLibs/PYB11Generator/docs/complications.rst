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

.. Note::

   It is critical here in our bindings for ``Bmodule`` that we use ``import Amodule_bindings``, and do *not* import ``A`` into our scope using ``from Amodule_bindings import A``!  If we put ``A`` in the top-level scope of our bindings for ``B``, the binding code for ``A`` will be generated redundantly in the new bindings, and cause a conflict when we try to import the two modules together.

The ``@PYB11module`` decoration on ``A`` tells PYB11Generator how to generate the pybind11 code to correctly import ``A`` rather than generate ``A`` locally, as described in the `pybind11 documentation <https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules>`_.

