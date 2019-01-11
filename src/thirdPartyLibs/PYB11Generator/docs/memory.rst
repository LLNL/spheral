.. _memory-policies:

=================
Memory management
=================

Generally memory management "just works" when binding C++ and Python with pybind11 without worrying about how the memory/lifetime of the objects is handled.  However, since C++ allows memory and objects to be allocated/deallocated in a variety of ways, at times it is necessary to pay attention to this issue.  In this section we discuss the pybind11 methods of handling memory and object lifetimes PYB11Generator provides wrappers for.  In order to understand this section it is advisable to read the `pybind11 documentation on the use of smart pointers <https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html#smart-pointers>`_, pybind11 :ref:`pybind11:return_value_policies`, and :ref:`pybind11:call_policies`.

.. _class-holder:

Class holder types
==================

When pybind11 creates a new instance of a bound C++ class, it uses a smart pointer type to hold and manage that instance.  The default type used by pybind11 for this purpose is ``std::unique_ptr``, which means new objects created in this manner by Python will be deallocated when their reference count goes to zero.  In most circumstances this is fine, but some C++ applications may have a smart pointer type they already are working with.  In such cases it might be preferable to make pybind11 manage these object using the same sort of smart pointer.  In PYB11 we specify this by decorating class declarations with ``@PYB11holder``.  For instance, to make pybind11 use ``std::shared_ptr`` to hold a class type ``A``::

  @PYB11holder("std::shared_ptr")
  class A:

      def pyinit(self):
          "Default constructor"

This tells pybind11 any new instance of ``A`` created from python should be managed by ``std::shared_ptr``.  pybind11 supports ``std::unique_ptr`` and ``std::shared_ptr`` without further work.  It is possible to use any reference counting smart pointer for this purpose, but types other than ``std::unique_ptr`` and ``std::shared_ptr`` require more information be specified to pybind11.  PYB11 does not provide any convenience methods for adding that additional information, but it can be done directly with pybind11  as described in the `pybind11 documentation <https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html#custom-smart-pointers>`_.

.. Note::

   Overriding the holder smart pointer type can result in subtleties that lead to hard to understand memory errors.  If using this capability, read the `pybind11 description <https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html#std-shared-ptr>`_ carefully!

.. _return-policies:

Return value policies
=====================

Return value policies are important but at times tricky for C++ types returned by reference or pointer.  By default pybind11 assumes Python will take ownership of returned values, implying Python can delete those objects when the Python reference count falls to zero.  If a C++ library returns a pointer to memory it expects to manage, however, the result of this conflict over who can manage (delete) that memory is an error.  For this reason pybind11 provides :ref:`pybind11:return_value_policies`, allowing the developer to explicitly state how memory returned from C++ should be handled.  Before using these policies it is *critical* to read and understand these policies from pybind11.  These return value policies are applied (for functions or methods) using the ``@PYB11returnpolicy`` decorator, with allowed values ``take_ownership``, ``copy``, ``move``, ``reference``, ``reference_internal``, ``automatic``, and ``automatic_reference``.  The default policy is ``automatic``.

Consider the example C++ function from the pybind11 documentation:

.. code-block:: cpp

    /* Function declaration */
    Data *get_data() { return _data; /* (pointer to a static data structure) */ }

If we want to tell pybind11 the C++ side will manage the memory for the returned ``Data*`` from this method, the ``reference`` return policy is appropriate.  We can express this by decorating the function binding as::

  @PYB11returnpolicy("reference")
  def get_data():
      return "Data*"

Decorating return values from class methods is identical to functions: simply use ``@PYB11returnpolicy`` to decorate the method declaration.

.. _call-policies:

Call policies
=============

While :ref:`return-policies` are specific to return types from functions or methods, call policies allow the user to tie together the lifetimes of return values and/or arguments.  This is discussed in depth in the pybind11 documentation :ref:`pybind11:call_policies`.  PYB11 provides the decorator ``@PYB11keepalive(a, b)`` for direct access to the pybind11 command ``py::keep_alive<a, b>``.  The arguments to the decorator ``a`` and ``b`` are integers indicating arguments in the call signature by position index, with the convention:

   * 0 denotes the return value of the function/method.

   * If decorating a class method, index 1 is the ``this`` (or ``self``) pointer.

To recreate the example from the pybind11 documentation, if we have a custom ``List`` class which we are binding in Python, we might want to decorate the ``append`` method like::

  class List:

      @PYB11keepalive(1, 2)
      def append(self, x):
          return "void"

This tells pybind11 to keep the the second argument (``x``, the element being appended) alive so long as the first argument (``self``, our container class) is alive.

.. _call_guards:

Call guards
===========

Another variation on wrapping functions/methods is to provide ``call_guard`` as discussed in the `pybind11 call_guard documentation <https://pybind11.readthedocs.io/en/stable/advanced/functions.html#call-guard>`_.  Call guards must be default constructable C++ types, and will be built in scope just before calling the wrapped method.  One typical usage would be to release the Python Global Interpreter Lock (GIL) before calling a wrapped method, so something like::

  @PYB11call_guard("py::gil_scoped_release")
  def my_wrapped_method():
      "Some C++ method that needs to have the GIL released."
      return
