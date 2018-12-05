.. _memory-policies:

=================
Memory management
=================

It is often possible to blindly bind C++ and Python via pybind11 without worrying about how the memory/lifetime of the objects is handled.  However, since C++ allows memory and objects to be allocated/deallocated in a variety of ways, at times it is necessary to pay attention to this issue.  In this section we discuss the pybind11 methods of handling memory and object lifetimes PYB11Generator provides wrappers for.  In order to understand this section it is advisable to read the `pybind11 documentation on the use of smart pointers <https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html#smart-pointers>`_, pybind11 :ref:`pybind11:return_value_policies`, and :ref:`pybind11:call_policies`.

Class holder types
==================

When pybind11 creates a new instance of a bound C++ class, it uses a smart pointer type to hold and manage that instance.  The default type used by pybind11 for this purpose is ``std::unique_ptr``, which means new objects created in this manner by Python will be deallocated when their reference count goes to zero.  In most circumstances this is fine, but some C++ applications may have a smart pointer type they already are working with.  In such cases we can take control of the smart pointer type used to allocate and manage classes from Python.  In PYB11 we specify this by decorating class declarations with ``@PYB11holder``.  For instance, to make pybind11 use ``std::shared_ptr`` to hold a class type ``A``::

  @PYB11holder("std::shared_ptr")
  class A:

      def pyinit(self):
          "Default constructor"

      def some_func(self, x="int"):
          return "int"

pybind11 supports ``std::unique_ptr`` and ``std::shared_ptr`` without further work.  It is possible to use any reference counting smart pointer for this purpose, but other types require more information specified to pybind11.  See the `pybind11 documentation for further information <https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html#custom-smart-pointers>`_.

.. Note::

   Overriding the holder smart pointer type can result in subtleties that lead to hard to understand memory errors.  If using this capability, read the `pybind11 description <https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html#std-shared-ptr>`_ carefully!

