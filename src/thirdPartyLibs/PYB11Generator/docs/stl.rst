.. _stl:

==============
STL containers
==============

In many cases STL containers (such as ``vector``, ``deque``, ``map``, etc.) can be used transparently with pybind11: python lists will automatically convert to ``std::vector`` and vice versa for example with no extra work or notation needed.  The main caveat to this convenience, however, is that this is accomplished by copying the data in the container, which has two potential drawbacks: for large containers this may not be practical; second, any attempt to change data on either the C++ or Python side will be lost due to the fact you would be operating on a copy.  Even if your function/method specification is passing STL containers by reference, this silent copying will make modifying them across the Python/C++ interface impossible without further work.

If you need to pass an STL container without all this magic copying, it becomes necessary to directly bind such containers and define the behavior you want.  PYB11Generator provides some interfaces to pybind11's machinery for such bindings, but it is essential to first read and understand what `pybind11 is doing for STL types <https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html#stl-containers>`_.

PYB11Generator currently only provides Python convenience methods for handling two STL containers: ``std::vector`` via :func:`PYB11_bind_vector`, and ``std::map`` with :func:`PYB11_bind_map`.  It is possible to use the pybind11 semantics directly in C++ in combination with PYB11Generator to handle cases beyond ``std::vector`` and ``std::map`` of course, it simply involves using the `pybind11 C++ interface <https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html#stl-containers>`_ directly.

std::vector
===========

Suppose we want to bind ``std::vector<int>`` and ``std::vector<Aclass>`` in our module such that they will be modifiable through the interface (no copying).  We can accomplish this by adding two lines to our Python module definition::

  vector_of_int = PYB11_bind_vector("int", opaque=True)
  vector_of_Aclass = PYB11_bind_vector("Aclass", opaque=True)

When we import the resulting compiled module it will now include the types ``vector_of_int`` and ``vector_of_Aclass`` explicitly, and we need to deal in those types rather than the more convenient Python lists for arguments of those types.  The ``opaque`` argument here is what makes pybind11 treat these vector's as mutable references through the Python/C++ interface.  We also have the option of making these bindings local to the module or global: see :func:`PYB11_bind_vector` for the full description.

std::map
========

Making ``std::map`` instances mutable through the Python/C++ interface (opaque as described in pybind11 terms) is similar to our treatment of ``std::vector``.  If in our module we need to use ``std::map<std::string, Aclass>`` as an opaque argument we simply add a line::

  map_of_int_to_Aclass = PYB11_bind_map("int", "Aclass", opaque=True)

Just as in our prior ``std::vector`` examples, our module will now include a type ``map_of_int_to_Aclass`` which we can use explicitly to pass this container between Python and C++ mutably.
