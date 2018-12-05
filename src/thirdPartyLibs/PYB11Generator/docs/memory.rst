.. _memory-policies:

=================
Memory management
=================

It is often possible to blindly bind C++ and Python via pybind11 without worrying about how the memory/lifetime of the objects is handled.  However, since C++ allows memory and objects to be allocated/deallocated in a variety of ways, at times it is necessary to pay attention to this issue.
