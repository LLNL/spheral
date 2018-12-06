.. _variables:

========================
PYB11 reserved variables
========================

For the most part Python variables declared in module bindings are ignored by PYB11Generator.  There are a few exceptions to this rule though -- some variables are used to communicate information to PYB11Generator as described below.

PYB11preamble = "..."
  PYB11preamble is used to specify a string of arbitrary C++ code that will be inserted near the top of the generated pybind11 source file.



