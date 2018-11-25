Introduction to PYB11Generator
==============================

Using PYB11Generator to write the interface for a C++ binding is intended to emulate writing that same interface in Python, so if you're familiar with Python it should be easy to get started with PYB11Generator.  As an example, if we have a header for a C++ function that looks like::
   // A really important function!
   int func();

We can simply define this method for binding with PYB11Generator by writing a python file including the spec::
  def func():
      "A really important function!"
      return "int"

