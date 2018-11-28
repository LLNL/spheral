.. _PYB11-functions:

PYB11 special functions and classes
===================================

This section describes the special functions and classes defined in PYB11Generator for use in createing python bindings.  Note we use the convention that PYB11 internals always start with the ``PYB11`` prefix.

``PYB11generateModule(pymodule, basename=None)``
  Inspect the function and class definitions in ``pymodule``, and write a C++ file containing pybind11 statements to bind those interfaces.

  * ``pymodule``: the module to be introspected for the interface

  * ``"basename"``: a basename for the generated C++ file.  If specified, the output is written to ``basename.cc``, otherwise output will be written to ``mymodule.cc``

``PYB11atr(value=None, pyname=None)``
  Create an attribute in a module; corresponds to the pybind11 command ``attr``.

  * ``value``: define the C++ name this variable corresponds to.  If ``None``, defaults to the name of the local python variable.

  * ``pyname``: define the generated python attribte name.  If ``None``, defaults to the name of the local python variable.

``PYB11readwrite(static=False, pyname=False, cppname=False, doc=None)``
  Define a readwrite class attribute; corresponds to pybind11 ``def_readwrite``.

  * ``static``: If ``True``, specifies the bound attribute is static.

  * ``pyname``: Optionally specify the Python name of the attribute.  If ``None``, assumes the Python name is the name of Python variable instance.

  * ``cppname``: Optionally specify the C++ name of the attribute.  If ``None``, assumes the C++ name is the name of Python variable instance.

  * ``doc``: Optionally give a docstring.
