.. _PYB11-functions:

PYB11 special functions and classes
===================================

This section describes the special functions and classes defined in PYB11Generator for use in createing python bindings.  Note we use the convention that PYB11 internals always start with the ``PYB11`` prefix.

.. #############################################################################
.. py:function:: PYB11generateModule(pymodule[, basename=None])

  Inspect the function and class definitions in ``pymodule``, and write a C++ file containing pybind11 statements to bind those interfaces.

  * ``pymodule``: the module to be introspected for the interface

  * ``"basename"``: a basename for the generated C++ file.  If specified, the output is written to ``basename.cc``, otherwise output will be written to ``mymodule.cc``

.. #############################################################################
.. py:function:: PYB11attr([value=None, pyname=None])

  Create an attribute in a module; corresponds to the pybind11 command ``attr``.

  * ``value``: define the C++ name this variable corresponds to.  If ``None``, defaults to the name of the local python variable.

  * ``pyname``: define the generated python attribte name.  If ``None``, defaults to the name of the local python variable.

.. #############################################################################
.. py:function:: PYB11readwrite([static=False, pyname=False, cppname=False, doc=None])

  Define a readwrite class attribute; corresponds to pybind11 ``def_readwrite``.

  * ``static``: If ``True``, specifies the bound attribute is static.

  * ``pyname``: Optionally specify the Python name of the attribute.  If ``None``, assumes the Python name is the name of Python variable instance.

  * ``cppname``: Optionally specify the C++ name of the attribute.  If ``None``, assumes the C++ name is the name of Python variable instance.

  * ``doc``: Optionally give a docstring.

.. #############################################################################
.. py:function:: PYB11readonly([static=False, pyname=False, cppname=False, doc=None])

  Define a readonly class attribute; corresponds to pybind11 ``def_readonly``.

  * ``static``: If ``True``, specifies the bound attribute is static.

  * ``pyname``: Optionally specify the Python name of the attribute.  If ``None``, assumes the Python name is the name of Python variable instance.

  * ``cppname``: Optionally specify the C++ name of the attribute.  If ``None``, assumes the C++ name is the name of Python variable instance.

  * ``doc``: Optionally give a docstring.

.. #############################################################################
.. py:function:: PYB11TemplateClass(klass_template, template_parameters[, cppname = None, pyname = None, docext = ""])

  Instantiate a class template (``klass_template``) that was decorated by ``@PYB11template``.

  * ``klass_template``: The template class definition

  * ``template_parameters``: A single string (for a single template parameter class) or tuple of strings (for multiple template parameters), one for each template parameter defined by ``@PYB11template`` on ``klass_template``.

  * ``cppname``: The name of the C++ class template, if different from that used for ``klass_template``.

  * ``pyname``: The name of the resulting Python class; defaults to the name of the instance created for this invocation of ``PYB11TemplateClass``.

  * ``docext``: An optional string extension to be applied to the docstring associated with ``klass_template``.

.. #############################################################################
.. py:function:: PYB11_bind_vector(element[, opaque=False, local=None])

   Bind an STL::vector explicitly.  This is essentially a thin wrapper around the pybind11 ``py::bind_vector`` function (see :ref:`pybind11::stl_bind`).

   * ``element``: the C++ element type of the ``std::vector``

   * ``opaque``: if ``True``, causes the bound STL vector to be "opaque", so elements can be changed in place rather than accessed as copies.  See :ref:`pybind11:stl_bind`.

   * ``local``: determines whether the binding of the STL vector should be module local or not; once again, see :ref:`pybind11:stl_bind`.

