.. _PYB11-functions:

PYB11 special functions and classes
===================================

This section describes the special functions and classes defined in PYB11Generator for use in createing python bindings.  Note we use the convention that PYB11 internals always start with the ``PYB11`` prefix.

.. #############################################################################
.. py:function:: PYB11generateModule(pymodule[, basename=None])

  Inspect the function and class definitions in ``pymodule``, and write a C++ file containing pybind11 statements to bind those interfaces.

  * ``pymodule``: the module to be introspected for the interface

  * ``"basename"``: a basename for the generated C++ file.  If specified, the output is written to ``basename.cc``, otherwise output will be written to ``pymodule.cc``

.. #############################################################################
.. py:function:: PYB11TemplateFunction(func_template, template_parameters[, cppname = None, pyname = None, docext = ""])

  Instantiate a function template (``func_template``) that was decorated by ``@PYB11template``.

  * ``func_template``: The template function definition

  * ``template_parameters``: A single string (for a single template parameter function) or tuple of strings (for multiple template parameters), one for each template parameter defined by ``@PYB11template`` on ``func_template``.

  * ``cppname``: The name of the C++ function template, if different from that used for ``func_template``.

  * ``pyname``: The name of the resulting Python function; defaults to the name of the instance created for this invocation of ``PYB11TemplateFunction``.

  * ``docext``: An optional string extension to be applied to the docstring associated with ``func_template``.

.. #############################################################################
.. py:function:: PYB11attr([value=None, pyname=None])

  Create an attribute in a module; corresponds to the pybind11 command ``attr``.

  * ``value``: define the C++ name this variable corresponds to.  If ``None``, defaults to the name of the local python variable.

  * ``pyname``: define the generated python attribte name.  If ``None``, defaults to the name of the local python variable.

.. #############################################################################
.. py:function:: PYB11readwrite([static=False, pyname=None, cppname=None, returnpolicy=None, doc=None])

  Define a readwrite class attribute; corresponds to pybind11 ``def_readwrite``.

  * ``static``: If ``True``, specifies the bound attribute is static.

  * ``pyname``: Optionally specify the Python name of the attribute.  If ``None``, assumes the Python name is the name of Python variable instance.

  * ``cppname``: Optionally specify the C++ name of the attribute.  If ``None``, assumes the C++ name is the name of Python variable instance.

 * ``returnpolicy``: Specify a special return policy for how to handle the memory of the return value.  Read pybind11 documentation at :ref:`pybind11:return_value_policies`.

  * ``doc``: Optionally give a docstring.

.. #############################################################################
.. py:function:: PYB11readonly([static=False, pyname=None, cppname=None, returnpolicy=None, doc=None])

  Define a readonly class attribute; corresponds to pybind11 ``def_readonly``.

  * ``static``: If ``True``, specifies the bound attribute is static.

  * ``pyname``: Optionally specify the Python name of the attribute.  If ``None``, assumes the Python name is the name of Python variable instance.

  * ``cppname``: Optionally specify the C++ name of the attribute.  If ``None``, assumes the C++ name is the name of Python variable instance.

  * ``returnpolicy``: Specify a special return policy for how to handle the memory of the return value.  Read pybind11 documentation at :ref:`pybind11:return_value_policies`.

  * ``doc``: Optionally give a docstring.

.. #############################################################################
.. py:function:: PYB11property([returnType = None, getter = None, setter = None, doc = None, getterraw = None, setterraw = None,  getterconst = True, setterconst = False, static = None, returnpolicy = None])
                 
   Helper to setup a class property.

   * ``returnType``: Specify the C++ type of the property

   * ``getter``: A string with the name of the getter method.  If ``None``, assumes the getter C++ specification looks like ``returnType (klass::*)() const``.

   * ``setter``: A string with the name of the setter method.  If ``None``, assumes the setter C++ specification looks like ``void (klass::*)(returnType& val)``.

   * ``doc``: Specify a document string for the property.

   * ``getterraw``: Optionally specify raw coding for the getter method.  Generally this is used to insert a C++ lambda function.  Only one of ``getter`` or ``getterraw`` may be specified.

   * ``setterraw``: Optionally specify raw coding for the setter method.  Generally this is used to insert a C++ lambda function.  Only one of ``setter`` or ``setterraw`` may be specified.

   * ``getterconst``: Specify if ``getter`` is a const method.

   * ``setterconst``: Specify if ``setter`` is a const method.

   * ``static``: If ``True``, make this a static property.

   * ``returnpolicy``: Specify a special return policy for how to handle the memory of the return value.  Read pybind11 documentation at :ref:`pybind11:return_value_policies`.

.. #############################################################################
.. py:function:: PYB11TemplateMethod(func_template, template_parameters[, cppname = None, pyname = None, docext = ""])

  Instantiate a class method (``func_template``) that was decorated by ``@PYB11template``.

  * ``func_template``: The template method definition

  * ``template_parameters``: A single string (for a single template parameter method) or tuple of strings (for multiple template parameters), one for each template parameter defined by ``@PYB11template`` on ``func_template``.

  * ``cppname``: The name of the C++ method template, if different from that used for ``func_template``.

  * ``pyname``: The name of the resulting Python method; defaults to the name of the instance created for this invocation of ``PYB11TemplateMethod``.

  * ``docext``: An optional string extension to be applied to the docstring associated with ``func_template``.

.. #############################################################################
.. py:function:: PYB11TemplateClass(klass_template, template_parameters[, cppname = None, pyname = None, docext = ""])

  Instantiate a class template (``klass_template``) that was decorated by ``@PYB11template``.

  * ``klass_template``: The template class definition

  * ``template_parameters``: A single string (for a single template parameter class) or tuple of strings (for multiple template parameters), one for each template parameter defined by ``@PYB11template`` on ``klass_template``.

  * ``cppname``: The name of the C++ class template, if different from that used for ``klass_template``.

  * ``pyname``: The name of the resulting Python class; defaults to the name of the instance created for this invocation of ``PYB11TemplateClass``.

  * ``docext``: An optional string extension to be applied to the docstring associated with ``klass_template``.

.. #############################################################################
.. py:function:: PYB11enum(values[, name=None, namespace="", cppname=None, export_values=False, doc=None])

   Declare a C++ enum for wrapping in pybind11 -- see `pybind11 docs <https://pybind11.readthedocs.io/en/stable/classes.html#enumerations-and-internal-types>`_.

   * ``values``: a tuple of strings listing the possible values for the enum

   * name: set the name of enum type in Python.  ``None`` defaults to the name of the instance given this enum declaration instance.

   * namespace: an optional C++ namespace the enum lives in.

   * cppname: the C++ name of the enum.  ``None`` defaults to the same as ``name``.

   * export_values: if ``True``, causes the enum values to be exported into the enclosing scope (like an old-style C enum).

   * doc: an optional document string.

.. #############################################################################
.. py:function:: PYB11_bind_vector(element[, opaque=False, local=None])

   Bind an STL::vector explicitly.  This is essentially a thin wrapper around the pybind11 ``py::bind_vector`` function (see :ref:`pybind11:stl_bind`).

   * ``element``: the C++ element type of the ``std::vector``

   * ``opaque``: if ``True``, causes the bound STL vector to be "opaque", so elements can be changed in place rather than accessed as copies.  See :ref:`pybind11:stl_bind`.

   * ``local``: determines whether the binding of the STL vector should be module local or not; once again, see :ref:`pybind11:stl_bind`.

.. #############################################################################
.. py:function:: PYB11_bind_map(key, value[, opaque=False, local=None])

   Bind an STL::map explicitly.  This is a thin wrapper around the pybind11 ``py::bind_map`` function (see :ref:`pybind11:stl_bind`).

   * ``key``: the C++ key type

   * ``value``: the C++ value type

   * ``opaque``: if ``True``, causes the bound STL map to be "opaque", so elements can be changed in place rather than accessed as copies.  See :ref:`pybind11:stl_bind`.

   * ``local``: determines whether the binding of the STL map should be module local or not; once again, see :ref:`pybind11:stl_bind`.

.. #############################################################################
.. py:function:: PYB11_inject(fromcls, tocls[, virtual=None, pure_virtual=None])

   Convenience method to inject methods from class ``fromcls`` into ``tocls``.  This is intended as a utility to help avoiding writing redundant methods common to many classes over and over again.  Instead a convenience class can be defined containing the shared methods (typically screened from generation by ``@PYB11ignore``), and then ``PYB11_inject`` is used to copy those methods into the target classes.

   * ``fromcls``: Python class with methods we want to copy from.

   * ``tocls``: Python class we're copying methods to.

   * ``virtual``: if ``True``, force all methods we're copying to be treated as virtual.

   * ``pure_virtual``: if ``True``, force all methods we're copying to be treated as pure virtual.
