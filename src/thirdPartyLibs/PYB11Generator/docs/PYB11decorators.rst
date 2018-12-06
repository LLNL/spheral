.. _decorators:

================
PYB11 decorators
================

.. #############################################################################
.. decorator:: PYB11ignore

  Specifies that the decorated object should be ignored by PYB11Generator, i.e., not processed to produce any pybind11 binding output.

.. #############################################################################
.. decorator:: PYB11template("type1", "type2", ...)

  Indicates the object should be treated as a C++ template.  Accepts any number of strings which represent the names of the template arguments.

  The succeeding python class or function can use the specified template argument strings as patterns for substitation with python dictionary string replacement.  So if we are binding a C++ templated function::

    template<typename T>
    T manipulate(const T& val);

  The corresponding PYB11 template specfication would look like::

    @PYB11template("T")
    def manipulate(val = "const %(T)s"):
        return "%(T)s"

.. #############################################################################
.. decorator:: PYB11template_dict

  Explicitly specifies the dictionary of template args to values for use with ``@PYB11template`` types.

  *NOTE: this is a highly unusual pattern to need/use.  It is preferable to use the ordinary PYB11 template instantion methods* ``PYB11TemplateClass``, ``PYB11TemplateMethod``, *or* ``PYB11TemplateFunction``.

.. #############################################################################
.. decorator:: PYB11singleton

  Specifies that the decorated object should be treated as a C++ singleton.

.. #############################################################################
.. decorator:: PYB11holder(holder_type)

  Specify a special C++ holder for the generated type in ``pybind``, rather than the usual default ``std::unique_ptr``.  See pybind11 documentation on using `shared_ptr as a holder type <https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html#std-shared-ptr>`_.

.. #############################################################################
.. decorator:: PYB11dynamic_attr

  Make the wrapped class modifiable, i.e., allow attributes to be added dynamically to an instance of the class in python.  See pybind11 documentation about `dynamic attributes <https://pybind11.readthedocs.io/en/stable/classes.html?highlight=dynamic_attr#dynamic-attributes>`_.

.. #############################################################################
.. decorator:: PYB11namespace("val")

  Set the namespace the C++ type should be found in.

.. #############################################################################
.. decorator:: PYB11cppname("val")

  Give a value for the C++ name of the decorated function, class, or method.  Overrides the default assumption that the C++ name is the same as that given for the object in the PYB11 python binding file.

.. #############################################################################
.. decorator PYB11pyname("val")

  Give a value for the generated Python name of the decorated function, class, or method.  Overrides the default assumption that the Python name is the same as that given for the object in the PYB11 python binding file.

.. #############################################################################
.. decorator:: PYB11pycppname("val")

  Simultaneously set the Python and C++ name of the decorated function, class, or method.  Shorthand for specifying both ``@PYB11pyname`` and ``@PYB11cppname`` to the given ``"val"``.

.. #############################################################################
.. decorator:: PYB11virtual

  Mark a class method as virtual.

.. #############################################################################
.. decorator:: PYB11pure_virtual

  Mark a class method as pure virtual, making the class abstract.

.. #############################################################################
.. decorator:: PYB11protected

  Mark a class method as protected.

.. #############################################################################
.. decorator:: PYB11const

  Mark a class method as const.

.. #############################################################################
.. decorator:: PYB11static

  Mark a class method as static.

.. #############################################################################
.. decorator:: PYB11implementation("val")

  Give an implementation for the bound function or method.  This is typically used to specify lambda function implementations, or explicitly call a helper method.

.. #############################################################################
.. _returnpolicy:
.. decorator:: PYB11returnpolicy("val")

  Specify a pybind11 return policy for the return value of a function or method.  This is a tricky topic that if misused can create memory errors, but is at times absolutely necessary to get the expected behavior from the underlying C++ code and types.  Before using this method carefully read the pybind11 discussion about :ref:`pybind11:return_value_policies`.

.. #############################################################################
.. decorator:: PYB11keepalive(a, b)

  Tie the lifetime of objects in the return value/argument spec together, where the arguments (``a``, ``b``) are integers indicating the order of the arguments to tie together (0 refers to the return value).  This is another way of specifying memory policies, similar to returnpolicy_.  Carefully read the pybind11 discussion of the ``keep_alive`` directive in :ref:`pybind11:call_policies`.

.. #############################################################################
.. decorator:: PYB11module("val")

  Indicate the object should be imported from the specified python module.  This is useful for classes wrapped in one module which are needed in another, such as for inheritance.
