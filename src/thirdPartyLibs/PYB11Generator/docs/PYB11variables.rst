.. _variables:

========================
PYB11 reserved variables
========================

For the most part Python variables declared in module bindings are ignored by PYB11Generator.  There are a few exceptions to this rule though -- some variables are used to communicate information to PYB11Generator as described below.

PYB11includes = [...]
  A list of strings, each of which represents a file that should be ``#include``-ed in the final C++ generated file.  For instance, if we needed the C++ to have the following include statements:

  .. code-block:: cpp

    #include "A.hh"
    #include <vector>

  We would put this in our file for PYB11Generator::

    PYB11includes = ['"A.hh"', '<vector>']

PYB11namespaces = [...]
  A list of strings for C++ namespaces we wish to use in the generated C++ -- as an example, the following statement::

    PYB11namespaces = ["extreme", "measures"]

  results in the following in the generated pybind11 C++ source:

  .. code-block:: cpp

    using namespace extreme;
    using namespace measures;

PYB11scopenames = [...]
  A list of C++ types we want to directly declare with ``use`` (more focused than importing an entire namespace).  In order to decare we want to use ``std::vector`` and ``MyNameSpace::A``, we could use::

    PYB11scopenames = ["std::vector", "MyNameSpace::A"]

  which simply inserts the following in the generated C++:

  .. code-block:: cpp

    using std::vector
    using MyNameSpace::A

PYB11preamble = "..."
  PYB11preamble is used to specify a string of arbitrary C++ code that will be inserted near the top of the generated pybind11 source file.  ``PYB11preamble`` is a bit of catch-all, we could for instance directly perform the tasks of ``PYB11includes`` and ``PYB11namespace`` using ``PYB11preamble`` by simply typing the final C++ code in here.  One typical usage of this preamble variable is to insert small inline utility methods directly in the final C++ code.  For instance, if we had need of a simple function we want to use in the subsequent bindings, we could do something like::
  
    PYB11preamble = """
    namespace JustForBindings {
      inline int square(int x) { return x*x; }
    }"""

  and now our generated code will include this function.

PYB11opaque = [...]
  A list of C++ types we want to be treated as opaque: typically STL types declared as opaque and global in another module.  See :ref:`stl` for further information.  As an example, to declare that ``std::vector<int>`` and ``std::vector<std::string>`` are declared as opaque in another module::

    PYB11opaque = ["std::vector<int>", "std::vector<std::string>"]

  which results in the following inserted into the generated C++:

  .. code-block:: cpp

    PYBIND11_MAKE_OPAQUE(std::vector<int>)
    PYBIND11_MAKE_OPAQUE(std::vector<std::string>)

Note all of these reserved variables affect the start of the generated pybind11 C++ code, coming before any of the function, class, module, or other pybind11 declarations that are subsequently generated.  The order that these methods are executed in is the same as they are listed above: first any ``PYB11includes``, then ``PYB11namespaces``, ``PYB11scopenames``, ``PYB11preamble``, and finally ``PYB11opaque``.  If we were to include all of the above examples (in any order) in a single source code for instance like so::

  PYB11includes = ['"A.hh"', '<vector>']
  PYB11namespaces = ["extreme", "measures"]
  PYB11scopenames = ["std::vector", "MyNameSpace::A"]
  PYB11preamble = """
  namespace JustForBindings {
    inline int square(int x) { return x*x; }
  }
  """
  PYB11opaque = ["std::vector<int>", "std::vector<std::string>"]

the generated pybind11 code would look like:

.. code-block:: cpp
  
  ...
  #include "A.hh"
  #include <vector>

  using namespace extreme;
  using namespace measures;

  using std::vector
  using MyNameSpace::A


  namespace JustForBindings {
    inline int square(int x) { return x*x; }
  }


  PYBIND11_MAKE_OPAQUE(std::vector<int>)
  PYBIND11_MAKE_OPAQUE(std::vector<std::string>)

  //------------------------------------------------------------------------------
  // Make the module
  //------------------------------------------------------------------------------
  ...

.. _class-variables:

---------------------------
PYB11 variables for classes
---------------------------

There is also one reserved variable for class scope in PYB11:

PYB11typedefs = "...."
  A string of C++ to be inserted at the beginning of a class declaration.  This is a bit of a misnomer; the string can be any valid C++ for use in the class declaration scope.  Suppose for instance we are defining a class ``A``, and we want to declare some typedefs only for use in the scope of ``A``::

    class A:

        PYB11typedefs = """
    typedef int    IntType;
    typedef double ScalarType;
    """

  which results in the following C++ code:

  .. code-block:: cpp

    //............................................................................
    // Class A
    {

      typedef int    IntType;
      typedef double ScalarType;

      class_<A> obj(m, "A");
      ...
    }

  so now ``IntType`` and ``ScalarType`` are avaiable as types in the scope where we are defining A.
