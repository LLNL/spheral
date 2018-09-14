#-------------------------------------------------------------------------------
# PYB11Generator
#-------------------------------------------------------------------------------
import inspect
import sys
from PYB11utils import *
from PYB11Decorators import *
from PYB11STLmethods import *
from PYB11class import *
from PYB11property import *

#-------------------------------------------------------------------------------
# PYB11generateModule
#-------------------------------------------------------------------------------
def PYB11generateModule(modobj):
    name = modobj.__name__
    with open(name + ".cc", "w") as f:
        ss = f.write
        PYB11generateModuleStart(modobj, ss, name)

        # Bind methods.
        PYB11generateModuleFunctions(modobj, ss)

        # Bind classes.
        PYB11generateModuleClasses(modobj, ss)

        # Bind STL types.
        PYB11generateModuleSTL(modobj, ss)

        # Closing
        ss("}\n")
        f.close()

    return

#-------------------------------------------------------------------------------
# PYB11generateModuleStart
#
# All the stuff up to the methods.
#-------------------------------------------------------------------------------
def PYB11generateModuleStart(modobj, ss, name):

    # Compiler guard.
    ss("""//------------------------------------------------------------------------------
// Module %(name)s
//------------------------------------------------------------------------------
// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

namespace py = pybind11;
using namespace pybind11::literals;

""" % {"name" : name})

    # Includes
    if hasattr(modobj, "includes"):
        for inc in modobj.includes:
            ss('#include %s\n' % inc)
        ss("\n")

    # Use  namespaces
    if hasattr(modobj, "namespaces"):
        for ns in modobj.namespaces:
            ss("using namespace " + ns + ";\n")
        ss("\n")

    # Use objects from scopes
    if hasattr(modobj, "scopenames"):
        for scopename in modobj.scopenames:
            ss("using " + scopename + "\n")
        ss("\n")

    # Preamble
    if hasattr(modobj, "preamble"):
        if modobj.preamble:
            ss(modobj.preamble + "\n")
        ss("\n")

    # Some pybind11 types need their own preamble.
    for objname, obj in PYB11STLobjs(modobj):
        obj.preamble(modobj, ss, objname)
    ss("\n")

    # Declare the module
    ss("""
//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(%(name)s, m) {

""" % {"name"     : name,
      })

    if inspect.getdoc(modobj):
        ss("  m.doc() = ")
        for i, line in enumerate(inspect.getdoc(modobj).split('\n')):
            if i > 0:
                ss("            ")
            ss('"%s"\n' % line);
        ss("  ;\n")
    ss("\n")

    return

#-------------------------------------------------------------------------------
# PYB11generateModuleFunctions
#
# Bind the methods in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleFunctions(modobj, ss):
    methods = PYB11functions(modobj)
    if methods:
        ss("  // Methods\n")
    for name, meth in methods:
        methattrs = PYB11attrs(meth)
        if not methattrs["ignore"]:

            # Arguments
            stuff = inspect.getargspec(meth)
            nargs = len(stuff.args)

            # Return type
            returnType = meth(*tuple(stuff.args))
            methattrs["returnType"] = returnType

            # Write the binding
            ss('  m.def("%(pyname)s", ' % methattrs)
            if returnType:
                assert not stuff.args is None
                assert not stuff.defaults is None
                assert len(stuff.args) == len(stuff.defaults)
                argNames = stuff.args
                argTypes, argDefaults = [], []
                for thing in stuff.defaults:
                    if isinstance(thing, tuple):
                        assert len(thing) == 2
                        argTypes.append(thing[0])
                        argDefaults.append(thing[1])
                    else:
                        argTypes.append(thing)
                        argDefaults.append(None)
                assert len(argNames) == nargs
                assert len(argTypes) == nargs
                assert len(argDefaults) == nargs
                ss("(%s (*)(" % returnType)
                for i, argType in enumerate(argTypes):
                    ss(argType)
                    if i < nargs - 1:
                        ss(", ")
                ss(")) &%(namespace)s%(cppname)s" % methattrs)
                for argType, argName, default in zip(argTypes, argNames, argDefaults):
                    ss(', "%s"_a' % argName)
                    if not default is None:
                        ss("=" + default)
            else:
                ss("&%(namespace)s%(cppname)s" % methattrs)

            # Write the doc string
            if inspect.getdoc(meth):
                doc = inspect.getdoc(meth)
                ss(',\n        "%s"' % inspect.getdoc(meth))
            ss(");\n")

#-------------------------------------------------------------------------------
# PYB11functions
#
# Get the functions to bind from a module
#-------------------------------------------------------------------------------
def PYB11functions(modobj):
    return [(name, meth) for (name, meth) in inspect.getmembers(modobj, predicate=inspect.isfunction)
            if name[:5] != "PYB11"]

