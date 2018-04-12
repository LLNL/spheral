#-------------------------------------------------------------------------------
# PYB11Generator
#-------------------------------------------------------------------------------
import inspect
import sys
from PYB11ClassDecorators import *
from PYB11FunctionDecorators import *

#-------------------------------------------------------------------------------
# generateModule
#-------------------------------------------------------------------------------
def generateModule(modobj):
    ss = sys.stdout.write
    name = modobj.__name__
    generateModuleStart(modobj, ss, name)

    # Bind methods.
    generateModuleFunctions(modobj, ss)

    # Closing
    ss("}\n")
    return

#-------------------------------------------------------------------------------
# generateModuleStart
#
# All the stuff up to the methods.
#-------------------------------------------------------------------------------
def generateModuleStart(modobj, ss, name):

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
            ss('#include "%s"\n' % inc)
        ss("\n")

    # Use  namespaces
    if hasattr(modobj, "namespaces"):
        for ns in modobj.namespaces:
            ss("using namespace " + ns + "\n")
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
# generateModuleFunctions
#
# Bind the methods in the module
#-------------------------------------------------------------------------------
def generateModuleFunctions(modobj, ss):
    methods = [(name, meth) for (name, meth) in inspect.getmembers(modobj, predicate=inspect.isfunction)
               if name[:2] != "__"]
    if methods:
        ss("  // Methods\n")
    for name, meth in methods:

        # Get the return type and arguments.
        returnType = meth()
        stuff = inspect.getargspec(meth)
        args = None
        if "args" in stuff.args:
            args = stuff.defaults[stuff.args.index("args") - 1]
            nargs = len(args)

        # Because python does not have function overloading, we provide the ability
        # to rename the c++ method.
        if "name" in stuff.args:
            name = stuff.defaults[stuff.args.index("name") - 1]

        # Write the binding
        dvals = {"name" : name, "returnType" : returnType}
        ss('  m.def("%s", ' % name)
        if returnType:
            assert not args is None
            ss("(%s (*)(" % returnType)
            for i, (argType, argName, default) in enumerate(__parseArgs(args)):
                ss(argType)
                if i < nargs - 1:
                    ss(", ")
            ss(")) &%s" % name)
            for argType, argName, default in __parseArgs(args):
                ss(', "%s"_a' % argName)
                if default:
                    ss("=" + default)
        else:
            ss("&%s" % name)

        # Write the doc string
        if inspect.getdoc(meth):
            doc = inspect.getdoc(meth)
            ss(',\n        "%s"' % inspect.getdoc(meth))
        ss(");\n")

#-------------------------------------------------------------------------------
# generateModuleClasses
#
# Bind the classes in the module
#-------------------------------------------------------------------------------
def generateModuleClasses(modobj, ss):
    classes = [(name, klas) for (name, klas) in inspect.getmembers(modobj, predicate=inspect.isclass)
               if name[:2] != "__"]
    for name, klass in classes:
        ss("  //............................................................................\n")
        ss("  // Class %s\n" % klass.__name__)

        # # Get the return type and arguments.
        # returnType = meth()
        # stuff = inspect.getargspec(meth)
        # args = None
        # if "args" in stuff.args:
        #     args = stuff.defaults[stuff.args.index("args") - 1]
        #     nargs = len(args)

        # # Because python does not have function overloading, we provide the ability
        # # to rename the c++ method.
        # if "name" in stuff.args:
        #     name = stuff.defaults[stuff.args.index("name") - 1]

        # # Write the binding
        # dvals = {"name" : name, "returnType" : returnType}
        # ss('  m.def("%s", ' % name)
        # if returnType:
        #     assert not args is None
        #     ss("(%s (*)(" % returnType)
        #     for i, (argType, argName, default) in enumerate(__parseArgs(args)):
        #         ss(argType)
        #         if i < nargs - 1:
        #             ss(", ")
        #     ss(")) &%s" % name)
        #     for argType, argName, default in __parseArgs(args):
        #         ss(', "%s"_a' % argName)
        #         if default:
        #             ss("=" + default)
        # else:
        #     ss("&%s" % name)

        # # Write the doc string
        # if inspect.getdoc(meth):
        #     doc = inspect.getdoc(meth)
        #     ss(',\n        "%s"' % inspect.getdoc(meth))
        # ss(");\n")

#-------------------------------------------------------------------------------
# __parseArgs
#
# Return (argType, argName, default_value (optional)
#-------------------------------------------------------------------------------
def __parseArgs(args):
    result = []
    for tup in args:
        if len(tup) == 2:
            result.append((tup[0], tup[1], None))
        else:
            result.append(tup)
    return result
