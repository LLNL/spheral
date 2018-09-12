#-------------------------------------------------------------------------------
# PYB11Generator
#-------------------------------------------------------------------------------
import inspect
import sys
from PYB11ClassDecorators import *
from PYB11FunctionDecorators import *
from PYB11STLmethods import *

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
            for i, (argType, argName, default) in enumerate(PYB11parseArgs(args)):
                ss(argType)
                if i < nargs - 1:
                    ss(", ")
            ss(")) &%s" % name)
            for argType, argName, default in PYB11parseArgs(args):
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
# PYB11generateModuleClasses
#
# Bind the classes in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleClasses(modobj, ss):

    #...........................................................................
    # Generate a generic class method spec.
    def generic_class_method(meth, methattrs, args):
        ss('    obj.def("%(pyname)s", ' % methattrs)
        if methattrs["returnType"] is None:
            ss(("&%(cppname)s::" % klassattrs) + methattrs["cppname"])
        else:
            argString = ""
            ss(("(%(returnType)s " % methattrs) + ("(%(cppname)s::*)(" % klassattrs))
            for i, (argType, argName, default) in enumerate(args):
                ss(argType)
                if i < len(args) - 1:
                    ss(", ")
                argString += ', "%s"_a' % argName
                if default:
                    argString += "=%s" % default
            if methattrs["const"]:
                ss((") const) &%(cppname)s::" % klassattrs) + methattrs["cppname"] + argString)
            else:
                ss((")) &%(cppname)s::" % klassattrs) + methattrs["cppname"] + argString)
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")

    #...........................................................................
    # Ignore a method
    def ignore(mesh, methattrs, args):
        pass

    #...........................................................................
    # readwrite attribute
    def readwrite_class_attribute(aname, attrs, args):
        ss('    obj.def_readwrite("%(pyname)s", ' % methattrs)
        ss(("&%(cppname)s::" % klassattrs) + methattrs["cppname"])
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")

    #...........................................................................
    # Property
    def class_property(meth, methattrs, args):

        ss('    obj.def("%(pyname)s", ' % methattrs)
        if methattrs["returnType"] is None:
            ss(("&%(cppname)s::" % klassattrs) + methattrs["cppname"])
        else:
            argString = ""
            ss(("(%(returnType)s " % methattrs) + ("(%(cppname)s::*)(" % klassattrs))
            for i, (argType, argName, default) in enumerate(args):
                ss(argType)
                if i < len(args) - 1:
                    ss(", ")
                argString += ', "%s"_a' % argName
                if default:
                    argString += "=%s" % default
            if methattrs["const"]:
                ss((") const) &%(cppname)s::" % klassattrs) + methattrs["cppname"] + argString)
            else:
                ss((")) &%(cppname)s::" % klassattrs) + methattrs["cppname"] + argString)
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")

    #...........................................................................
    # pyinit<>
    def pyinit(meth, methattrs, args):
        ss("    obj.def(py::init<")
        argString = ""
        for i, (argType, argName, default) in enumerate(args):
            if i < len(args) - 1:
                ss("%s, " % argType)
            else:
                ss("%s" % argType)
            argString += ', "%s"_a' % argName
            if default:
                argString += "=%s" % default
        ss(">()%s);\n" % argString)
        return

    # Some methods we want to skip.
    ignores = ["__init__"]

    # Iterate over the module classes.
    classes = PYB11classes(modobj)
    for kname, klass in classes:
        objinst = klass()
        klassattrs = PYB11attrs(klass)
        ss("""
  //............................................................................
  // Class %(pyname)s
  {
    py::class_<%(cppname)s""" % klassattrs)
        if klassattrs["singleton"]:
            ss(", std::unique_ptr<RestartRegistrar, py::nodelete>")
        ss('> obj(m, "%(pyname)s");\n' % klassattrs)

        # Is there a doc string?
        if inspect.getdoc(klass):
            ss("    obj.doc() = ")
            for i, line in enumerate(inspect.getdoc(klass).split('\n')):
                if i > 0:
                    ss("            ")
                ss('"%s"\n' % line);
            ss("  ;\n")

        # Bind methods of the class.
        for mname, meth in PYB11methods(klass):
            if mname not in ignores:
                methattrs = PYB11attrs(meth)
                methattrs["returnType"] = eval("objinst." + mname + "()")
                args = PYB11parseArgs(meth)
                if mname[:6] == "pyinit":
                    pyinit(meth, methattrs, args)
                elif methattrs["readwrite"]:
                    readwrite_class_attribute(meth, methattrs, args)
                else:
                    generic_class_method(meth, methattrs, args)

        ss("  }\n")

#-------------------------------------------------------------------------------
# PYB11generateModuleSTL
#
# Bind the STL containers in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleSTL(modobj, ss):
    stuff = PYB11STLobjs(modobj)
    for (name, obj) in stuff:
        ss("  ")
        obj(modobj, ss, name)
    ss("\n")
    return

#-------------------------------------------------------------------------------
# PYB11parseArgs
#
# Return (argType, argName, <default_value>)
#-------------------------------------------------------------------------------
def PYB11parseArgs(meth):
    stuff = inspect.getargspec(meth)
    result = []
    if stuff.defaults:
        nargs = len(stuff.defaults)
        for argName, val in zip(stuff.args[-nargs:], stuff.defaults):
            if isinstance(val, tuple):
                assert len(val) == 2
                argType, default = val
            else:
                argType, default = val, None
            result.append((argType, argName, default))
    return result

#-------------------------------------------------------------------------------
# PYB11functions
#
# Get the functions to bind from a module
#-------------------------------------------------------------------------------
def PYB11functions(modobj):
    return [(name, meth) for (name, meth) in inspect.getmembers(modobj, predicate=inspect.isfunction)
            if name[:5] != "PYB11"]

#-------------------------------------------------------------------------------
# PYB11classes
#
# Get the classes to bind from a module
#-------------------------------------------------------------------------------
def PYB11classes(modobj):
    return [(name, cls) for (name, cls) in inspect.getmembers(modobj, predicate=inspect.isclass)
            if name[:5] != "PYB11"]

#-------------------------------------------------------------------------------
# PYB11STLobjs
#
# Get the STL objects to bind from a module
#-------------------------------------------------------------------------------
def PYB11STLobjs(modobj):
    return [(name, obj) for (name, obj) in inspect.getmembers(modobj)
            if name[:5] != "PYB11" and
            (isinstance(obj, PYB11_bind_vector) or
             isinstance(obj, PYB11_bind_map))]

#-------------------------------------------------------------------------------
# PYB11methods
#
# Get the methods to bind from a class
#-------------------------------------------------------------------------------
def PYB11methods(obj):
    return inspect.getmembers(obj, predicate=inspect.ismethod)

#-------------------------------------------------------------------------------
# PYB11attrs
#
# Read the possible PYB11 generation attributes from the obj
#-------------------------------------------------------------------------------
def PYB11attrs(obj):
    d = {"pyname"       : obj.__name__,
         "cppname"      : obj.__name__,
         "namespace"    : None,
         "singleton"    : False,
         "virtual"      : False,
         "pure_virtual" : False,
         "const"        : False,
         "static"       : False,
         "property"     : None,           # Property
         "getter"       : None,           # Property
         "setter"       : None,           # Property
         "readwrite"    : False,          # Attribute
         "readonly"     : False}          # Attribute
    for key in d:
        if hasattr(obj, "PYB11" + key):
            d[key] = eval("obj.PYB11%s" % key)
    return d
