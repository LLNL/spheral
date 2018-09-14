#-------------------------------------------------------------------------------
# PYB11Generator
#-------------------------------------------------------------------------------
import inspect
import sys
from PYB11utils import *
from PYB11Decorators import *
from PYB11STLmethods import *
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
            ss(("&%(namespace)s%(cppname)s::" % klassattrs) + methattrs["cppname"])
        else:
            argString = ""
            ss(("(%(returnType)s " % methattrs) + ("(%(namespace)s%(cppname)s::*)(" % klassattrs))
            for i, (argType, argName, default) in enumerate(args):
                ss(argType)
                if i < len(args) - 1:
                    ss(", ")
                argString += ', "%s"_a' % argName
                if default:
                    argString += "=%s" % default
            if methattrs["const"]:
                ss((") const) &%(namespace)s%(cppname)s::" % klassattrs) + methattrs["cppname"] + argString)
            else:
                ss((")) &%(namespace)s%(cppname)s::" % klassattrs) + methattrs["cppname"] + argString)
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")

    #...........................................................................
    # Ignore a method
    def ignore(mesh, methattrs, args):
        pass

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
        ss(">()%s" % argString)
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")


        return

    #...........................................................................
    # readonly attribute
    def readonly_class_attribute(aname, attrs, args):
        if attrs["static"]:
            ss('    obj.def_readonly_static("%(pyname)s", ' % methattrs)
        else:
            ss('    obj.def_readonly("%(pyname)s", ' % methattrs)
        ss(("&%(namespace)s%(cppname)s::" % klassattrs) + methattrs["cppname"])
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")

    #...........................................................................
    # readwrite attribute
    def readwrite_class_attribute(aname, attrs, args):
        if attrs["static"]:
            ss('    obj.def_readwrite_static("%(pyname)s", ' % methattrs)
        else:
            ss('    obj.def_readwrite("%(pyname)s", ' % methattrs)
        ss(("&%(namespace)s%(cppname)s::" % klassattrs) + methattrs["cppname"])
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")

    #...........................................................................
    # Property
    def class_property(getmeth, getmethattrs, getargs,
                       setmeth, setmethattrs, setargs):
        ss('    obj.def_property("%(pyname)s", ' % getmethattrs)

        # Get method
        if not getmethattrs["returnType"] is None:
            ss("(" + getmethattrs["returnType"] + " (%(namespace)s%(cppname)s::*)()" % klassattrs)
            if getmethattrs["const"]:
                ss(" const) ")
            else:
                ss(") ")
        ss(("&%(namespace)s%(cppname)s::" % klassattrs) + getmethattrs["cppname"])

        # Set method
        if setmeth:
            if not setmethattrs["returnType"] is None:
                ss(", (" + setmethattrs["returnType"] + " (%(namespace)s%(cppname)s::*)()" % klassattrs)
                if setmethattrs["const"]:
                    ss(" const) ")
                else:
                    ss(") ")
            ss(("&%(namespace)s%(cppname)s::" % klassattrs) + setmethattrs["cppname"])
        doc = inspect.getdoc(meth)
        if doc:
            ss(',\n            "%s"' % doc)
        ss(");\n")

    #...........................................................................
    # Binary operators
    def binary_operator(meth, methattrs, args, op):
        assert len(args) == 1
        argType, argName, default = args[0]
        ss('    obj.def(py::self %s %s);\n' % (op, argType))

    #...........................................................................
    # Reverse binary operators
    def reverse_binary_operator(meth, methattrs, args, op):
        assert len(args) == 1
        argType, argName, default = args[0]
        ss('    obj.def(%s %s py::self);\n' % (argType, op))

    #...........................................................................
    # Unary operators
    def unary_operator(meth, methattrs, args, op):
        assert len(args) == 0
        ss('    obj.def(%s py::self);\n' % op)

    #...........................................................................
    # Tabulate the dispatch for special operations.
    special_operators =  {"__init__": (ignore, ""),

                          "__add__" : (binary_operator, "+"),
                          "__sub__" : (binary_operator, "-"),
                          "__mul__" : (binary_operator, "*"),
                          "__div__" : (binary_operator, "/"),
                          "__mod__" : (binary_operator, "%"),
                          "__and__" : (binary_operator, "&"),
                          "__xor__" : (binary_operator, "^"),
                          "__or__"  : (binary_operator, "|"),
                          
                          "__radd__" : (reverse_binary_operator, "+"),
                          "__rsub__" : (reverse_binary_operator, "-"),
                          "__rmul__" : (reverse_binary_operator, "*"),
                          "__rdiv__" : (reverse_binary_operator, "/"),
                          "__rmod__" : (reverse_binary_operator, "%"),
                          "__rand__" : (reverse_binary_operator, "&"),
                          "__rxor__" : (reverse_binary_operator, "^"),
                          "__ror__"  : (reverse_binary_operator, "|"),
                          
                          "__iadd__" : (binary_operator, "+="),
                          "__isub__" : (binary_operator, "-="),
                          "__imul__" : (binary_operator, "*="),
                          "__idiv__" : (binary_operator, "/="),
                          "__imod__" : (binary_operator, "%="),
                          "__iand__" : (binary_operator, "&="),
                          "__ixor__" : (binary_operator, "^="),
                          "__ior__"  : (binary_operator, "|="),

                          "__neg__"    : (unary_operator, "-"),
                          "__invert__" : (unary_operator, "~"),

                          "__lt__" : (binary_operator, "<"),
                          "__le__" : (binary_operator, "<="),
                          "__eq__" : (binary_operator, "=="),
                          "__ne__" : (binary_operator, "!="),
                          "__gt__" : (binary_operator, ">"),
                          "__ge__" : (binary_operator, ">=")}

    #...........................................................................
    # Iterate over the module classes.
    classes = PYB11classes(modobj)
    for kname, klass in classes:
        objinst = klass()
        klassattrs = PYB11attrs(klass)
        ss("""
  //............................................................................
  // Class %(pyname)s
  {
    py::class_<%(namespace)s%(cppname)s""" % klassattrs)
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

        # Grab all the methods
        allmethods = [(mname, meth) for (mname, meth) in PYB11methods(klass)
                      if not PYB11attrs(meth)["ignore"]]

        # Bind constructors of the class.
        ss("\n    // Constructors\n")
        kills = []
        for i, (mname, meth) in enumerate(allmethods):
            if mname[:6] == "pyinit":
                methattrs = PYB11attrs(meth)
                args = PYB11parseArgs(meth)
                pyinit(meth, methattrs, args)
                kills.append(i)
        for i in reversed(kills):
            del allmethods[i]

        # Bind special operators.
        ss("\n    // Operators\n")
        kills = []
        for i, (mname, meth) in enumerate(allmethods):
            methattrs = PYB11attrs(meth)
            methattrs["returnType"] = eval("objinst." + mname + "()")
            args = PYB11parseArgs(meth)
            if mname in special_operators:
                func, op = special_operators[mname]
                func(meth, methattrs, args, op)
                kills.append(i)
        for i in reversed(kills):
            del allmethods[i]

        # Bind attributes
        ss("\n    // Attributes\n")
        kills = []
        for i, (mname, meth) in enumerate(allmethods):
            methattrs = PYB11attrs(meth)
            methattrs["returnType"] = eval("objinst." + mname + "()")
            args = PYB11parseArgs(meth)
            if methattrs["readonly"]:
                readonly_class_attribute(meth, methattrs, args)
                kills.append(i)
            elif methattrs["readwrite"]:
                readwrite_class_attribute(meth, methattrs, args)
                kills.append(i)
        for i in reversed(kills):
            del allmethods[i]

        # Bind the remaining methods of the class.
        ss("\n    // Methods\n")
        for i, (mname, meth) in enumerate(allmethods):
            methattrs = PYB11attrs(meth)
            methattrs["returnType"] = eval("objinst." + mname + "()")
            args = PYB11parseArgs(meth)
            generic_class_method(meth, methattrs, args)

        # Bind properties
        PYB11GenerateClassProperties(objinst, klassattrs, ss)

        ss("  }\n\n")

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
# PYB11methods
#
# Get the methods to bind from a class
#-------------------------------------------------------------------------------
def PYB11methods(obj):
    return inspect.getmembers(obj, predicate=inspect.ismethod)

