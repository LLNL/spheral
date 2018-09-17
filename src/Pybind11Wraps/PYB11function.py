#-------------------------------------------------------------------------------
# PYB11function
#
# Handle functions in pybind11
#-------------------------------------------------------------------------------
from PYB11utils import *
from PYB11property import *
import copy, StringIO

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
                PYB11generateFunction(meth, methattrs, ss)

#-------------------------------------------------------------------------------
# Generate a function definition
#-------------------------------------------------------------------------------
def PYB11generateFunction(meth, methattrs, ss):
    fs = StringIO.StringIO()
    ss = fs.write

    # Arguments
    stuff = inspect.getargspec(meth)
    nargs = len(stuff.args)

    # Return type
    returnType = meth(*tuple(stuff.args))
    methattrs["returnType"] = returnType

    # Write the binding
    ss('  m.def("%(pyname)s", ' % methattrs)

    # If there is an implementation, short-circuit the rest of our work.
    if methattrs["implementation"]:
        ss(methattrs["implementation"])

    elif returnType:
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
    return

#-------------------------------------------------------------------------------
# PYB11functions
#
# Get the functions to bind from a module
#-------------------------------------------------------------------------------
def PYB11functions(modobj):
    return [(name, meth) for (name, meth) in inspect.getmembers(modobj, predicate=inspect.isfunction)
            if name[:5] != "PYB11"]

