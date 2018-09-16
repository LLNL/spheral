#--------------------------------------------------------------------------------
# PYB11property
#
# Stuff for handling properties in pybind11
#--------------------------------------------------------------------------------
from PYB11utils import *

#--------------------------------------------------------------------------------
# Make a property
#--------------------------------------------------------------------------------
def PYB11property(propname,
                  prop,
                  klassinst,
                  klassattrs,
                  ss):

    getter = prop.fget
    setter = prop.fset
    doc = prop.__doc__
    assert not getter is None

    # Write the getter property formula.
    returnType = eval("klassinst." + getter.__name__ + "()")
    getterattrs = PYB11attrs(getter)
    getterattrs["returnType"] = returnType
    ss('    obj.def_property("%s", (%s ' % (propname, returnType))
    ss('(%(namespace)s%(cppname)s::*)()' % klassattrs)
    if getterattrs["const"]:
        ss(' const)')
    else:
        ss(')')
    ss(' &%(namespace)s%(cppname)s::' % klassattrs + getterattrs["cppname"])

    # setter, if any
    if setter:
        setterattrs = PYB11attrs(setter)
        args = PYB11parseArgs(setter)
        assert len(args) == 1
        ss(', (void (%(namespace)s%(cppname)s::*)' % klassattrs)
        ss('(%s) ' % args[0][0])
        if setterattrs["const"]:
            ss(' const)')
        else:
            ss(')')
        ss(' &%(namespace)s%(cppname)s::' % klassattrs + setterattrs["cppname"])

    # Is there a docstring?
    if doc:
        ss(',\n                     "%s");\n' % doc)
    else:
        ss(');\n')
    return

#-------------------------------------------------------------------------------
# PYB11generateClassProperties
#
# Bind the class properties
#-------------------------------------------------------------------------------
def PYB11GenerateClassProperties(klass, klassinst, klassattrs, ss):
    props = [x for x in dir(klass) if isinstance(eval("klass.%s" % x), property)]
    if props:
        ss("\n    // Properties\n")
        for propname in props:
            prop = eval("klass.%s" % propname)
            PYB11property(propname, prop, klassinst, klassattrs, ss)
        ss("\n")
    return
