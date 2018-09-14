#--------------------------------------------------------------------------------
# PYB11property
#
# Stuff for handling classes in pybind11
#--------------------------------------------------------------------------------
from PYB11utils import *

#--------------------------------------------------------------------------------
# Make a property
#--------------------------------------------------------------------------------
class PYB11property:

    def __init__(self,
                 getter,
                 setter = None,
                 doc = None):
        self.getter = getter
        self.setter = setter
        self.doc = doc
        return

    def __call__(self,
                 propname,
                 klassinst,
                 klassattrs,
                 ss):

        # C++ type of the property
        returnType = eval("klassinst." + self.getter.__name__ + "()")

        # getter attributes
        getterattrs = PYB11attrs(self.getter)
        getterattrs["returnType"] = returnType

        # Write the getter property formula.
        ss('    obj.def_property("%s", (%s ' % (propname, returnType))
        ss('(%(namespace)s%(cppname)s::*)()' % klassattrs)
        if getterattrs["const"]:
            ss(' const)')
        else:
            ss(')')
        ss(' &%(namespace)s%(cppname)s::' % klassattrs + getterattrs["cppname"])

        # setter, if any
        if self.setter:
            setterattrs = PYB11attrs(self.setter)
            setterattrs["returnType"] = returnType
            args = PYB11parseArgs(self.setter)
            assert len(args) == 1
            ss(', (void (%(namespace)s%(cppname)s::*)' % klassattrs)
            ss('(%s) ' % args[0][0])
            if setterattrs["const"]:
                ss(' const)')
            else:
                ss(')')
            ss(' &%(namespace)s%(cppname)s::' % klassattrs + setterattrs["cppname"])

        # Is there a docstring?
        if self.doc:
            ss(',\n                     "%s");\n' % self.doc)
        else:
            ss(');\n')
        return

#-------------------------------------------------------------------------------
# PYB11generateClassProperties
#
# Bind the class properties
#-------------------------------------------------------------------------------
def PYB11GenerateClassProperties(klassinst, klassattrs, ss):
    props = [attr for attr in dir(klassinst) if isinstance(eval("klassinst.%s" % attr), PYB11property)]
    if props:
        ss("\n    // Properties\n")
        for propname in props:
            prop = eval("klassinst.%s" % propname)
            prop(propname, klassinst, klassattrs, ss)
        ss("\n")
    return
