#--------------------------------------------------------------------------------
# PYB11property
#
# Stuff for handling properties in pybind11
#--------------------------------------------------------------------------------
from PYB11Decorators import *
from PYB11utils import *
import inspect, types

#--------------------------------------------------------------------------------
# Class to help making a property, in cases where python's ordinary "property"
# method isn't what we want.
#--------------------------------------------------------------------------------
class PYB11property:

    def __init__(self,
                 returnType = None,
                 getter = None,
                 setter = None,
                 doc = None,
                 getterraw = None,
                 setterraw = None,
                 getterconst = True,
                 setterconst = False,
                 static = None,
                 returnpolicy = None):
        self.returnType = returnType
        self.getter = getter
        self.setter = setter
        self.doc = doc
        self.getterraw = getterraw
        self.setterraw = setterraw
        self.getterconst = getterconst
        self.setterconst = setterconst
        self.static = static
        self.returnpolicy = returnpolicy

        assert self.getter is None or self.getterraw is None, "PYB11property: cannot specify both getter and getterraw"
        assert self.setter is None or self.setterraw is None, "PYB11property: cannot specify both setter and setterraw"
        return

    def __call__(self, propname, klassattrs, ss):
        if not self.getterraw and self.getter is None:
            self.getter = propname

        if self.static:
            if self.setter:
                proptype = "_readwrite_static"
            else:
                proptype = "_readonly_static"
        else:
            if self.setter or self.setterraw:
                proptype = ""
            else:
                proptype = "_readonly"

        ss('    obj.def_property%s("%s", ' % (proptype, propname))

        # getter code
        if self.getterraw:
            ss(self.getterraw)
        else:
            if self.returnType:
                ss('(%s ' % self.returnType)
                if self.static:
                    ss('(%(namespace)s*)()' % klassattrs)
                else:
                    ss('(%(namespace)s%(cppname)s::*)()' % klassattrs)
                    if self.getterconst:
                        ss(' const')
                ss(') ')
            ss('&%(namespace)s%(cppname)s::' % klassattrs + self.getter)

        # setter, if any
        if self.setterraw:
            ss(', %s' % self.setterraw)
        else:
            if self.setter:
                ss(', ')
                if self.returnType:
                    if self.static:
                        ss('(void (%(namespace)s*)' % klassattrs)
                    else:
                        ss('(void (%(namespace)s%(cppname)s::*)' % klassattrs)
                    ss('(%s) ' % self.returnType)
                    if self.setterconst:
                        ss(' const')
                    ss(')')
                ss(' &%(namespace)s%(cppname)s::' % klassattrs + self.setter)

        # Is there a return policy?
        if self.returnpolicy:
            ss(", py::return_value_policy::" + self.returnpolicy)

        # Is there a docstring?
        if self.doc:
            ss(", ")
            PYB11docstring(self.doc, ss)
        ss(');\n')
        return

#--------------------------------------------------------------------------------
# Make a property
#--------------------------------------------------------------------------------
def PYB11GenerateProperty(propname,
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
    getterattrs["namespace"] = "%(namespace)s" % klassattrs
    getterattrs["classcppname"] = "%(cppname)s" % klassattrs
    if getterattrs["protected"]:
        getterattrs["classcppname"] = "PYB11Publicist" + getterattrs["classcppname"]

    # What kind of property do we have (readwrite, readonly, static, etc.)?
    if getterattrs["static"]:
        if setter:
            proptype = "_readwrite_static"
        else:
            proptype = "_readonly_static"
    else:
        if setter:
            proptype = ""
        else:
            proptype = "_readonly"

    ss('    obj.def_property%s("%s", (%s ' % (proptype, propname, returnType))

    if getterattrs["static"]:
        ss('(%(namespace)s*)()' % klassattrs)
    else:
        ss('(%(namespace)s%(cppname)s::*)()' % klassattrs)
    if getterattrs["const"]:
        ss('const)')
    else:
        ss(')')
    ss(' &%(namespace)s%(classcppname)s::%(cppname)s' % getterattrs)

    # setter, if any
    if setter:
        setterattrs = PYB11attrs(setter)
        setterattrs["namespace"] = "%(namespace)s" % klassattrs
        setterattrs["classcppname"] = "%(cppname)s" % klassattrs
        if setterattrs["protected"]:
            setterattrs["classcppname"] = "PYB11Publicist" + setterattrs["classcppname"]
        args = PYB11parseArgs(setter)
        assert len(args) == 1, "Bad number of arguments to property %s" % propname
        if setterattrs["static"]:
            ss(', (void (%(namespace)s*)' % klassattrs)
        else:
            ss(', (void (%(namespace)s%(cppname)s::*)' % klassattrs)
        ss('(%s) ' % args[0][0])
        if setterattrs["const"]:
            ss(' const)')
        else:
            ss(')')
        ss(' &%(namespace)s%(classcppname)s::%(cppname)s' % setterattrs)

    # Is there a return policy?
    if getterattrs["returnpolicy"]:
        ss(", py::return_value_policy::%(returnpolicy)s" % getterattrs)

    # Is there a docstring?
    if doc:
        ss(", ")
        PYB11docstring(doc, ss)
    ss(');\n')
    return

#-------------------------------------------------------------------------------
# PYB11generateClassProperties
#
# Bind the class properties
#-------------------------------------------------------------------------------
def PYB11GenerateClassProperties(klass, klassinst, klassattrs, ss):

    props = [x for x in dir(klass) if isinstance(eval("klass.%s" % x), property) and x in klass.__dict__]
    PYB11props = [x for x in dir(klass) if isinstance(eval("klass.%s" % x), PYB11property) and x in klass.__dict__]
    if props or PYB11props:
        ss("\n    // Properties\n")

    for propname in props:
        prop = eval("klass.%s" % propname)
        PYB11GenerateProperty(propname, prop, klassinst, klassattrs, ss)
    for propname in PYB11props:
        exec('klassinst.%s("%s", klassattrs, ss)' % (propname, propname))
    ss("\n")

    return
