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
                 returnType,
                 getter,
                 setter = None,
                 doc = None,
                 getterraw = None,
                 setterraw = None,
                 getterconst = True,
                 setterconst = False,
                 static = None):
        self.returnType = returnType
        self.getter = getter
        self.setter = setter
        self.doc = doc
        self.getterraw = getterraw
        self.setterraw = setterraw
        self.getterconst = getterconst
        self.setterconst = setterconst
        self.static = static

        assert self.getter or self.getterraw, "PYB11property: must specify getter or getterraw"
        assert self.getter is None or self.getterraw is None, "PYB11property: cannot specify both getter and getterraw"
        assert self.setter is None or self.setterraw is None, "PYB11property: cannot specify both setter and setterraw"
        return

    def __call__(self, propname, klassinst, klassattrs, ss):
        if self.static:
            if self.setter:
                proptype = "_readwrite_static"
            else:
                proptype = "_readonly_static"
        else:
            if self.setter:
                proptype = ""
            else:
                proptype = "_readonly"

        ss('    obj.def_property%s("%s", (%s ' % (proptype, propname, self.returnType))

        # getter code
        if self.getterraw:
            ss(self.getterraw)
        else:
            if self.static:
                ss('(%(namespace)s*)()' % klassattrs)
            else:
                ss('(%(namespace)s%(cppname)s::*)()' % klassattrs)
            if self.getterconst:
                ss(' const)')
            else:
                ss(')')
            ss(' &%(namespace)s%(cppname)s::' % klassattrs + self.getter)

        # setter, if any
        if self.setterraw:
            ss(', %s' % self.setterraw)
        else:
            if self.setter:
                if self.static:
                    ss(', (void (%(namespace)s*)' % klassattrs)
                else:
                    ss(', (void (%(namespace)s%(cppname)s::*)' % klassattrs)
                ss('(%s) ' % self.returnType)
                if self.setterconst:
                    ss(' const)')
                else:
                    ss(')')
                ss(' &%(namespace)s%(cppname)s::' % klassattrs + self.setter)

        # Is there a docstring?
        if self.doc:
            ss(',\n                     "%s");\n' % self.doc)
        else:
            ss(');\n')
        return

#         class Prop:
#             def __init__(self):
#                 pass
#         prop = Prop()
#         prop.__doc__ = self.doc
        
#         # getter
#         if inspect.ismethod(self.getter):
#             prop.fget = self.getter
#         else:
#             getter = '@PYB11cppname("%s")' % self.getter
#             if self.static:
#                 getter += '@PYB11static'
#             getter += '''
# def get_PYB11_%(propname)s(self):
#     return "%(returnType)s"
# ''' % {"propname"   : self.propname,
#        "returnType" : self.returnType}
#             print getter
#             exec(getter)
#             exec("klass.get_PYB11_%s = types.MethodType(get_PYB11_%s, klassinst, klass)" % (self.propname, self.propname))
#             prop.fget = eval("klass.get_PYB11_%s" % self.propname)

#         # setter
#         if inspect.ismethod(self.setter) or self.setter is None:
#             prop.fset = self.setter
#         PYB11GenerateProperty(self.propname, prop, klassinst, klassattrs, ss)
#         return

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
    ss(' &%(namespace)s%(cppname)s::' % klassattrs + getterattrs["cppname"])

    # setter, if any
    if setter:
        setterattrs = PYB11attrs(setter)
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
    PYB11props = [x for x in dir(klass) if isinstance(eval("klass.%s" % x), PYB11property)]
    if props or PYB11props:
        ss("\n    // Properties\n")

    for propname in props:
        prop = eval("klass.%s" % propname)
        PYB11GenerateProperty(propname, prop, klassinst, klassattrs, ss)
    for propname in PYB11props:
        exec('klassinst.%s("%s", klassinst, klassattrs, ss)' % (propname, propname))
    ss("\n")

    return
