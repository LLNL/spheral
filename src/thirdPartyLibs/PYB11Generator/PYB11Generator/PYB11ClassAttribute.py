from PYB11utils import *
import sys, inspect

#-------------------------------------------------------------------------------
# PYB11generateClassAttributes
#
# Bind the class attributes
#-------------------------------------------------------------------------------
def PYB11GenerateClassAttributes(klass, klassinst, klassattrs, ss):

    PYB11attrs = [x for x in dir(klass) if isinstance(eval("klass.%s" % x), PYB11ClassAttribute) and x in klass.__dict__]
    if PYB11attrs:
        ss("\n    // Properties\n")

    for attrname in PYB11attrs:
        exec('klassinst.%s("%s", klassattrs, ss)' % (attrname, attrname))
    ss("\n")

    return

#-------------------------------------------------------------------------------
# The base class for attributes, most of the implementation
#-------------------------------------------------------------------------------
class PYB11ClassAttribute:

    def __init__(self,
                 static,
                 pyname,
                 cppname,
                 returnpolicy,
                 doc,
                 deftype):
        self.static = static
        self.pyname = pyname
        self.cppname = cppname
        self.returnpolicy = returnpolicy
        self.doc = doc
        self.deftype = deftype
        return

    def __call__(self,
                 name,
                 klassattrs,
                 ss):
        if self.pyname:
            pyname = self.pyname
        else:
            pyname = name
        if self.cppname:
            cppname = self.cppname
        else:
            cppname = name
        if self.static:
            ss('    obj.def_%s_static("%s", ' % (self.deftype, pyname))
        else:
            ss('    obj.def_%s("%s", ' % (self.deftype, pyname))
        ss(("&%(namespace)s%(cppname)s::" % klassattrs) + cppname)
        if self.returnpolicy:
            ss(', py::return_value_policy::%s' % self.returnpolicy)
        if self.doc:
            ss(', ')
            PYB11docstring(self.doc, ss)
        ss(");\n")
        return

#-------------------------------------------------------------------------------
# readwrite
#-------------------------------------------------------------------------------
class PYB11readwrite(PYB11ClassAttribute):

    def __init__(self,
                 static = False,
                 pyname = None,
                 cppname = None,
                 returnpolicy = None,
                 doc = None):
        PYB11ClassAttribute.__init__(self, static, pyname, cppname, returnpolicy, doc, "readwrite")

#-------------------------------------------------------------------------------
# readonly
#-------------------------------------------------------------------------------
class PYB11readonly(PYB11ClassAttribute):

    def __init__(self,
                 static = False,
                 pyname = None,
                 cppname = None,
                 returnpolicy = None,
                 doc = None):
        PYB11ClassAttribute.__init__(self, static, pyname, cppname, returnpolicy, doc, "readonly")

