from PYB11utils import *

import sys, inspect

#-------------------------------------------------------------------------------
# PYB11generateModuleEnums
#
# Bind the classes in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleEnums(modobj, ss):
    enums = [x for x in dir(modobj) if isinstance(eval("modobj.%s" % x), PYB11enum)]
    for name in enums:
        inst = eval("modobj.%s" % name)
        inst(modobj, ss)
    return

#-------------------------------------------------------------------------------
# Generate an enum
#-------------------------------------------------------------------------------
class PYB11enum:

    def __init__(self,
                 values,
                 name = None,
                 namespace = "",
                 cppname = None,
                 export_values = False,
                 doc = None):
        self.values = values
        self.name = name
        self.namespace = namespace
        if self.namespace and self.namespace[:-2] != "::":
            self.namespace += "::"
        self.cppname = cppname
        self.export_values = export_values
        self.doc = doc
        return

    def __call__(self,
                 modobj,
                 ss,
                 klass = None):
        if self.name:
            self.__name__ = self.name
        else:
            self.__name__ = self.getInstanceName(modobj)
        if self.cppname is None:
            self.cppname = self.__name__
        enumattrs = PYB11attrs(self)
        enumattrs["namespace"] = self.namespace
        enumattrs["cppname"] = self.cppname
        if klass:
            klassattrs = PYB11attrs(klass)

        ss("  py::enum_<%(namespace)s%(cppname)s>(" % enumattrs)
        if klass:
            ss('%(pyname)s, ' % klassattrs + '"%(pyname)s"' % enumattrs)
        else:
            ss('m, "%(pyname)s"' % enumattrs)

        if self.doc:
            ss(', "%s")\n' % self.doc)
        else:
            ss(')\n')

        for value in self.values:
            ss('    .value("%s", ' % value)
            if klass:
                ss('%(namespace)s%(cppname)s::' % klassattrs)
            ss('%(namespace)s%(cppname)s::' % enumattrs + value + ')\n')

        if self.export_values:
            ss('    .export_values();\n')
        else:
            ss('    ;')

        return

    def getInstanceName(self, modobj):
        for name in dir(modobj):
            thing = eval("modobj.%s" % name)
            if isinstance(thing, PYB11enum) and thing == self:
                return name
        raise RuntimeError, "PYB11enum: unable to find myself!"
