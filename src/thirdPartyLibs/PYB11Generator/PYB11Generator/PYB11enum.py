from PYB11utils import *

import sys, inspect

#-------------------------------------------------------------------------------
# PYB11generateModuleEnums
#
# Bind the classes in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleEnums(modobj, ss):
    
    # Module enums
    enums = [x for x in dir(modobj) if isinstance(eval("modobj.%s" % x), PYB11enum)]
    if enums:
        ss("  //..............................................................................\n")
        ss("  // enum types\n")
    for name in enums:
        inst = eval("modobj.%s" % name)
        inst(modobj, ss)
    if enums:
        ss("\n")

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
                 scope,
                 ss,
                 scopeattrs=None):
        if self.name:
            self.__name__ = self.name
        else:
            self.__name__ = self.getInstanceName(scope)
        if self.cppname is None:
            self.cppname = self.__name__
        enumattrs = PYB11attrs(self)
        enumattrs["namespace"] = self.namespace
        enumattrs["cppname"] = self.cppname
        if inspect.isclass(scope):
            klass = scope
            klassattrs = scopeattrs
        else:
            klass = False

        if klass:
            ss('  py::enum_<%(namespace)s%(cppname)s::' % klassattrs)
            ss('%(cppname)s>(obj, "%(pyname)s"' % enumattrs)
        else:
            ss('  py::enum_<%(namespace)s%(cppname)s>(m, "%(pyname)s"' % enumattrs)

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
            ss('    .export_values();\n\n')
        else:
            ss('    ;\n\n')

        return

    def getInstanceName(self, scope):
        for name in dir(scope):
            thing = eval("scope.%s" % name)
            if isinstance(thing, PYB11enum) and thing == self:
                return name
        raise RuntimeError, "PYB11enum: unable to find myself!"
