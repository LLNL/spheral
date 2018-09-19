from PYB11utils import *

import inspect

#-------------------------------------------------------------------------------
# PYB11generateModuleEnums
#
# Bind the classes in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleEnums(modobj, ss):
    enums = [x for x in dir(modobj) if isinstance(eval("modobj.%s" % x), PYB11enumType)]
    for name in enums:
        inst = eval("modobj.%s" % name)
        inst(modobj, ss)
    return

#-------------------------------------------------------------------------------
# Generate an enum
#-------------------------------------------------------------------------------
class PYB11enumType:

    def __init__(self,
                 values,
                 name,
                 namespace = "",
                 cppname = None,
                 export_values = False,
                 doc = None):
        self.values = values
        self.__name__ = name
        self.namespace = namespace
        if self.namespace and self.namespace[:-2] != "::":
            self.namespace += "::"
        if cppname:
            self.cppname = cppname
        else:
            self.cppname = name
        self.export_values = export_values
        self.doc = doc
        return

    def __call__(self,
                 modobj,
                 ss,
                 klass = None):
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
