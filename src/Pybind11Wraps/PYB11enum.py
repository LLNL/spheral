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
        inst(ss)
    return

#-------------------------------------------------------------------------------
# Generate an enum
#-------------------------------------------------------------------------------
class PYB11enumType:

    def __init__(self,
                 values,
                 export_values = False):
        self.values = values
        self.export_values = export_values
        print type(self).__name__
        return

    def __call__(self,
                 ssout,
                 klass = None):
        print dir(self)
        print dir(self.__class__)
        print type(self).__name__
        enumattrs = PYB11attrs(self)
        if klass:
            klassattrs = PYB11attrs(klass)

        ss("  py::enum_<%(namespace)s%(cppname)s>(" % enumattrs)
        if klass:
            ss('%(pyname)s, ' % klassattrs + '"%(pyname)s")\n' % enumattrs)
        else:
            ss('"%(pyname)s")\n' % enumattrs)

        for value in self.values:
            ss('    .value("%s", ' % value)
            if klass:
                ss('%(namespace)s%(cppname)s::' % klassattrs)
            ss('%(cppname)s)\n' % enumattrs)

        if self.export_values:
            ss('    .export_values();\n')
        else:
            ss('    ;')

        return
