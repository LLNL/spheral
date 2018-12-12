from PYB11utils import *

import sys, inspect

#-------------------------------------------------------------------------------
# PYB11generateModuleEnums
#
# Bind the classes in the module
#-------------------------------------------------------------------------------
def PYB11generateModuleAttrs(modobj, ss):
    
    # Module attrs
    stuff = [x for x in dir(modobj) if isinstance(eval("modobj.%s" % x), PYB11attr)]
    if stuff:
        ss('\n  // module attributes\n')
        for pyname in stuff:
            inst = eval("modobj.%s" % pyname)
            inst(pyname, ss)

    return

#-------------------------------------------------------------------------------
# Generate an attr
#-------------------------------------------------------------------------------
class PYB11attr:

    def __init__(self,
                 value = None,
                 pyname = None):
        self.value = value
        self.pyname = pyname
        return

    def __call__(self,
                 pyname,
                 ss):
        if self.pyname:
            self.__name__ = self.pyname
        else:
            self.__name__ = pyname
        attrattrs = PYB11attrs(self)
        if self.value:
            attrattrs["cppname"] = self.value
        else:
            attrattrs["cppname"] = pyname
        ss('  m.attr("%(pyname)s") = %(cppname)s;\n' % attrattrs)
        return
