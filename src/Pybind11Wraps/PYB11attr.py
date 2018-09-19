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
        ss('  // module attributes\n')
        for name in stuff:
            inst = eval("modobj.%s" % name)
            inst(modobj, ss)

    return

#-------------------------------------------------------------------------------
# Generate an attr
#-------------------------------------------------------------------------------
class PYB11attr:

    def __init__(self,
                 value,
                 name = None):
        self.value = value
        self.name = name
        return

    def __call__(self,
                 scope,
                 ss):
        if self.name:
            self.__name__ = self.name
        else:
            self.__name__ = self.getInstanceName(scope)
        attrattrs = PYB11attrs(self)
        attrattrs["cppname"] = self.value
        ss('  m.attr("%(pyname)s") = %(cppname)s;\n' % attrattrs)
        return

    def getInstanceName(self, scope):
        for name in dir(scope):
            thing = eval("scope.%s" % name)
            if isinstance(thing, PYB11attr) and thing == self:
                return name
        raise RuntimeError, "PYB11attr: unable to find myself!"
