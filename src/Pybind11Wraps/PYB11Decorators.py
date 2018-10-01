#-------------------------------------------------------------------------------
# Decorators for PYB11 generation.
#-------------------------------------------------------------------------------
from functools import wraps as PYB11wraps    # Add PYB11 to screen out in generation
import decorator as PYB11decorator           # To preserve wrapped functions args
import types

#-------------------------------------------------------------------------------
# ignore
#-------------------------------------------------------------------------------
def PYB11ignore(thing):
    thing.PYB11ignore = True
    return thing

#-------------------------------------------------------------------------------
# template (for now just cause the method/class to be ignored)
#-------------------------------------------------------------------------------
class PYB11template:
    def __init__(self, *args):
        self.template = tuple(args)
        return
    def __call__(self, thing):
        if self.template:
            thing.PYB11ignore = True
        else:
            thing.PYB11ignore = False
        thing.PYB11template = self.template
        return thing

#-------------------------------------------------------------------------------
# Singleton (class)
#-------------------------------------------------------------------------------
def PYB11singleton(cls):
    cls.PYB11singleton = True
    return cls

#-------------------------------------------------------------------------------
# namespace (class or method)
#-------------------------------------------------------------------------------
class PYB11namespace:
    def __init__(self, x):
        self.namespace = x
        if self.namespace[:-2] != "::":
            self.namespace += "::"
        return
    def __call__(self, thing):
        thing.PYB11namespace = self.namespace
        return thing

#-------------------------------------------------------------------------------
# pycppname (class or method)
#-------------------------------------------------------------------------------
class PYB11pycppname:
    def __init__(self, x):
        self.x = x
        return
    def __call__(self, thing):
        thing.PYB11cppname = self.x
        thing.PYB11pyname = self.x
        return thing

#-------------------------------------------------------------------------------
# cppname (class or method)
#-------------------------------------------------------------------------------
class PYB11cppname:
    def __init__(self, x):
        self.cppname = x
        return
    def __call__(self, thing):
        thing.PYB11cppname = self.cppname
        return thing

#-------------------------------------------------------------------------------
# pyname (class or method)
#-------------------------------------------------------------------------------
class PYB11pyname:
    def __init__(self, x):
        self.pyname = x
        return
    def __call__(self, thing):
        thing.PYB11pyname = self.pyname
        return thing

#-------------------------------------------------------------------------------
# Virtual (method)
#-------------------------------------------------------------------------------
def PYB11virtual(f):
    f.PYB11virtual = True
    return f

#-------------------------------------------------------------------------------
# Pure virtual (method)
#-------------------------------------------------------------------------------
def PYB11pure_virtual(f):
    f.PYB11pure_virtual = True
    return f

#-------------------------------------------------------------------------------
# const (method)
#-------------------------------------------------------------------------------
def PYB11const(f):
    f.PYB11const = True
    return f

#-------------------------------------------------------------------------------
# static (method)
#-------------------------------------------------------------------------------
def PYB11static(f):
    f.PYB11static = True
    return f

#-------------------------------------------------------------------------------
# implementation -- provide an inline implementation in C++ (only for experts!)
#-------------------------------------------------------------------------------
class PYB11implementation:
    def __init__(self, x):
        self.val = x
        return
    def __call__(self, thing):
        thing.PYB11implementation = self.val
        return thing

#-------------------------------------------------------------------------------
# returnpolicy
#-------------------------------------------------------------------------------
class PYB11returnpolicy:
    def __init__(self, x):
        self.val = x
        return
    def __call__(self, thing):
        thing.PYB11returnpolicy = self.val
        return thing
