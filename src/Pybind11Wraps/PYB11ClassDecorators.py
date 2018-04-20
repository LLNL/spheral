#-------------------------------------------------------------------------------
# Decorators for attributes to classes.
#-------------------------------------------------------------------------------
from functools import wraps as PYB11wraps    # Add PYB11 to screen out in generation
import decorator as PYB11decorator           # To preserve wrapped functions args

#-------------------------------------------------------------------------------
# Singleton
#-------------------------------------------------------------------------------
def PYB11singleton(Cls):
    @PYB11wraps(Cls)
    class Wrapper:
        PYB11singleton = True
        def __init__(self, *args, **kwargs):
            self.instance = Cls(*args, **kwargs)
            return
    return Wrapper

#-------------------------------------------------------------------------------
# Virtual method
#-------------------------------------------------------------------------------
def PYB11virtual(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11virtual = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# Pure virtual method
#-------------------------------------------------------------------------------
def PYB11pure_virtual(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11pure_virtual = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# const method
#-------------------------------------------------------------------------------
def PYB11const(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11const = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# property
#-------------------------------------------------------------------------------
def PYB11property(f, x):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11property = x
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# getter (for property)
#-------------------------------------------------------------------------------
def PYB11getter(f, x):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11getter = x
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# setter (for property)
#-------------------------------------------------------------------------------
def PYB11setter(f, x):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11setter = x
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# static method
#-------------------------------------------------------------------------------
def PYB11static(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11static = True
    return PYB11decorator.decorate(f, wrapper)
