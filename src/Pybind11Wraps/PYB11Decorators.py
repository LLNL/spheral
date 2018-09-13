#-------------------------------------------------------------------------------
# Decorators for PYB11 generation.
#-------------------------------------------------------------------------------
from functools import wraps as PYB11wraps    # Add PYB11 to screen out in generation
import decorator as PYB11decorator           # To preserve wrapped functions args
import types

#-------------------------------------------------------------------------------
# ignore
#-------------------------------------------------------------------------------
class PYB11ignore:
    def __init__(self, func):
        self.func = func
        func.PYB11ignore = True
        return
    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

#-------------------------------------------------------------------------------
# Singleton (class)
#-------------------------------------------------------------------------------
def PYB11singleton(cls):
    def __init__(self):
        return
    def __call__(self, cls):
        @PYB11wraps(cls)
        class WrappedCls(cls):
            PYB11singleton = True
        return WrappedCls

#-------------------------------------------------------------------------------
# namespace (class or method)
#-------------------------------------------------------------------------------
def PYB11namespace(cls, x):
    def __init__(self, x):
        self.namespace = x
        return
    def __call__(self, thing):
        if type(thing) == types.ClassType:
            @PYB11wraps(thing)
            class WrappedCls(thing):
                PYB11namespace = self.namespace
            return WrappedCls
        else:
            def wrapper(f, *args, **kwargs):
                return f(*args, **kwargs)
            thing.PYB11namespace = self.namespace
            return PYB11decorator.decorate(thing, wrapper)

#-------------------------------------------------------------------------------
# cppname (class or method)
#-------------------------------------------------------------------------------
class PYB11cppname:
    def __init__(self, x):
        self.cppname = x
        return
    def __call__(self, thing):
        if type(thing) == types.ClassType:
            @PYB11wraps(thing)
            class WrappedCls(thing):
                PYB11cppname = self.cppname
            return WrappedCls
        else:
            def wrapper(f, *args, **kwargs):
                return f(*args, **kwargs)
            thing.PYB11cppname = self.cppname
            return PYB11decorator.decorate(thing, wrapper)

#-------------------------------------------------------------------------------
# Virtual (method)
#-------------------------------------------------------------------------------
def PYB11virtual(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11virtual = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# Pure virtual (method)
#-------------------------------------------------------------------------------
def PYB11pure_virtual(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11pure_virtual = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# const (method)
#-------------------------------------------------------------------------------
def PYB11const(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11const = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# static (method)
#-------------------------------------------------------------------------------
def PYB11static(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11static = True
    return PYB11decorator.decorate(f, wrapper)

#-------------------------------------------------------------------------------
# attribute
#-------------------------------------------------------------------------------
def PYB11readwrite(f):
    def wrapper(f, *args, **kwargs):
        return f(*args, **kwargs)
    f.PYB11readwrite = True
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

