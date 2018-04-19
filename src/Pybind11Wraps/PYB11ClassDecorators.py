#-------------------------------------------------------------------------------
# Decorators for attributes to classes.
#-------------------------------------------------------------------------------
from functools import wraps as PYB11wraps    # Add PYB11 to screen out in generation
import decorator                             # To preserve wrapped functions args

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
    @PYB11wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)
    wrapper.PYB11virtual = True
    return wrapper

#-------------------------------------------------------------------------------
# Pure virtual method
#-------------------------------------------------------------------------------
def PYB11pure_virtual(f):
    @PYB11wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)
    wrapper.PYB11pure_virtual = True
    return wrapper

#-------------------------------------------------------------------------------
# const method
#-------------------------------------------------------------------------------
def PYB11const(f):
    @PYB11wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)
    wrapper.PYB11const = True
    return wrapper

# @decorator.decorator
# def PYB11const(f, *args, **kwargs):
#     f.PYB11const = True
#     return f(*args, **kwargs)
